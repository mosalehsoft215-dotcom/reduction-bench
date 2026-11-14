"""Coordinated final-two-acquisition reconstruction from full physical beliefs.

No policy file is read. This changes two comparisons together, escaping a
single-comparison local optimum while checking the complete affected clock.
"""
import copy
import numpy as np
from scipy.special import xlogy
from policy_exchange import motion_seconds

def rebuild_tails(tree,actions,K0,K1,catalog,cfg,limits,shortlist=6):
    S=K0.shape[1];truth=np.repeat(np.arange(12),S//12)
    ids={int(row[1]):i for i,row in enumerate(catalog)};lookup={tuple(a):i for i,a in enumerate(actions.tolist())}
    positions=np.asarray([ids[int(a[0])] for a in actions]);movement=np.array([[motion_seconds(a,b,cfg) for b in catalog] for a in catalog])
    costs=actions[:,2]*cfg['exposure_unit_s']+cfg['readout_s_per_spectrum'];row0=K0.sum(axis=-1);row1=K1.sum(axis=-1);lost=K0[...,2]+K1[...,2]
    bounds=np.array([limits[k] for k in ('worst_error_max','mean_error_max','worst_memory_loss_max','mean_memory_loss_max')])
    budget=limits['time_budget_s'];maximum=limits['max_policy_depth'];changed=0
    def score(r,l):
        return (np.maximum(r/bounds[0]-.95,0)**2+np.maximum(l/bounds[2]-.95,0)**2).sum(axis=-1)+S*np.maximum(r.mean(axis=-1)/bounds[1]-.95,0)**2+S*np.maximum(l.mean(axis=-1)/bounds[3]-.95,0)**2+1e-7*(r.mean(axis=-1)/bounds[1]+l.mean(axis=-1)/bounds[3])
    def terminal_options(w):
        a=np.einsum('asf,sf->as',row0,w);b=np.einsum('asf,sf->as',row1,w)
        da=a.reshape(len(actions),12,-1).sum(axis=-1).argmax(axis=-1);db=b.reshape(len(actions),12,-1).sum(axis=-1).argmax(axis=-1)
        r=a*(da[:,None]!=truth)+b*(db[:,None]!=truth);l=np.einsum('asf,sf->as',lost,w)
        return a,b,da,db,r,l
    def reconstruct():
        info={};nodes=[]
        def back(node):
            if 'domain' in node:
                e=np.broadcast_to((truth!=node['domain'])[:,None],(S,3)).astype(float);l=np.broadcast_to([0.,0.,1.],(S,3));return e,l
            e0,l0=back(node['children'][0]);e1,l1=back(node['children'][1]);k=lookup[tuple(node['action'])]
            e=np.einsum('sfg,sg->sf',K0[k],e0)+np.einsum('sfg,sg->sf',K1[k],e1)
            l=np.einsum('sfg,sg->sf',K0[k],l0)+np.einsum('sfg,sg->sf',K1[k],l1);info[id(node)]=(e,l);return e,l
        e,l=back(tree);risk=e[:,0];loss=l[:,0]
        def forward(node,w,depth,previous,clock):
            if 'domain' in node:return
            if depth==maximum-2:nodes.append((node,w,previous,clock));return
            k=lookup[tuple(node['action'])];pos=positions[k];end=clock+costs[k]+(movement[previous,pos] if previous is not None else 0.)
            for child,K in zip(node['children'],(K0[k],K1[k])):forward(child,np.einsum('sf,sfg->sg',w,K),depth+1,pos,end)
        w=np.zeros((S,3));w[:,0]=1.;forward(tree,w,0,None,0.)
        return info,nodes,risk,loss
    info,nodes,risk,loss=reconstruct()
    for number,(node,w,previous,clock) in enumerate(nodes):
        eold,lold=info[id(node)];outside_r=risk-np.einsum('sf,sf->s',w,eold);outside_l=loss-np.einsum('sf,sf->s',w,lold)
        a,b,da,db,rr,ll=terminal_options(w)
        def entropy(p):
            d=p.reshape(len(actions),12,-1).sum(axis=-1);m=d.sum(axis=-1);return xlogy(m,m)-xlogy(d,d).sum(axis=-1)
        first=entropy(a)+entropy(b);ends=clock+costs+(movement[previous,positions] if previous is not None else 0.)
        first[ends>budget-np.min(costs)]=np.inf
        # Neighboring thresholds at one view should not consume the whole
        # planning horizon. Keep distinct optical comparisons in the shortlist.
        candidate=[];seen=set();order=np.argsort(first)
        for j in order:
            if not np.isfinite(first[j]):break
            key=tuple(actions[j,:2])
            if key not in seen:candidate.append(int(j));seen.add(key)
            if len(candidate)>=shortlist:break
        leading=set(tuple(actions[j,:2]) for j in candidate[:2]);dwells=set()
        for j in order:
            if not np.isfinite(first[j]):break
            key=tuple(actions[j,:3])
            if tuple(actions[j,:2]) in leading and key not in dwells:
                candidate.append(int(j));dwells.add(key)
            if len(dwells)>=4:break
        candidate.append(lookup[tuple(node['action'])])
        for laser in np.unique(catalog[:,2]):
            subset=np.flatnonzero(catalog[positions,2]==laser)
            candidate.extend(subset[np.argsort(first[subset])[:3]])
        best=float(score(risk,loss));chosen=None
        for k in dict.fromkeys(candidate):
            if not np.isfinite(first[k]):continue
            branches=[];parts=[]
            for kernel in (K0[k],K1[k]):
                wb=np.einsum('sf,sfg->sg',w,kernel);x,y,dx,dy,r,l=terminal_options(wb)
                values=score(outside_r[None,:]+r,outside_l[None,:]+l)
                final=ends[k]+costs+movement[positions[k],positions];values[final>budget+1e-9]=np.inf
                j=int(values.argmin())
                masses=wb.sum(axis=1);label=int(masses.reshape(12,-1).sum(axis=1).argmax());stop_r=masses*(truth!=label);stop_l=wb[:,2]
                if score(outside_r+stop_r,outside_l+stop_l)<=values[j]:branches.append({'domain':label});parts.append((stop_r,stop_l))
                else:branches.append({'action':actions[j].tolist(),'children':[{'domain':int(dx[j])},{'domain':int(dy[j])}]});parts.append((r[j],l[j]))
            rnew=outside_r+parts[0][0]+parts[1][0];lnew=outside_l+parts[0][1]+parts[1][1];value=float(score(rnew,lnew))
            if value<best-1e-10:best=value;chosen={'action':actions[k].tolist(),'children':branches};best_r=rnew;best_l=lnew
        if chosen is not None:
            node.clear();node.update(chosen);risk=best_r;loss=best_l;changed+=1
        if number%32==0:print('Tail reconstruction',number,changed,float(risk.mean()),float(risk.max()),flush=True)
    print('Tail reconstruction finished',changed,float(risk.mean()),float(risk.max()),flush=True)
    return tree
