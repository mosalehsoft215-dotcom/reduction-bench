"""Full-belief measurement planning with a two-comparison horizon.

Small posterior entries are ignored only in candidate ranking. Every child
belief, final risk, material-loss probability and clock is evaluated in full.
No policy or selected-setting list is loaded.
"""
import numpy as np
from scipy.special import xlogy
from adaptive_rates import leg

def lookahead_construct(actions,K0,K1,catalog,cfg,limits,penalty=6.,beta=.1,prefix=2,shortlist=8,calibration_tolerance=.0005):
    S=K0.shape[1];A=len(actions);D=limits['max_policy_depth'];budget=limits['time_budget_s']
    p0=K0.sum(-1);p1=K1.sum(-1);lost=K0[...,2]+K1[...,2];lost[:,:,2]-=1
    move=np.array([[leg(a,b,cfg) for b in catalog] for a in catalog]);dt=actions[:,2]*cfg['exposure_unit_s']+cfg['readout_s_per_spectrum']
    domains=np.arange(S)//(S//12);registrations=(np.arange(S)//cfg['kinetic_population_count'])%3;nodes=0
    def grouped(x,ids,count):return np.column_stack([x[:,ids==j].sum(1) for j in range(count)])
    def ent(x,ids,cal=False,joint_weight=None):
        q=grouped(x,ids,3 if cal else 12);m=q.sum(1)
        value=xlogy(m,m)-xlogy(q,q).sum(1)
        weight=beta if joint_weight is None else joint_weight
        if weight and not cal:value+=weight*(xlogy(m,m)-xlogy(x,x).sum(1))
        return value
    def rank(w,depth,previous,clock,cal=False):
        mass=w.sum();active=w.sum(1)>mass*1e-10;wa=w[active]
        a=np.maximum(np.einsum('asf,sf->as',p0[:,active],wa),0)
        b=np.maximum(np.einsum('asf,sf->as',p1[:,active],wa),0)
        loss=np.einsum('asf,sf->a',lost[:,active],wa)
        end=clock+dt+(move[previous,actions[:,0]] if previous is not None else 0)
        cost=penalty*loss+.0005*(end-clock)*mass
        if cal:value=ent(a,registrations[active],True)+ent(b,registrations[active],True)+.0003*(end-clock)*mass
        elif depth==D-1:value=-grouped(a,domains[active],12).max(1)-grouped(b,domains[active],12).max(1)+cost
        else:
            weight=beta if depth<5 else 0.
            value=ent(a,domains[active],joint_weight=weight)+ent(b,domains[active],joint_weight=weight)+cost
        valid=end<=budget-2.4*(D-depth-1)+1e-9
        if cal:valid&=catalog[actions[:,0],0]==0
        value[~valid]=np.inf
        return value,end,cost
    def grow(w,depth,previous,clock):
        nonlocal nodes
        m=w.sum(1);mass=float(m.sum());domain=np.bincount(domains,weights=m,minlength=12);label=int(domain.argmax())
        if depth==D or mass<1e-13 or mass-domain[label]<1e-7:return {'domain':label}
        q=np.bincount(registrations,weights=m,minlength=3);cal=depth<prefix and mass-q.max()>calibration_tolerance*mass
        values,ends,cost=rank(w,depth,previous,clock,cal);choice=int(values.argmin())
        if not np.isfinite(values[choice]):return {'domain':label}
        if not cal and depth==D-1 and -domain.max()<=values[choice]:return {'domain':label}
        if shortlist and not cal and mass>1e-3 and w[:,2].sum()<.95*mass and (depth<5 or depth==D-2):
            candidates=[];seen=set()
            for k in np.argsort(values):
                if not np.isfinite(values[k]):break
                key=tuple(actions[k,:2])
                if key not in seen:candidates.append(int(k));seen.add(key)
                if len(candidates)>=shortlist:break
            for laser in (532,633):
                subset=np.flatnonzero(catalog[actions[:,0],2]==laser)
                k=int(subset[np.argmin(values[subset])])
                if np.isfinite(values[k]):candidates.append(k)
            best=np.inf
            for k in dict.fromkeys(candidates):
                estimate=cost[k]
                for kernel in (K0[k],K1[k]):
                    child=np.einsum('sf,sfg->sg',w,kernel)
                    if child.sum()<1e-13:continue
                    future,_,_=rank(child,depth+1,int(actions[k,0]),float(ends[k]))
                    score=float(future.min())
                    if not np.isfinite(score):
                        cm=child.sum(1)
                        score=-np.bincount(domains,weights=cm,minlength=12).max() if depth==D-2 else float(ent(cm[None,:],domains)[0])
                    estimate+=score
                if estimate<best:best=estimate;choice=k
        nodes+=1
        if nodes%32==0:print('Physical lookahead',nodes,depth,flush=True)
        i=int(actions[choice,0]);end=float(ends[choice])
        return dict(action=actions[choice].tolist(),children=[grow(np.einsum('sf,sfg->sg',w,k),depth+1,i,end) for k in (K0[choice],K1[choice])])
    initial=np.zeros((S,3));initial[:,0]=1/S
    return grow(initial,0,None,0.)
