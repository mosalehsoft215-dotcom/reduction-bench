import numpy as np
from heterogeneity import predict,KEYS,SCALE

def forecast_all(x,case,sens=False):
    rows=[];final=[];grad=[]
    for run in case['forecasts']:
        p=predict(x,run,sens=sens,full=True,rtol=3e-10,atol=2e-13);snapshots=[]
        for j,t in enumerate(run['times']):
            snapshots.append(dict(time=t,sample_temperature=float(p['sample_temperature'][j]),coverages=p['coverages'][j].tolist(),gas_yields=p['yields'][j].tolist(),gas_rates=p['flux'][j].tolist(),pair_fractions=p['pair_fractions'][j].tolist()))
        rows.append(dict(id=run['id'],snapshots=snapshots));final.append(p['yields'][-1,1]);grad.append(p['dyields'][-1,1])
    gradient=np.asarray(grad)*SCALE[None,:] if sens else np.zeros((len(final),len(KEYS)))
    return rows,np.asarray(final),gradient

def forecast(x,case,sens=False):
    rows,values,grad=forecast_all(x,case,sens)
    assert len(values)==2
    return rows,float(values[0]-values[1]),grad[0]-grad[1]
