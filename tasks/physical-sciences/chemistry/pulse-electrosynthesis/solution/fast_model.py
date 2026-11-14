# Conservative splitting solver for reference optimization.
import json, time
import numpy as np
from numba import njit

@njit
def thomas(a,b,c,r):
    n=len(b); cc=c.copy(); rr=r.copy()
    cc[0]/=b[0]; rr[0]/=b[0]
    for i in range(1,n):
        den=b[i]-a[i]*cc[i-1]
        cc[i]/=den; rr[i]=(rr[i]-a[i]*rr[i-1])/den
    for i in range(n-2,-1,-1): rr[i]-=cc[i]*rr[i+1]
    return rr

@njit
def decay(r,h,kd,kh,dx):
    # Concentrations in mol/L; kd in L/mol/s; volume integrals normalized by C0*L.
    ex=np.exp(-kh*h)
    new=r*ex/(1+2*kd*r*(1-ex)/kh)
    integral=np.log1p(2*kd*r*(1-ex)/kh)/(2*kd)
    side=kh*np.sum(integral)*dx
    desired=(np.sum(r-new)*dx-side)/2
    return new, desired, side

@njit
def trajectory(schedule, pars, n=48, dtmax=.003):
    L,c0,ds,dr,kdim,kh,k1,k2,shift,ads_scale,des_scale=pars
    dx=L/n
    s=np.full(n,c0); r=np.zeros(n)
    p=0.; w=0.; q=0.; theta=0.
    clock=0.; edge=.05; bin_index=0; last_q=0.; bins=np.zeros(120)
    for duration,e in schedule:
        end=clock+duration
        u0=k1*np.exp(-.5*38.9217*(e-shift+1.15))
        v0=k2*np.exp(-.5*38.9217*(e-shift+1.05))
        ka=6*ads_scale/(1+np.exp((e-shift+1.08)/.018))
        kb=des_scale*(.02+10/(1+np.exp(-(e-shift+1.03)/.020)))
        rate=ka+kb; equilibrium=ka/rate
        while clock < end-1e-12:
            dt=min(dtmax,end-clock,edge-clock)
            middle=equilibrium+(theta-equilibrium)*np.exp(-rate*dt/2)
            theta=equilibrium+(theta-equilibrium)*np.exp(-rate*dt)
            u=u0*(1-middle); v=v0*(1-middle)
            # Robin boundary, eliminate interface concentrations over half a cell.
            us=u/(1+u*dx/(2*ds)); vr=v/(1+v*dx/(2*dr))
            r,dp,dw=decay(r,dt/2,kdim,kh,dx); p+=dp; w+=dw
            lam=ds*dt/dx**2
            a=np.full(n,-lam); b=np.full(n,1+2*lam); c=np.full(n,-lam)
            a[0]=0.; c[-1]=0.; b[-1]=1+lam; b[0]=1+lam+dt*us/dx
            s=thomas(a,b,c,s)
            j1=us*s[0]
            lam=dr*dt/dx**2
            a[:]=-lam; b[:]=1+2*lam; c[:]=-lam
            a[0]=0.; c[-1]=0.; b[-1]=1+lam; b[0]=1+lam+dt*vr/dx
            # Surface R = (cell R + j1*dx/(2Dr))/(1+v*dx/(2Dr)).
            # Only the fraction of newly produced R that escapes immediate overreduction enters the cell.
            escape=1/(1+v*dx/(2*dr))
            r[0]+=dt*j1*escape/dx
            r=thomas(a,b,c,r)
            j2=vr*r[0]+j1*(1-escape)
            w+=dt*j2; q+=dt*(j1+j2)
            r,dp,dw=decay(r,dt/2,kdim,kh,dx); p+=dp; w+=dw
            clock+=dt
            if clock >= edge-1e-11:
                bins[bin_index]=(q-last_q)/(c0*L)
                last_q=q; bin_index+=1; edge=.05*(bin_index+1)
    total=c0*L
    consumed=total-np.sum(s)*dx
    metrics=np.array([2*p/total,w/total,q/total,np.sum(s)*dx/total,np.sum(r)*dx/total,
                     2*p/max(consumed,1e-20),theta])

    return metrics,bins

@njit
def simulate(schedule, pars, n=48, dtmax=.003):
    return trajectory(schedule,pars,n,dtmax)[0]
