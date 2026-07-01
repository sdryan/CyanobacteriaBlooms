"""
Local (one-at-a-time) normalized sensitivity analysis of the non-dim bloom model.
Elasticity S_p = (dY/dp)(p/Y) via central differences (+/-5%) for outputs:
  B_peak  : peak total bloom concentration
  t_dec   : decay time (bloom falls below 1 ug/L threshold)
  t_ot    : time Ana overtakes Mike
Parameters: delta*, gamma*, qmax*, eta1*(=zeta1*), zeta2*, f0*, Tc.
Goal: check robustness AND confirm the decay-time theorem predictions
(t_dec ~ 1/delta  => elasticity ~ -1; t_dec increases with 1/gamma => positive small).
"""
import numpy as np

base = dict(a1n=1.9929e5, a2n=4.286e5, dnd=0.08929, qnd=0.25, gnd=0.08929,
            e1=1.2217, z1=1.2217, e2=0.0, z2=0.3062, f0n=20.22, Tc=20.0)
rhoC=3.76; t0=1/0.028
L=4e10; M=400; dx=L/M; x=np.linspace(0,L,M+1); dt_=0.008; Tend=80.0
ns=int(Tend/dt_)
a1,b1=L/4.0,3.0*L/8.0
fbump=np.zeros(M+1); m=(x>a1)&(x<b1)
fbump[m]=np.exp(-1.0/(1.0-(-1.0+(2.0/(b1-a1))*(x[m]-a1))**2))  # shape; scaled by f0n

def lap(y):
    d=np.zeros_like(y); d[1:-1]=(y[2:]-2*y[1:-1]+y[:-2])/(dx*dx); return d

def run(p):
    u=np.full(M+1,0.1); w=np.full(M+1,0.1); N=np.full(M+1,0.01); P=np.full(M+1,0.01)
    f=p['f0n']*fbump
    T=[]; pU=[]; pW=[]
    for k in range(1,ns+1):
        t=k*dt_; src=f if t<=p['Tc'] else 0.0
        rN=N/(1+N); rP=p['qnd']*P/(1+P)
        u=u+dt_*(p['a1n']*lap(u)+rN*u*(1-u)-p['dnd']*u)
        w=w+dt_*(p['a1n']*lap(w)+rP*w*(1-w)-p['dnd']*w)
        N=N+dt_*(p['a2n']*lap(N)-p['gnd']*N-p['e1']*N*u-p['z1']*N*w+(src if np.ndim(src) else 0))
        P=P+dt_*(p['a2n']*lap(P)-p['gnd']*P-p['e2']*P*u-p['z2']*P*w+(src if np.ndim(src) else 0))
        u[0]=u[1];u[-1]=u[-2];w[0]=w[1];w[-1]=w[-2];N[0]=N[1];N[-1]=N[-2];P[0]=P[1];P[-1]=P[-2]
        if k%20==0:
            T.append(t); pU.append(u.max()); pW.append(w.max())
    T=np.array(T); pU=np.array(pU); pW=np.array(pW); days=T*t0/24
    tot=(pU+pW)*rhoC
    Bpeak=tot.max(); ip=np.argmax(tot)
    post=np.where((days>days[ip])&(tot<1.0))[0]
    tdec=days[post[0]] if len(post) else np.nan
    ot=np.where(pW>pU)[0]; tot_=days[ot[0]] if len(ot) else np.nan
    return Bpeak, tdec, tot_

Y0=run(base)
names=['B_peak','t_dec','t_ot']
print("baseline:", dict(zip(names,[round(v,3) if v==v else v for v in Y0])))

params=['dnd','gnd','qnd','e1','z2','f0n','Tc']
plabel={'dnd':'delta','gnd':'gamma','qnd':'qmax*','e1':'uptake eta1','z2':'uptake zeta2','f0n':'f0','Tc':'Tc'}
# tie e1 and z1 together
elast={n:{} for n in names}
for pr in params:
    Sp={}
    for sign in (+1,-1):
        p=dict(base); p[pr]=base[pr]*(1+sign*0.05)
        if pr=='e1': p['z1']=p['e1']  # keep zeta1=eta1
        Sp[sign]=run(p)
    for j,n in enumerate(names):
        yp,ym=Sp[+1][j],Sp[-1][j]
        if (yp==yp) and (ym==ym) and Y0[j]==Y0[j] and Y0[j]!=0:
            dY=(yp-ym)/(2*0.05*base[pr]); elast[n][pr]=dY*base[pr]/Y0[j]
        else:
            elast[n][pr]=np.nan

print("\nNormalized sensitivities (elasticities):")
hdr="param".ljust(14)+"".join(n.rjust(10) for n in names)
print(hdr)
for pr in params:
    row=plabel[pr].ljust(14)+"".join((f"{elast[n][pr]:+.2f}" if elast[n][pr]==elast[n][pr] else "   n/a").rjust(10) for n in names)
    print(row)

# figure: tornado of t_dec elasticities + B_peak
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(10,3.8))
for a,out,ttl in [(ax[0],'t_dec','Decay time  $t_{dec}$'),(ax[1],'B_peak','Peak bloom  $B_{peak}$')]:
    vals=[elast[out][pr] for pr in params]; labs=[plabel[pr] for pr in params]
    order=np.argsort([abs(v) if v==v else 0 for v in vals])
    vals=[vals[i] for i in order]; labs=[labs[i] for i in order]
    cols=['tab:red' if (v==v and v<0) else 'tab:blue' for v in vals]
    a.barh(range(len(vals)),[v if v==v else 0 for v in vals],color=cols)
    a.set_yticks(range(len(vals))); a.set_yticklabels(labs,fontsize=8)
    a.axvline(0,color='k',lw=0.6); a.set_title(ttl,fontsize=10); a.set_xlabel('elasticity')
plt.tight_layout(); plt.savefig('sensitivity_analysis.png',dpi=140)
print("\nsaved sensitivity_analysis.png")
