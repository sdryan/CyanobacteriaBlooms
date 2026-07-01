"""
Calendar-aligned comparison of the model's total-bloom trajectory against
documented western Lake Erie bloom phenology (citable public sources).

Observational anchors (all citable):
  - Bloom peak timing: early-to-mid August. Peak surface Chl-a 6784 ug/L on
    10 Aug 2015 (WE08); 532 ug/L on 4 Aug 2017 (WE09); 593 ug/L on 5 Aug 2019
    (WE09). Source: Boegehold et al., ESSD 15:3853 (2023); NCEI 10.25921/11da-3x54.
  - 2011 succession: Microcystis dominant through mid-July-late Aug, then a
    secondary nitrogen-fixing Anabaena bloom as N depleted (late Aug-Sept).
    Source: Michalak et al., PNAS 110:6448 (2013); IJC Cyanobacterial Blooms report.
This is a PHENOLOGICAL alignment (timing/sequence), not a least-squares calibration.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import datetime as dt

# ---- run the model (base case, same nondim params as algae_sim.py) ----
kmax=0.028; qmax=0.007; l0=10.0; rhoC=3.76
alpha1=155.0; alpha2=333.33; kn=530.0
delta=0.0025; gamma=0.0025; eta1=0.0091; zeta1=eta1; eta2=0.0; zeta2=0.01*0.228
f0=300.0; t0=1.0/kmax
a1n=alpha1/(l0*l0*kmax); a2n=alpha2/(l0*l0*kmax); dnd=delta/kmax; qnd=qmax/kmax
gnd=gamma/kmax; e1=eta1*rhoC/kmax; z1=zeta1*rhoC/kmax; e2=eta2*rhoC/kmax
z2=zeta2*rhoC/kmax; f0n=f0/(kmax*kn)

L=4e10; M=800; dx=L/M; x=np.linspace(0,L,M+1); dt_=0.004; Tc=20.0; Tend=80.0
ns=int(Tend/dt_)
a1,b1=L/4.0,3.0*L/8.0; f=np.zeros(M+1); m=(x>a1)&(x<b1)
f[m]=f0n*np.exp(-1.0/(1.0-(-1.0+(2.0/(b1-a1))*(x[m]-a1))**2))
u=np.full(M+1,0.1); w=np.full(M+1,0.1); N=np.full(M+1,0.01); P=np.full(M+1,0.01)
def lap(y):
    d=np.zeros_like(y); d[1:-1]=(y[2:]-2*y[1:-1]+y[:-2])/(dx*dx); return d
T=[]; pU=[]; pW=[]
for k in range(1,ns+1):
    t=k*dt_; src=f if t<=Tc else 0.0*f
    rN=N/(1+N); rP=qnd*P/(1+P)
    u2=u+dt_*(a1n*lap(u)+rN*u*(1-u)-dnd*u)
    w2=w+dt_*(a1n*lap(w)+rP*w*(1-w)-dnd*w)
    N2=N+dt_*(a2n*lap(N)-gnd*N-e1*N*u-z1*N*w+src)
    P2=P+dt_*(a2n*lap(P)-gnd*P-e2*P*u-z2*P*w+src)
    u,w,N,P=u2,w2,N2,P2
    u[0]=u[1];u[-1]=u[-2];w[0]=w[1];w[-1]=w[-2];N[0]=N[1];N[-1]=N[-2];P[0]=P[1];P[-1]=P[-2]
    if k%25==0:
        T.append(t); pU.append(u.max()); pW.append(w.max())
T=np.array(T); pU=np.array(pU); pW=np.array(pW)
days=T*t0/24.0
totB=(pU+pW)*rhoC  # ug/L peak total

# characteristic times (days)
def first_cross(arr, thr, days):
    idx=np.where(arr>thr)[0]
    return days[idx[0]] if len(idx) else np.nan
t_on = first_cross(totB,1.0,days)
i_peak=np.argmax(totB); t_peak=days[i_peak]
# decay: last day above 1 ug/L after peak
post=np.where((days>t_peak)&(totB<1.0))[0]
t_dec=days[post[0]] if len(post) else np.nan
# Ana overtakes Mike
ot=np.where(pW>pU)[0]
t_ot=days[ot[0]] if len(ot) else np.nan
print(f"MODEL (days from runoff onset): t_on={t_on:.1f}, t_peak={t_peak:.1f}, "
      f"t_ot(Ana>Mike)={t_ot}, t_dec={t_dec:.1f}")

# ---- calendar alignment: runoff onset -> July 1 ----
anchor = dt.date(2011,7,1)
def to_date(d): return anchor + dt.timedelta(days=float(d))
model_dates = [to_date(d) for d in days]

# observed phenology anchors (date, label)
obs_peak = dt.date(2011,8,10)      # documented peak window (10 Aug 2015 exemplar)
obs_onset= dt.date(2011,7,15)      # blooms typically detectable mid-July
obs_ana  = dt.date(2011,8,25)      # 2011 Anabaena secondary bloom, late Aug
obs_end  = dt.date(2011,9,30)      # decline through Sept-Oct

# ---- figure ----
fig,ax=plt.subplots(figsize=(9,4.6))
ax.plot(model_dates, totB, 'm-', lw=2.2, label='Model: peak total bloom')
ax.plot(model_dates, pU*rhoC, color='tab:red', lw=1, ls='--', label='Model: Mike')
ax.plot(model_dates, pW*rhoC, color='tab:blue', lw=1, ls='--', label='Model: Ana')
ax.set_ylim(0,6.6)
ax.axhline(1.0,color='gray',ls=':',lw=0.8)
for d,lab,c in [(obs_onset,'obs. onset\n(mid-Jul)','green'),
                (obs_peak,'obs. peak\n(early-mid Aug)','black'),
                (obs_ana,'2011 Anabaena\nsecondary (late Aug)','blue'),
                (obs_end,'obs. decline\n(Sep-Oct)','gray')]:
    ax.axvline(d,color=c,ls='-',lw=1,alpha=0.5)
    ax.text(d,5.9,lab,rotation=90,va='top',ha='right',fontsize=7,color=c)
ax.set_ylabel('peak total bloom (ug/L)'); ax.set_xlabel('calendar date (runoff onset set to 1 Jul)')
ax.set_title('Model bloom phenology vs documented western Lake Erie timing')
ax.legend(fontsize=8,loc='upper right'); fig.autofmt_xdate()
plt.tight_layout(); plt.savefig('model_vs_lakeerie.png',dpi=140)
print("saved model_vs_lakeerie.png")

# ---- comparison table ----
print("\nPHENOLOGY COMPARISON (calendar)")
print(f"  onset : model {to_date(t_on)}   | obs ~ {obs_onset} (mid-Jul)")
print(f"  peak  : model {to_date(t_peak)} | obs ~ {obs_peak} (early-mid Aug)")
if not np.isnan(t_ot):
    print(f"  Ana>Mike: model {to_date(t_ot)} | obs ~ {obs_ana} (2011 late Aug)")
else:
    print(f"  Ana>Mike: model NONE in season | obs ~ {obs_ana} (2011 only)")
print(f"  decay : model {to_date(t_dec)}  | obs ~ {obs_end} (Sep-Oct decline)")
