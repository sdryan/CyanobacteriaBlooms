"""
Python port of algae_1d_ext_nondim_paper_2025.m
Purpose: reproduce the non-dimensional 4-field reaction-diffusion bloom model
and numerically test the analytical decay-time bound and KPP spreading speed.

Fields: u (Mike), w (Ana), N (nitrogen), P (phosphorus) on x in [0,L], no-flux BCs.
Non-dimensional parameters are computed exactly as in the MATLAB script.
"""
import numpy as np

# ---- dimensional params (from MATLAB) ----
kmax = 0.028      # 1/hr max N absorption by Mike
qmax = 0.007      # 1/hr max P absorption by Ana
l0   = 10.0       # micron
rhoC = 3.76       # ug/L char cyano conc
rhoN0= 1874.4
rhoP0= 84.4
alpha1 = 155.0    # algae diffusion (um^2/s)  -- NOTE: per-second while rates per-hour
alpha2 = 333.33   # nutrient diffusion
kn = 530.0; qn = 530.0
delta = 0.0025    # 1/hr cyano mortality
gamma = 0.0025    # 1/hr nutrient degradation
eta1 = 0.0091     # N consumption by Mike
zeta1= eta1       # N consumption by Ana
eta2 = 0.0        # P by Mike (=0)
zeta2= 0.01*0.228 # P by Ana
f0   = 300.0
t0   = 1.0/kmax   # 35.7 hr

# ---- non-dimensionalization (exactly as MATLAB does it) ----
alpha1_nd = alpha1/(l0*l0*kmax)
alpha2_nd = alpha2/(l0*l0*kmax)
delta_nd  = delta/kmax
qmax_nd   = qmax/kmax
gamma_nd  = gamma/kmax
eta1_nd   = eta1*rhoC/kmax
zeta1_nd  = zeta1*rhoC/kmax
eta2_nd   = eta2*rhoC/kmax
zeta2_nd  = zeta2*rhoC/kmax
f0_nd     = f0/(kmax*kn)

print("Non-dim params:")
print(f"  alpha1_nd={alpha1_nd:.4g}, alpha2_nd={alpha2_nd:.4g}")
print(f"  delta_nd ={delta_nd:.4g}, gamma_nd={gamma_nd:.4g}, qmax_nd={qmax_nd:.4g}")
print(f"  eta1_nd={eta1_nd:.4g}, zeta1_nd={zeta1_nd:.4g}, zeta2_nd={zeta2_nd:.4g}, f0_nd={f0_nd:.4g}")

# ---- grid ----
L = 4e10            # nondim domain length (=400 km / l0)
M = 800             # spatial cells (coarser than MATLAB 1e5; diffusion negligible)
dx = L/M
x = np.linspace(0, L, M+1)
dt = 0.004          # nondim time step
Tc = 20.0           # source shut-off (nondim) ~ 30 days
Tend = 80.0         # ~120 days
nsteps = int(Tend/dt)

# source bump on [L/4, 3L/8]
a1, b1 = L/4.0, 3.0*L/8.0
f = np.zeros(M+1)
mask = (x > a1) & (x < b1)
arg = -1.0/(1.0 - (-1.0 + (2.0/(b1-a1))*(x[mask]-a1))**2)
f[mask] = f0_nd*np.exp(arg)

# initial conditions
u = np.full(M+1, 0.1)
w = np.full(M+1, 0.1)
N = np.full(M+1, 0.01)
P = np.full(M+1, 0.01)

beta = 0.0  # competition off (base case)

def lap(y):
    d = np.zeros_like(y)
    d[1:-1] = (y[2:]-2*y[1:-1]+y[:-2])/(dx*dx)
    return d

# diagnostics
times=[]; totU=[]; totW=[]; maxU=[]; maxW=[]; maxN=[]; maxP=[]
N_at_Tc = None; P_at_Tc=None; uTc=None; wTc=None

for k in range(1, nsteps+1):
    t = k*dt
    src = f if t <= Tc else 0.0*f
    rN = N/(1.0+N)
    rP = qmax_nd*P/(1.0+P)
    un = u + dt*(alpha1_nd*lap(u) + rN*u*(1-u) - beta*u*w - delta_nd*u)
    wn = w + dt*(alpha1_nd*lap(w) + rP*w*(1-w) - beta*u*w - delta_nd*w)
    Nn = N + dt*(alpha2_nd*lap(N) - gamma_nd*N - eta1_nd*N*u - zeta1_nd*N*w + src)
    Pn = P + dt*(alpha2_nd*lap(P) - gamma_nd*P - eta2_nd*P*u - zeta2_nd*P*w + src)
    u,w,N,P = un,wn,Nn,Pn
    # no-flux
    u[0]=u[1]; u[-1]=u[-2]; w[0]=w[1]; w[-1]=w[-2]
    N[0]=N[1]; N[-1]=N[-2]; P[0]=P[1]; P[-1]=P[-2]
    if abs(t-Tc) < dt/2:
        N_at_Tc=N.copy(); P_at_Tc=P.copy(); uTc=u.copy(); wTc=w.copy()
    if k % 50 == 0:
        times.append(t); totU.append(u.sum()*dx); totW.append(w.sum()*dx)
        maxU.append(u.max()); maxW.append(w.max()); maxN.append(N.max()); maxP.append(P.max())

times=np.array(times); maxU=np.array(maxU); maxW=np.array(maxW)
maxN=np.array(maxN); maxP=np.array(maxP)
totU=np.array(totU); totW=np.array(totW)

# total bloom in ug/L (peak concentration) -> threshold 1 ug/L => nondim u_th = 1/rhoC
u_th = 1.0/rhoC
peakB = (maxU+maxW)   # nondim peak total concentration
peakB_ugL = peakB*rhoC

# ---- empirical decay time: last time peak total conc > 1 ug/L ----
above = np.where(peakB_ugL > 1.0)[0]
t_dec_nd = times[above[-1]] if len(above)>0 else np.nan
t_dec_days = t_dec_nd*t0/24.0

# ---- analytical bound ----
# ||u||_inf(t) <= ||u||_inf(Tc) * exp(C_N) * exp(-delta_nd (t-Tc)),  C_N = max_x N(Tc)/gamma_nd
Nmax_Tc = N_at_Tc.max(); Pmax_Tc = P_at_Tc.max()
uinf_Tc = uTc.max();    winf_Tc = wTc.max()
C_N = Nmax_Tc/gamma_nd
C_P = qmax_nd*Pmax_Tc/gamma_nd
# refined constant using rho=N/(1+N)<=min(N,1): Phi = (1+ln Nmax)/gamma if Nmax>1
Phi_N = (1.0+np.log(Nmax_Tc))/gamma_nd if Nmax_Tc>1 else Nmax_Tc/gamma_nd
# bound for time to fall below threshold (Mike branch, the dominant one near decay)
tdec_bound_u  = Tc + (1.0/delta_nd)*(np.log(uinf_Tc/u_th) + C_N)
tdec_bound_ref= Tc + (1.0/delta_nd)*(np.log(uinf_Tc/u_th) + Phi_N)
tdec_bound_days = tdec_bound_u*t0/24.0
tdec_ref_days   = tdec_bound_ref*t0/24.0

print("\n--- DECAY VERIFICATION ---")
print(f"Tc (nondim)={Tc}  = {Tc*t0/24:.1f} days")
print(f"N_max(Tc)={Nmax_Tc:.4g}, u_inf(Tc)={uinf_Tc:.4g}, delta_nd={delta_nd:.4g}, gamma_nd={gamma_nd:.4g}")
print(f"C_N = Nmax/gamma_nd = {C_N:.3f}")
print(f"empirical t_dec  = {t_dec_nd:.2f} nondim = {t_dec_days:.1f} days")
print(f"analytic  t_dec <= {tdec_bound_u:.2f} nondim = {tdec_bound_days:.1f} days (crude C_N={C_N:.1f})")
print(f"analytic  t_dec <= {tdec_bound_ref:.2f} nondim = {tdec_ref_days:.1f} days (refined Phi_N={Phi_N:.1f})")
print(f"bound holds (empirical <= analytic)? {t_dec_nd <= tdec_bound_ref}")

# check asymptotic exponential decay rate of peak after Tc matches ~ delta_nd
post = times > (Tc+5)
if post.sum() > 5:
    lp = np.log(np.maximum(peakB[post],1e-12))
    slope = np.polyfit(times[post], lp, 1)[0]
    print(f"\nfitted post-shutoff decay rate of peak total = {-slope:.4f} (nondim);  delta_nd={delta_nd:.4f}")

# ---- KPP spreading speed (dimensional) ----
# c* = 2 sqrt(alpha1 * r_eff), r_eff = max growth - mortality (dimensional, per hour)
# dimensional max growth rate = 1/t0 = kmax (per hr). r_eff = kmax - delta.
alpha1_um2_per_hr = alpha1*3600.0
c_star_um_per_hr = 2.0*np.sqrt(alpha1_um2_per_hr*(kmax-delta))
c_star_m_per_day = c_star_um_per_hr*1e-6*24.0
print("\n--- KPP SPREADING SPEED ---")
print(f"r_eff = kmax-delta = {kmax-delta:.4g} /hr")
print(f"c* = 2 sqrt(D*r) = {c_star_um_per_hr:.3g} um/hr = {c_star_m_per_day:.4g} m/day")
print(f"over 120-day season: front advances ~ {c_star_m_per_day*120:.3g} m (domain is 4e5 m)")

# ---- verification figure ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
days = times*t0/24.0
fig, ax = plt.subplots(1,2, figsize=(11,4.2))
ax[0].semilogy(days, peakB_ugL, 'm-', lw=2, label='sim: peak total bloom')
# envelope bound after Tc (Mike-dominated decay), in ug/L
tt = times[times>=Tc]
env = uinf_Tc*np.exp(Phi_N)*np.exp(-delta_nd*(tt-Tc))*rhoC
ax[0].semilogy(tt*t0/24, env, 'k--', lw=1.5, label=r'bound $\propto e^{-\delta(t-T_c)}$')
ax[0].axhline(1.0, color='gray', ls=':'); ax[0].axvline(Tc*t0/24, color='r', ls=':')
ax[0].set_xlabel('days'); ax[0].set_ylabel('peak conc (ug/L)')
ax[0].set_title('Bloom decay after source shutoff'); ax[0].legend(); ax[0].set_ylim(1e-2,1e2)
ax[1].plot(days, maxN*rhoN0, 'g-', label='peak N')
ax[1].plot(days, maxP*rhoP0, 'c-', label='peak P')
ax[1].axvline(Tc*t0/24, color='r', ls=':')
ax[1].set_xlabel('days'); ax[1].set_ylabel('peak nutrient (ug/L)')
ax[1].set_title('Nutrient flush-out (rate >= gamma)'); ax[1].legend()
plt.tight_layout(); plt.savefig('decay_verification.png', dpi=130)
print("\nsaved decay_verification.png")
