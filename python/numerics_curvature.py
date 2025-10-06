# numerical solution for curvature perturbations in the cold regime [numerics_curvature.py]
#
# import basic tools and integration routines
import numpy as np
from scipy.integrate import solve_ivp 

# parameters [mpl]
fa = 1.25          # inflaton decay constant
m = 1.09e-6        # inflaton mass
koa_ref = 7.35e52  # momentum of perturbations at initial time
delta = 1.e-15     # regulator for circumventing poles [mpl^2]
Pi = np.pi

# potential of inflaton (phi) [mpl^4]
def V(phi): return m*m*fa*fa*( 1 - np.cos(phi/fa) )

# derivative of inflaton potential by phi [mpl^3]
def Vd(phi): return m*m*fa*np.sin(phi/fa)

# total energy density of inflaton field [mpl^4]
def etot(phi,phid): return phid*phid/2 + V(phi)

# Hubble rate [mpl]
def H(energy_density): return np.sqrt(8*Pi*np.abs(energy_density)/3)

# time derivative of Hubble rate [mpl^2]
def dotH(phid): return - 4*Pi*phid*phid

# solve for [dot phi, ddot phi, number of efolds] for the background evolution
# derivatives are taken with respect to t = Href*time
def derivatives1(t,y):
    Phi, Phid, Nfolds = y
    Hubble = H(etot(Phi, Phid));        dN_dt = Hubble/Href
    dphi_dt = Phid/Href;                ddphi_ddt = ( -3*Hubble*Phid - Vd(Phi) )/Href
    dy_dt = [dphi_dt, ddphi_ddt, dN_dt]
    return dy_dt

# solve for [dot phi, ddot phi, number of efolds, dot Rk, ddot Rk] including curvature perturbations
def derivatives2(t,y):
    Phi, Phid, efolds, reRk, imRk, reRkd, imRkd = y
    Hubble = H(etot(Phi, Phid));        dHubble = dotH(Phid);             dN_dt = Hubble/Href
    dphi_dt = Phid/Href;                ddphi_ddt = ( -3*Hubble*Phid - Vd(Phi) )/Href
    recalF = ddphi_ddt*Href*np.real(1/(Phid + 1j* delta)) - dHubble/Hubble
    imcalF = ddphi_ddt*Href*np.imag(1/(Phid + 1j* delta))
    dreRk_dt = reRkd/Href;              dimRk_dt = imRkd/Href
    ddreRk_ddt = - ( 2*(recalF*reRkd-imcalF*imRkd) + 3*Hubble*reRkd 
                     + koa_ref*koa_ref*np.exp(-2*efolds)*reRk )/Href
    ddimRk_ddt = - ( 2*(recalF*imRkd+imcalF*reRkd) + 3*Hubble*imRkd 
                     + koa_ref*koa_ref*np.exp(-2*efolds)*imRk )/Href
    dy_dt = [dphi_dt, ddphi_ddt, dN_dt, dreRk_dt, dimRk_dt, ddreRk_ddt, ddimRk_ddt]
    return dy_dt

# solve for [e_phi, number of efolds, dot Rk, ddot Rk] after averaging over oscillations
def derivatives3(t,y):
    e_phi, efolds, reRk, imRk, reRkd, imRkd = y
    Hubble = H(e_phi);                  dHubble = -4*Pi*e_phi;            dN_dt = Hubble/Href
    recalF =  - dHubble/Hubble;         imcalF = 0
    dreRk_dt = reRkd/Href;              dimRk_dt = imRkd/Href
    ddreRk_ddt = - ( 2*(recalF*reRkd-imcalF*imRkd) + 3*Hubble*reRkd 
                     + koa_ref*koa_ref*np.exp(-2*efolds)*reRk )/Href
    ddimRk_ddt = - ( 2*(recalF*imRkd+imcalF*reRkd) + 3*Hubble*imRkd 
                     + koa_ref*koa_ref*np.exp(-2*efolds)*imRk )/Href
    dy_dt = [-3*Hubble*e_phi/Href, dN_dt, dreRk_dt, dimRk_dt, ddreRk_ddt, ddimRk_ddt]
    return dy_dt

# initial conditions and reference values
phi_0 = 3.5                             # mpl
Href = np.sqrt(4*Pi/3)*m*phi_0          # mpl
phid_0 = - Vd(phi_0)/(3*Href)           # approximate slow-roll value for initial derivative [mpl^2]
N_0 = 0                                 # initial e-folds
koaH_start = 1e3                        # when to start solving for perturbations
time_end = 250                          # choose large enough that desired koaH is reached

# event function monitoring k/aH (here _not_ the correct H, rather the initial approximation)
def koaH(t, y):  return koa_ref*np.exp(-y[2])/Href - koaH_start

# integrate without curvature perturbations until k/aH obtains prescribed value koaH_start
sol1 = solve_ivp(derivatives1, [1, time_end], [phi_0, phid_0, N_0], events=koaH, rtol=1e-10)
t_match1 = sol1.t_events[0][0]          # point at which desired value was reached
time1 = sol1.t[sol1.t <= t_match1];     n_match1 = len(time1)
phi1 = sol1.y[0][:n_match1];            phid1 = sol1.y[1][:n_match1]
ephi1 = etot(phi1, phid1);              efolds1 = sol1.y[2][:n_match1]
reRk1 = np.ones(n_match1);              imRk1 = np.zeros(n_match1) 
slowroll1 = H(ephi1)**4/(2*Pi*phid1)**2

# define initial conditions for curvature perturbations at moment t_match1
time_0 = t_match1
phi_0 = sol1.y_events[0][0][0]; phid_0 = sol1.y_events[0][0][1]; N_0 = sol1.y_events[0][0][2]
reRk_0 = -H(etot(phi_0,phid_0))/(2*Pi*phid_0)*koa_ref*np.exp(-N_0);  imRk_0 = 0.
dreRk_0 = 0.;                           dimRk_0 = - reRk_0* koa_ref*np.exp(-N_0)
print("start Rk-evolution at t/tref = ",time_0,", with k/aH approx ",koaH_start)
koaH_start = 10                         # re-adjust target k/aH
time_end = 250                          # choose large enough that desired koaH is reached

# integrate with curvature perturbations until k/aH obtains prescribed value koaH_start
sol2 = solve_ivp(derivatives2, [time_0,time_end],
                 [phi_0, phid_0, N_0, reRk_0, imRk_0, dreRk_0, dimRk_0],
                 events=koaH,t_eval=np.linspace(time_0,time_end,1000),
                 atol=1e-12,rtol=1e-13,method='DOP853')  
t_match2 = sol2.t_events[0][0]          # point at which desired value was reached
time2 = sol2.t[sol2.t <= t_match2];     n_match2 = len(time2)
phi2 = sol2.y[0][:n_match2];            phid2 = sol2.y[1][:n_match2]
ephi2 = etot(phi2,phid2);               efolds2 = sol2.y[2][:n_match2]
reRk2 = sol2.y[3][:n_match2];           imRk2 = sol2.y[4][:n_match2]
slowroll2 = H(ephi2)**4/(2*Pi*phid2)**2

# define initial conditions for next range (crossing outside of Hubble horizon)
time_0 = t_match2
phi_0 = sol2.y_events[0][0][0] ;  phid_0 = sol2.y_events[0][0][1];  N_0 = sol2.y_events[0][0][2]
reRk_0 = sol2.y_events[0][0][3];        imRk_0 = sol2.y_events[0][0][4]
dreRk_0 = sol2.y_events[0][0][5];       dimRk_0 = sol2.y_events[0][0][6]
print("start horizon crossing at t/tref = ",time_0,", with k/aH approx ",koaH_start)
time_end = 1e3                          # choose large enough that 10 full oscillations take place

# define the event function counting the zeros of phi
def phi_crossings(t, y):  return y[0]

# integrate up to 10 oscillations (21 crossings of zero)
sol3 = solve_ivp(derivatives2, [time_0,time_end],
                 [phi_0, phid_0, N_0, reRk_0, imRk_0, dreRk_0, dimRk_0],
                 events=phi_crossings,t_eval=np.linspace(time_0,time_end,1000),
                 atol=1e-12,rtol=1e-13,method='DOP853')  
t_match3 = sol3.t_events[0][20]         # point at which desired value was reached
time3 = sol3.t[sol3.t <= t_match3];     n_match3 = len(time3)
phi3 = sol3.y[0][:n_match3];            phid3 = sol3.y[1][:n_match3]
ephi3 = etot(phi3,phid3);               efolds3 = sol3.y[2][:n_match3]
reRk3 = sol3.y[3][:n_match3];           imRk3 = sol3.y[4][:n_match3]
slowroll3 = H(ephi3)**4/(2*Pi*phid3)**2

# define initial conditions for averaged regime
time_0 = t_match3;                      time_end = 1e8
phi_0 = sol3.y_events[0][20][0];        phid_0 = sol3.y_events[0][20][1]
etot_0 = etot(phi_0,phid_0);            N_0 = sol3.y_events[0][20][2]
reRk_0 = sol3.y_events[0][20][3];       imRk_0 = sol3.y_events[0][20][4] 
dreRk_0 = sol3.y_events[0][20][5];      dimRk_0 = sol3.y_events[0][20][6]
print("start averaged regime at t/tref = ",time_0)

# integrate in averaged regime (matter-dominated era)
sol4 = solve_ivp(derivatives3, [time_0, time_end],
                 [etot_0, N_0, reRk_0, imRk_0, dreRk_0, dimRk_0],
                 atol=1e-12,rtol=1e-13,method='DOP853')
time4 = sol4.t;                         n_match4 = len(time4)
ephi4 = sol4.y[0];                      efolds4 = sol4.y[1]
reRk4 = sol4.y[2];                      imRk4 = sol4.y[3];            slowroll4 = np.zeros(n_match4)

# assemble together the complete solution
time = np.concatenate((time1,time2,time3,time4)); e_phi = np.concatenate((ephi1,ephi2,ephi3,ephi4))
efolds = np.concatenate((efolds1,efolds2,efolds3,efolds4)); koa = koa_ref*np.exp(-efolds)
reRk = np.concatenate((reRk1,reRk2,reRk3,reRk4)); imRk = np.concatenate((imRk1,imRk2,imRk3,imRk4))
slowroll = np.concatenate((slowroll1,slowroll2,slowroll3,slowroll4))

# print to file
np.savetxt('numerics_curvature.dat',
           np.c_[time, e_phi, H(e_phi), efolds, koa, reRk**2+imRk**2, slowroll],fmt='%.6e',newline='\n',
           header='columns: t*H_ref, e_phi/mpl**4, H/mpl, efolds, k_ref/a/mpl, P_R, slowroll')
