# numerical solution for reheating and/or warm inflation [numerics_bg_warm.py]
#
# import basic tools and integration routines
import numpy as np
from scipy.integrate import solve_ivp

# parameters [mpl]
fa = 1.25;       m = 1.09e-6;    Pi = np.pi     # inflaton decay constant and mass, shorthand for pi
gstar = 106.75;  hstar = gstar;  istar = gstar  # energy density, entropy density, heat capacity

# parameters for the cases considered
time_end = 1e20        # time sufficiently large that the system has reheated
case = 'thin' # 'thick'
if case == 'thin': 
    kappaT = 1e0;      kappam = 1e0;     T_0 = 1.e-9        # initial seed temperature [mpl]
else:
    kappaT = 1e7;      kappam = 1e9;     T_0 = 1.e-8        

# potential of inflaton (phi) [mpl^4]
def V(phi): return m*m*fa*fa*( 1 - np.cos(phi/fa) ) 

# derivative of inflaton potential by phi [mpl^3]
def Vd(phi): return m*m*fa*np.sin(phi/fa)

# friction coefficient [mpl]
def Ups(T): return (kappaT*Pi*Pi*Pi*T*T*T + kappam*m*m*m)/(4*4*4*Pi*Pi*Pi)/(fa*fa)

# radiation energy density [mpl^4]
def er(T): return gstar*Pi*Pi*T*T*T*T/30

# radiation entropy density [mpl^3]
def sr(T): return 2*hstar*Pi*Pi*T*T*T/45

# radiation heat capacity [mpl^3]
def cr(T): return 2*istar*Pi*Pi*T*T*T/15

# energy density of inflaton field [mpl^4]
def efield(phi,phid): return phid*phid/2 + V(phi)

# Hubble rate [mpl]
def H(energy_density): return np.sqrt(8*Pi*np.abs(energy_density)/3)

# evolution at early times: solve for [phi, dot phi, efolds, T]
# derivatives are taken with respect to t=Href*time
def derivatives1(t,y):
    Phi, Phid, Nfolds, T = y
    Hubble = H(efield(Phi, Phid)+er(T))
    dphi_dt = Phid/Href;      ddphi_ddt = (-(3*Hubble+Ups(T))*Phid - Vd(Phi))/Href
    dN_dt = Hubble/Href;      dT_dt = (Ups(T)*Phid*Phid-3*Hubble*T*sr(T))/cr(T)/Href
    dy_dt = [dphi_dt, ddphi_ddt, dN_dt, dT_dt]; return dy_dt

# evolution in averaged regime: solve for [e_phi, efolds, T]
def derivatives2(t,y):
    e_phi, Nfolds, T = y
    Hubble = H(e_phi+er(T))
    dephi_dt = (-(3*Hubble+Ups(T))*e_phi)/Href
    dN_dt = Hubble/Href;      dT_dt = (Ups(T)*e_phi-3*Hubble*T*sr(T))/cr(T)/Href
    dy_dt = [dephi_dt, dN_dt, dT_dt];  return dy_dt   

# evolution after e_phi has become insignificant: solve for [efolds, T]
def derivatives3(t,y):
    Nfolds, T = y
    dN_dt = H(er(T))/Href;    dT_dt = (-3*H(er(T))*T*sr(T))/cr(T)/Href
    dy_dt = [dN_dt, dT_dt];  return dy_dt   

# initial conditions and reference values [mpl]
phi_0 = 3.5;   Href = np.sqrt(4*Pi/3)*m*phi_0;   phid_0 = - Vd(phi_0)/(3*Href);   N_0 = 0

# define the event function counting the zeros of phi
def phi_crossings(t, y): return y[0]  

# integrate up to 10 oscillations (21 crossings of zero)
sol1 = solve_ivp(derivatives1, [1, 1e4], [phi_0, phid_0, N_0, T_0],
                 events=phi_crossings,method='DOP853',atol=1e-12,rtol=1e-13) 
t_match = sol1.t_events[0][20]         # time at which condition was met
time1 = sol1.t[sol1.t <= t_match];     n_match1 = len(time1)  
ephi1 = efield(sol1.y[0][:n_match1], sol1.y[1][:n_match1]);
efolds1 = sol1.y[2][:n_match1];        T1 = sol1.y[3][:n_match1]

# initial conditions for averaged regime
time_0 = t_match;   ephi_0 = efield(sol1.y_events[0][20][0], sol1.y_events[0][20][1])
efolds_0 = sol1.y_events[0][20][2];   T_0 = sol1.y_events[0][20][3] 
print("go over to averaged evolution at t/t_ref = ",time_0)

# define the event function locating when e_field = e_r/500
def phi_decays(t, y): return y[0]-er(y[2])/500
setattr(phi_decays,'terminal',True)    # for stopping immediately when condition is met

# integrate in averaged regime until phi decays
if(ephi_0 > er(T_0)/500):
    sol2 = solve_ivp(derivatives2, [time_0, time_end], [ephi_0, efolds_0, T_0],
                     events=phi_decays,method='DOP853',atol=1e-12,rtol=1e-13)  
    t_match = sol2.t_events[0][0]      # time at which condition was met
    time2 = sol2.t[sol2.t <= t_match]; n_match2 = len(time2)  
    ephi2 = sol2.y[0][:n_match2];      efolds2 = sol2.y[1][:n_match2];    T2 = sol2.y[2][:n_match2]
    time_0 = t_match;                  efolds_0 = sol2.y_events[0][0][1]; T_0 = sol2.y_events[0][0][2] 
else:
    t_match = time_0;                  time2 = np.array([])
    ephi2 = np.array([]);              efolds2 = np.array([]);            T2 = np.array([])
print("go over to evolution without e_phi at t/t_ref = ",time_0)

# parameters for matching to Standard Model
mplMpc = 1.9092e57                     # conversion between microscopic and macroscopic scales 
Tswitchmpl = 1.e-12                    # temperature at which we switch to Standard Model [mpl]
efolds_after_switch = 46.5             # efolds from Tswitch until today

# define the event function locating the switching temperature
def Tswitch(t, y): return y[1]-Tswitchmpl  

# integrate until temperature is reached
sol3 = solve_ivp(derivatives3, [time_0, time_end], [efolds_0, T_0],
                 events=Tswitch,method='DOP853',atol=1e-12,rtol=1e-13)  
t_match = sol3.t_events[0][0]          # time at which condition was met
time3 = sol3.t[sol3.t <= t_match];     n_match3 = len(time3)  
ephi3 = np.zeros(n_match3);            efolds3 = sol3.y[0][:n_match3];    T3 = sol3.y[1][:n_match3]

# collect information about switching point 
efolds_till_switch = sol3.y_events[0][0][0]
print("Tswitch reached at t/tref = ",t_match," after ",efolds_till_switch," e-folds")
coeff = np.exp(-efolds_till_switch)*np.exp(-efolds_after_switch)*mplMpc

# initial conditions for the last domain
time_0 = t_match;                      efolds_0 = sol3.y_events[0][0][0]; T_0 = sol3.y_events[0][0][1] 

# integrate till the very end
sol4 = solve_ivp(derivatives3, [time_0, time_end], [efolds_0, T_0],
                 method='DOP853',atol=1e-12,rtol=1e-13)  
time4 = sol4.t;                        n_match4 = len(time4)  
t_match = time4[n_match4-1]            # time reached at the end
ephi4 = np.zeros(n_match4);            efolds4 = sol4.y[0];               T4 = sol4.y[1]
print("integration terminated at t/tref = ",t_match," after ",efolds4[n_match4-1]," e-folds")

# assemble together the complete solution
time = np.concatenate((time1,time2,time3,time4)); e_phi = np.concatenate((ephi1,ephi2,ephi3,ephi4))
efolds = np.concatenate((efolds1,efolds2,efolds3,efolds4)); T = np.concatenate((T1,T2,T3,T4))
e_r = er(T)

# print out background solution
np.savetxt('numerics_bg_warm.dat',
           np.c_[time, e_phi, H(e_phi+e_r), T, Ups(T), efolds],fmt='%.6e',newline='\n',
           header='columns: t*H_ref, e_phi/mpl**4, H/mpl, T/mpl, Ups/mpl, efolds')

# print out relation of microscopic and macroscopic momentum scales
koampl = np.geomspace(1e-10,1e60,100); koaMpc = coeff*koampl
np.savetxt('numerics_redshift.dat',
           np.c_[koaMpc,koampl],fmt='%.6e',newline='\n',
           header='columns: k/a_0*Mpc, k/a_i/mpl')
