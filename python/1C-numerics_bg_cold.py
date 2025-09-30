# numerical solution for inflation and end of inflation [numerics_bg_cold.py]
#
# import basic tools and integration routines
import numpy as np
from scipy.integrate import solve_ivp

# parameters [mpl]
fa = 1.25          # inflaton decay constant
m = 1.09e-6        # inflaton mass
koa_ref = 7.35e52  # momentum of perturbations at initial time
Pi = np.pi

# potential of inflaton (phi) [mpl^4]
def V(phi): return m*m*fa*fa*( 1 - np.cos(phi/fa) )

# derivative of inflaton potential by phi [mpl^3]
def Vd(phi): return m*m*fa*np.sin(phi/fa)

# total energy density of scalar inflaton field [mpl^4]
def etot(phi,phid): return phid*phid/2 + V(phi)

# Hubble rate [mpl]
def H(energy_density): return np.sqrt(8*Pi*np.abs(energy_density)/3)

# 2nd-order ODE as two 1st-order ODEs: solve for [dot phi, ddot phi, number of efolds]
# derivatives are taken with respect to t = Href*time
def derivatives1(t,y):
    Phi, Phid, Nfolds = y
    dphi_dt = Phid/Href
    ddphi_ddt = ( -3*H(etot(Phi, Phid))*Phid - Vd(Phi) )/Href
    dN_dt = H(etot(Phi, Phid))/Href
    dy_dt = [dphi_dt, ddphi_ddt, dN_dt]
    return dy_dt

# average over fast oscillations at late times => one 1st-order ODE
def derivatives2(t,y):
    e_phi, Nfolds = y
    dy_dt = [-3*H(e_phi)*e_phi/Href,  H(e_phi)/Href]
    return dy_dt

# initial conditions and reference values
phi_0 = 3.5                     # mpl
Href = np.sqrt(4*Pi/3)*m*phi_0  # mpl
phid_0 = - Vd(phi_0)/(3*Href)   # mpl^2

# define the event function counting the zeros of phi
def phi_crossings(t, y):  return y[0]

# integrate up to 50 oscillations (101 crossings of zero)
sol1 = solve_ivp(derivatives1, [1, 1e4], [phi_0, phid_0, 0], events=phi_crossings)
t_match = sol1.t_events[0][100]
time1 = sol1.t[sol1.t <= t_match]
n_match = len(time1)
ephi1 = etot(sol1.y[0], sol1.y[1])
efolds1 = sol1.y[2]

# initial conditions for oscillatory regime
time_0 = time1[n_match-1]
etot_0 = ephi1[n_match-1]
efolds_0 = efolds1[n_match-1]

# integrate in oscillatory regime (matter-dominated era)
sol2 = solve_ivp(derivatives2, [time_0, 1e8], [etot_0, efolds_0], method='DOP853')
ephi2 = sol2.y[0]
efolds2 = sol2.y[1]

# append solutions
time = np.append(time1,sol2.t)
e_phi = np.append(ephi1[:n_match], ephi2)
efolds = np.append(efolds1[:n_match], efolds2)
koa = koa_ref*np.exp(-efolds)

# print to file
np.savetxt('numerics_bg_cold.dat', np.c_[time, e_phi, H(e_phi), efolds, koa],
           fmt='%.6e',newline='\n',
           header='columns: t*H_ref, e_phi/mpl**4, H/mpl, efolds, k_ref/a/mpl')

