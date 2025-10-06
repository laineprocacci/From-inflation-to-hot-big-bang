# acoustic oscillations and Jeans instability after re-entry [numerics_acoustic_jeans.py]
#
# import basic tools and integration routines
import numpy as np
from scipy.integrate import solve_ivp

# parameters 
koaH = 0.2   # initial momentum over Hubble rate
alpha = 0.5  # rate of background expansion (radiation dominated)

# sources for [Rv, RT, Rvdm, Redm]
# derivatives are taken with respect to x = ln(t/t_out)
def derivatives(x,y):
    Rv, RT, Rvdm, Redm = y
    dRv_dx = alpha*(Rv - RT)
    dRT_dx = Rv - RT + alpha/3*np.exp(2*(1-alpha)*x)*koaH*koaH*Rv
    dRvdm_dx = Rv - Rvdm
    dRedm_dx = Rv - RT + alpha/3*np.exp(2*(1-alpha)*x)*koaH*koaH*Rvdm
    dy_dx = [dRv_dx, dRT_dx, dRvdm_dx, dRedm_dx]
    return dy_dx

# initial conditions 
Rv_0 = 1.0; RT_0 = 1.0; Rvdm_0 = 1.0; Redm_0 = 1.0

# integrate up to chosen point 
sol = solve_ivp(derivatives,
                [0, 15/2/(1-alpha)],
                [Rv_0, RT_0, Rvdm_0, Redm_0],method='DOP853',
                t_eval=np.linspace(0, 15/2/(1-alpha),3000)) 
x_grid = sol.t
Rv_sol = sol.y[0]; RT_sol = sol.y[1]; Rvdm_sol = sol.y[2]; Redm_sol = sol.y[3]

# print to file
np.savetxt('numerics_acoustic_jeans.dat',
           np.c_[np.exp(x_grid),Rv_sol**2, RT_sol**2, Rvdm_sol**2, Redm_sol**2], 
           fmt='%.6e',newline='\n',
           header='columns: t/t_out, Rv^2, RT^2, Rvdm^2, Redm^2')

