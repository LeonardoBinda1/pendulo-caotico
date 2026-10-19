"""Independent short-time integration and exact variational contraction check."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad
from reproduce import orbit,TD,OM,DATA

def rhs(t,y):
    theta,omega=y[:2]
    J=np.array([[0.,1.],[-np.cos(theta),-.5]])
    return np.r_[omega,-np.sin(theta)-.5*omega+1.2*np.cos(OM*t),
                 (J@y[2:].reshape(2,2)).ravel()]

reference=solve_ivp(rhs,[0,TD],np.r_[.2,0.,np.eye(2).ravel()],
                    method='DOP853',rtol=1e-11,atol=1e-13)
assert reference.success
state=reference.y[:2,-1].copy()
state[0]=(state[0]+np.pi)%(2*np.pi)-np.pi
points,_=orbit(1.2,240,0,1)
error=float(np.max(np.abs(state-points[0])))
determinant=float(np.linalg.det(reference.y[2:,-1].reshape(2,2)))
exact=float(np.exp(-.5*TD))
integral=quad(lambda t:np.cos(OM*t)/np.cosh(t),-40,40,epsabs=1e-11)[0]
transform=np.pi/np.cosh(np.pi*OM/2)
assert error<1e-5,(error,'RK4 differs from the independent one-period integration')
assert abs(determinant-exact)<1e-9,(determinant,exact)
assert abs(integral-transform)<1e-10,(integral,transform)
result={'one_period_RK4_DOP853_max_error':error,'one_period_det_Phi':determinant,
        'exact_exp_minus_q_T':exact,'sech_integral_error':float(abs(integral-transform))}
(DATA/'independent_checks.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
