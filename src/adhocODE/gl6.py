import numpy as np
_eps = np.finfo(float).eps
from scipy.optimize import fsolve

### constants for 6th-order Gauss-Legendre method
_AG = np.empty((3,3))

_AG[0,0] = 5/36
_AG[0,1] = 2/9 - np.sqrt(15)/15
_AG[0,2] = 5/36 - np.sqrt(15)/30

_AG[1,0] = 5/36 + np.sqrt(15)/24
_AG[1,1] = 2/9
_AG[1,2] = 5/36 - np.sqrt(15)/24

_AG[2,0] = 5/36 + np.sqrt(15)/30
_AG[2,1] = 2/9 + np.sqrt(15)/15
_AG[2,2] = 5/36

_CG = 0.5 + np.sqrt(15)*0.1*np.array([-1, 0, 1])
_BG = np.array([2.5, 4.0, 2.5])/9

class Solver:
  def __init__(self, dydt, Ndim, atol):
    self.Ndim = Ndim
    self._dydt = dydt
    self.tol = 2**-51
    self.tol1 = 10.0*np.max([atol, 10**-14])

  def imp1(self, yg, y0, t0, dt):
    dy = self._dydt(t0+dt, yg)*dt
    dyg = yg - y0
    return dy - dyg

  def imp3(self, yg, y0, dt, t0):
    df0 = self._dydt(t0 + dt*_CG[0], yg[:self.Ndim])
    df1 = self._dydt(t0 + dt*_CG[1], yg[self.Ndim:2*self.Ndim])
    df2 = self._dydt(t0 + dt*_CG[2], yg[2*self.Ndim:])

    dy = np.empty((3*self.Ndim))
    dy[:self.Ndim] = df0*_AG[0,0] + df1*_AG[0,1] + df2*_AG[0,2]
    dy[self.Ndim:2*self.Ndim] = df0*_AG[1,0] + df1*_AG[1,1] + df2*_AG[1,2]
    dy[2*self.Ndim:] = df0*_AG[2,0] + df1*_AG[2,1] + df2*_AG[2,2]

    dyg = np.copy(yg)
    dyg[:self.Ndim] -= y0
    dyg[self.Ndim:2*self.Ndim] -= y0
    dyg[2*self.Ndim:] -= y0
    return dy*dt - dyg

  def update(self, t0, x0, dt):
    x0g, info, ier, mesg = fsolve(self.imp1, x0, args = (x0, t0, dt*_CG[0]), full_output=1, xtol=self.tol1)
    x1g, info, ier, mesg = fsolve(self.imp1, x0g, args = (x0, t0 + dt*_CG[0], dt*(_CG[1]-_CG[0]) ), full_output=1, xtol=self.tol1)
    x2g, info, ier, mesg = fsolve(self.imp1, x1g, args = (x0, t0 + dt*_CG[1], dt*(_CG[2]-_CG[1]) ), full_output=1, xtol=self.tol1)

    XG0 = np.append(x0g, x1g)
    XG0 = np.append(XG0, x2g)
    yarr, info, ierr, mesg = fsolve(self.imp3, x0=XG0, args=(x0, dt, t0), full_output=1, xtol=self.tol)

    df0 = self._dydt(t0 + dt*_CG[0], yarr[:self.Ndim])
    df1 = self._dydt(t0 + dt*_CG[1], yarr[self.Ndim:2*self.Ndim])
    df2 = self._dydt(t0 + dt*_CG[2], yarr[2*self.Ndim:])

    update = dt*(df0*_BG[0] + df1*_BG[1] + df2*_BG[2])
    # N.B. embedded method is only 2nd order
    yerr = (-df0*2.5 + 8*df1 - 2.5*df2)*dt/3
    EE = update-yerr

    return update, EE*dt**3

  def getDtNorm(self, EE, ynow, atol, rtol):
    arg = ((rtol*np.abs(ynow) + atol)/(np.abs(EE)+_eps))**2
    return np.sqrt(np.mean(arg))**(1/6)

