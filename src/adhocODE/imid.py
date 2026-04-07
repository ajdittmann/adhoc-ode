import numpy as np
_eps = np.finfo(float).eps
from scipy.optimize import fsolve

class Solver:
  def __init__(self, dydt, Ndim, atol):
    self.Ndim = Ndim
    self._dydt = dydt
    self.tol = np.max([atol*0.1, 10**-14])

  def imp1(self, yg, y0, t0, dt):
    dy = self._dydt(t0+dt, yg)*dt
    dyg = yg - y0
    return dy - dyg

  def imp2(self, yg, y0, dydt0, t0, dt):
    dydtg = self._dydt(t0+dt, yg)
    dy = 0.5*(dydtg + dydt0)*dt
    dyg = yg - y0
    return dy - dyg

  def update(self, t0, x0, dt):
    xg = np.copy(x0)

    y1, info, ier, mesg = fsolve(self.imp1, x0 = xg, args=(x0, t0, dt), full_output=1, xtol=self.tol)

    dy1 = self._dydt(t0, x0)
    y2, info, ier, mesg = fsolve(self.imp2, x0 = y1, args=(x0, dy1, t0, dt), full_output=1, xtol=self.tol)
    
    EE = y2-y1

    return y2, EE

  def getDt(self, dt0, EE, ynow, atol, rtol):
    arg = ((rtol*np.abs(ynow) + atol)/(np.abs(EE)+_eps))**2
    return dt0*np.sqrt(np.min(arg))**(0.5)


