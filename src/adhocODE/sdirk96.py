import numpy as np
_eps = np.finfo(float).eps
from scipy.optimize import fsolve

### constants for SDIRK(9,6) method ###
#_a2 = np.empty(2)
#_a2[0] =  0.5
#_a2[1] =  0.5
#_c2 = np.sum(_a2)

_a2 = np.empty(2)
_a2[0] = -9.035148561194185e-02
_a2[1] =  2.181277819449076e-01
_c2 = np.sum(_a2)

_a3 = np.empty(3)
_a3[0] =  1.729520391389366e-01
_a3[1] = -3.536550103628203e-01
_a3[2] =  2.181277819449076e-01
_c3 = np.sum(_a3)

_a4 = np.empty(4)
_a4[0] =  5.119998759191926e-01
_a4[1] =  2.896403322019248e-02
_a4[2] = -1.440309456570937e-02
_a4[3] =  2.181277819449076e-01
_c4 = np.sum(_a4)

_a5 = np.empty(5)
_a5[0] =  4.653034955067823e-03
_a5[1] = -7.563581876659697e-02
_a5[2] =  2.172730307867122e-01
_a5[3] = -2.065194287254723e-02
_a5[4] =  2.181277819449076e-01
_c5 = np.sum(_a5)

_a6 = np.empty(6)
_a6[0] =  8.961455017624717e-01
_a6[1] =  1.392673277004985e-01
_a6[2] = -1.869209797528052e-01
_a6[3] =  6.729710123717235e-02
_a6[4] = -3.508919634421756e-01
_a6[5] =  2.181277819449076e-01
_c6 = np.sum(_a6)

_a7 = np.empty(7)
_a7[0] =  5.529597018857514e-01
_a7[1] = -4.393605797936621e-01
_a7[2] =  3.337040023250907e-01
_a7[3] = -3.394265207784165e-02
_a7[4] = -1.519474459125954e-01
_a7[5] =  2.138256610269428e-02
_a7[6] =  2.181277819449076e-01
_c7 = np.sum(_a7)

_a8 = np.empty(8)
_a8[0] =  6.313603740364756e-01
_a8[1] =  7.247336196414658e-01
_a8[2] = -4.321706254252584e-01
_a8[3] =  5.986113821824766e-01
_a8[4] = -7.090871970343450e-01
_a8[5] = -4.839866856969341e-01
_a8[6] =  3.783915629051305e-01
_a8[7] =  2.181277819449076e-01
_c8 = np.sum(_a8)

_a9 = np.empty(9)
_a9[0] =  0.0
_a9[1] = -1.550445253086903e-01
_a9[2] =  1.945184786607890e-01
_a9[3] =  6.351564027920301e-01
_a9[4] =  8.117227866417299e-01
_a9[5] =  1.107361086915850e-01
_a9[6] = -4.953046924144789e-01
_a9[7] = -3.199123410078724e-01 
_a9[8] =  2.181277819449076e-01
_c9 = np.sum(_a9)

_bl = np.empty(9)
_bl[0] =  0.0
_bl[1] =  7.366155582789420e-02 
_bl[2] =  1.035273972622287e-01
_bl[3] =  1.002474819354989e+00
_bl[4] =  3.613772892500572e-01
_bl[5] = -7.854259299613646e-01
_bl[6] = -1.704990479607844e-02
_bl[7] =  2.963212522147690e-01
_bl[8] = -3.488647915249531e-02

class Solver:
  def __init__(self, dydt, Ndim):
    self.Ndim = Ndim
    self._dydt = dydt

    self._A2 = np.tile(_a2, (Ndim,1)).T
    self._A3 = np.tile(_a3, (Ndim,1)).T
    self._A4 = np.tile(_a4, (Ndim,1)).T
    self._A5 = np.tile(_a5, (Ndim,1)).T
    self._A6 = np.tile(_a6, (Ndim,1)).T
    self._A7 = np.tile(_a7, (Ndim,1)).T
    self._A8 = np.tile(_a8, (Ndim,1)).T
    self._A9 = np.tile(_a9, (Ndim,1)).T
    self._BL = np.tile(_bl, (Ndim,1)).T

  def imp1(self, yg, y0, t0, dt):
    dy = self._dydt(t0+dt, yg)*dt
    dyg = yg - y0
    return dy - dyg

  #implicit update for stage 2
  def imp2(self, yg, y0, dydt0, t0, dt):
    dydtg = self._dydt(t0 + dt*_c2, yg)
    dyg = yg - y0
    dy = dt*np.sum(np.array([dydt0, dydtg])*self._A2, axis=0)
    return dy - dyg

  def impN(self, yg, y0, dys, _AN, _cn, t0, dt):
    dydtg = self._dydt(t0 + dt*_cn, yg)
    dyg = yg - y0
    dy = dt*np.sum(np.array(dys + [dydtg])*_AN, axis=0)
    return dy - dyg

  def update(self, t0, x0, dt):
    xg = np.copy(x0)

    #y1, info, ier, mesg = fsolve(self.imp1, x0 = xg, args=(x0, t0, dt), full_output=1, xtol=1.e-14)
    dy1 = self._dydt(t0, x0)

    y2, info, ier, mesg = fsolve(self.imp2, x0 = x0, args=(x0, dy1, t0, dt), full_output=1, xtol=1.e-14)
    dys = [dy1]
    dys.append(self._dydt(t0 + dt*_c2, y2))

    y3, info, ier, mesg = fsolve(self.impN, x0 = y2, args=(x0, dys, self._A3, _c3, t0, dt), full_output=1, xtol=1.e-14)
    dys.append(self._dydt(t0 + dt*_c3, y3))

    y4, info, ier, mesg = fsolve(self.impN, x0 = y3, args=(x0, dys, self._A4, _c4, t0, dt), full_output=1, xtol=1.e-14)
    dys.append(self._dydt(t0 + dt*_c4, y4))

    y5, info, ier, mesg = fsolve(self.impN, x0 = y4, args=(x0, dys, self._A5, _c5, t0, dt), full_output=1, xtol=1.e-14)
    dys.append(self._dydt(t0 + dt*_c5, y5))

    y6, info, ier, mesg = fsolve(self.impN, x0 = y5, args=(x0, dys, self._A6, _c6, t0, dt), full_output=1, xtol=1.e-14)
    dys.append(self._dydt(t0 + dt*_c6, y6))

    y7, info, ier, mesg = fsolve(self.impN, x0 = y6, args=(x0, dys, self._A7, _c7, t0, dt), full_output=1, xtol=1.e-14)
    dys.append(self._dydt(t0 + dt*_c7, y7))

    y8, info, ier, mesg = fsolve(self.impN, x0 = y7, args=(x0, dys, self._A8, _c8, t0, dt), full_output=1, xtol=1.e-14)
    dys.append(self._dydt(t0 + dt*_c8, y8))

    y9, info, ier, mesg = fsolve(self.impN, x0 = y8, args=(x0, dys, self._A9, _c9, t0, dt), full_output=1, xtol=1.e-14)
    dys.append(self._dydt(t0 + dt*_c9, y9))

    yL = x0 + dt*np.sum(np.array(dys)*self._BL, axis=0)

    EE = (y9-yL)
    
    return y9, EE

  def getDt(self, dt0, EE, ynow, target):
    arg = ((target*np.abs(ynow) + target)/(EE+_eps))**2
    return dt0*np.sqrt(np.mean(arg))**(1/6)




