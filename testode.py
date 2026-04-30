import numpy as np
import adhocODE as adhoc
import matplotlib.pyplot as plt

def fixedDT(t, y, k=1, dtout=0.1):
  return dtout

def dydt(t, y, k=1):
  y0 = y[0]
  y1 = y[1]

  df = np.array([y1, -k*y0])
  return df

y0 = np.array([0.0, 1.0]).T
en0 = 0.5*y0[1]**2

Ns = np.geomspace(1.0, 10**5, 100)
errs = []

for N in Ns:
  dt = 2.0/N
  res = adhoc.solve_ivp(dydt, [0,2], y0, args=[3.0], t_eval = np.linspace(0,2,2), method='rk87', dtfunc=fixedDT, dtfunc_args=[dt], rtol=10**10, atol=10**10)
  en = 0.5*3*res.y[0,-1]**2 + 0.5*res.y[1,-1]**2
  errs.append( np.abs((en-en0)/en0) )

plt.scatter(Ns, errs)
plt.xscale('log')
plt.yscale('log')
plt.show()
