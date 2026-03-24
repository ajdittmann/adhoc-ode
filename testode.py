import numpy as np
import adhocODE as adhoc
import matplotlib.pyplot as plt

def customDT1(t, y, k=1):
  return 0.3

def customDT2(t, y, k=1):
  return 0.1

def customDT3(t, y, k=1):
  return 0.03

def dydt(t, y, k=1):
  y0 = y[0]
  y1 = y[1]
  #y2 = y[2]
  #y3 = y[3]

  #df = np.array([y1, y2, y3, y2*(12*y0*y0 + 8.0)])
  df = np.array([y1, -k*y0])
  return df

y0 = np.array([0.0, 1.0]).T #, 0.0, 2.0]).T

res = adhoc.solve_ivp(dydt, [0,-10], y0, args=[3.0], tol=1e-2, t_eval = np.linspace(0,-10,30), dtfunc=customDT1)
plt.plot(res.t, res.y[0,:])

print('first finished')

res = adhoc.solve_ivp(dydt, [0,-10], y0, args=[3.0], tol=1e-2, t_eval = np.linspace(0,-10,30), dtfunc=customDT2)
plt.plot(res.t, res.y[0,:])

print('second finished')

res = adhoc.solve_ivp(dydt, [0,-10], y0, args=[3.0], tol=1e-16, t_eval = np.linspace(0,-10,30), dtfunc=customDT3 )
plt.plot(res.t, res.y[0,:])

plt.show()
