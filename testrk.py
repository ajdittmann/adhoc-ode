import numpy as np
import rk8

def dydt(y):
  y0 = y[0]
  y1 = y[1]
  y2 = y[2]
  y3 = y[3]

  df = np.array([y1, y2, y3, y2*(12*y0*y0 + 8.0)])
  return df

tend = 1.0

ns = np.unique(np.geomspace(1, 2**11, 80, dtype='int'))

errs = []
dts = []

truth = np.tan(tend)
print(truth)

terrs = np.geomspace(10**-3, 10**-12, 20)
y0 = np.array([0.0, 1.0, 0.0, 2.0]).T

solver = rk8.Solver(dydt, len(y0) )

for terr in terrs:
  yn = np.copy(y0)
  dt = tend*0.1
  f0, ee = solver.update(yn, dt)  
  for i in range(2):
    dt = dt*(terr/np.max(np.abs(ee)))**(1/8)
    f0, ee = solver.update(yn, dt)  
  print(ee, terr)



#print(errs)
'''
import matplotlib.pyplot as plt
plt.scatter(dts, errs)
plt.plot(dts, errs[0]*10.0*(dts/dts[0])**8, color='gray', ls ='--' )
plt.xscale('log')
plt.yscale('log')
plt.show()

'''
