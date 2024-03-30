# 22 degree or 46 degree halo

import matplotlib.pyplot as plt
from scipy.constants import degree
from halo import halo

zeta = 30*degree

plt.figure(figsize=(5,3.75))

x,y = halo(zeta) # 22 deg
plt.plot(x/degree, y/degree, '.', ms=1)

x,y = halo(zeta, io=(1,2), N=30000) # 46 deg
plt.plot(x[0]/degree, y[0]/degree, '.', ms=1)
plt.plot(x[1]/degree, y[1]/degree, '.', ms=1)

plt.axis('equal')
plt.tight_layout()
plt.savefig('fig2.png')
plt.show()
