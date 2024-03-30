# contact arcs to the 46 degree halo

import matplotlib.pyplot as plt
from scipy.constants import degree
from halo import *
from math import pi

zeta = [0,5,10,20,30,40,50,60,70,80] # sun elevation / deg
psi0 = 0 # max spin angle

r22 = min_deflection()
r46 = min_deflection(pi/2)

plt.figure(figsize=(5,10))

for i,zt in enumerate(zeta):
    plt.subplot(5,2,i+1)
    z = zt*degree; r = pi/2 - z
    x,y = halo(z, h=0.1, spin_max=psi0, io=(1,2), N=30000)
    plt.plot(x[0],y[0], '.', ms=1)
    plt.plot(x[1],y[1], '.', ms=1)
    x,y = circle(z, r22)
    plt.plot(x,y, '--', lw=1)
    x,y = circle(pi/2, pi/2, -r)
    plt.plot(x,y, ':', lw=1)
    x,y = circle(z, r46)
    plt.plot(x,y, '--', lw=1)
    x,y = circle(pi/2, r, -r)
    plt.plot(x,y, ':', lw=1)
    plt.xticks(())
    plt.yticks(())
    plt.axis('equal')
    plt.axis([-1, 1, -1, 1])
    plt.text(-1.3, -1, r'$\zeta = %d^\circ$'%zt)

plt.tight_layout()
plt.savefig('fig10.png')
plt.show()
