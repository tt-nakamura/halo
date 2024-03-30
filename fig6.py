# Lowitz arc

import matplotlib.pyplot as plt
from scipy.constants import degree
from halo import *
from math import pi

zeta = [0,10,20,30,40,50,60,70,75,80] # sun elevation / deg
psi0 = 0 # max spin angle

r22 = min_deflection()

plt.figure(figsize=(5,10))

for i,zt in enumerate(zeta):
    plt.subplot(5,2,i+1)
    z = zt*degree; r = pi/2 - z
    x,y = halo(z, h=0.1, spin_max=psi0)
    plt.plot(x,y, '.', ms=1)
    x,y = circle(z, r22)
    plt.plot(x,y, ':', lw=1)
    x,y = circle(pi/2, pi/2, -r)
    plt.plot(x,y, '--', lw=1)
    x,y = circle(pi/2, r, -r)
    plt.plot(x,y, ':', lw=1)
    plt.xticks(())
    plt.yticks(())
    plt.axis('equal')
    plt.axis([-0.7, 0.7, -0.7, 0.7])
    plt.text(-0.9, -0.7, r'$\zeta = %d^\circ$'%zt)

plt.tight_layout()
plt.savefig('fig6.png')
plt.show()
