# circumscribed halo

import matplotlib.pyplot as plt
from scipy.constants import degree
from halo import *
from math import pi

zeta = [0,5,10,15,20,25,30,40,50,70] # sun elevation / deg
th0 = 0 # max tilt of crystals

r22 = min_deflection()

plt.figure(figsize=(5,10))

for i,zt in enumerate(zeta):
    plt.subplot(5,2,i+1)
    z = zt*degree; r = pi/2 - z
    x,y = halo(z, tilt_max=th0)
    plt.plot(x,y, '.', ms=1)
    x,y = circle(z, r22)
    plt.plot(x,y, '--', lw=1)
    x,y = circle(pi/2, pi/2, -r)
    plt.plot(x,y, ':', lw=1)
    plt.xticks(())
    plt.yticks(())
    plt.axis('equal')
    plt.axis([-0.7, 0.7, -0.7, 0.7])
    plt.text(-0.9, -0.7, r'$\zeta = %d^\circ$'%zt)

plt.tight_layout()
plt.savefig('fig4.png')
plt.show()
