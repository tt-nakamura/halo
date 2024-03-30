# sun dogs

import matplotlib.pyplot as plt
from scipy.constants import degree
from halo import *
from math import pi

zeta = [0, 15, 30, 45] # sun elevation / deg
th0 = 3*degree # max tilt of crystals

r22 = min_deflection()

plt.figure(figsize=(5,3.75))

for i,zt in enumerate(zeta):
    plt.subplot(2,2,i+1)
    z = zt*degree; r = pi/2 - z
    x,y = halo(z, h=0.1, tilt_max=th0)
    plt.plot(x,y, '.', ms=1)
    x,y = circle(z, r22)
    plt.plot(x,y, '--', lw=1)
    x,y = circle(pi/2, r, -r)
    plt.plot(x,y, ':', lw=1)
    plt.xticks(())
    plt.yticks(())
    plt.axis('equal')
    plt.axis([-0.6, 0.6, -0.6, 0.6])
    plt.text(-0.8, -0.55, r'$\zeta = %d^\circ$'%zt)

plt.tight_layout()
plt.savefig('fig3.png')
plt.show()
