# circumhorizontal arc

import matplotlib.pyplot as plt
from scipy.constants import degree
from halo import *
from math import pi

zeta = [60,63,70,80] # sun elevation / deg
th0 = 3*degree # max tilt of crystals

r46 = min_deflection(pi/2)

plt.figure(figsize=(5,3.75))

for i,zt in enumerate(zeta):
    plt.subplot(2,2,i+1)
    z = zt*degree
    x,y = halo(z, h=0.1, tilt_max=th0, io=1, LOS=-r46)
    plt.plot(x,y, '.', ms=1)
    x,y = circle(z, r46, -r46)
    plt.plot(x,y, '--', lw=1)
    plt.xticks(())
    plt.yticks(())
    plt.axis('equal')
    plt.axis([-0.9, 0.9, -0.9, 0.9])
    plt.text(-1.2, -0.85, r'$\zeta = %d^\circ$'%zt)

plt.tight_layout()
plt.savefig('fig8.png')
plt.show()
