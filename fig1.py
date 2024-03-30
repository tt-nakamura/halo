import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import degree
from halo import halo_

N = 100
alp = [np.linspace(10, 85, N), # angle of incidence
       np.linspace(55, 89, N)]
tau = np.r_[0, 10, 20, 30] # skew angle
A = np.eye(3)

plt.figure(figsize=(5,5))

for i in range(2):
    plt.subplot(2,1,i+1)
    for t in tau:
        td = t*degree
        dlt = []
        ct,st = np.cos(td), np.sin(td)
        for a in alp[i]:
            a *= degree
            ca,sa = np.cos(a), np.sin(a)
            if i: s = np.r_[ca*ct, st, sa*ct]
            else: s = np.r_[ca*ct, sa*ct, st]
            e = halo_(s, A, io=i, outI=False)
            d = np.arccos(np.dot(e[:,0],s))
            dlt.append(d/degree) # deflection angle

        plt.plot(alp[i], dlt, label=r'$\tau=%d^\circ$'%t)

    plt.ylabel(r'$\delta$ = deflection angle / deg')
    plt.legend()

plt.xlabel(r'$\alpha$ = angle of incidence / deg')
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
