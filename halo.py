# referece: R. Greenler "Rainbows, Halos and Glories" chapter 2

import numpy as np
from numpy.random import random

sqrt_3 = np.sqrt(3)
R = np.array([[1,-sqrt_3,0], [sqrt_3,1,0], [0,0,2]])/2
R = np.array([np.eye(3), R, R.T]) # three side faces

def halo_(s, A, io=0, h=10, m=1.31, outI=True):
    """
    s: unit vector to sun's direction (shape(3,))
    A: crystal orientation matrix (shape(3,3,N))
       where N is number of (random) samples
    A[2]: unit vector along prism axis (shape(3,N))
    A[0]: unit normal to one of side face (shape(3,N))
    A[1]: A[2] cross A[0] (so that det(A) == +1)
    io: 0,1,2 specifies incident and outgoing surface
       if io==0, (in,out) = (side,side),
       if io==1, (in,out) = (side,bottom),
       if io==2, (in,out) = (top,side).
    h: (prism height)/(hexagonal edge) (scalar or shape(N,))
    m: index of refraction (scalar)
    outI: output intensity or not
    return e,I where
       e = -(unit vector of outgoing ray),
       I = intensity of outgoing ray (0<=I<=1),
       e.shape=(3,3,N), I.shape=(3,N),
       e[:,i,:] and I[i,:] corresponds to i-th side face.
    For totally reflected rays, e=nan, I=-inf.
    If outI is False, return only e.
    Assume A is numpy array
    """
    s = np.broadcast_to(s, A.T.shape).T
    a1,a2,a3 = np.einsum('lik,kj...->ijl...', R, A)
    if io>=2: a1,a2,a3 = a3,-a2,a1
    clp = np.sum(a1*s, axis=0) # \pm cos(incident angle)
    b1,b2 = np.where(clp>0, (a1,a2), (-a1,-a2))

    c = np.cross(s, b1, axis=0)
    c /= np.linalg.norm(c, axis=0)
    if io: a3,b2 = -b2,a3

    ctu = np.sum(c*a3, axis=0) # \pm cos(skew angle)
    d = np.where(ctu>0, -b2, b2)
    if io==0: d = (b1 + sqrt_3*d)/2

    clp = np.abs(clp)
    alp = np.arccos(clp) # incident angle
    bet1 = np.arcsin(np.sin(alp)/m) # snell's law
    dlt1 = alp - bet1 # deflection angle
    s1 = np.cross(c, s, axis=0)
    r = s*np.cos(dlt1) + s1*np.sin(dlt1) # refracted ray

    bet2 = np.arccos(np.sum(r*d, axis=0))
    msb = m*np.sin(bet2)
    msb[msb>1] = np.nan # total reflection
    gam = np.arcsin(msb) # snell's law
    dlt2 = gam - bet2 # deflection angle

    g = np.cross(d, r, axis=0)
    g /= np.linalg.norm(g, axis=0)
    r1 = np.cross(g, r, axis=0)
    e = r*np.cos(dlt2) + r1*np.sin(dlt2) # exiting ray
    if not outI: return e

    # intensity
    if   io==0: a,b = sqrt_3, 1/h
    elif io==1: a,b = sqrt_3/h, h
    else: a,b = h/sqrt_3, sqrt_3
    ctu = np.abs(ctu)
    s = ctu*np.tan(bet1)*a
    t = np.tan(np.arccos(ctu))*b
    st = s*t; st2 = s*(1-st/2)
    if io==0: # geometrical factor
        I = (3-s)/2*(1-(t+st)/4)
        J = np.where(st<t+2, (1-(st-t)/2)**2/t/2, 0)
        J = np.where(s<1, 2/t, J)
    else: # approximate hexagon by rectangle of 1 x sqrt3
        I,J = 1-t/2, 1/2/t
    J = np.where(st<1, st2, J)
    I = np.where(s<1, st2, I)
    I = np.where(t<1, I, J)
    if io>=2: I /= a
    I *= clp # cross section
    # Fresnel reflection formula
    R1 = (np.sin(alp-bet1)/np.sin(alp+bet1))**2
    R1+= (np.tan(alp-bet1)/np.tan(alp+bet1))**2
    R2 = (np.sin(gam-bet2)/np.sin(gam+bet2))**2
    R2+= (np.tan(gam-bet2)/np.tan(gam+bet2))**2
    I *= (1-R1/2)*(1-R2/2) # neglect polarization of refracted light
    I[np.isnan(I)] = -np.inf # total reflection

    return e,I

def halo(sun_elev, h=10, tilt_max=np.pi/2, spin_max=np.pi,
         io=0, N=10000, m=1.31, LOS=(0,0), out3=False, Imax=1):
    """
    sun_elev: elevation of sun / radian (scalar)
    h: (prism height)/(hexagonal edge) (scalar)
       h>1 for pencil crystal, h<1 for plate crystal
    tilt_max: maximum tilt angle (scalar)
       if h>1, tilt==0 when prism axis is horizontal
       else,   tilt==0 when prism axis is vertical
    spin_max: maximum spin angle (scalar)
       e.g., spin_max==0 for lowitz arc
    io: 0,1,2 (int or tuple of int)
      specifies incident and outgoing surface
       if io==0, (in,out) = (side,side),
       if io==1, (in,out) = (side,bottom),
       if io==2, (in,out) = (top,side).
    N: number of random samples (scalar)
    m: index of refraction (scalar)
    LOS: observer's line of sight (th,ph) where
       th: elevation of LOS from sun_elev / rad,
       ph: azimuthal angle from sun-zenith plane / rad,
       default (0,0) is to sun's direction.
       If LOS is scalar, ph is set to zero
    out3: specifies shape of output arrays (boolean)
    Imax: maximum intensity so that if ray intensity is
       larger than random value from [0,Imax), then plot.
       If Imax is None, Imax is set to max(ray intensity)
    return x,y = r*cos(th), r*sin(th) where,
       r: angle from LOS to crystal / rad
       th: atan2 angle from rightward horizontal
    For invisible rays (dim or totally reflected), x=y=nan.
    If out3 is False, x.shape=y.shape=(3*N)
    else x.shape=y.shape(3,N) where
       x[i,:] and y[i,:] corresponds to i-th side face.
    If io is tuple, newaxis is inserted at axis==0
       of x,y and thier shape[0] becomes len(io)
    """
    czt, szt = np.cos(sun_elev), np.sin(sun_elev)
    s = np.r_[czt, 0, szt] # direction to sun

    phi = random(N)*2*np.pi # random azimuth
    psi = spin_max*(2*random(N) - 1) # random spin
    cph,sph = np.cos(phi), np.sin(phi)
    cps,sps = np.cos(psi), np.sin(psi)
    # random tilt
    if h>1: cth = np.sin(tilt_max)*random(N) # pencil
    else: cth = 1 - (1 - np.cos(tilt_max))*random(N) # plate
    sth = np.sqrt(1 - cth**2)
    # ZYZ rotation
    a3 = [cph*sth, sph*sth, cth]
    a1 = [cph*cth*cps - sph*sps,  sph*cth*cps + cph*sps, -sth*cps]
    a2 = np.cross(a3, a1, axis=0)
    A = np.asarray([a1,a2,a3])
    
    if np.isscalar(io): # simulate
        e,I = halo_(s, A, io, h, m)
    else:
        eI = [halo_(s,A,i,h,m) for i in io]
        e = np.moveaxis([e[0] for e in eI], 0, -1)
        I = np.asarray([I[1] for I in eI])

    th,ph = (LOS,0) if np.isscalar(LOS) else LOS
    ct,cp = np.cos(sun_elev + th), np.cos(ph)
    st,sp = np.sin(sun_elev + th), np.sin(ph)
    o3 = np.r_[ct*cp, ct*sp, st] # observer's LOS
    o1 = np.r_[sp, -cp, 0] # righward horizontal
    o2 = np.cross(o1,o3)

    e = np.einsum('ki,i...->k...', [o1,o2,o3], e)
    r,xy = np.arccos(e[2]), e[:2]
    xy *= r/np.linalg.norm(xy, axis=0)

    if not np.isscalar(io):
        xy = np.moveaxis(xy, -1, 1)

    if Imax is None: Imax = np.max(I)
    dim = (I < Imax*random(I.shape))
    xy[:,dim] = np.nan # invisibly dim

    if out3: return xy
    else: return xy.reshape(xy.shape[:-2] + (-1,))

def min_deflection(wedge=np.pi/3, m=1.31):
    """ minimum deflection angle
    wedge: angle between incident and exit surfaces
    m: index of refraction
    """
    return 2*np.arcsin(m*np.sin(wedge/2)) - wedge

def circle(c, r, LOS=(0,0), N=100):
    """ circle on sphere (for plot utility)
    c: elevation of center of circle / rad (scalar)
    r: radius of circle / rad (scalar)
    LOS: observer's line of sight (th,ph) where
       th: elevation of LOS from center / rad,
       ph: azimuthal angle from center-zenith plane / rad,
       default (0,0) is to center of circle.
       If LOS is scalar, ph is set to zero
    N: length of ouput arrays
    return x,y on the same plane as output of halo()
       x.shape = y.shape = (N,)
    """
    cr,sr = np.cos(r), np.sin(r)
    cc,sc = np.cos(c), np.sin(c)
    e3 = np.r_[cc, 0, sc][:,np.newaxis]
    e1 = np.r_[0, -1, 0][:,np.newaxis]
    e2 = np.cross(e1, e3, axis=0)
    t = np.linspace(0, 2*np.pi, N)
    ct,st = np.cos(t), np.sin(t)
    r = e3*cr + (e1*ct + e2*st)*sr

    th,ph = (LOS,0) if np.isscalar(LOS) else LOS
    ct,cp = np.cos(c+th), np.cos(ph)
    st,sp = np.sin(c+th), np.sin(ph)
    e3 = np.r_[ct*cp, ct*sp, st]
    e1 = np.r_[sp, -cp, 0]
    e2 = np.cross(e1,e3)
    r = np.dot([e1,e2,e3], r)
    r,xy = np.arccos(r[2]), r[:2]
    xy *= r/np.linalg.norm(xy, axis=0)
    return xy
