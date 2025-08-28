#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright 2021 anonimous <shkolnick-kun@gmail.com> and contributors.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    See the License for the specific language governing permissions
"""
import matplotlib.pyplot as plt
import numpy as np
import time

"""
import pyximport
import sys

sys.path.insert(0,'../../src/yaflpy')

pyximport.install(
    build_dir='../projects/obj',
    pyimport=True,
    reload_support=True
    )
"""

#from yaflpy import MerweSigmaPoints as SP
from yaflpy import JulierSigmaPoints as SP

#from yaflpy import UnscentedRobustBierman as KF
from yaflpy import UnscentedAdaptiveRobustBierman as KF


def _fx(x, dt, **fx_args):
    x = x.copy()
    x[0] += x[1] * dt
    x[2] += x[3] * dt
    return x


def _hx(x, **hx_args):
    if hx_args:
        print(hx_args)
    return np.array([x[0], x[2]])

def _gz(beta):
    # +- 3*sigma
    if 3.0 >= np.abs(beta):
        return float(beta)

    # +- 6*sigma - uncertain measurements
    if 6.0 >= np.abs(beta):
        return float(beta/3.0)

    # outliers
    return float(np.sign(beta))

def _gdotz(beta):
    # +- 3*sigma
    if 3.0 >= np.abs(beta):
        return 1.0

    # +- 6*sigma - uncertain measurements
    if 6.0 >= np.abs(beta):
        return 1.0/3.0

    # outliers
    return 0.0

def _zrf(a,b):
    #print(a - b)
    return a - b

STD = 100.

#sp = SP(4, 0.1, 2., 0)
sp = SP(4, 0.0)
kf = KF(4, 2, 1., _hx, _fx, sp, gz=_gz, gdotz=_gdotz)
#kf = KF(4, 2, 1., _hx, _fx, sp, residual_z=_zrf)

kf.x[0] = 0.
kf.x[1] = 0.3
kf.Dp *= .00001
kf.Dq *= 1.0e-8
#This is robust filter, so no square here
kf.Dr *= STD*STD*0.001

kf.Dr[0] *= .75
kf.Ur += .5

# sp = MerweScaledSigmaPoints(4, alpha=.1, beta=2., kappa=0)
# kf = UnscentedKalmanFilter(4, 2, dt=1., fx=_fx, hx=_hx, points=sp)

# kf.P *= .1
# kf.Q *= 1e-8
# kf.R *= STD*STD

kf.rff = 0.001
kf.qff = 0.0001

print(kf.x)

N = 10000

clean = np.zeros((N, 2))
noisy = np.zeros((N, 2))
t     = np.zeros((N,), dtype=np.float)
# for i in range(1, len(clean)//2):
#     clean[i] = clean[i-1] + np.array([1.5,1.])
#     noisy[i] = clean[i]   + np.random.normal(scale=STD, size=2)
#     t[i] = i

# for i in range(i, len(clean)):
#     clean[i] = clean[i-1] + np.array([1.,1.5])
#     noisy[i] = clean[i]   + np.random.normal(scale=STD, size=2)
#     t[i] = i

for i in range(1, len(clean)):
    clean[i] = clean[i-1] + np.array([1.,1.])
    noisy[i] = clean[i]   + np.random.normal(scale=STD, size=2)
    t[i] = i

kf_out = np.zeros((N, 2))



start = time.time()
for i, z in enumerate(noisy):
    kf.predict()
    kf.update(z)
    kf_out[i] = kf.zp
end = time.time()
print(end - start)



plt.plot(t, noisy - kf_out)
plt.show()

plt.plot(t, clean - kf_out)
plt.show()

plt.plot(clean[:,0], clean[:,1], kf_out[:,0], kf_out[:,1])
plt.show()

plt.plot(noisy[:,0], noisy[:,1],  kf_out[:,0], kf_out[:,1])
plt.show()

plt.plot(t, noisy[:,1], t, kf_out[:,1], t, clean[:,1])
plt.show()

plt.plot(t, noisy[:,0], t, kf_out[:,0], t, clean[:,0])
plt.show()


print('Done!')
