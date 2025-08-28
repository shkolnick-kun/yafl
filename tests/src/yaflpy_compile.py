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

from yaflpy import Bierman as KF
#from yaflpy import Joseph as KF
#from yaflpy import AdaptiveBierman as KF
#from yaflpy import AdaptiveJoseph as KF
#from yaflpy import DoNotUseThisFilter as KF

def _fx(x, dt, **fx_args):
    x = x.copy()
    x[0] += x[1] * dt
    x[2] += x[3] * dt
    return x

def _jfx(x, dt, **fx_args):
    F = np.array([
        [1., dt, 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., 1., dt],
        [0., 0., 0., 1.],
        ])
    return F

def _hx(x, **hx_args):
    if hx_args:
        print(hx_args)
    return np.array([x[0], x[2]])

def _jhx(x, **hx_args):
    H = np.array([
        [1., 0., 0., 0.],
        [0., 0., 1., 0.],
        ])
    return H

def _zrf(a,b):
    return a - b

STD = 100.

#kf = KF(4, 2, 1., _fx, _jfx, _hx, _jhx, residual_z=_zrf)
kf = KF(4, 2, 1., _fx, _jfx, _hx, _jhx)
kf.x[0] = 0.
kf.x[1] = 0.3
kf.Dp *= .00001
kf.Dq *= 1.0e-6
#This is robust filter, so no square here
kf.Dr *= STD*STD*0.000001
kf.Dr[0] *= .75
kf.Ur += 0.5

kf.rff = 0.001
kf.qff = 0.0001

N = 20000

#kf.chi2 = 8.807468393511947

clean = np.zeros((N, 2))
noisy = np.zeros((N, 2))
t     = np.zeros((N,), dtype=np.float)
# for i in range(1, len(clean)//2):
#     clean[i] = clean[i-1] + np.array([1.5,1.])
#     noisy[i] = clean[i]   + np.random.normal(scale=STD, size=2)
#     t[i] = i

# for i in range(i, len(clean)):
#     clean[i] = clean[i-1] + np.array([1.,10.])
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
    kf_out[i] = kf.x[::2]
end = time.time()
print(end - start)

plt.plot(t, noisy - kf_out)
plt.show()

plt.plot(t, clean - kf_out)
plt.show()

plt.plot(clean[:,0], clean[:,1], kf_out[:,0], kf_out[:,1])
plt.show()

plt.plot(noisy[:,0], noisy[:,1], kf_out[:,0], kf_out[:,1])
plt.show()

plt.plot(t, noisy[:,1], t, kf_out[:,1], t, clean[:,1])
plt.show()

plt.plot(t, noisy[:,0], t, kf_out[:,0], t, clean[:,0])
plt.show()

print('Done!')
