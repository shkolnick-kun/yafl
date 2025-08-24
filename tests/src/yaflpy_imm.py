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
import pyximport
import scipy.stats
import sys
import time

from yaflpy import Bierman as KF
from yaflpy import IMMEstimator

def _ca(x, dt, **fx_args):
    x = x.copy()
    x[0] += (x[1] + 0.5 * x[2] * dt) * dt
    x[1] += x[2] * dt
    return x    

def _jca(x, dt, **fx_args):
    F = np.array([
        [1., dt, dt*dt*0.5],
        [0., 1.,    dt    ],
        [0., 0.,    1.    ],
        ])
    return F


def _cv(x, dt, **fx_args):
    x = x.copy()
    x[0] += x[1] * dt
    return x    

def _jcv(x, dt, **fx_args):
    F = np.array([
        [1., dt, 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        ])
    return F

def _hx(x, **hx_args):
    return np.array([x[0]])

def _jhx(x, **hx_args):
    H = np.array([[1., 0., 0.]])
    return H


_dt = 0.01
STD = 0.01

N = 1000
time  = []
clean = []
noisy = []

#CA
t = 0.
tgt = np.array([0., 0., 0.1])
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)
#CV
tgt[2] = 0.
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)
#CA
tgt[2] = -0.1
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

plt.plot(time, noisy, time, clean)
plt.show()

kf = KF(3, 1, _dt, _cv, _jcv, _hx, _jhx)
kf.Dp *= .0001
kf.Dq *= .0000001
kf.Dr = STD*STD
kf.x[0] = 0.
kf.x[1] = 0.
kf.x[2] = 0.

out = []
for i,z in enumerate(clean):
     kf.predict()
     kf.update(z)
     out.append(kf.x[0])    

plt.plot(time, out, time, clean)
plt.show()

cv = KF(3, 1, _dt, _cv, _jcv, _hx, _jhx)
cv.Dp *= .0001
cv.Dq *= .0001
cv.Dr = STD*STD
cv.x[0] = 0.
cv.x[1] = 0.
cv.x[2] = 0.

ca = KF(3, 1, _dt, _ca, _jca, _hx, _jhx)
ca.Dp *= .0001
ca.Dq *= .0001
ca.Dr = STD*STD
ca.x[0] = 0.
ca.x[1] = 0.
ca.x[2] = 0.

mu = np.array([0.5, 0.5])
M  = np.array([[0.95, 0.05],
               [0.05, 0.95]])

imm = IMMEstimator([cv, ca], mu, M, _dt)

out = []
mu  = []
cvl  = []
cal  = []


for i,z in enumerate(clean):
     imm.predict()
     imm.update(np.array([z]))
     out.append(imm.x[0])
     mu.append(imm.mu.copy())
     cvl.append(cv.l)
     cal.append(ca.l)

plt.plot(time, out, time, clean)
plt.show()

plt.plot(mu)
plt.show()

# plt.plot(t, noisy - kf_out)
# plt.show()

# plt.plot(t, clean - kf_out)
# plt.show()

# plt.plot(clean[:,0], clean[:,1], kf_out[:,0], kf_out[:,1])
# plt.show()

# plt.plot(noisy[:,0], noisy[:,1], kf_out[:,0], kf_out[:,1])
# plt.show()

# plt.plot(t, noisy[:,1], t, kf_out[:,1], t, clean[:,1])
# plt.show()

# plt.plot(t, noisy[:,0], t, kf_out[:,0], t, clean[:,0])
# plt.show()



print('Done!')
