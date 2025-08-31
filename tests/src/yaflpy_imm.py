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


_dt = 0.1
STD = 1.

N = 100
time  = []
clean = []
noisy = []

#CV
t = 0.
tgt = np.array([-0., 0., 0.])
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CV
tgt[2] = 0.1
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
tgt[2] = -0.2
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
tgt[2] = 0.1
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

plt.plot(time, noisy, time, clean)
plt.show()

kf = KF(3, 1, _dt, _ca, _jca, _hx, _jhx)
kf.Dp *= .1
kf.Dq *= .000001
kf.Dr = STD*STD
kf.x[0] = 0.
kf.x[1] = 0.
kf.x[2] = 0.

out = []
for i,z in enumerate(noisy):
     kf.predict()
     kf.update(z)
     out.append(kf.x[0])    

plt.plot(time, out, time, clean)
plt.show()

def _zrf(a,b):
    return a - b

cv = KF(3, 1, _dt, _cv, _jcv, _hx, _jhx, residual_z=_zrf)
cv.Dp *= .1
cv.Dq *= .000001
cv.Dr = STD*STD
cv.x[0] = 0.
cv.x[1] = 0.
cv.x[2] = 0.

ca = KF(3, 1, _dt, _ca, _jca, _hx, _jhx, residual_z=_zrf)
ca.Dp *= .1
ca.Dq *= .000001
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
cvl = []
cal = []
cvx = []
cax = []


for i,z in enumerate(noisy):
     imm.predict(_dt)
     imm.update(np.array([z]))
     out.append(imm.x[0])
     mu.append(imm.mu.copy())
     cvl.append(cv.l)
     cal.append(ca.l)
     cvx.append(cv.x.copy())
     cax.append(ca.x.copy())

plt.plot(time, out, time, clean)
plt.show()

plt.plot(time, mu)
plt.show()

plt.plot(time, cvl)
plt.show()

plt.plot(time, cal)
plt.show()

print('Done!')
