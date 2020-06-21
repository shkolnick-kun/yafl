#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright 2020 anonimous <shkolnick-kun@gmail.com> and contributors.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    See the License for the specific language governing permissions
"""
import numpy as np
import pyximport

pyximport.install(
    build_dir='../tests/projects/obj', 
    pyimport=True,
    reload_support=True, 
    language_level=3,
    setup_args={
        'include_dirs': ['./', '../tests/src'],
        }
    )

#from yakf_py import Bierman as KF
from yakf_py import Joseph as KF

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

#kf = FK(4, 2, 1., _fx, _jfx, _hx, _jhx, residual_z=_zrf)
kf = KF(4, 2, 1., _fx, _jfx, _hx, _jhx)
kf.x[0] = 0
kf.x[1] = 10
kf.Dp *= .1
kf.Dq *= 1.0e-6
kf.Dr *= 400
kf.Ur += 0.5

N = 1000

clean = np.zeros((N, 2))
noisy = np.zeros((N, 2))
t     = np.zeros((N,), dtype=np.float)
for i in range(1, len(clean)):
    clean[i] = clean[i-1] + np.array([1.,1.])
    noisy[i] = clean[i]   + np.random.normal(scale=20., size=2)
    t[i] = i

kf_out = np.zeros((N, 2))

import time

start = time.time()
for i, z in enumerate(noisy):
    kf.predict()
    kf.update(z)
    kf_out[i] = kf.x[::2]
end = time.time()
print(end - start)


import matplotlib.pyplot as plt
plt.plot(t, clean - kf_out)
plt.show()

plt.plot(clean[:,0], clean[:,1], kf_out[:,0], kf_out[:,1])
plt.show()
print('Done!')
