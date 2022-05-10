#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright 2022 anonimous <shkolnick-kun@gmail.com> and contributors.

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
from yaflpy import JulierSigmaPoints as SP
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
def case_data(n, std):
    clean = np.zeros((n, 2))
    noisy = np.zeros((n, 2))
    t     = np.zeros((n,), dtype=np.float)

    for i in range(1, len(clean)):
        clean[i] = clean[i-1] + np.array([1.,1.])
        noisy[i] = clean[i]   + np.random.normal(scale=std, size=2)
        t[i] = i

    return clean, noisy, t

#------------------------------------------------------------------------------
def _kf_init(kf, std):
    kf.x[0] = 0.
    kf.x[1] = 0.3
    kf.Dp *= .00001
    kf.Dq *= 1.0e-6

    kf.Dr *= std * std
    kf.Dr[0] *= .75
    kf.Ur += 0.5

#------------------------------------------------------------------------------
def case_ekf(clsEKF, std):
    kf = clsEKF(4, 2, 1., _fx, _jfx, _hx, _jhx)
    _kf_init(kf, std)
    return kf

#------------------------------------------------------------------------------
def case_ukf(clsUKF, std, sp=None):

    if not sp:
        sp = SP(4, 0.0)

    kf = clsUKF(4, 2, 1., _hx, _fx, sp)
    _kf_init(kf, std)
    return kf
