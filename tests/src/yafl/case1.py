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
import yaflpy
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

#Default noice model is Normal with Poisson outliers
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

#------------------------------------------------------------------------------
def case_data(n, std):
    clear = np.zeros((n, 2))
    noisy = np.zeros((n, 2))
    t     = np.zeros((n,), dtype=np.float)

    for i in range(1, len(clear)):
        clear[i] = clear[i-1] + np.array([2.,1.])
        noisy[i] = clear[i]   + np.random.normal(scale=std, size=2)
        t[i] = i

    return clear, noisy, t

#------------------------------------------------------------------------------
def _kf_init(kf, std):

    kf.x[0] = 0.
    kf.x[1] = 0.3

    kf.Up = 1e-8
    kf.Dp *= .00001

    kf.Uq = 1e-8
    kf.Dq *= 1.0e-6

    kf.Dr *= std * std
    kf.Dr[0] *= .75
    kf.Ur += 0.5

#------------------------------------------------------------------------------
ROBUST_EKF = [yaflpy.RobustBierman,
              yaflpy.RobustJoseph,
              yaflpy.AdaptiveRobustBierman,
              yaflpy.AdaptiveRobustJoseph]

def case_ekf(clsEKF, std):
    if clsEKF in ROBUST_EKF:
        kf = clsEKF(4, 2, 1., _fx, _jfx, _hx, _jhx, gz=_gz, gdotz=_gdotz)
    else:
        kf = clsEKF(4, 2, 1., _fx, _jfx, _hx, _jhx)
    _kf_init(kf, std)
    return kf

#------------------------------------------------------------------------------
ROBUST_UKF = [yaflpy.UnscentedRobustBierman,
              yaflpy.UnscentedAdaptiveRobustBierman]

def case_ukf(clsUKF, std, sp=None):

    if not sp:
        sp = SP(4, 0.0)

    if clsUKF in ROBUST_UKF:
        kf = clsUKF(4, 2, 1., _hx, _fx, sp, gz=_gz, gdotz=_gdotz)
    else:
        kf = clsUKF(4, 2, 1., _hx, _fx, sp)

    _kf_init(kf, std)
    return kf
