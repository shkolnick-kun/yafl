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
import sys
sys.path.append("../yafl")

import numpy as np

from filterpy.kalman import ExtendedKalmanFilter
from filterpy.kalman import UnscentedKalmanFilter
from filterpy.kalman import JulierSigmaPoints

from ab_tests import *
from case1    import *
from case1    import _fx, _jfx, _hx, _jhx

from yaflpy import _mwgs, _ruv, _rum, _set_u

class SeqEKF(ExtendedKalmanFilter):

    def shx(self, x):
        return _ruv(self.Ur, _hx(x))[1]

    def sjhx(self, x):
        sjhx = _jhx(x)
        _rum(sjhx, self.Ur)
        return sjhx

    def update(self, z):
        # Decorrelate noice z
        zd = _ruv(self.Ur, z)[1]
        super().update(zd, self.sjhx, self.shx)
        _, self.Up, self.Dp = _mwgs(np.linalg.cholesky(self.P), np.ones(self.dim_x))
        _, self.Uq, self.Dq = _mwgs(np.linalg.cholesky(self.Q), np.ones(self.dim_x))

def new_seq_ekf(std):
    sz = 4

    ekf = SeqEKF(sz, 2)
    ekf.x = np.array([0, 0.3, 0, 0])
    ekf.F = _jfx(ekf.x, 1.0)

    up = _set_u(1e-8 * np.ones(((sz*(sz-1))//2,)))[1]
    dp = 0.00001 * np.eye(sz)
    ekf.P = up.dot(dp.dot(up.T))

    uq = up.copy()
    dq = 1e-6 * np.eye(sz)
    ekf.Q = uq.dot(dq.dot(uq.T))

    ekf.Ur = np.array([0.5])
    ekf.Dr = std * std * np.array([0.75, 1.0])

    #ur = _set_u(ekf.Ur)[1]
    #dr = np.diag(ekf.Dr)

    ekf.R = np.diag(ekf.Dr) #ur.dot(dr.dot(ur.T))

    return ekf
