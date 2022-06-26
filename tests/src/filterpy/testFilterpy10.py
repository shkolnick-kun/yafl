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
import time
import matplotlib.pyplot as plt
import numpy as np

import init_tests

from ab_tests import *
from case1    import *
from case1    import _fx, _jfx, _hx, _jhx

from yaflpy import _mwgs, _ruv, _rum, _set_u

from filterpy.kalman import ExtendedKalmanFilter
from yaflpy          import Joseph as B
#------------------------------------------------------------------------------
N = 10000
STD = 100.

#------------------------------------------------------------------------------
clean, noisy, t = case_data(N, STD)

#------------------------------------------------------------------------------
class A(ExtendedKalmanFilter):
    def update(self, z):
        super().update(z, self.jhx, self.hx)
        _, self.Up, self.Dp = _mwgs(np.linalg.cholesky(self.P), np.ones(self.dim_x))
        _, self.Uq, self.Dq = _mwgs(np.linalg.cholesky(self.Q), np.ones(self.dim_x))
        _, self.Ur, self.Dr = _mwgs(np.linalg.cholesky(self.R), np.ones(self.dim_z))
        _, self.y = _ruv(self.Ur, self.y)

a = A(4,2)
a.x = np.array([0, 0.3, 0, 0])
a.F = _jfx(a.x, 1.0)

aup = np.array([[1, 1e-8, 1e-8, 1e-8],
                [0, 1,    1e-8, 1e-8],
                [0, 0,    1,    1e-8],
                [0, 0,    0,    1]])
adp = 0.00001 * np.eye(4)
a.P = aup.dot(adp.dot(aup.T))

auq = np.array([[1, 1e-8, 1e-8, 1e-8],
                [0, 1,    1e-8, 1e-8],
                [0, 0,    1,    1e-8],
                [0, 0,    0,    1]])
adq = 1e-6 * np.eye(4)
a.Q = auq.dot(adq.dot(auq.T))


#aur = np.array([[1, 0.5],
#                [0, 1   ]])

#adr = STD * STD * np.eye(2)
#adr[0,0] *= 0.75
#a.R = aur.dot(adr.dot(aur.T))

a.R = STD * STD * np.eye(2)

a.hx  = _hx
a.jhx = _jhx

#------------------------------------------------------------------------------
b = case_ekf(B, STD)

b.Ur *= 0.
b.Dr[0] *= 4./3.
#------------------------------------------------------------------------------
start = time.time()

dxp = []
dxu = []
dy = []

for z in noisy:
    #b.x = a.x

    #Не влияет на расхождение
    #_, b.Up, b.Dp = _mwgs(np.linalg.cholesky(a.P), np.ones(a.dim_x))

    #Тут всё 1 в 1
    a.predict()
    b.predict()
    dxp.append(2. * norm(a.x - b.x) / norm(a.x + b.x))

    #Расхождение тут!
    a.update(z)
    b.update(z)
    dxu.append(2. * norm(a.x - b.x) / norm(a.x + b.x))
    dy.append(2. * norm(a.y - b.y) / norm(a.y + b.y))

end = time.time()
print(end - start)

plt.plot(dxp)
plt.show()

plt.plot(dxu)
plt.show()

plt.plot(dy)
plt.show()

#------------------------------------------------------------------------------

