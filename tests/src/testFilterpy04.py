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

from ab_tests import *
from case1    import *
from case1    import _fx, _jfx, _hx, _jhx

from yaflpy import MWGS, RUV

from filterpy.kalman import UnscentedKalmanFilter
from filterpy.kalman import JulierSigmaPoints

from yaflpy          import Unscented as B
#------------------------------------------------------------------------------
N = 10000
STD = 100.

#------------------------------------------------------------------------------
clean, noisy, t = case_data(N, STD)

#------------------------------------------------------------------------------
class A(UnscentedKalmanFilter):
    def update(self, z):
        super().update(z)
        _, self.Up, self.Dp = MWGS(np.linalg.cholesky(self.P), np.ones(self.dim_x))
        _, self.Uq, self.Dq = MWGS(np.linalg.cholesky(self.Q), np.ones(self.dim_x))
        _, self.Ur, self.Dr = MWGS(np.linalg.cholesky(self.R), np.ones(self.dim_z))
        _, us, ds = MWGS(np.linalg.cholesky(self.S), np.ones(self.dim_z))
        _, self.y = RUV(us, self.y)

sp = JulierSigmaPoints(4, 0.0)
a = A(4, 2, 1.0, _hx, _fx, sp)
a.dim_x = 4
a.dim_z = 2

a.x = np.array([0, 0.3, 0, 0])

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


aur = np.array([[1, 0.5],
                [0, 1   ]])

adr = STD * STD * np.eye(2)
adr[0,0] *= 0.75
a.R = aur.dot(adr.dot(aur.T))

#------------------------------------------------------------------------------
b = case_ukf(B, STD)

#------------------------------------------------------------------------------
start = time.time()

rpa,rua,xa, rpb,rub,xb, nup,ndp, nuq,ndq, nur,ndr, nx, ny = yafl_ab_test(a, b, noisy)

end = time.time()
print(end - start)

#------------------------------------------------------------------------------
plt.plot(nup)
plt.show()
plt.plot(ndp)
plt.show()

plt.plot(nuq)
plt.show()
plt.plot(ndq)
plt.show()

plt.plot(nur)
plt.show()
plt.plot(ndr)
plt.show()

plt.plot(nx)
plt.show()
plt.plot(ny)
plt.show()

plt.plot(noisy[:,0], noisy[:,1], xa[:,0], xa[:,2])
plt.show()
plt.plot(clean[:,0], clean[:,1], xa[:,0], xa[:,2])
plt.show()
