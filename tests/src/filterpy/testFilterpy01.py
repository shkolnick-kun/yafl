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

from UDEKF import UDExtendedKalmanFilter

from yaflpy import Bierman as B
#------------------------------------------------------------------------------
N = 10000
STD = 100.

#------------------------------------------------------------------------------
clean, noisy, t = case_data(N, STD)

#------------------------------------------------------------------------------
class A(UDExtendedKalmanFilter):
    def update(self, z):
        super().update(z, self.jhx, self.hx)
        self.y = np.dot(self.Dm, self.y)

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


aur = np.array([[1, 0.5],
                [0, 1   ]])

adr = STD * STD * np.eye(2)
adr[0,0] *= 0.75
a.R = aur.dot(adr.dot(aur.T))

a.hx  = _hx
a.jhx = _jhx

#------------------------------------------------------------------------------
b = case_ekf(B, STD)

#------------------------------------------------------------------------------
start = time.time()

xa,xb, nP,nQ,nR, nx,ny = filterpy_ab_test(a, b, noisy)

end = time.time()
print(end - start)

#------------------------------------------------------------------------------
plt.plot(nP)
plt.show()

plt.plot(nQ)
plt.show()

plt.plot(nR)
plt.show()

plt.plot(nx)
plt.show()

plt.plot(ny)
plt.show()

plt.plot(noisy[:,0], noisy[:,1], xa[:,0], xa[:,2])
plt.show()
plt.plot(clean[:,0], clean[:,1], xa[:,0], xa[:,2])
plt.show()
