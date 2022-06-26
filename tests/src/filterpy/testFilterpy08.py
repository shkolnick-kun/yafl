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

from filterpy.kalman import ExtendedKalmanFilter
from UDEKF import UDExtendedKalmanFilter
#------------------------------------------------------------------------------
N = 10000
STD = 100.

#------------------------------------------------------------------------------
clean, noisy, t = case_data(N, STD)

#------------------------------------------------------------------------------
class A(ExtendedKalmanFilter):
    def update(self, z):
        super().update(z, self.jhx, self.hx)

a = A(4,2)
a.x = np.array([0, 0.3, 0, 0])
a.F = _jfx(a.x, 1.0)
a.P *= 0.0001
a.Q *= 1e-6
a.R *= STD * STD
a.hx  = _hx
a.jhx = _jhx

#------------------------------------------------------------------------------
class B(UDExtendedKalmanFilter):
    def update(self, z):
        super().update(z, self.jhx, self.hx)


b = B(4,2)
b.x = np.array([0, 0.3, 0, 0])
b.F = _jfx(a.x, 1.0)
b.P *= 0.0001
b.Q *= 1e-6
b.R *= STD * STD
b.hx  = _hx
b.jhx = _jhx

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
