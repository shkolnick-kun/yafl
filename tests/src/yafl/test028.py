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

from yaflpy import Joseph  as A
from yaflpy import Bierman as B

#------------------------------------------------------------------------------
N = 10000
STD = 100.

#------------------------------------------------------------------------------
clean, noisy, t = case_data(N, STD)

#------------------------------------------------------------------------------
a = case_ekf(A, STD)
a.rff = 0.0001
a.qff = 0.0001

b = case_ekf(B, STD)
b.rff = 0.0001
b.qff = 0.0001

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
