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

N = 1000

aa = []
dd = []
ss = []

for i in range(1,N):
    rnr = int(np.random.uniform(low=2, high=i+2))
    rnc = int(np.random.uniform(low=2, high=i+2))
    sz  = int(np.random.uniform(low=2, high=min(rnr, rnc)))
    r  = int(np.random.uniform(low=0, high=rnr-sz-1))
    c  = int(np.random.uniform(low=0, high=rnc-sz-1))

    res = np.random.uniform(low=1e-15, high=1.0, size=(rnr,rnc))
    a = np.random.uniform(low=1e-15, high=1.0, size=(sz,))
    b = np.random.uniform(low=1e-15, high=1.0, size=(sz,))

    p = res.copy()
    p[r:r+sz, c:c+sz] -= np.outer(a,b)
    s = yaflpy._bsub_vvt(res, a, b, r, c)

    aa.append(np.linalg.norm(p))
    ss.append(s)
    dd.append(np.linalg.norm(p - res)/np.linalg.norm(p))

print(np.min(aa))
print(np.max(aa))
print(np.average(aa))
print(np.median(aa))

print(np.min(dd))
print(np.max(dd))
print(np.average(dd))
print(np.median(dd))

