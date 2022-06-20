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
    r = np.random.uniform(low=1e-15, high=1.0, size=(i,))
    a = np.random.uniform(low=1e-15, high=1.0, size=(i,))
    b = np.random.uniform(low=1e-15, high=1.0, size=(i,))

    p = r + a/b
    s = yaflpy._add_vrv(r, a, b)

    aa.append(np.linalg.norm(p))
    ss.append(s)
    dd.append(np.linalg.norm(p - r))

print(np.min(aa))
print(np.max(aa))
print(np.average(aa))
print(np.median(aa))

print(np.min(dd))
print(np.max(dd))
print(np.average(dd))
print(np.median(dd))
