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

from numpy.linalg import norm

#------------------------------------------------------------------------------
def yafl_ab_test(a, b, measurements):

    rpa = []
    rua = []
    xa  = []

    rpb = []
    rub = []
    xb  = []

    nup = []
    ndp = []
    nuq = []
    ndq = []
    nur = []
    ndr = []
    nx  = []
    ny  = []

    for i,z in enumerate(measurements):
        # Run A
        rpa.append(a.predict())
        rua.append(a.update(z))
        xa.append(a.x.copy())

        # Run B
        rpb.append(b.predict())
        rub.append(b.update(z))
        xb.append(b.x.copy())

        nx.append(2. * norm(a.x - b.x) / norm(a.x + b.x))
        ny.append(2. * norm(a.y - b.y) / norm(a.y + b.y))

        nup.append(2. * norm(a.Up - b.Up) / norm(a.Up + b.Up))
        ndp.append(2. * norm(a.Dp - b.Dp) / norm(a.Dp + b.Dp))

        nuq.append(2. * norm(a.Uq - b.Uq) / norm(a.Uq + b.Uq))
        ndq.append(2. * norm(a.Dq - b.Dq) / norm(a.Dq + b.Dq))

        nur.append(2. * norm(a.Ur - b.Ur) / norm(a.Ur + b.Ur))
        ndr.append(2. * norm(a.Dr - b.Dr) / norm(a.Dr + b.Dr))

    rpa = np.array(rpa)
    rua = np.array(rua)
    xa  = np.array(xa)

    rpb = np.array(rpb)
    rub = np.array(rub)
    xb  = np.array(xb)

    nup = np.array(nup)
    ndp = np.array(ndp)

    nuq = np.array(nuq)
    ndq = np.array(ndq)

    nur = np.array(nur)
    ndr = np.array(ndr)

    nx = np.array(nx)
    ny = np.array(ny)

    return rpa,rua,xa, rpb,rub,xb, nup,ndp, nuq,ndq, nur,ndr, nx, ny

#------------------------------------------------------------------------------

