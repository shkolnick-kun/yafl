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
import h5py
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

    for z in measurements:
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
def yafl_record_test_data(a, clear, noisy, t, fname):
    rp = []
    ru = []
    x  = []
    y  = []
    up = []
    dp = []
    uq = []
    dq = []
    ur = []
    dr = []
    for z in noisy:
        rp.append(a.predict())
        ru.append(a.update(z))

        x.append(a.x.copy())
        y.append(a.y.copy())

        up.append(a.Up.copy())
        dp.append(a.Dp.copy())

        uq.append(a.Uq.copy())
        dq.append(a.Dq.copy())

        ur.append(a.Ur.copy())
        dr.append(a.Dr.copy())

    rp = np.array(rp)
    ru = np.array(ru)
    x  = np.array(x)
    y  = np.array(y)
    up = np.array(up)
    dp = np.array(dp)
    uq = np.array(uq)
    dq = np.array(dq)
    ur = np.array(ur)
    dr = np.array(dr)

    with h5py.File(fname, 'w') as h5f:
        h5f.create_dataset('clear', data=clear)
        h5f.create_dataset('noisy', data=noisy)
        h5f.create_dataset('t', data=noisy)

        h5f.create_dataset('rp', data=rp)
        h5f.create_dataset('ru', data=ru)

        h5f.create_dataset('x', data=x)
        h5f.create_dataset('y', data=y)

        h5f.create_dataset('up', data=up)
        h5f.create_dataset('dp', data=dp)

        h5f.create_dataset('uq', data=uq)
        h5f.create_dataset('dq', data=dq)

        h5f.create_dataset('ur', data=ur)
        h5f.create_dataset('dr', data=dr)

#------------------------------------------------------------------------------
def yafl_file_test(b, fname):

    clear = None
    noisy = None
    t     = None

    rpa = None
    rua = None

    xa  = None
    ya  = None

    aup = None
    adp = None

    auq = None
    adq = None

    aur = None
    adr = None

    with h5py.File(fname, 'r') as h5f:
        print(h5f.keys())
        clear = h5f['clear'][:]
        noisy = h5f['noisy'][:]
        t     = h5f['t'][:]

        rpa = h5f['rp'][:]
        rua = h5f['ru'][:]

        xa  = h5f['x'][:]
        ya  = h5f['y'][:]

        aup = h5f['up'][:]
        adp = h5f['dp'][:]

        auq = h5f['uq'][:]
        adq = h5f['dq'][:]

        aur = h5f['ur'][:]
        adr = h5f['dr'][:]

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

    for i,z in enumerate(noisy):
        # Run B
        rpb.append(b.predict())
        rub.append(b.update(z))
        xb.append(b.x.copy())

        nx.append(2. * norm(xa[i] - b.x) / norm(xa[i] + b.x))
        ny.append(2. * norm(ya[i] - b.y) / norm(ya[i] + b.y))

        nup.append(2. * norm(aup[i] - b.Up) / norm(aup[i] + b.Up))
        ndp.append(2. * norm(adp[i] - b.Dp) / norm(adp[i] + b.Dp))

        nuq.append(2. * norm(auq[i] - b.Uq) / norm(auq[i] + b.Uq))
        ndq.append(2. * norm(adq[i] - b.Dq) / norm(adq[i] + b.Dq))

        nur.append(2. * norm(aur[i] - b.Ur) / norm(aur[i] + b.Ur))
        ndr.append(2. * norm(adr[i] - b.Dr) / norm(adr[i] + b.Dr))

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

    return clear,noisy,t, rpa,rua,xa, rpb,rub,xb, nup,ndp, nuq,ndq, nur,ndr, nx, ny

#------------------------------------------------------------------------------
def filterpy_ab_test(a, b, measurements):

    xa  = []
    xb  = []

    nP = []
    nQ = []
    nR = []

    nx  = []
    ny  = []

    for z in measurements:
        # Run A
        a.predict()
        a.update(z)
        xa.append(a.x.copy())

        # Run B
        b.predict()
        b.update(z)
        xb.append(b.x.copy())

        nx.append(2. * norm(a.x - b.x) / norm(a.x + b.x))
        ny.append(2. * norm(a.y - b.y) / norm(a.y + b.y))

        nP.append(2. * norm(a.P - b.P) / norm(a.P + b.P))
        nQ.append(2. * norm(a.Q - b.Q) / norm(a.Q + b.Q))
        nR.append(2. * norm(a.R - b.R) / norm(a.R + b.R))

    xa  = np.array(xa)
    xb  = np.array(xb)

    nP = np.array(nP)
    nQ = np.array(nQ)
    nR = np.array(nR)

    nx = np.array(nx)
    ny = np.array(ny)

    return xa,xb, nP,nQ,nR, nx,ny
