#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright 2021 anonimous <shkolnick-kun@gmail.com> and contributors.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    See the License for the specific language governing permissions
"""
import matplotlib.pyplot as plt
import numpy as np

import numpy as np

from numpy.linalg import norm

##############################################################################
from yaflpy import NP_DTYPE
from yaflpy import Bierman as KFB
from yaflpy import IMMEstimator as IMMEstimatorB

###############################################################################
from filterpy.kalman import ExtendedKalmanFilter as KFA

class IMMEstimatorA(object):
    def __init__(self, filters, mu, M, HJacobian, Hx):
        self.filters   = filters
        self.mu        = mu
        self.M         = M
        self.HJacobian = HJacobian
        self.Hx        = Hx 
        self.cbar      = np.zeros_like(mu)
        self.omega     = np.zeros_like(M)
        self.P         = np.zeros_like(filters[0].P)
        self.x         = np.zeros_like(filters[0].x)
        self.N         = len(filters)

    def predict(self):
        xs, Ps = [], []
        
        self.cbar = np.dot(self.mu, self.M)
        for j in range(self.N):
            for i in range(self.N):
                self.omega[j,i] = self.M[j,i] * self.mu[j] / self.cbar[i]
                
        for i, (f, w) in enumerate(zip(self.filters, self.omega.T)):
            x = np.zeros(self.x.shape)
            for kf, wj in zip(self.filters, w):
                x += kf.x * wj
            xs.append(x)

            P = np.zeros(self.P.shape)
            for kf, wj in zip(self.filters, w):
                y = kf.x - x
                P += wj * (np.outer(y, y) + kf.P)
            Ps.append(P)
            
        for i, f in enumerate(self.filters):
            f.x = xs[i]
            f.P = Ps[i]
            f.predict()
            
    def update(self, z):

        s = 0.0
        
        for i, f in enumerate(self.filters):
            f.update(z, self.HJacobian, self.Hx)    
            self.mu[i] = f.likelihood * self.cbar[i]
            s += self.mu[i]
            
        for i in range(self.N):
            self.mu[i] /= s
        
        self.mu /= np.sum(self.mu)
        
        self.x.fill(0.)
        self.P.fill(0.)
        
        for f, w in zip(self.filters, self.mu):
            self.x += f.x * w

        for f, w in zip(self.filters, self.mu):
            y = f.x - self.x
            self.P += w * (np.outer(y, y) + f.P)

###############################################################################
def imm_ab_test(a, b, measurements):

    xa = []
    xb = []
    
    ua = [] 
    ub = []

    dP = []
    dx = []
    du = []
    
    nx = len(a.x)
    nb = len(a.mu)

    for z in measurements:
        # Run A
        a.predict()
        a.update(np.array([z]))
        xa.append(a.x.copy().reshape((nx,)))
        ua.append(a.mu.copy().reshape((nb,)))
        pa = a.P.copy().reshape((nx,nx))

        # Run B
        b.predict()
        b.update(np.array([z]))
        xb.append(b.x.copy().reshape((nx,)))
        ub.append(b.mu.copy().reshape((nb,)))
        pb = b.P.copy().reshape((nx,nx))

        dP.append(2. * norm(pa - pb) / norm(pa + pb))
        dx.append(2. * norm(xa[-1] - xb[-1]) / norm(xa[-1] + xb[-1]))
        du.append(2. * norm(ua[-1] - ub[-1]) / norm(ua[-1] + ub[-1]))

    xa  = np.array(xa)
    xb  = np.array(xb)

    dP = np.array(dP)
    dx = np.array(dx)
    du = np.array(du)

    return xa,xb, dP,dx,du

###############################################################################
# Common IMM things
def _ca(x, dt, **fx_args):
    x = x.copy()
    x[0] += (x[1] + 0.5 * x[2] * dt) * dt
    x[1] += x[2] * dt
    return x    

def _jca(x, dt, **fx_args):
    F = np.array([
        [1., dt, dt*dt*0.5],
        [0., 1.,    dt    ],
        [0., 0.,    1.    ],
        ])
    return F


def _cv(x, dt, **fx_args):
    x = x.copy()
    x[0] += x[1] * dt
    return x    

def _jcv(x, dt, **fx_args):
    F = np.array([
        [1., dt, 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        ])
    return F

def _hx_imm(x, **hx_args):
    return np.array([x[0]])

def _jhx_imm(x, **hx_args):
    H = np.array([[1., 0., 0.]])
    return H

#IMM estimator data
mu = np.array([0.5, 0.5], dtype=NP_DTYPE)
M  = np.array([[0.95, 0.05],
               [0.05, 0.95]], dtype=NP_DTYPE)

###############################################################################
_dt = 0.1
STD = 1.

N = 100
time  = []
clean = []
noisy = []

#CV
t = 0.
tgt = np.array([-0., 0., 0.])
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CV
tgt[2] = 0.1
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)
#CV
tgt[2] = 0.
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)
#CA
tgt[2] = -0.2
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CV
tgt[2] = 0.
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CA
tgt[2] = 0.1
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CV
tgt[2] = 0.
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

plt.plot(time, noisy, time, clean)
plt.show()

###############################################################################
# filterpy based estimator
acv = KFA(3, 1)
acv.F = _jcv(acv.x, _dt)
acv.P *= .1
acv.Q *= .000001
acv.R *= STD*STD
acv.x[0] = 0.
acv.x[1] = 0.
acv.x[2] = 0.

aca   = KFA(3, 1)
aca.F = _jca(aca.x, _dt)
aca.P *= .1
aca.Q *= .000001
aca.R *= STD*STD
aca.x[0] = 0.
aca.x[1] = 0.
aca.x[2] = 0.

a = IMMEstimatorA([aca, acv], mu, M, _jhx_imm, _hx_imm)

###############################################################################
# yaflpy based estimator
bcv = KFB(3, 1, _dt, _cv, _jcv, _hx_imm, _jhx_imm)
bcv.Dp *= .1
bcv.Dq *= .000001
bcv.Dr = STD*STD
bcv.x[0] = 0.
bcv.x[1] = 0.
bcv.x[2] = 0.

bca = KFB(3, 1, _dt, _ca, _jca, _hx_imm, _jhx_imm)
bca.Dp *= .1
bca.Dq *= .000001
bca.Dr = STD*STD
bca.x[0] = 0.
bca.x[1] = 0.
bca.x[2] = 0.

b = IMMEstimatorB([bca, bcv], mu, M, _dt)

###############################################################################
# Data processing
xa,xb, dP,dx,du = imm_ab_test(a,b, noisy)

plt.plot(time, dP)
plt.show()

plt.plot(time, dx)
plt.show()

plt.plot(time, du)
plt.show()
