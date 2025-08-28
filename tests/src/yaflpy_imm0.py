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

from filterpy.kalman import ExtendedKalmanFilter as KF

class IMMEstimator0(object):
    def __init__(self, filters, mu, M):
        self.filters = filters
        self.mu      = mu
        self.M       = M
        self.cbar    = np.zeros_like(mu)
        self.omega   = np.zeros_like(M)
        self.P       = np.zeros_like(filters[0].P)
        self.x       = np.zeros_like(filters[0].x)
        self.N       = len(filters)

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
            
    def update(self, z, HJacobian, Hx):

        s = 0.0
        
        for i, f in enumerate(self.filters):
            f.update(z, HJacobian, Hx)    
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

def _hx(x, **hx_args):
    return np.array([x[0]])

def _jhx(x, **hx_args):
    H = np.array([[1., 0., 0.]])
    return H


_dt = 0.1
STD = 10.

N = 1000
time  = []
clean = []
noisy = []

#CA
t = 0.
tgt = np.array([0., 0., 0.1])
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
tgt[2] = -0.1
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

plt.plot(time, noisy, time, clean)
plt.show()

kf = KF(3, 1)
kf.F = _jcv(kf.x, _dt)
kf.P *= .0001
kf.Q *= .0000001
kf.R = STD*STD
kf.x[0] = 0.
kf.x[1] = 0.
kf.x[2] = 0.

out = []
for i,z in enumerate(clean):
     kf.predict()
     kf.update(z, _jhx, _hx)
     out.append(kf.x[0])    

plt.plot(time, out, time, clean)
plt.show()

cv = KF(3, 1)
cv.F = _jcv(cv.x, _dt)
cv.P *= .0001
cv.Q *= .0000001
cv.R *= STD*STD
cv.x[0] = 0.
cv.x[1] = 0.
cv.x[2] = 0.

ca   = KF(3, 1)
ca.F = _jca(ca.x, _dt)
ca.P *= .0001
ca.Q *= .0000001
ca.R *= STD*STD
ca.x[0] = 0.
ca.x[1] = 0.
ca.x[2] = 0.

mu = np.array([0.5, 0.5])
M  = np.array([[0.95, 0.05],
               [0.05, 0.95]])

imm = IMMEstimator0([ca, cv], mu, M)

out = []
mu  = []
cvl = []
cal = []
cvx = []
cax = []


for i,z in enumerate(clean):
     imm.predict()
     imm.update(np.array([z]), _jhx, _hx)
     out.append(imm.x[0].copy())
     mu.append(imm.mu.copy())
     cvl.append(cv.log_likelihood)
     cal.append(ca.log_likelihood)
     cvx.append(cv.x.copy())
     cax.append(ca.x.copy())

plt.plot(time, out, time, clean)
plt.show()

plt.plot(time, mu)
plt.show()

plt.plot(time, cvl)
plt.show()

plt.plot(time, cal)
plt.show()

print('Done!')
