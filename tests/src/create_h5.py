#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 17:08:55 2020

@author: anon
"""
import subprocess
import numpy as np
import h5py

import matplotlib.pyplot as plt

N = 10000

clear = np.zeros((N, 2))
noisy = np.zeros((N, 2))
for i in range(1, len(clear)):
    clear[i] = clear[i-1] + np.array([1.,1.])
    noisy[i] = clear[i]   + np.random.normal(scale=20., size=2)

with h5py.File('../data/input.h5', 'w') as h5f:
    h5f.create_dataset('clear', data=clear)
    h5f.create_dataset('noisy', data=noisy)

#Run kf aplication here
subprocess.call(['../projects/ada_ukf/bin/Debug/hdf5_test'])

with h5py.File('../data/output.h5', 'r') as h5f:
    kf_out = h5f['kf_out'][:]

plt.plot(clear - kf_out)
plt.show()
