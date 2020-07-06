#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 19:49:39 2020

@author: anon
"""
import numpy as np
import pyximport

pyximport.install(
    build_dir='../projects/obj', 
    pyimport=True,
    reload_support=True, 
    language_level=3,
    setup_args={
        'include_dirs': ['./'],
        }
    )

import cyexp

cyexp.hello()



a = np.array([[1,2,3],[4,5,6.0]])
cyexp.cprint_mat(a)

def hello_cb(o):
    print('Will print some object:\n', o)
    
cyexp.pycallback(hello_cb, a)

cyexp.czero_mat(a)
print(a)

cyexp.set_view(a)

cyexp.print_view(0,1)

a[0,1] = 7

cyexp.print_view(0,1)

a[0,1] = 34

cyexp.print_view(0,1)
cyexp.print_view(0,2)

f = cyexp.Fail(3)

cyexp.print_view_static()
