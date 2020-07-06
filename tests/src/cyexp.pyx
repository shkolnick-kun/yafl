# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 19:49:19 2020

@author: anon
"""

#from cython.view cimport array as cvarray

import  numpy as np
cimport numpy as np

cdef extern from "cyexp_c.c":
    
    ctypedef _cHelloSt cHelloSt
    
    ctypedef void (*cyWrapperP)(cHelloSt *)
    
    ctypedef struct _cHelloSt:
        cyWrapperP cy_wrapper
        void * py_func
        void * py_object
    
    cdef void do_callback(cHelloSt * arg)
    
    cpdef void hello()
    
    cdef void print_mat(int nr, int nc, double * m)
    cdef void zero_mat(int nr, int nc, double * m)
    
cdef void callback(cHelloSt * arg):
    pyfunc   = <object>(arg.py_func)
    pyobject = <object>(arg.py_object)
    print(type(pyfunc))
    print(callable(pyfunc))
    print(type(pyobject))
    print(type(pyobject) == np.ndarray)
    pyfunc(pyobject)
    
def pycallback(f, o):
    cdef cHelloSt arg
    arg.cy_wrapper = callback
    arg.py_func    = <void *>f
    arg.py_object  = <void *>o
    do_callback(&arg)
    
#Для матриц Си в последнем инексе указывается ::1
def cprint_mat(m):
    cdef double [:, ::1] m_view
    
    m_view = m
    print(m_view[0,1])
    print_mat(m.shape[0], m.shape[1], &m_view[0,0])
    
def czero_mat(m):
    cdef double [:, ::1] m_view = m
    zero_mat(m.shape[0], m.shape[1], &m_view[0,0])
    
cdef double [:, ::1] view_static

def set_view(m):
    global view_static
    view_static = m
    
def print_view(int i, int j):
    global view_static
    print('view_static:', view_static[i,j])

cdef class Fail:
    cdef double * _zp;
    cdef double [:, ::1] _v
    cdef object _o
    cdef dict   _d
    cdef np.ndarray _a
    cdef dict __dict__
    
    
    def __init__(self, int n):
        self.z = np.zeros((n*2, n), dtype=np.float64)
        self._v = self.z
        self._zp = <double *>&self._v[0,0]
        self._a = self.z
        
        self._o = self
        self._d = {'a':1, 'b':2}
        #self._c = cprint_mat
        
    def print_a(self):
        print(self._a)
        
    def print_zp(self, i):
        print(self._zp[i])

cdef void print_as_array(np.float64_t * x, int a, int b):
    print(np.asarray(<np.float64_t[:a, :b]> x))
    print(np.asarray(<np.float64_t[:a*b]> x))

def print_view_static():
    print('view_static:\n')
    print_as_array(&view_static[0, 0], view_static.shape[0], view_static.shape[1])
