# -*- coding: utf-8 -*-
"""
    Copyright 2020 anonimous <shkolnick-kun@gmail.com> and contributors.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    See the License for the specific language governing permissions
    and limitations under the License.
"""
#==============================================================================
#                                YAKF C API
#==============================================================================
# distutils: language = c
from libc cimport stdint

cdef extern from "yakf_config.h":
    ctypedef double         yakfFloat
    ctypedef stdint.int32_t yakfInt

cdef extern from "yakf.c":

    ctypedef _yakfBaseSt yakfBaseSt

    ctypedef void (* yakfFuncP)(yakfBaseSt *)
    ctypedef void (* yakfResFuncP)(yakfBaseSt *, yakfFloat *)
    ctypedef void (* yakfScalarUpdateP)(yakfBaseSt *, yakfInt)

    ctypedef struct _yakfBaseSt:        
        yakfFuncP f      #
        yakfFuncP jf     #

        yakfFuncP h      #
        yakfFuncP jh     #
        yakfResFuncP zrf #

        yakfFloat * x    #
        yakfFloat * y    #
        yakfFloat * H    #

        yakfFloat * Up   #
        yakfFloat * Dp   #

        yakfFloat * Uq   #
        yakfFloat * Dq   #

        yakfFloat * Ur   #
        yakfFloat * Dr   #

        yakfFloat * W    #
        yakfFloat * D    #

        yakfInt   Nx     #
        yakfInt   Nz     #

    cdef void yakf_base_predict(yakfBaseSt * self)
    cdef void yakf_base_update(yakfBaseSt * self, yakfFloat * z, yakfScalarUpdateP scalar_update)

    cdef void yakf_bierman_update(yakfBaseSt * self, yakfFloat * z)

    cdef void yakf_joseph_update(yakfBaseSt * self, yakfFloat * z)

    ctypedef struct yakfAdaptiveSt:
        yakfBaseSt base
        yakfFloat  chi2

    cdef void yakf_adaptive_bierman_update(yakfAdaptiveSt * self, yakfFloat * z)

    cdef void yakf_adaptive_joseph_update(yakfAdaptiveSt * self, yakfFloat * z)

    # For demonstration purposes only
    cdef void yakf_do_not_use_this_update(yakfAdaptiveSt * self, yakfFloat * z)

cdef extern from "yakf_math.c":
    cdef void yakfm_set_u(yakfInt sz, yakfFloat *res, yakfFloat *u)

#==============================================================================
#Extension API
#==============================================================================
# We need numpy for Pythonic interfaces
cimport numpy as np
import  numpy as np

# Kalman filter C-structure with Python callback info
ctypedef struct yakfPySt:
    # Kalman filter base structure
    yakfAdaptiveSt base

    # Python/Cython self
    void * py_self

cdef int _U_sz(int dim_u):
    return max(1, (dim_u * (dim_u - 1))//2)

# Bierman filter class
cdef class yakfBase:
    # Kalman filter C-self
    cdef yakfPySt c_self
    
    # Kalman filter memory views
    cdef yakfFloat [::1]    v_x
    cdef yakfFloat [::1]    v_y
    cdef yakfFloat [::1]    v_z
    cdef yakfFloat [:, ::1] v_H

    cdef yakfFloat [::1]    v_Up
    cdef yakfFloat [::1]    v_Dp

    cdef yakfFloat [::1]    v_Uq
    cdef yakfFloat [::1]    v_Dq

    cdef yakfFloat [::1]    v_Ur
    cdef yakfFloat [::1]    v_Dr

    cdef yakfFloat [:, ::1] v_W
    cdef yakfFloat [::1]    v_D
    
    # Kalman filter numpy arrays
    cdef np.ndarray  _x
    cdef np.ndarray  _y
    cdef np.ndarray  _z
    cdef np.ndarray  _H

    cdef np.ndarray  _Up
    cdef np.ndarray  _Dp

    cdef np.ndarray _Uq
    cdef np.ndarray _Dq

    cdef np.ndarray _Ur
    cdef np.ndarray _Dr

    cdef np.ndarray _W
    cdef np.ndarray _D
    
    # Callback info
    cdef yakfFloat _dt
    cdef dict      _fx_args
    cdef object    _fx
    cdef object    _jfx
    
    cdef dict      _hx_args
    cdef object    _hx
    cdef object    _jhx
    cdef object    _residual_z

    #The object will be Extensible
    cdef dict __dict__

    def __init__(self, int dim_x, int dim_z, yakfFloat dt, \
                 fx, jfx, hx, jhx, residual_z = None):

        #Init chi2
        self.c_self.base.chi2 = 10.8275662

        #Store dimensions
        self.c_self.base.base.Nx = dim_x
        self.c_self.base.base.Nz = dim_z

        #Setup callbacks
        self.c_self.py_self = <void *>self

        self._dt = dt
        self._fx_args = {}
        self._hx_args = {}

        if not callable(fx):
            raise ValueError('fx must be callable!')
        self.c_self.base.base.f = <yakfFuncP> yakf_py_fx
        self._fx = fx

        if not callable(jfx):
            raise ValueError('jfx must be callable!')
        self.c_self.base.base.jf = <yakfFuncP> yakf_py_jfx
        self._jfx = jfx

        if not callable(hx):
            raise ValueError('hx must be callable!')
        self.c_self.base.base.h = <yakfFuncP> yakf_py_hx
        self._hx = hx

        if not callable(jhx):
            raise ValueError('jhx must be callable!')
        self.c_self.base.base.jh = <yakfFuncP> yakf_py_jhx
        self._jhx = jhx

        if residual_z:   
            if not callable(residual_z):
                raise ValueError('residual_z must be callable!')
            self.c_self.base.base.zrf = <yakfResFuncP> yakf_py_zrf
            self._residual_z = residual_z
        else:
            self.c_self.base.base.zrf = <yakfResFuncP> 0
            self._residual_z = None

        # Allocate memories and setup the rest of c_self
        self._z  = np.zeros((dim_z,), dtype=np.float64)
        self.v_z = self._z

        self._x  = np.zeros((dim_x,), dtype=np.float64)
        self.v_x = self._x
        self.c_self.base.base.x = &self.v_x[0]

        self._y  = np.zeros((dim_z,), dtype=np.float64)
        self.v_y = self._y
        self.c_self.base.base.y = &self.v_y[0]

        self._H  = np.zeros((dim_z, dim_x), dtype=np.float64)
        self.v_H = self._H
        self.c_self.base.base.H = &self.v_H[0, 0]

        self._Up  = np.zeros((_U_sz(dim_x),), dtype=np.float64)
        self.v_Up = self._Up
        self.c_self.base.base.Up = &self.v_Up[0]

        self._Dp  = np.ones((dim_x,), dtype=np.float64)
        self.v_Dp = self._Dp
        self.c_self.base.base.Dp = &self.v_Dp[0]

        self._Uq  = np.zeros((_U_sz(dim_x),), dtype=np.float64)
        self.v_Uq = self._Uq
        self.c_self.base.base.Uq = &self.v_Uq[0]

        self._Dq  = np.ones((dim_x,), dtype=np.float64)
        self.v_Dq = self._Dq
        self.c_self.base.base.Dq = &self.v_Dq[0]

        self._Ur  = np.zeros((_U_sz(dim_z),), dtype=np.float64)
        self.v_Ur = self._Ur
        self.c_self.base.base.Ur = &self.v_Ur[0]

        self._Dr  = np.ones((dim_z,), dtype=np.float64)
        self.v_Dr = self._Dr
        self.c_self.base.base.Dr = &self.v_Dr[0]

        self._W  = np.zeros((dim_x, 2 * dim_x), dtype=np.float64)
        self.v_W = self._W
        self.c_self.base.base.W = &self.v_W[0,0]

        self._D  = np.ones((2 * dim_x,), dtype=np.float64)
        self.v_D = self._D
        self.c_self.base.base.D = &self.v_D[0]

    #==========================================================================
    #Decorators
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x[:] = value
    #--------------------------------------------------------------------------
    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        raise AttributeError('yakfBase does not support this!')
    #--------------------------------------------------------------------------
    @property
    def Up(self):
        return self._Up

    @Up.setter
    def Up(self, value):
        self._Up[:] = value
    #--------------------------------------------------------------------------
    @property
    def Dp(self):
        return self._Dp

    @Dp.setter
    def Dp(self, value):
        self._Dp[:] = value
    #--------------------------------------------------------------------------
    @property
    def Uq(self):
        return self._Uq

    @Uq.setter
    def Uq(self, value):
        self._Uq[:] = value
    #--------------------------------------------------------------------------
    @property
    def Dq(self):
        return self._Dq

    @Dq.setter
    def Dq(self, value):
        self._Dq[:] = value
    #--------------------------------------------------------------------------
    @property
    def Ur(self):
        return self._Ur

    @Ur.setter
    def Ur(self, value):
        self._Ur[:] = value
    #--------------------------------------------------------------------------
    @property
    def Dr(self):
        return self._Dr

    @Dr.setter
    def Dr(self, value):
        self._Dr[:] = value
   
    #==========================================================================
    def _predict(self):
        yakf_base_predict(&(self.c_self.base.base))

    def predict(self, dt=None, **fx_args):
        old_dt = self._dt

        if dt:
            if np.isnan(dt):
                raise ValueError('Invalid dt value (nan)!')
            self._dt = <yakfFloat>dt
            
        self._fx_args = fx_args
        self._predict()

        self._dt = old_dt
    
    #==========================================================================
    def _update(self):
        raise NotImplementedError('yakfBase is the base class!')
    
    def update(self, z, **hx_args):

        self._z[:] = z           
        self._hx_args = hx_args
        self._update()
        
#==============================================================================
# State transition function 
cdef void yakf_py_fx(yakfPySt * self):
    
    py_self = <yakfBase>(self.py_self)
    if not isinstance(py_self, yakfBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfBase)!')
    
    fx = py_self._fx
    if not callable(fx):
        raise ValueError('fx must be callable!')
    
    dt = py_self._dt
    if np.isnan(dt):
        raise ValueError('Invalid dt value (nan)!')

    fx_args = py_self._fx_args
    if not isinstance(fx_args, dict):
        raise ValueError('Invalid fx_args type (must be dict)!')
    
    #How about handling exceptions here???
    py_self._x[:] = fx(py_self._x, dt, **fx_args)

#==============================================================================
# State transition function Jacobian
cdef void yakf_py_jfx(yakfPySt * self):
    
    py_self = <yakfBase>(self.py_self)
    if not isinstance(py_self, yakfBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfBase)!')
       
    jfx = py_self._jfx
    if not callable(jfx):
        raise ValueError('jfx must be callable!')
    
    dt = py_self._dt
    if np.isnan(dt):
        raise ValueError('Invalid dt value (nan)!')

    fx_args = py_self._fx_args
    if not isinstance(fx_args, dict):
        raise ValueError('Invalid fx_args type (must be dict)!')
    
    #How about handling exceptions here???
    py_self._W[:, :self.base.base.Nx] = jfx(py_self._x, dt, **fx_args)

#==============================================================================
# State transition function 
cdef void yakf_py_hx(yakfPySt * self):
    
    py_self = <yakfBase>(self.py_self)
    if not isinstance(py_self, yakfBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfBase)!')
    
    hx = py_self._hx
    if not callable(hx):
        raise ValueError('hx must be callable!')

    hx_args = py_self._hx_args
    if not isinstance(hx_args, dict):
        raise ValueError('Invalid hx_args type (must be dict)!')
    
    #How about handling exceptions here???
    py_self._y[:] = hx(py_self._x, **hx_args)

#==============================================================================
# State transition function Jacobian
cdef void yakf_py_jhx(yakfPySt * self):
    
    py_self = <yakfBase>(self.py_self)
    if not isinstance(py_self, yakfBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfBase)!')
    
    jhx = py_self._jhx
    if not callable(jhx):
        raise ValueError('jhx must be callable!')

    hx_args = py_self._hx_args
    if not isinstance(hx_args, dict):
        raise ValueError('Invalid hx_args type (must be dict)!')
    
    #How about handling exceptions here???
    py_self._H[:,:] = jhx(py_self._x, **hx_args)

#==============================================================================    
# Measurement residual function
cdef void yakf_py_zrf(yakfPySt * self, yakfFloat * zp):
    
    py_self = <yakfBase>(self.py_self)
    if not isinstance(py_self, yakfBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfBase)!')
    
    zrf = py_self._residual_z
    if not callable(zrf):
        raise ValueError('jhx must be callable!')
    
    #How about handling exceptions here???
    py_self._y[:] = zrf(py_self._z, py_self._y)
    
#==============================================================================
cdef class Bierman(yakfBase):
    def _update(self):
        yakf_bierman_update(&self.c_self.base.base, &self.v_z[0])

#==============================================================================
cdef class Joseph(yakfBase):
    def _update(self):
        yakf_joseph_update(&self.c_self.base.base, &self.v_z[0])

#==============================================================================
cdef class AdaptiveBierman(yakfBase):
    def _update(self):
        yakf_adaptive_bierman_update(&self.c_self.base, &self.v_z[0])

#==============================================================================
cdef class AdaptiveJoseph(yakfBase):
    def _update(self):
        yakf_adaptive_joseph_update(&self.c_self.base, &self.v_z[0])

#==============================================================================
cdef class DoNotUseThisFilter(yakfBase):
    """
    WARNING!!!
    DO NOT USE THIS variant of Adaptive Joseph filter !!!
    It was implemented to show some flaws of the corresponding algorithm!
    """
    def _update(self):
        yakf_do_not_use_this_update(&self.c_self.base, &self.v_z[0])
