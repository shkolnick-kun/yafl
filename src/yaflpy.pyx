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
#                                YAFL C API
#==============================================================================
#cython: language_level=3
#distutils: language=c
from libc cimport stdint

#------------------------------------------------------------------------------
cdef extern from "yafl_config.h":
    ctypedef double         yaflFloat
    ctypedef stdint.int32_t yaflInt

#------------------------------------------------------------------------------
cdef extern from "yafl_math.c":
    ctypedef enum yaflStatusEn:
        #Warning flag masks
        YAFL_ST_MSK_REGULARIZED  = 0x01 #YAFL_ST_R
        YAFL_ST_MSK_GLITCH_SMALL = 0x02 #YAFL_ST_S
        YAFL_ST_MSK_GLITCH_LARGE = 0x04 #YAFL_ST_L
        YAFL_ST_MSK_ANOMALY      = 0x08 #YAFL_ST_A
        #Everthing is OK
        YAFL_ST_OK           = 0x00
        #Wagnings
        YAFL_ST_R            = 0x01
        YAFL_ST_S            = 0x02
        YAFL_ST_SR           = 0x03
        YAFL_ST_L            = 0x04
        YAFL_ST_LR           = 0x05
        YAFL_ST_SL           = 0x06
        YAFL_ST_SLR          = 0x07
        YAFL_ST_A            = 0x08
        YAFL_ST_AR           = 0x09
        YAFL_ST_SA           = 0x0a
        YAFL_ST_SAR          = 0x0b
        YAFL_ST_LA           = 0x0c
        YAFL_ST_LAR          = 0x0d
        YAFL_ST_SLA          = 0x0e
        YAFL_ST_SLAR         = 0x0f
        # Error threshold value (greater values are errors)
        YAFL_ST_ERR_THR      = 0x010
        # Invalid argument numer
        YAFL_ST_INV_ARG_1    = 0x100
        YAFL_ST_INV_ARG_2    = 0x110
        YAFL_ST_INV_ARG_3    = 0x120
        YAFL_ST_INV_ARG_4    = 0x130
        YAFL_ST_INV_ARG_5    = 0x140
        YAFL_ST_INV_ARG_6    = 0x150
        YAFL_ST_INV_ARG_7    = 0x160
        YAFL_ST_INV_ARG_8    = 0x170
        YAFL_ST_INV_ARG_9    = 0x180
        YAFL_ST_INV_ARG_10   = 0x190
        YAFL_ST_INV_ARG_11   = 0x1a0

    #--------------------------------------------------------------------------
    #cdef yaflStatusEn yafl_math_set_u(yaflInt sz, yaflFloat *res, yaflFloat *u)

#------------------------------------------------------------------------------
cdef extern from "yafl.c":
    ctypedef _yaflKalmanBaseSt yaflKalmanBaseSt

    ctypedef yaflStatusEn (* yaflKalmanFuncP)(yaflKalmanBaseSt *, \
                                              yaflFloat *, yaflFloat *)

    ctypedef yaflStatusEn (* yaflKalmanResFuncP)(yaflKalmanBaseSt *, \
                                                 yaflFloat *, yaflFloat *, \
                                                     yaflFloat *)

    ctypedef yaflStatusEn (* yaflKalmanScalarUpdateP)(yaflKalmanBaseSt *, \
                                                      yaflInt)

    ctypedef yaflFloat (* yaflKalmanRobFuncP)(yaflKalmanBaseSt *, yaflFloat)

    ctypedef struct _yaflKalmanBaseSt:
        yaflKalmanFuncP f      #
        yaflKalmanFuncP h      #
        yaflKalmanResFuncP zrf #

        yaflFloat * x    #
        yaflFloat * y    #

        yaflFloat * Up   #
        yaflFloat * Dp   #

        yaflFloat * Uq   #
        yaflFloat * Dq   #

        yaflFloat * Ur   #
        yaflFloat * Dr   #

        yaflInt   Nx     #
        yaflInt   Nz     #

    #==========================================================================
    #                     UD-factorized EKF definitions
    #==========================================================================
    ctypedef struct yaflEKFBaseSt:
        yaflKalmanBaseSt base

        yaflKalmanFuncP jf #
        yaflKalmanFuncP jh #

        yaflFloat * H      #
        yaflFloat * W      #
        yaflFloat * D      #

    cdef yaflStatusEn yafl_ekf_base_predict(yaflKalmanBaseSt * self)

    cdef yaflStatusEn \
        yafl_ekf_base_update(yaflKalmanBaseSt * self, yaflFloat * z, \
                             yaflKalmanScalarUpdateP scalar_update)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_bierman_update_scalar(yaflKalmanBaseSt * self, yaflInt i)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_joseph_update_scalar(yaflKalmanBaseSt * self, yaflInt i)


    #==========================================================================
    ctypedef struct yaflEKFAdaptiveSt:
        yaflEKFBaseSt base
        yaflFloat  chi2

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_adaptive_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                yaflInt i)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_adaptive_joseph_update_scalar(yaflKalmanBaseSt * self, \
                                               yaflInt i)

    #--------------------------------------------------------------------------
    # For demonstration purposes only
    #cdef yaflStatusEn \
    #    yafl_ekf_do_not_use_this_update(yaflEKFAdaptiveSt * self, \
    #                                    yaflFloat * z)

    #==========================================================================
    ctypedef struct yaflEKFRobustSt:
        yaflEKFBaseSt   base
        yaflKalmanRobFuncP g
        yaflKalmanRobFuncP gdot

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                              yaflInt i)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_robust_joseph_update_scalar(yaflKalmanBaseSt * self, \
                                             yaflInt i)

    #==========================================================================
    ctypedef struct yaflEKFAdaptiveRobustSt:
        yaflEKFRobustSt base
        yaflFloat  chi2

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_adaptive_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                       yaflInt i)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ekf_adaptive_robust_joseph_update_scalar(yaflKalmanBaseSt * self, \
                                                      yaflInt i)

    #==========================================================================
    #                     UD-factorized UKF definitions
    #==========================================================================
    ctypedef _yaflUKFBaseSt yaflUKFBaseSt

    ctypedef yaflStatusEn (* yaflUKFSigmaAddP)(yaflUKFBaseSt *, yaflFloat *, \
                                               yaflFloat *, yaflFloat)

    ctypedef struct yaflUKFSigmaSt:
        yaflInt     np
        yaflUKFSigmaAddP addf

    #--------------------------------------------------------------------------
    ctypedef yaflStatusEn (* yaflUKFSigmaGenWeigthsP)(yaflUKFBaseSt *)
    ctypedef yaflStatusEn (* yaflUKFSigmaGenSigmasP)(yaflUKFBaseSt *)

    ctypedef struct yaflUKFSigmaMethodsSt:
        yaflUKFSigmaGenWeigthsP   wf
        yaflUKFSigmaGenSigmasP  spgf

    #--------------------------------------------------------------------------

    ctypedef struct _yaflUKFBaseSt:
        yaflKalmanBaseSt base

        yaflUKFSigmaSt              * sp_info
        const yaflUKFSigmaMethodsSt * sp_meth

        yaflKalmanFuncP    xmf
        yaflKalmanResFuncP xrf

        yaflKalmanFuncP    zmf
        yaflFloat * zp

        yaflFloat * Sx
        yaflFloat * Pzx

        yaflFloat * sigmas_x
        yaflFloat * sigmas_z
        yaflFloat * wm
        yaflFloat * wc

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_ukf_post_init(yaflUKFBaseSt * self)  #static inline

    cdef yaflStatusEn yafl_ukf_gen_sigmas(yaflUKFBaseSt * self) #static inline

    cdef yaflStatusEn yafl_ukf_base_predict(yaflUKFBaseSt * self)

    cdef yaflStatusEn \
        yafl_ukf_base_update(yaflUKFBaseSt * self, yaflFloat * z, \
                             yaflKalmanScalarUpdateP scalar_update)

    #==========================================================================
    cdef yaflStatusEn \
        yafl_ukf_bierman_update_scalar(yaflKalmanBaseSt * self, yaflInt i)

    #==========================================================================
    ctypedef struct yaflUKFAdaptivedSt:
        yaflUKFBaseSt base
        yaflFloat  chi2

    cdef yaflStatusEn \
        yafl_ukf_adaptive_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                yaflInt i)

    #==========================================================================
    ctypedef struct yaflUKFRobustSt:
        yaflUKFBaseSt   base
        yaflKalmanRobFuncP g
        yaflKalmanRobFuncP gdot

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ukf_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                              yaflInt i)

    #==========================================================================
    ctypedef struct yaflUKFAdaptiveRobustSt:
        yaflUKFRobustSt   base
        yaflFloat chi2

    #--------------------------------------------------------------------------
    cdef yaflStatusEn \
        yafl_ukf_adaptive_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                              yaflInt i)

    #==========================================================================
    ctypedef struct yaflUKFSt:
        yaflUKFBaseSt base

        yaflFloat * Us
        yaflFloat * Ds

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_ukf_update(yaflUKFBaseSt * self, yaflFloat * z)

    #==========================================================================
    ctypedef struct yaflUKFFullAdapiveSt:
        yaflUKFSt base
        yaflFloat chi2

    yaflStatusEn yafl_ukf_adaptive_update(yaflUKFBaseSt * self, yaflFloat * z)

    #==========================================================================
    #                  Van der Merwe sigma point generator
    #==========================================================================
    ctypedef struct yaflUKFMerweSt:
        yaflUKFSigmaSt base
        yaflFloat alpha
        yaflFloat beta
        yaflFloat kappa

    cdef const yaflUKFSigmaMethodsSt yafl_ukf_merwe_spm

#==============================================================================
#Extension API
#==============================================================================
# We need numpy for Pythonic interfaces
cimport numpy as np
import  numpy as np#WTF??

import scipy.stats as st

# We need traceback to print pythonic callback exceptions
import traceback as tb

#==============================================================================
#Status masks
ST_MSK_REGULARIZED  = YAFL_ST_MSK_REGULARIZED
ST_MSK_GLITCH_SMALL = YAFL_ST_MSK_GLITCH_SMALL
ST_MSK_GLITCH_LARGE = YAFL_ST_MSK_GLITCH_LARGE
ST_MSK_ANOMALY      = YAFL_ST_MSK_ANOMALY
ST_MSK_ERROR        = 0xFF0
ST_MSK_WARNING      = 0xF

#Everthing is OK
ST_OK = YAFL_ST_OK

#Error threshold
ST_ERR_THR = YAFL_ST_ERR_THR
#Errors:
# Invalid argument numer
ST_INV_ARG_1  = YAFL_ST_INV_ARG_1
ST_INV_ARG_2  = YAFL_ST_INV_ARG_2
ST_INV_ARG_3  = YAFL_ST_INV_ARG_3
ST_INV_ARG_4  = YAFL_ST_INV_ARG_4
ST_INV_ARG_5  = YAFL_ST_INV_ARG_5
ST_INV_ARG_6  = YAFL_ST_INV_ARG_6
ST_INV_ARG_7  = YAFL_ST_INV_ARG_7
ST_INV_ARG_8  = YAFL_ST_INV_ARG_8
ST_INV_ARG_9  = YAFL_ST_INV_ARG_9
ST_INV_ARG_10 = YAFL_ST_INV_ARG_10
ST_INV_ARG_11 = YAFL_ST_INV_ARG_11

#==============================================================================
#                          UD-factorized EKF API
#==============================================================================
#------------------------------------------------------------------------------
#                       Kalman filter basic union
#------------------------------------------------------------------------------
ctypedef union yaflPyKalmanBaseUn:
    yaflKalmanBaseSt        base
    yaflEKFBaseSt           ekf
    yaflEKFAdaptiveSt       ekf_adaptive
    yaflEKFRobustSt         ekf_robust
    yaflEKFAdaptiveRobustSt ekf_ada_rob

    yaflUKFBaseSt           ukf
    yaflUKFAdaptivedSt      ukf_adaptive
    yaflUKFRobustSt         ukf_robust
    yaflUKFAdaptiveRobustSt ukf_ada_rob
    yaflUKFSt               ukf_full
    yaflUKFFullAdapiveSt    ukf_full_adaptive

#------------------------------------------------------------------------------
# Kalman filter C-structure with Python callback
#------------------------------------------------------------------------------
ctypedef struct yaflPyKalmanBaseSt:
    # Kalman filter base union
    yaflPyKalmanBaseUn base

    # Python/Cython self
    void * py_self

#------------------------------------------------------------------------------
cdef int _U_sz(int dim_u):
    return max(1, (dim_u * (dim_u - 1))//2)

#------------------------------------------------------------------------------
#                             Basic Filter class
#------------------------------------------------------------------------------
cdef class yaflKalmanBase:
    # Kalman filter C-self
    cdef yaflPyKalmanBaseSt c_self

    # Kalman filter memory views
    cdef yaflFloat [::1]    v_x
    cdef yaflFloat [::1]    v_y
    cdef yaflFloat [::1]    v_z

    cdef yaflFloat [::1]    v_Up
    cdef yaflFloat [::1]    v_Dp

    cdef yaflFloat [::1]    v_Uq
    cdef yaflFloat [::1]    v_Dq

    cdef yaflFloat [::1]    v_Ur
    cdef yaflFloat [::1]    v_Dr

    # Kalman filter numpy arrays
    cdef np.ndarray  _x
    cdef np.ndarray  _y
    cdef np.ndarray  _z

    cdef np.ndarray  _Up
    cdef np.ndarray  _Dp

    cdef np.ndarray _Uq
    cdef np.ndarray _Dq

    cdef np.ndarray _Ur
    cdef np.ndarray _Dr

    # Callback info
    cdef yaflFloat _dt
    cdef dict      _fx_args
    cdef object    _fx

    cdef dict      _hx_args
    cdef object    _hx

    cdef object    _residual_z

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, hx, residual_z = None):

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
        self.c_self.base.base.f = <yaflKalmanFuncP>yafl_py_kalman_fx
        self._fx = fx

        if not callable(hx):
            raise ValueError('hx must be callable!')
        self.c_self.base.base.h = <yaflKalmanFuncP>yafl_py_kalman_hx
        self._hx = hx

        if residual_z:
            if not callable(residual_z):
                raise ValueError('residual_z must be callable!')
            self.c_self.base.base.zrf = <yaflKalmanResFuncP>yafl_py_kalman_zrf
            self._residual_z = residual_z
        else:
            self.c_self.base.base.zrf = <yaflKalmanResFuncP>0
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
        raise AttributeError('yaflKalmanBase does not support this!')
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
        raise NotImplementedError('yaflKalmanBase is the base class!')

    def predict(self, dt=None, **fx_args):
        old_dt = self._dt

        if dt:
            if np.isnan(dt):
                raise ValueError('Invalid dt value (nan)!')
            self._dt = <yaflFloat>dt

        self._fx_args = fx_args

        res = self._predict()
        if res > YAFL_ST_ERR_THR:
            raise ValueError('Bad return value on yaflKalmanBase.predict!')

        self._dt = old_dt
        return res

    #==========================================================================
    def _update(self):
        raise NotImplementedError('yaflKalmanBase is the base class!')

    def update(self, z, **hx_args):

        self._z[:] = z
        self._hx_args = hx_args
        res = self._update()
        if res > YAFL_ST_ERR_THR:
            raise ValueError('Bad return value on yaflKalmanBase.update!')
        return res

#------------------------------------------------------------------------------
cdef yaflStatusEn yafl_py_kalman_fx(yaflPyKalmanBaseSt * self, \
                                    yaflFloat * new_x, yaflFloat * old_x):
    try:
        if not isinstance(<object>(self.py_self), yaflKalmanBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflKalmanBase)!')

        py_self = <yaflKalmanBase>(self.py_self)

        fx = py_self._fx
        if not callable(fx):
            raise ValueError('fx must be callable!')

        dt = py_self._dt
        if np.isnan(dt):
            raise ValueError('Invalid dt value (nan)!')

        fx_args = py_self._fx_args
        if not isinstance(fx_args, dict):
            raise ValueError('Invalid fx_args type (must be dict)!')

        nx = self.base.base.Nx
        if nx <= 0:
            raise ValueError('nx must be > 0!')

        _new_x = np.asarray(<yaflFloat[:nx]> new_x) #np.float64_t
        _old_x = np.asarray(<yaflFloat[:nx]> old_x) #np.float64_t

        _new_x[:] = fx(_old_x, dt, **fx_args)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#------------------------------------------------------------------------------
cdef yaflStatusEn yafl_py_kalman_hx(yaflPyKalmanBaseSt * self, \
                         yaflFloat * z, yaflFloat * x):
    try:
        if not isinstance(<object>(self.py_self), yaflKalmanBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflKalmanBase)!')

        py_self = <yaflKalmanBase>(self.py_self)

        hx = py_self._hx
        if not callable(hx):
            raise ValueError('hx must be callable!')

        hx_args = py_self._hx_args
        if not isinstance(hx_args, dict):
            raise ValueError('Invalid hx_args type (must be dict)!')

        nx = self.base.base.Nx
        if nx <= 0:
            raise ValueError('nx must be > 0!')

        nz = self.base.base.Nz
        if nx <= 0:
            raise ValueError('nz must be > 0!')

        _x = np.asarray(<yaflFloat[:nx]> x) #np.float64_t
        _z = np.asarray(<yaflFloat[:nz]> z) #np.float64_t

        _z[:] = hx(_x, **hx_args)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#------------------------------------------------------------------------------
cdef yaflStatusEn yafl_py_kalman_zrf(yaflPyKalmanBaseSt * self, yaflFloat * res, \
                             yaflFloat * sigma, yaflFloat * pivot):
    try:
        if not isinstance(<object>(self.py_self), yaflKalmanBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflKalmanBase)!')

        py_self = <yaflKalmanBase>(self.py_self)

        residual_z = py_self._residual_z
        if not callable(residual_z):
            raise ValueError('residual_z must be callable!')

        nz = self.base.base.Nz
        if nz <= 0:
            raise ValueError('nx must be > 0!')

        _res   = np.asarray(<yaflFloat[:nz]> res)   #np.float64_t
        _sigma = np.asarray(<yaflFloat[:nz]> sigma) #np.float64_t
        _pivot = np.asarray(<yaflFloat[:nz]> pivot) #np.float64_t

        _res[:] = residual_z(_sigma, _pivot)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#==============================================================================
#                             Basic Filter class
#==============================================================================
cdef class yaflExtendedBase(yaflKalmanBase):

    # Kalman filter memory views
    cdef yaflFloat [:, ::1] v_H
    cdef yaflFloat [:, ::1] v_W
    cdef yaflFloat [::1]    v_D

    # Kalman filter numpy arrays
    cdef np.ndarray _H
    cdef np.ndarray _W
    cdef np.ndarray _D

    # Callback info
    cdef object    _jfx
    cdef object    _jhx

    #The object will be Extensible
    cdef dict __dict__

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, jfx, hx, jhx, residual_z = None):

        super().__init__(dim_x, dim_z, dt, fx, hx, residual_z)

        if not callable(jfx):
            raise ValueError('jfx must be callable!')
        self.c_self.base.ekf.jf = <yaflKalmanFuncP>yafl_py_ekf_jfx
        self._jfx = jfx

        if not callable(jhx):
            raise ValueError('jhx must be callable!')
        self.c_self.base.ekf.jh = <yaflKalmanFuncP>yafl_py_ekf_jhx
        self._jhx = jhx


        # Allocate memories and setup the rest of c_self
        self._H  = np.zeros((dim_z, dim_x), dtype=np.float64)
        self.v_H = self._H
        self.c_self.base.ekf.H = &self.v_H[0, 0]

        self._W  = np.zeros((dim_x, 2 * dim_x), dtype=np.float64)
        self.v_W = self._W
        self.c_self.base.ekf.W = &self.v_W[0,0]

        self._D  = np.ones((2 * dim_x,), dtype=np.float64)
        self.v_D = self._D
        self.c_self.base.ekf.D = &self.v_D[0]

    #==========================================================================
    def _predict(self):
        return yafl_ekf_base_predict(&(self.c_self.base.base))

    #==========================================================================
    def _update(self):
        raise NotImplementedError('yaflExtendedBase is the base class!')

#------------------------------------------------------------------------------
# State transition function Jacobian
cdef yaflStatusEn yafl_py_ekf_jfx(yaflPyKalmanBaseSt * self, yaflFloat * w, \
                                  yaflFloat * x):
    try:
        if not isinstance(<object>(self.py_self), yaflExtendedBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')

        py_self = <yaflExtendedBase>(self.py_self)

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

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#------------------------------------------------------------------------------
# State transition function Jacobian
cdef yaflStatusEn yafl_py_ekf_jhx(yaflPyKalmanBaseSt * self, yaflFloat * h, \
                                  yaflFloat * x):
    try:
        if not isinstance(<object>(self.py_self), yaflExtendedBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')

        py_self = <yaflExtendedBase>(self.py_self)

        jhx = py_self._jhx
        if not callable(jhx):
            raise ValueError('jhx must be callable!')

        hx_args = py_self._hx_args
        if not isinstance(hx_args, dict):
            raise ValueError('Invalid hx_args type (must be dict)!')

        #How about handling exceptions here???
        py_self._H[:,:] = jhx(py_self._x, **hx_args)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#==============================================================================
cdef class Bierman(yaflExtendedBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_bierman_update_scalar)

#------------------------------------------------------------------------------
cdef class Joseph(yaflExtendedBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_joseph_update_scalar)

#==============================================================================
#                        Adaptive filter basic class
#==============================================================================
cdef class yaflAdaptiveBase(yaflExtendedBase):
    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, jfx, hx, jhx, **kwargs):

        super().__init__(dim_x, dim_z, dt, fx, jfx, hx, jhx, **kwargs)

        #Init chi2 with scipy.stats.chi2.ppf(0.999, 1)
        self.c_self.base.ekf_adaptive.chi2 = 10.8275662
    #==========================================================================
    #Decorators
    @property
    def chi2(self):
        return self.c_self.base.ekf_adaptive.chi2

    @chi2.setter
    def chi2(self, value):
        self.c_self.base.ekf_adaptive.chi2 = <yaflFloat>value
#==============================================================================
cdef class AdaptiveBierman(yaflAdaptiveBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_adaptive_bierman_update_scalar)

#------------------------------------------------------------------------------
cdef class AdaptiveJoseph(yaflAdaptiveBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_adaptive_joseph_update_scalar)

#------------------------------------------------------------------------------
# cdef class DoNotUseThisFilter(yaflAdaptiveBase):
#     """
#     WARNING!!!
#     DO NOT USE THIS variant of Adaptive Joseph filter !!!
#     It was implemented to show some flaws of the corresponding algorithm!
#     """
#     def _update(self):
#         return yafl_ekf_do_not_use_this_update(&self.c_self.base.adaptive, \
#                                                &self.v_z[0])



#==============================================================================
#                        Robust filter basic class
#==============================================================================
cdef class yaflRobustBase(yaflExtendedBase):

    cdef object _gz
    cdef object _gdotz

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, jfx, hx, jhx, gz=None, gdotz=None, **kwargs):

        super().__init__(dim_x, dim_z, dt, fx, jfx, hx, jhx, **kwargs)

        if gz:
            if not callable(gz):
                raise ValueError('gz must be callable!')

            if not gdotz:
                raise ValueError('gdotz must be passed!')

            if not callable(gdotz):
                raise ValueError('gdotz must be callable!')

            self._gz = gz
            self._gdotz = gdotz

            self.c_self.base.ekf_robust.g    = <yaflKalmanRobFuncP>yafl_py_ekf_rob_gz
            self.c_self.base.ekf_robust.gdot = <yaflKalmanRobFuncP>yafl_py_ekf_rob_gdotz

        else:
            self.c_self.base.ekf_robust.g    = <yaflKalmanRobFuncP>0
            self.c_self.base.ekf_robust.gdot = <yaflKalmanRobFuncP>0

#------------------------------------------------------------------------------
# Influence limiting function
cdef yaflFloat yafl_py_ekf_rob_gz(yaflPyKalmanBaseSt * self, yaflFloat nu):
    try:
        if not isinstance(<object>(self.py_self), yaflRobustBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflRobustBase)!')

        py_self = <yaflRobustBase>(self.py_self)

        gz = py_self._gz
        if not callable(gz):
            raise ValueError('gz must be callable!')

        hx_args = py_self._hx_args
        if not isinstance(hx_args, dict):
            raise ValueError('Invalid hx_args type (must be dict)!')

        #How about handling exceptions here???
        ret = gz(nu, **hx_args)
        if type(ret) != float:
            raise ValueError('gz must return float!')

        return <yaflFloat>ret

    except Exception as e:
        print(tb.format_exc())
        return <yaflFloat>0.0
#------------------------------------------------------------------------------
# Influence limiting function derivative
cdef yaflFloat yafl_py_ekf_rob_gdotz(yaflPyKalmanBaseSt * self, yaflFloat nu):
    try:
        if not isinstance(<object>(self.py_self), yaflRobustBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflRobustBase)!')

        py_self = <yaflRobustBase>(self.py_self)

        gdotz = py_self._gdotz
        if not callable(gdotz):
            raise ValueError('gdotz must be callable!')

        hx_args = py_self._hx_args
        if not isinstance(hx_args, dict):
            raise ValueError('Invalid hx_args type (must be dict)!')

        #How about handling exceptions here???
        ret = gdotz(nu, **hx_args)
        if type(ret) != float:
            raise ValueError('gdotz must return float!')

        return <yaflFloat>ret

    except Exception as e:
        print(tb.format_exc())
        return <yaflFloat>0.0
#==============================================================================
cdef class RobustBierman(yaflRobustBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_robust_bierman_update_scalar)

#------------------------------------------------------------------------------
cdef class RobustJoseph(yaflRobustBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_robust_joseph_update_scalar)

#==============================================================================
#                   Adaptive robust filter basic class
#==============================================================================
cdef class yaflAdaptiveRobustBase(yaflRobustBase):

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, jfx, hx, jhx, **kwargs):

        super().__init__(dim_x, dim_z, dt, fx, jfx, hx, jhx, **kwargs)

        #Init chi2 with scipy.stats.chi2.ppf(0.997, 1)
        self.c_self.base.ekf_ada_rob.chi2 = 8.807468393511947

    #==========================================================================
    #Decorators
    @property
    def chi2(self):
        return self.c_self.base.ekf_ada_rob.chi2

    @chi2.setter
    def chi2(self, value):
        self.c_self.base.ekf_ada_rob.chi2 = <yaflFloat>value

#==============================================================================
cdef class AdaptiveRobustBierman(yaflAdaptiveRobustBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_adaptive_robust_bierman_update_scalar)

#------------------------------------------------------------------------------
cdef class AdaptiveRobustJoseph(yaflAdaptiveRobustBase):
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_adaptive_robust_joseph_update_scalar)

#==============================================================================
#                          UD-factorized UKF API
#==============================================================================
#                   Sigma points generator basic definitions
#------------------------------------------------------------------------------
ctypedef union yaflPySigmaBaseUn:
    yaflUKFSigmaSt base
    yaflUKFMerweSt merwe

#------------------------------------------------------------------------------
ctypedef struct yaflPySigmaSt:
    # Sigma point generator base structure
    yaflPySigmaBaseUn base

    # Python/Cython self
    void * py_self

#------------------------------------------------------------------------------
cdef class yaflSigmaBase:
    # Sigma point generator C-self
    cdef yaflPySigmaSt c_self

    # Callback info
    cdef object _addf

    #The object will be Extensible
    cdef dict __dict__

    def __init__(self, yaflInt dim_x, addf=None):
        #Setup callbacks
        self.c_self.py_self = <void *>self

        if addf:
            if not callable(addf):
                raise ValueError('addf must be callable!')
            self.c_self.base.base.addf = <yaflUKFSigmaAddP>yafl_py_sigma_addf
            self._addf = addf
        else:
            self.c_self.base.base.addf = <yaflUKFSigmaAddP>0

        pnum = self.get_np(dim_x)

        self.c_self.base.base.np = pnum

    cdef yaflInt get_np(self, int dim_x):
        raise NotImplementedError('yaflSigmaBase is the base class!')

    cdef const yaflUKFSigmaMethodsSt * get_spm(self):
        raise NotImplementedError('yaflSigmaBase is the base class!')

    @property
    def pnum(self):
        return self.c_self.base.base.np

    @pnum.setter
    def pnum(self, value):
        raise AttributeError('yaflSigmaBase does not support this!')

#------------------------------------------------------------------------------
#                         UD-factorized UKF definitions
#------------------------------------------------------------------------------
cdef class yaflUnscentedBase(yaflKalmanBase):

    # Kalman filter memory views
    cdef yaflFloat [::1]    v_zp

    cdef yaflFloat [::1]    v_Sx
    cdef yaflFloat [:, ::1] v_Pzx

    cdef yaflFloat [:, ::1] v_sigmas_x
    cdef yaflFloat [:, ::1] v_sigmas_z
    cdef yaflFloat [::1]    v_wm
    cdef yaflFloat [::1]    v_wc

    # Kalman filter numpy arrays
    cdef np.ndarray  _zp

    cdef np.ndarray  _Sx
    cdef np.ndarray  _Pzx

    cdef np.ndarray  _sigmas_x
    cdef np.ndarray  _sigmas_z
    cdef np.ndarray  _wm
    cdef np.ndarray  _wc

    # Callback info
    cdef object    _points

    cdef object    _mean_x
    cdef object    _resudual_x

    cdef object    _mean_z

    #The object will be Extensible
    cdef dict __dict__

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, hx, fx, points,\
                 x_mean_fn=None, z_mean_fn=None, \
                 residual_x=None, residual_z=None):

        super().__init__(dim_x, dim_z, dt, fx, hx, residual_z)


        if x_mean_fn:
            if not callable(x_mean_fn):
                raise ValueError('x_mean_fn must be callable!')
            self.c_self.base.ukf.xmf = <yaflKalmanFuncP>yafl_py_ukf_xmf
            self._mean_x = x_mean_fn
        else:
            self.c_self.base.ukf.xmf = <yaflKalmanFuncP>0
            self._mean_x = None

        if residual_x:
            if not callable(residual_x):
                raise ValueError('residual_x must be callable!')
            self.c_self.base.ukf.xrf = <yaflKalmanResFuncP>yafl_py_ukf_xrf
            self._residual_x = residual_x
        else:
            self.c_self.base.ukf.xrf = <yaflKalmanResFuncP>0
            self._residual_x = None

        if z_mean_fn:
            if not callable(z_mean_fn):
                raise ValueError('z_mean_fn must be callable!')
            self.c_self.base.ukf.zmf = <yaflKalmanFuncP>yafl_py_ukf_zmf
            self._mean_z = z_mean_fn
        else:
            self.c_self.base.ukf.zmf = <yaflKalmanFuncP>0
            self._mean_z = None

        #Setup sigma points generator
        if not isinstance(points, yaflSigmaBase):
            raise ValueError('Invalid points type (must be subclass of yaflSigmaBase)!')

        _points = <yaflSigmaBase>points

        self._points = _points
        self.c_self.base.ukf.sp_info = &_points.c_self.base.base
        self.c_self.base.ukf.sp_meth = _points.get_spm()

        # Allocate memories and setup the rest of c_self
        # Sigma points and weights
        pnum = _points.pnum

        self._wm  = np.zeros((pnum,), dtype=np.float64)
        self.v_wm = self._wm
        self.c_self.base.ukf.wm = &self.v_wm[0]

        self._wc  = np.zeros((pnum,), dtype=np.float64)
        self.v_wc = self._wc
        self.c_self.base.ukf.wc = &self.v_wc[0]

        self._sigmas_x  = np.zeros((pnum, dim_x), dtype=np.float64)
        self.v_sigmas_x = self._sigmas_x
        self.c_self.base.ukf.sigmas_x = &self.v_sigmas_x[0, 0]

        self._sigmas_z  = np.zeros((pnum, dim_z), dtype=np.float64)
        self.v_sigmas_z = self._sigmas_z
        self.c_self.base.ukf.sigmas_z = &self.v_sigmas_z[0, 0]

        #Rest of rhe UKF
        self._zp  = np.zeros((dim_z,), dtype=np.float64)
        self.v_zp = self._zp
        self.c_self.base.ukf.zp = &self.v_zp[0]

        self._Pzx  = np.zeros((dim_z, dim_x), dtype=np.float64)
        self.v_Pzx = self._Pzx
        self.c_self.base.ukf.Pzx = &self.v_Pzx[0, 0]

        self._Sx  = np.zeros((dim_x,), dtype=np.float64)
        self.v_Sx = self._Sx
        self.c_self.base.ukf.Sx = &self.v_Sx[0]

        #Call C-post init
        yafl_ukf_post_init(&self.c_self.base.ukf)

    #==========================================================================
    #Decorators
    @property
    def Pzx(self):
        return self._Pzx

    #--------------------------------------------------------------------------
    @property
    def zp(self):
        return self._zp

    #--------------------------------------------------------------------------
    @property
    def sigmas_x(self):
        return self._sigmas_x

    @sigmas_x.setter
    def sigmas_x(self, value):
        raise AttributeError('yaflUnscentedBase does not support this!')

    #--------------------------------------------------------------------------
    @property
    def sigmas_z(self):
        return self._sigmas_z

    @sigmas_z.setter
    def sigmas_z(self, value):
        raise AttributeError('yaflUnscentedBase does not support this!')
    #--------------------------------------------------------------------------
    @property
    def wc(self):
        return self._wc

    @wc.setter
    def wc(self, value):
        raise AttributeError('yaflUnscentedBase does not support this!')
    #--------------------------------------------------------------------------
    @property
    def wm(self):
        return self._wm

    @wm.setter
    def wm(self, value):
        raise AttributeError('yaflUnscentedBase does not support this!')

    #==========================================================================
    def _predict(self):
        return yafl_ukf_base_predict(&self.c_self.base.ukf)

    #==========================================================================
    def _update(self):
        raise NotImplementedError('yaflUnscentedBase is the base class!')

#==============================================================================
cdef yaflStatusEn yafl_py_sigma_addf(yaflPyKalmanBaseSt * self, yaflFloat * delta, \
                             yaflFloat * pivot, yaflFloat mult):
    try:
        if not isinstance(<object>(self.py_self), yaflUnscentedBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflUnscentedBase)!')

        py_self = <yaflUnscentedBase>(self.py_self)

        _addf = py_self._points._addf
        if not callable(_addf):
            raise ValueError('_addf must be callable!')

        nx = self.base.base.Nx
        if nx <= 0:
            raise ValueError('nx must be > 0!')

        _delta = np.asarray(<yaflFloat[:nx]> delta) #np.float64_t
        _pivot = np.asarray(<yaflFloat[:nx]> pivot) #np.float64_t

        _delta[:] = _addf(_delta, _pivot, mult)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#------------------------------------------------------------------------------
cdef yaflStatusEn yafl_py_ukf_xmf(yaflPyKalmanBaseSt * self, \
                          yaflFloat * res, yaflFloat * sigmas):
    try:
        if not isinstance(<object>(self.py_self), yaflUnscentedBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflUnscentedBase)!')

        py_self = <yaflUnscentedBase>(self.py_self)

        mean_x = py_self._mean_x
        if not callable(mean_x):
            raise ValueError('mean_x must be callable!')

        nx = self.base.base.Nx
        if nx <= 0:
            raise ValueError('nx must be > 0!')

        _points = py_self._points
        if not isinstance(py_self, yaflSigmaBase):
            raise ValueError('Invalid _points type (must be subclass of yaflSigmaBase)!')

        pnum = _points.base.base.np
        if pnum <= 0:
            raise ValueError('pnum must be > 0!')

        _sigmas = np.asarray(<yaflFloat[:pnum, :nx]> sigmas) #np.float64_t
        _res    = np.asarray(<yaflFloat[:nx]> res)           #np.float64_t

        _res[:] = mean_x(_sigmas, _points._wm)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#------------------------------------------------------------------------------
cdef yaflStatusEn yafl_py_ukf_xrf(yaflPyKalmanBaseSt * self, yaflFloat * res, \
                             yaflFloat * sigma, yaflFloat * pivot):
    try:
        if not isinstance(<object>(self.py_self), yaflUnscentedBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflUnscentedBase)!')

        py_self = <yaflUnscentedBase>(self.py_self)

        residual_x = py_self._residual_x
        if not callable(residual_x):
            raise ValueError('residual_x must be callable!')

        nx = self.base.base.Nx
        if nx <= 0:
            raise ValueError('nx must be > 0!')

        _res   = np.asarray(<yaflFloat[:nx]> res)   #np.float64_t
        _sigma = np.asarray(<yaflFloat[:nx]> sigma) #np.float64_t
        _pivot = np.asarray(<yaflFloat[:nx]> pivot) #np.float64_t

        _res[:] = residual_x(_sigma, _pivot)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#------------------------------------------------------------------------------
cdef yaflStatusEn yafl_py_ukf_zmf(yaflPyKalmanBaseSt * self, \
                          yaflFloat * res, yaflFloat * sigmas):
    try:
        if not isinstance(<object>(self.py_self), yaflUnscentedBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflUnscentedBase)!')

        py_self = <yaflUnscentedBase>(self.py_self)

        mean_z = py_self._mean_z
        if not callable(mean_z):
            raise ValueError('mean_z must be callable!')

        nz = self.base.base.Nz
        if nz <= 0:
            raise ValueError('nz must be > 0!')

        _points = py_self._points
        if not isinstance(py_self, yaflSigmaBase):
            raise ValueError('Invalid _points type (must be subclass of yaflSigmaBase)!')

        pnum = _points.base.base.np
        if pnum <= 0:
            raise ValueError('pnum must be > 0!')

        _sigmas = np.asarray(<yaflFloat[:pnum, :nz]> sigmas) #np.float64_t
        _res    = np.asarray(<yaflFloat[:nz]> res)           #np.float64_t

        _res[:] = mean_z(_sigmas, _points._wm)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#==============================================================================
cdef class UnscentedBierman(yaflUnscentedBase):
    """
    UD-factorized UKF implementation
    """
    def _update(self):
        return yafl_ukf_base_update(&self.c_self.base.ukf, &self.v_z[0], \
                                yafl_ukf_bierman_update_scalar)

#==============================================================================
cdef class UnscentedAdaptiveBierman(yaflUnscentedBase):
    def __init__(self, int dim_x, int dim_z, yaflFloat dt, hx, fx, points, \
                 **kwargs):

        super().__init__(dim_x, dim_z, dt, hx, fx, points, **kwargs)

        #Init chi2 with scipy.stats.chi2.ppf(0.99999, 1)
        self.c_self.base.ukf_adaptive.chi2 = 19.511420964666268

    #==========================================================================
    #Decorators
    @property
    def chi2(self):
        return self.c_self.base.ukf_adaptive.chi2

    @chi2.setter
    def chi2(self, value):
        self.c_self.base.ukf_adaptive.chi2 = <yaflFloat>value

    """
    UD-factorized UKF implementation
    """
    def _update(self):
        return yafl_ukf_base_update(&self.c_self.base.ukf, &self.v_z[0], \
                                yafl_ukf_adaptive_bierman_update_scalar)
#==============================================================================
#                        Robust filter basic class
#==============================================================================
cdef class yaflRobustUKFBase(yaflUnscentedBase):

    cdef object _gz
    cdef object _gdotz

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 hx, fx, points, gz=None, gdotz=None, **kwargs):

        super().__init__(dim_x, dim_z, dt, hx, fx, points, **kwargs)

        if gz:
            if not callable(gz):
                raise ValueError('gz must be callable!')

            if not gdotz:
                raise ValueError('gdotz must be passed!')

            if not callable(gdotz):
                raise ValueError('gdotz must be callable!')

            self._gz = gz
            self._gdotz = gdotz

            self.c_self.base.ukf_robust.g    = <yaflKalmanRobFuncP>yafl_py_ukf_rob_gz
            self.c_self.base.ukf_robust.gdot = <yaflKalmanRobFuncP>yafl_py_ukf_rob_gdotz

        else:
            self.c_self.base.ukf_robust.g    = <yaflKalmanRobFuncP>0
            self.c_self.base.ukf_robust.gdot = <yaflKalmanRobFuncP>0

#------------------------------------------------------------------------------
# Influence limiting function
cdef yaflFloat yafl_py_ukf_rob_gz(yaflPyKalmanBaseSt * self, yaflFloat nu):
    try:
        if not isinstance(<object>(self.py_self), yaflRobustUKFBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflRobustBase)!')

        py_self = <yaflRobustUKFBase>(self.py_self)

        gz = py_self._gz
        if not callable(gz):
            raise ValueError('gz must be callable!')

        hx_args = py_self._hx_args
        if not isinstance(hx_args, dict):
            raise ValueError('Invalid hx_args type (must be dict)!')

        #How about handling exceptions here???
        ret = gz(nu, **hx_args)
        if type(ret) != float:
            raise ValueError('gz must return float!')

        return <yaflFloat>ret

    except Exception as e:
        print(tb.format_exc())
        return <yaflFloat>0.0
#------------------------------------------------------------------------------
# Influence limiting function derivative
cdef yaflFloat yafl_py_ukf_rob_gdotz(yaflPyKalmanBaseSt * self, yaflFloat nu):
    try:
        if not isinstance(<object>(self.py_self), yaflRobustUKFBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflRobustBase)!')

        py_self = <yaflRobustUKFBase>(self.py_self)

        gdotz = py_self._gdotz
        if not callable(gdotz):
            raise ValueError('gdotz must be callable!')

        hx_args = py_self._hx_args
        if not isinstance(hx_args, dict):
            raise ValueError('Invalid hx_args type (must be dict)!')

        #How about handling exceptions here???
        ret = gdotz(nu, **hx_args)
        if type(ret) != float:
            raise ValueError('gdotz must return float!')

        return <yaflFloat>ret

    except Exception as e:
        print(tb.format_exc())
        return <yaflFloat>0.0
#==============================================================================
cdef class UnscentedRobustBierman(yaflRobustUKFBase):
    def _update(self):
        return yafl_ukf_base_update(&self.c_self.base.ukf, &self.v_z[0], \
                                    yafl_ukf_robust_bierman_update_scalar)

#==============================================================================
#                         Adaptive robust UKF base
#==============================================================================
cdef class yaflAdaptiveRobustUKFBase(yaflRobustUKFBase):

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 hx, fx, points, **kwargs):

        super().__init__(dim_x, dim_z, dt, hx, fx, points, **kwargs)

        #Init chi2 with scipy.stats.chi2.ppf(0.997, 1)
        self.c_self.base.ukf_ada_rob.chi2 = 8.807468393511947

    #==========================================================================
    #Decorators
    @property
    def chi2(self):
        return self.c_self.base.ukf_ada_rob.chi2

    @chi2.setter
    def chi2(self, value):
        self.c_self.base.ukf_ada_rob.chi2 = <yaflFloat>value

#==============================================================================
cdef class UnscentedAdaptiveRobustBierman(yaflAdaptiveRobustUKFBase):
    def _update(self):
        return yafl_ukf_base_update(&self.c_self.base.ukf, &self.v_z[0], \
                                    yafl_ukf_adaptive_robust_bierman_update_scalar)

#==============================================================================
#           Full UKF, not sequential square root version of UKF
#==============================================================================
cdef class Unscented(yaflUnscentedBase):

    cdef yaflFloat [::1]    v_Us
    cdef yaflFloat [::1]    v_Ds

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, hx, fx, points, \
                 **kwargs):

        super().__init__(dim_x, dim_z, dt, hx, fx, points, **kwargs)

        self._Us  = np.zeros((_U_sz(dim_z),), dtype=np.float64)
        self.v_Us = self._Us
        self.c_self.base.ukf_full.Us = &self.v_Us[0]

        self._Ds  = np.ones((dim_z,), dtype=np.float64)
        self.v_Ds = self._Ds
        self.c_self.base.ukf_full.Ds = &self.v_Ds[0]

    #==========================================================================
    #Decorators
    @property
    def Us(self):
        return self._Us

    @Us.setter
    def Us(self, value):
        raise AttributeError('Unscented does not support this!')

    #--------------------------------------------------------------------------
    @property
    def Ds(self):
        return self._Ds

    @Ds.setter
    def Ds(self, value):
        raise AttributeError('Unscented does not support this!')

    #==========================================================================
    """
    UD-factorized UKF implementation
    """
    def _update(self):
        return yafl_ukf_update(&self.c_self.base.ukf, &self.v_z[0])

#==============================================================================
#       Full adaptive UKF, not sequential square root version of UKF
#==============================================================================
cdef class UnscentedAdaptive(Unscented):

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, hx, fx, points, \
                 aplha = 0.000001, **kwargs):

        super().__init__(dim_x, dim_z, dt, hx, fx, points, **kwargs)

        self.c_self.base.ukf_full_adaptive.chi2 = \
            st.chi2.ppf(1.0 - aplha, dim_z)

    #==========================================================================
    #Decorators
    @property
    def chi2(self):
        return self.c_self.base.ukf_full_adaptive.chi2

    @chi2.setter
    def chi2(self, value):
        self.c_self.base.ukf_full_adaptive.chi2 = <yaflFloat>value

    #==========================================================================
    """
    UD-factorized UKF implementation
    """
    def _update(self):
        return yafl_ukf_adaptive_update(&self.c_self.base.ukf, &self.v_z[0])

#==============================================================================
cdef class MerweSigmaPoints(yaflSigmaBase):
    """
    Van der Merwe sigma point generator implementation
    """
    def __init__(self, yaflInt dim_x, yaflFloat alpha, yaflFloat beta, \
                 yaflFloat kappa=0.0, **kwargs):
        super().__init__(dim_x, **kwargs)
        self.c_self.base.merwe.alpha = alpha
        self.c_self.base.merwe.beta  = beta
        self.c_self.base.merwe.kappa = kappa

    cdef yaflInt get_np(self, int dim_x):
        return (2 * dim_x + 1)

    cdef const yaflUKFSigmaMethodsSt * get_spm(self):
        return &yafl_ukf_merwe_spm

    #==========================================================================
    #Decorators
    @property
    def alpha(self):
        return self.c_self.base.merwe.alpha

    @alpha.setter
    def alpha(self, value):
        raise AttributeError('MerweSigmaPoints does not support this!')

    @property
    def beta(self):
        return self.c_self.base.merwe.beta

    @beta.setter
    def beta(self, value):
        raise AttributeError('MerweSigmaPoints does not support this!')

    @property
    def kappa(self):
        return self.c_self.base.merwe.kappa

    @kappa.setter
    def kappa(self, value):
        raise AttributeError('MerweSigmaPoints does not support this!')
