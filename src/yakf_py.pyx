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

#------------------------------------------------------------------------------
cdef extern from "yafl_config.h":
    ctypedef double         yaflFloat
    ctypedef stdint.int32_t yaflInt

#------------------------------------------------------------------------------
cdef extern from "yafl_math.c":
    cdef enum yaflStatusEn:
        YAFL_ST_OK           = 0x00
        YAFL_ST_REGULARIZED  = 0x01
        YAFL_ST_GLITCH_SMALL = 0x02
        YAFL_ST_GSR          = 0x03
        YAFL_ST_GLITCH_LARGE = 0x04
        YAFL_ST_GLR          = 0x05
        YAFL_ST_GLITCH_ALL   = 0x06
        YAFL_ST_GAR          = 0x07
        # Error threshold value (greater values are errors)
        YAFL_ST_ERR_THR      = 0x100
        # Invalid argument numer
        YAFL_ST_INV_AGR_1    = 0x100
        YAFL_ST_INV_AGR_2    = 0x110
        YAFL_ST_INV_AGR_3    = 0x120
        YAFL_ST_INV_AGR_4    = 0x130
        YAFL_ST_INV_AGR_5    = 0x140
        YAFL_ST_INV_AGR_6    = 0x150
        YAFL_ST_INV_AGR_7    = 0x160

    #--------------------------------------------------------------------------
    #cdef yaflStatusEn yafl_math_set_u(yaflInt sz, yaflFloat *res, yaflFloat *u)

#------------------------------------------------------------------------------
cdef extern from "yafl.c":
    #==========================================================================
    #                     UD-factorized EKF definitions
    #==========================================================================
    ctypedef _yaflEKFBaseSt yaflEKFBaseSt

    ctypedef void (* yaflEKFFuncP)(yaflEKFBaseSt *)
    ctypedef void (* yaflEKFResFuncP)(yaflEKFBaseSt *, yaflFloat *)
    ctypedef yaflStatusEn (* yaflEKFScalarUpdateP)(yaflEKFBaseSt *, yaflInt)

    ctypedef struct _yaflEKFBaseSt:
        yaflEKFFuncP f      #
        yaflEKFFuncP jf     #

        yaflEKFFuncP h      #
        yaflEKFFuncP jh     #
        yaflEKFResFuncP zrf #

        yaflFloat * x    #
        yaflFloat * y    #
        yaflFloat * H    #

        yaflFloat * Up   #
        yaflFloat * Dp   #

        yaflFloat * Uq   #
        yaflFloat * Dq   #

        yaflFloat * Ur   #
        yaflFloat * Dr   #

        yaflFloat * W    #
        yaflFloat * D    #

        yaflInt   Nx     #
        yaflInt   Nz     #

    cdef void yafl_ekf_base_predict(yaflEKFBaseSt * self)
    cdef void yafl_ekf_base_update(yaflEKFBaseSt * self, yaflFloat * z, \
                               yaflEKFScalarUpdateP scalar_update)

    cdef void yafl_ekf_bierman_update(yaflEKFBaseSt * self, yaflFloat * z)

    cdef void yafl_ekf_joseph_update(yaflEKFBaseSt * self, yaflFloat * z)

    #--------------------------------------------------------------------------
    ctypedef struct yaflEKFAdaptiveSt:
        yaflEKFBaseSt base
        yaflFloat  chi2

    cdef void yafl_ekf_adaptive_bierman_update(yaflEKFAdaptiveSt * self, \
                                               yaflFloat * z)

    cdef void yafl_ekf_adaptive_joseph_update(yaflEKFAdaptiveSt * self, \
                                              yaflFloat * z)

    # For demonstration purposes only
    cdef void yafl_ekf_do_not_use_this_update(yaflEKFAdaptiveSt * self, \
                                              yaflFloat * z)

    #--------------------------------------------------------------------------
    ctypedef yaflFloat (* yaflEKFRobFuncP)(yaflEKFBaseSt *, yaflFloat)

    ctypedef struct yaflEKFRobustSt:
        yaflEKFBaseSt   base
        yaflEKFRobFuncP g
        yaflEKFRobFuncP gdot

    cdef void yafl_ekf_robust_bierman_update(yaflEKFRobustSt * self, \
                                             yaflFloat * z)

    cdef void yafl_ekf_robust_joseph_update(yaflEKFRobustSt * self, \
                                            yaflFloat * z)
    #--------------------------------------------------------------------------
    ctypedef struct yaflEKFAdaptiveRobustSt:
        yaflEKFRobustSt base
        yaflFloat  chi2

    cdef void yafl_ekf_adaptive_robust_bierman_update(yaflEKFAdaptiveRobustSt * self,\
                                                  yaflFloat * z)

    cdef void yafl_ekf_adaptive_robust_joseph_update(yaflEKFAdaptiveRobustSt * self, \
                                                 yaflFloat * z)

    #==========================================================================
    #                     UD-factorized UKF definitions
    #==========================================================================
    ctypedef _yaflUKFBaseSt yaflUKFBaseSt

    ctypedef void (* yaflUKFSigmaAddP)(yaflUKFBaseSt *, yaflFloat *, \
                                       yaflFloat *, yaflFloat)

    ctypedef struct yaflUKFSigmaSt:
        yaflInt     np
        yaflUKFSigmaAddP addf

    #--------------------------------------------------------------------------
    ctypedef yaflStatusEn (* yaflUKFScalarUpdateP)(yaflUKFBaseSt *, yaflInt)
    ctypedef void (* yaflUKFSigmaGenWeigthsP)(yaflUKFBaseSt *)
    ctypedef void (* yaflUKFSigmaGenSigmasP)(yaflUKFBaseSt *)

    ctypedef struct yaflUKFSigmaMethodsSt:
        yaflUKFSigmaGenWeigthsP   wf
        yaflUKFSigmaGenSigmasP  spgf

    #--------------------------------------------------------------------------
    ctypedef void (* yaflUKFFuncP)(yaflUKFBaseSt *, yaflFloat *, yaflFloat *)

    ctypedef void (* yaflUKFResFuncP)(yaflUKFBaseSt *, yaflFloat *, \
                                      yaflFloat *, yaflFloat *)

    ctypedef struct _yaflUKFBaseSt:
        yaflUKFSigmaSt              * sp_info
        const yaflUKFSigmaMethodsSt * sp_meth

        yaflUKFFuncP      f
        yaflUKFFuncP    xmf
        yaflUKFResFuncP xrf

        yaflUKFFuncP      h
        yaflUKFFuncP    zmf
        yaflUKFResFuncP zrf

        yaflFloat * x
        yaflFloat * zp
        yaflFloat * y

        yaflFloat * Up
        yaflFloat * Dp

        yaflFloat * Us
        yaflFloat * Ds

        yaflFloat * Pzx

        yaflFloat * Uq
        yaflFloat * Dq

        yaflFloat * Ur
        yaflFloat * Dr

        yaflFloat * sigmas_x
        yaflFloat * sigmas_z
        yaflFloat * wm
        yaflFloat * wc

        yaflFloat * Sx

        yaflInt   Nx
        yaflInt   Nz

    #--------------------------------------------------------------------------
    cdef void yafl_ukf_post_init(yaflUKFBaseSt * self)  #static inline

    cdef void yafl_ukf_gen_sigmas(yaflUKFBaseSt * self) #static inline

    cdef void yafl_ukf_predict(yaflUKFBaseSt * self)

    cdef void yafl_ukf_update(yaflUKFBaseSt * self, yaflFloat * z)

    #--------------------------------------------------------------------------
    cdef void yafl_ukf_base_update(yaflUKFBaseSt * self, yaflFloat * z, \
                                   yaflUKFScalarUpdateP scalar_update)

    #--------------------------------------------------------------------------
    cdef void yafl_ukf_bierman_update(yaflUKFBaseSt * self, yaflFloat * z)

    #==========================================================================
    ctypedef struct yaflUKFAdaptivedSt:
        yaflUKFBaseSt base
        yaflFloat  chi2

    cdef void yafl_ukf_adaptive_bierman_update(yaflUKFAdaptivedSt * self,\
                                               yaflFloat * z)
    #--------------------------------------------------------------------------
    #                  Van der Merwe sigma point generator
    #--------------------------------------------------------------------------
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
import  numpy as np

#==============================================================================
#                          UD-factorized EKF API
#==============================================================================
#------------------------------------------------------------------------------
#                       Kalman filter basic union
#------------------------------------------------------------------------------
ctypedef union yaflPyEkfBaseUn:
    yaflEKFBaseSt           base
    yaflEKFAdaptiveSt       adaptive
    yaflEKFRobustSt         robust
    yaflEKFAdaptiveRobustSt ada_rob

#------------------------------------------------------------------------------
# Kalman filter C-structure with Python callback
#------------------------------------------------------------------------------
ctypedef struct yaflPyEkfBaseSt:
    # Kalman filter base union
    yaflPyEkfBaseUn base

    # Python/Cython self
    void * py_self

#------------------------------------------------------------------------------
cdef int _U_sz(int dim_u):
    return max(1, (dim_u * (dim_u - 1))//2)

#------------------------------------------------------------------------------
#                             Basic Filter class
#------------------------------------------------------------------------------
cdef class yaflExtendedBase:
    # Kalman filter C-self
    cdef yaflPyEkfBaseSt c_self
    
    # Kalman filter memory views
    cdef yaflFloat [::1]    v_x
    cdef yaflFloat [::1]    v_y
    cdef yaflFloat [::1]    v_z
    cdef yaflFloat [:, ::1] v_H

    cdef yaflFloat [::1]    v_Up
    cdef yaflFloat [::1]    v_Dp

    cdef yaflFloat [::1]    v_Uq
    cdef yaflFloat [::1]    v_Dq

    cdef yaflFloat [::1]    v_Ur
    cdef yaflFloat [::1]    v_Dr

    cdef yaflFloat [:, ::1] v_W
    cdef yaflFloat [::1]    v_D
    
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
    cdef yaflFloat _dt
    cdef dict      _fx_args
    cdef object    _fx
    cdef object    _jfx
    
    cdef dict      _hx_args
    cdef object    _hx
    cdef object    _jhx
    cdef object    _residual_z

    #The object will be Extensible
    cdef dict __dict__

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, jfx, hx, jhx, residual_z = None):

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
        self.c_self.base.base.f = <yaflEKFFuncP> yafl_py_ekf_fx
        self._fx = fx

        if not callable(jfx):
            raise ValueError('jfx must be callable!')
        self.c_self.base.base.jf = <yaflEKFFuncP> yafl_py_ekf_jfx
        self._jfx = jfx

        if not callable(hx):
            raise ValueError('hx must be callable!')
        self.c_self.base.base.h = <yaflEKFFuncP> yafl_py_ekf_hx
        self._hx = hx

        if not callable(jhx):
            raise ValueError('jhx must be callable!')
        self.c_self.base.base.jh = <yaflEKFFuncP> yafl_py_ekf_jhx
        self._jhx = jhx

        if residual_z:   
            if not callable(residual_z):
                raise ValueError('residual_z must be callable!')
            self.c_self.base.base.zrf = <yaflEKFResFuncP> yafl_py_ekf_zrf
            self._residual_z = residual_z
        else:
            self.c_self.base.base.zrf = <yaflEKFResFuncP> 0
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
        raise AttributeError('yaflExtendedBase does not support this!')
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
        yafl_ekf_base_predict(&(self.c_self.base.base))

    def predict(self, dt=None, **fx_args):
        old_dt = self._dt

        if dt:
            if np.isnan(dt):
                raise ValueError('Invalid dt value (nan)!')
            self._dt = <yaflFloat>dt
            
        self._fx_args = fx_args
        self._predict()

        self._dt = old_dt
    
    #==========================================================================
    def _update(self):
        raise NotImplementedError('yaflExtendedBase is the base class!')
    
    def update(self, z, **hx_args):

        self._z[:] = z           
        self._hx_args = hx_args
        self._update()
        
#==============================================================================
#                             Basic C-callbacks
#==============================================================================
# State transition function 
cdef void yafl_py_ekf_fx(yaflPyEkfBaseSt * self):
    
    py_self = <yaflExtendedBase>(self.py_self)
    if not isinstance(py_self, yaflExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')
    
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

#------------------------------------------------------------------------------
# State transition function Jacobian
cdef void yafl_py_ekf_jfx(yaflPyEkfBaseSt * self):
    
    py_self = <yaflExtendedBase>(self.py_self)
    if not isinstance(py_self, yaflExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')
       
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

#------------------------------------------------------------------------------
# State transition function 
cdef void yafl_py_ekf_hx(yaflPyEkfBaseSt * self):
    
    py_self = <yaflExtendedBase>(self.py_self)
    if not isinstance(py_self, yaflExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')
    
    hx = py_self._hx
    if not callable(hx):
        raise ValueError('hx must be callable!')

    hx_args = py_self._hx_args
    if not isinstance(hx_args, dict):
        raise ValueError('Invalid hx_args type (must be dict)!')
    
    #How about handling exceptions here???
    py_self._y[:] = hx(py_self._x, **hx_args)

#------------------------------------------------------------------------------
# State transition function Jacobian
cdef void yafl_py_ekf_jhx(yaflPyEkfBaseSt * self):
    
    py_self = <yaflExtendedBase>(self.py_self)
    if not isinstance(py_self, yaflExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')
    
    jhx = py_self._jhx
    if not callable(jhx):
        raise ValueError('jhx must be callable!')

    hx_args = py_self._hx_args
    if not isinstance(hx_args, dict):
        raise ValueError('Invalid hx_args type (must be dict)!')
    
    #How about handling exceptions here???
    py_self._H[:,:] = jhx(py_self._x, **hx_args)

#------------------------------------------------------------------------------
# Measurement residual function
cdef void yafl_py_ekf_zrf(yaflPyEkfBaseSt * self, yaflFloat * zp):
    
    py_self = <yaflExtendedBase>(self.py_self)
    if not isinstance(py_self, yaflExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')
    
    zrf = py_self._residual_z
    if not callable(zrf):
        raise ValueError('jhx must be callable!')
    
    #How about handling exceptions here???
    py_self._y[:] = zrf(py_self._z, py_self._y)
    
#==============================================================================
cdef class Bierman(yaflExtendedBase):
    def _update(self):
        yafl_ekf_bierman_update(&self.c_self.base.base, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class Joseph(yaflExtendedBase):
    def _update(self):
        yafl_ekf_joseph_update(&self.c_self.base.base, &self.v_z[0])

#==============================================================================
#                        Adaptive filter basic class
#==============================================================================
cdef class yakfAdaptiveBase(yaflExtendedBase):
    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, jfx, hx, jhx, **kwargs):
        
        super().__init__(dim_x, dim_z, dt, fx, jfx, hx, jhx, **kwargs)
        
        #Init chi2 with scipy.stats.chi2.ppf(0.999, 1) 
        self.c_self.base.adaptive.chi2 = 10.8275662
    #==========================================================================
    #Decorators
    @property
    def chi2(self):
        return self.c_self.base.adaptive.chi2

    @chi2.setter
    def chi2(self, value):
        self.c_self.base.adaptive.chi2 = <yaflFloat>value
#==============================================================================
cdef class AdaptiveBierman(yakfAdaptiveBase):
    def _update(self):
        yafl_ekf_adaptive_bierman_update(&self.c_self.base.adaptive, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class AdaptiveJoseph(yakfAdaptiveBase):
    def _update(self):
        yafl_ekf_adaptive_joseph_update(&self.c_self.base.adaptive, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class DoNotUseThisFilter(yakfAdaptiveBase):
    """
    WARNING!!!
    DO NOT USE THIS variant of Adaptive Joseph filter !!!
    It was implemented to show some flaws of the corresponding algorithm!
    """
    def _update(self):
        yafl_ekf_do_not_use_this_update(&self.c_self.base.adaptive, &self.v_z[0])

#==============================================================================
#                        Robust filter basic class
#==============================================================================
cdef class yakfRobustBase(yaflExtendedBase):
    
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
            
            self.c_self.base.robust.g    = <yaflEKFRobFuncP>yafl_py_ekf_rob_gz
            self.c_self.base.robust.gdot = <yaflEKFRobFuncP>yafl_py_ekf_rob_gdotz
            
        else:
            self.c_self.base.robust.g    = <yaflEKFRobFuncP>0
            self.c_self.base.robust.gdot = <yaflEKFRobFuncP>0

#------------------------------------------------------------------------------
# State transition function 
cdef yaflFloat yafl_py_ekf_rob_gz(yaflPyEkfBaseSt * self, yaflFloat nu):

    py_self = <yakfRobustBase>(self.py_self)
    if not isinstance(py_self, yaflExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')

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

#------------------------------------------------------------------------------
# State transition function Jacobian
cdef yaflFloat yafl_py_ekf_rob_gdotz(yaflPyEkfBaseSt * self, yaflFloat nu):

    py_self = <yakfRobustBase>(self.py_self)
    if not isinstance(py_self, yaflExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yaflExtendedBase)!')

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

#==============================================================================
cdef class RobustBierman(yakfRobustBase):
    def _update(self):
        yafl_ekf_robust_bierman_update(&self.c_self.base.robust, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class RobustJoseph(yakfRobustBase):
    def _update(self):
        yafl_ekf_robust_joseph_update(&self.c_self.base.robust, &self.v_z[0])

#==============================================================================
#                        Robust filter basic class
#==============================================================================
cdef class yakfAdaptiveRobustBase(yakfRobustBase):

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, jfx, hx, jhx, **kwargs):
        
        super().__init__(dim_x, dim_z, dt, fx, jfx, hx, jhx, **kwargs)
        
        #Init chi2 with scipy.stats.chi2.ppf(0.997, 1)
        self.c_self.base.ada_rob.chi2 = 8.807468393511947
        
    #==========================================================================
    #Decorators
    @property
    def chi2(self):
        return self.c_self.base.ada_rob.chi2

    @chi2.setter
    def chi2(self, value):
        self.c_self.base.ada_rob.chi2 = <yaflFloat>value

#==============================================================================
cdef class AdaptiveRobustBierman(yakfAdaptiveRobustBase):
    def _update(self):
        yafl_ekf_adaptive_robust_bierman_update(&self.c_self.base.ada_rob, \
                                            &self.v_z[0])

#------------------------------------------------------------------------------
cdef class AdaptiveRobustJoseph(yakfAdaptiveRobustBase):
    def _update(self):
        yafl_ekf_adaptive_robust_joseph_update(&self.c_self.base.ada_rob, \
                                           &self.v_z[0])

#==============================================================================
#                          UD-factorized UKF API
#==============================================================================
#------------------------------------------------------------------------------
#                   Sigma points generator basic definitions
#------------------------------------------------------------------------------
ctypedef union yakfPySigmaBaseUn:
    yaflUKFSigmaSt base
    yaflUKFMerweSt merwe

#------------------------------------------------------------------------------
ctypedef struct yakfPySigmaSt:
    # Sigma point generator base structure
    yakfPySigmaBaseUn base

    # Python/Cython self
    void * py_self

#------------------------------------------------------------------------------
cdef class yakfSigmaBase:
    # Sigma point generator C-self
    cdef yakfPySigmaSt c_self

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
        raise NotImplementedError('yakfSigmaBase is the base class!')

    cdef const yaflUKFSigmaMethodsSt * get_spm(self):
        raise NotImplementedError('yakfSigmaBase is the base class!')

    @property
    def pnum(self):
        return self.c_self.base.base.np
    
    @pnum.setter
    def pnum(self, value):
        raise AttributeError('yakfSigmaBase does not support this!')
#------------------------------------------------------------------------------
#                         UD-factorized UKF definitions
#------------------------------------------------------------------------------
ctypedef union yakfPyUnscentedBaseUn:
    yaflUKFBaseSt base
    yaflUKFAdaptivedSt adaptive

#------------------------------------------------------------------------------
ctypedef struct yakfPyUnscentedSt:
    # Basic filter structure
    yakfPyUnscentedBaseUn base

    # Python/Cython self
    void * py_self

#------------------------------------------------------------------------------
cdef class yakfUnscentedBase:
    # Kalman filter C-self
    cdef yakfPyUnscentedSt c_self

    # Kalman filter memory views
    cdef yaflFloat [::1]    v_z

    cdef yaflFloat [::1]    v_x
    cdef yaflFloat [::1]    v_zp
    cdef yaflFloat [::1]    v_y

    cdef yaflFloat [::1]    v_Up
    cdef yaflFloat [::1]    v_Dp

    cdef yaflFloat [::1]    v_Us
    cdef yaflFloat [::1]    v_Ds

    cdef yaflFloat [:, ::1] v_Pzx

    cdef yaflFloat [::1]    v_Uq
    cdef yaflFloat [::1]    v_Dq

    cdef yaflFloat [::1]    v_Ur
    cdef yaflFloat [::1]    v_Dr

    cdef yaflFloat [:, ::1] v_sigmas_x
    cdef yaflFloat [:, ::1] v_sigmas_z
    cdef yaflFloat [::1]    v_wm
    cdef yaflFloat [::1]    v_wc

    cdef yaflFloat [::1]    v_Sx

    # Kalman filter numpy arrays
    cdef np.ndarray  _z

    cdef np.ndarray  _x
    cdef np.ndarray  _zp
    cdef np.ndarray  _y

    cdef np.ndarray  _Up
    cdef np.ndarray  _Dp

    cdef np.ndarray  _Us
    cdef np.ndarray  _Ds

    cdef np.ndarray  _Pzx

    cdef np.ndarray _Uq
    cdef np.ndarray _Dq

    cdef np.ndarray _Ur
    cdef np.ndarray _Dr

    cdef np.ndarray  _sigmas_x
    cdef np.ndarray  _sigmas_z
    cdef np.ndarray  _wm
    cdef np.ndarray  _wc

    cdef np.ndarray _Sx

    # Callback info
    cdef object    _points

    cdef yaflFloat _dt
    cdef dict      _fx_args
    cdef object    _fx
    cdef object    _mean_x
    cdef object    _resudual_x

    cdef dict      _hx_args
    cdef object    _hx
    cdef object    _mean_z
    cdef object    _residual_z

    #The object will be Extensible
    cdef dict __dict__

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, hx, fx, points,\
                 x_mean_fn=None, z_mean_fn=None, \
                 residual_x=None, residual_z=None):

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
        self.c_self.base.base.f = <yaflUKFFuncP> yafl_py_ukf_fx
        self._fx = fx

        if not callable(hx):
            raise ValueError('hx must be callable!')
        self.c_self.base.base.h = <yaflUKFFuncP> yafl_py_ukf_hx
        self._hx = hx

        if x_mean_fn:
            if not callable(x_mean_fn):
                raise ValueError('x_mean_fn must be callable!')
            self.c_self.base.base.xmf = <yaflUKFFuncP> yafl_py_ukf_xmf
            self._mean_x = x_mean_fn
        else:
            self.c_self.base.base.xmf = <yaflUKFFuncP> 0
            self._mean_x = None

        if residual_x:
            if not callable(residual_x):
                raise ValueError('residual_x must be callable!')
            self.c_self.base.base.xrf = <yaflUKFResFuncP> yafl_py_ukf_xrf
            self._residual_x = residual_x
        else:
            self.c_self.base.base.xrf = <yaflUKFResFuncP> 0
            self._residual_x = None

        if z_mean_fn:
            if not callable(z_mean_fn):
                raise ValueError('z_mean_fn must be callable!')
            self.c_self.base.base.zmf = <yaflUKFFuncP> yafl_py_ukf_zmf
            self._mean_z = z_mean_fn
        else:
            self.c_self.base.base.zmf = <yaflUKFFuncP> 0
            self._mean_z = None

        if residual_z:
            if not callable(residual_z):
                raise ValueError('residual_z must be callable!')
            self.c_self.base.base.zrf = <yaflUKFResFuncP> yafl_py_ukf_zrf
            self._residual_z = residual_z
        else:
            self.c_self.base.base.zrf = <yaflUKFResFuncP> 0
            self._residual_z = None

        #Setup sigma points generator
        if not isinstance(points, yakfSigmaBase):
            raise ValueError('Invalid points type (must be subclass of yakfSigmaBase)!')

        _points = <yakfSigmaBase>points

        self._points = _points
        self.c_self.base.base.sp_info = &_points.c_self.base.base
        self.c_self.base.base.sp_meth = _points.get_spm()

        # Allocate memories and setup the rest of c_self
        # Sigma points and weights
        pnum = _points.pnum

        self._wm  = np.zeros((pnum,), dtype=np.float64)
        self.v_wm = self._wm
        self.c_self.base.base.wm = &self.v_wm[0]

        self._wc  = np.zeros((pnum,), dtype=np.float64)
        self.v_wc = self._wc
        self.c_self.base.base.wc = &self.v_wc[0]

        self._sigmas_x  = np.zeros((pnum, dim_x), dtype=np.float64)
        self.v_sigmas_x = self._sigmas_x
        self.c_self.base.base.sigmas_x = &self.v_sigmas_x[0, 0]

        self._sigmas_z  = np.zeros((pnum, dim_z), dtype=np.float64)
        self.v_sigmas_z = self._sigmas_z
        self.c_self.base.base.sigmas_z = &self.v_sigmas_z[0, 0]

        #Rest of rhe UKF
        self._z  = np.zeros((dim_z,), dtype=np.float64)
        self.v_z = self._z

        self._x  = np.zeros((dim_x,), dtype=np.float64)
        self.v_x = self._x
        self.c_self.base.base.x = &self.v_x[0]

        self._zp  = np.zeros((dim_z,), dtype=np.float64)
        self.v_zp = self._zp
        self.c_self.base.base.zp = &self.v_zp[0]

        self._y  = np.zeros((dim_z,), dtype=np.float64)
        self.v_y = self._y
        self.c_self.base.base.y = &self.v_y[0]

        self._Up  = np.zeros((_U_sz(dim_x),), dtype=np.float64)
        self.v_Up = self._Up
        self.c_self.base.base.Up = &self.v_Up[0]

        self._Dp  = np.ones((dim_x,), dtype=np.float64)
        self.v_Dp = self._Dp
        self.c_self.base.base.Dp = &self.v_Dp[0]

        self._Us  = np.zeros((_U_sz(dim_z),), dtype=np.float64)
        self.v_Us = self._Us
        self.c_self.base.base.Us = &self.v_Us[0]

        self._Ds  = np.ones((dim_z,), dtype=np.float64)
        self.v_Ds = self._Ds
        self.c_self.base.base.Ds = &self.v_Ds[0]

        self._Pzx  = np.zeros((dim_z, dim_x), dtype=np.float64)
        self.v_Pzx = self._Pzx
        self.c_self.base.base.Pzx = &self.v_Pzx[0, 0]

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

        self._Sx  = np.zeros((dim_x,), dtype=np.float64)
        self.v_Sx = self._Sx
        self.c_self.base.base.Sx = &self.v_Sx[0]
        
        #Call C-post init
        yafl_ukf_post_init(&self.c_self.base.base)

    #==========================================================================
    #Decorators
    @property
    def dim_x(self):
        return self.c_self.base.base.Nx

    @property
    def dim_z(self):
        return self.c_self.base.base.Nz
        
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x[:] = value
        
    @property
    def Pzx(self):
        return self._Pzx
    
    @property
    def zp(self):
        return self._zp
    #--------------------------------------------------------------------------
    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        raise AttributeError('yakfUnscentedBase does not support this!')
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

    #--------------------------------------------------------------------------
    @property
    def sigmas_x(self):
        return self._sigmas_x

    @sigmas_x.setter
    def sigmas_x(self, value):
        raise AttributeError('yakfUnscentedBase does not support this!')

    #--------------------------------------------------------------------------
    @property
    def sigmas_z(self):
        return self._sigmas_z

    @sigmas_z.setter
    def sigmas_z(self, value):
        raise AttributeError('yakfUnscentedBase does not support this!')
    #--------------------------------------------------------------------------
    @property
    def wc(self):
        return self._wc

    @wc.setter
    def wc(self, value):
        raise AttributeError('yakfUnscentedBase does not support this!')
    #--------------------------------------------------------------------------
    @property
    def wm(self):
        return self._wm

    @wm.setter
    def wm(self, value):
        raise AttributeError('yakfUnscentedBase does not support this!')
    #==========================================================================
    def _predict(self):
        yafl_ukf_predict(&self.c_self.base.base)

    def predict(self, dt=None, **fx_args):
        old_dt = self._dt

        if dt:
            if np.isnan(dt):
                raise ValueError('Invalid dt value (nan)!')
            self._dt = <yaflFloat>dt

        self._fx_args = fx_args
        self._predict()

        self._dt = old_dt

    #==========================================================================
    def _update(self):
        raise NotImplementedError('yakfUnscentedBase is the base class!')

    def update(self, z, **hx_args):

        self._z[:] = z
        self._hx_args = hx_args
        self._update()

#==============================================================================
cdef void yafl_py_sigma_addf(yakfPyUnscentedSt * self, yaflFloat * delta, \
                             yaflFloat * pivot, yaflFloat mult):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

    _addf = py_self._points._addf
    if not callable(_addf):
        raise ValueError('_addf must be callable!')

    nx = self.base.base.Nx
    if nx <= 0:
        raise ValueError('nx must be > 0!')

    _delta = np.asarray(<yaflFloat[:nx]> delta) #np.float64_t
    _pivot = np.asarray(<yaflFloat[:nx]> pivot) #np.float64_t

    _delta[:] = _addf(_delta, _pivot, mult)

#==============================================================================
cdef void yafl_py_ukf_fx(yakfPyUnscentedSt * self, \
                         yaflFloat * new_x, yaflFloat * old_x):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

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

#------------------------------------------------------------------------------
cdef void yafl_py_ukf_xmf(yakfPyUnscentedSt * self, \
                          yaflFloat * res, yaflFloat * sigmas):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

    mean_x = py_self._mean_x
    if not callable(mean_x):
        raise ValueError('mean_x must be callable!')

    nx = self.base.base.Nx
    if nx <= 0:
        raise ValueError('nx must be > 0!')

    _points = py_self._points
    if not isinstance(py_self, yakfSigmaBase):
        raise ValueError('Invalid _points type (must be subclass of yakfSigmaBase)!')

    pnum = _points.base.base.np
    if pnum <= 0:
        raise ValueError('pnum must be > 0!')

    _sigmas = np.asarray(<yaflFloat[:pnum, :nx]> sigmas) #np.float64_t
    _res    = np.asarray(<yaflFloat[:nx]> res)           #np.float64_t

    _res[:] = mean_x(_sigmas, _points._wm)

#------------------------------------------------------------------------------
cdef void yafl_py_ukf_xrf(yakfPyUnscentedSt * self, yaflFloat * res, \
                             yaflFloat * sigma, yaflFloat * pivot):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

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

#==============================================================================
cdef void yafl_py_ukf_hx(yakfPyUnscentedSt * self, \
                         yaflFloat * z, yaflFloat * x):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

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
    #print(<np.int64_t>x - <np.int64_t>self.base.base.sigmas_x)
    #print(<np.int64_t>z - <np.int64_t>self.base.base.sigmas_z)

    _z[:] = hx(_x, **hx_args)

#------------------------------------------------------------------------------
cdef void yafl_py_ukf_zmf(yakfPyUnscentedSt * self, \
                          yaflFloat * res, yaflFloat * sigmas):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

    mean_z = py_self._mean_z
    if not callable(mean_z):
        raise ValueError('mean_z must be callable!')

    nz = self.base.base.Nz
    if nz <= 0:
        raise ValueError('nz must be > 0!')

    _points = py_self._points
    if not isinstance(py_self, yakfSigmaBase):
        raise ValueError('Invalid _points type (must be subclass of yakfSigmaBase)!')

    pnum = _points.base.base.np
    if pnum <= 0:
        raise ValueError('pnum must be > 0!')

    _sigmas = np.asarray(<yaflFloat[:pnum, :nz]> sigmas) #np.float64_t
    _res    = np.asarray(<yaflFloat[:nz]> res)           #np.float64_t

    _res[:] = mean_z(_sigmas, _points._wm)

#------------------------------------------------------------------------------
cdef void yafl_py_ukf_zrf(yakfPyUnscentedSt * self, yaflFloat * res, \
                             yaflFloat * sigma, yaflFloat * pivot):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

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

#==============================================================================
cdef class Unscented(yakfUnscentedBase):
    """
    UD-factorized UKF implementation
    """
    def _update(self):
        yafl_ukf_update(&self.c_self.base.base, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class UnscentedBierman(yakfUnscentedBase):
    """
    UD-factorized UKF implementation
    """
    def _update(self):
        yafl_ukf_bierman_update(&self.c_self.base.base, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class UnscentedAdaptiveBierman(yakfUnscentedBase):
    def __init__(self, int dim_x, int dim_z, yaflFloat dt, hx, fx, points, \
                 **kwargs):

        super().__init__(dim_x, dim_z, dt, hx, fx, points, **kwargs)

        #Init chi2 with scipy.stats.chi2.ppf(0.99999, 1)
        self.c_self.base.adaptive.chi2 = 19.511420964666268

    """
    UD-factorized UKF implementation
    """
    def _update(self):
        yafl_ukf_adaptive_bierman_update(&self.c_self.base.adaptive, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class MerweSigmaPoints(yakfSigmaBase):
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
