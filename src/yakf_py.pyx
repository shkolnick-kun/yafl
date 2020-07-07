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
cdef extern from "yakf_config.h":
    ctypedef double         yakfFloat
    ctypedef stdint.int32_t yakfInt

#------------------------------------------------------------------------------
cdef extern from "yakf.c":
    #==========================================================================
    #                     UD-factorized EKF definitions
    #==========================================================================
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
    cdef void yakf_base_update(yakfBaseSt * self, yakfFloat * z, \
                               yakfScalarUpdateP scalar_update)

    cdef void yakf_bierman_update(yakfBaseSt * self, yakfFloat * z)

    cdef void yakf_joseph_update(yakfBaseSt * self, yakfFloat * z)

    #--------------------------------------------------------------------------
    ctypedef struct yakfAdaptiveSt:
        yakfBaseSt base
        yakfFloat  chi2

    cdef void yakf_adaptive_bierman_update(yakfAdaptiveSt * self, yakfFloat * z)

    cdef void yakf_adaptive_joseph_update(yakfAdaptiveSt * self, yakfFloat * z)

    # For demonstration purposes only
    cdef void yakf_do_not_use_this_update(yakfAdaptiveSt * self, yakfFloat * z)

    #--------------------------------------------------------------------------
    ctypedef yakfFloat (* yakfRobFuncP)(yakfBaseSt *, yakfFloat)

    ctypedef struct yakfRobustSt:
        yakfBaseSt   base
        yakfRobFuncP g
        yakfRobFuncP gdot

    cdef void yakf_robust_bierman_update(yakfRobustSt * self, yakfFloat * z)

    cdef void yakf_robust_joseph_update(yakfRobustSt * self, yakfFloat * z)
    #--------------------------------------------------------------------------
    ctypedef struct yakfAdaptiveRobustSt:
        yakfRobustSt base
        yakfFloat  chi2

    cdef void yakf_adaptive_robust_bierman_update(yakfAdaptiveRobustSt * self,\
                                                  yakfFloat * z)

    cdef void yakf_adaptive_robust_joseph_update(yakfAdaptiveRobustSt * self, \
                                                 yakfFloat * z)

    #==========================================================================
    #                     UD-factorized UKF definitions
    #==========================================================================
    ctypedef _yakfUnscentedSt yakfUnscentedSt

    ctypedef void (* yakfSigmaAddP)(yakfUnscentedSt *, yakfFloat *, \
                                    yakfFloat *, yakfFloat)

    ctypedef struct yakfSigmaSt:
        yakfInt     np
        yakfFloat * wm
        yakfFloat * wc
        yakfSigmaAddP addf

    #--------------------------------------------------------------------------
    ctypedef void (* yakfSigmaGenWeigthsP)(yakfUnscentedSt *)
    ctypedef void (* yakfSigmaGenSigmasP)(yakfUnscentedSt *)

    ctypedef struct yakfSigmaMethodsSt:
        yakfSigmaGenWeigthsP   wf
        yakfSigmaGenSigmasP  spgf

    #--------------------------------------------------------------------------
    ctypedef void (* yakfUnscentedFuncP)(yakfUnscentedSt *, yakfFloat *, \
                                         yakfFloat *)

    ctypedef void (* yakfUnscentedResFuncP)(yakfUnscentedSt *, yakfFloat *, \
                                            yakfFloat *, yakfFloat *)

    ctypedef struct _yakfUnscentedSt:
        yakfSigmaSt           * points
        const yakfSigmaMethodsSt * spm

        yakfUnscentedFuncP      f
        yakfUnscentedFuncP    xmf
        yakfUnscentedResFuncP xrf

        yakfUnscentedFuncP      h
        yakfUnscentedFuncP    zmf
        yakfUnscentedResFuncP zrf

        yakfFloat * x
        yakfFloat * zp

        yakfFloat * Up
        yakfFloat * Dp

        yakfFloat * Us
        yakfFloat * Ds

        yakfFloat * Pzx

        yakfFloat * Uq
        yakfFloat * Dq

        yakfFloat * Ur
        yakfFloat * Dr

        yakfFloat * sigmas_x
        yakfFloat * sigmas_z

        yakfFloat * Sx
        yakfFloat * Sz

        yakfInt   Nx
        yakfInt   Nz

    #--------------------------------------------------------------------------
    cdef void yakf_unscented_post_init(yakfUnscentedSt * self)  #static inline

    cdef void yakf_unscented_gen_sigmas(yakfUnscentedSt * self) #static inline

    cdef void yakf_unscented_predict(yakfUnscentedSt * self)

    cdef void yakf_unscented_update(yakfUnscentedSt * self, yakfFloat * z)

    #--------------------------------------------------------------------------
    #                  Van der Merwe sigma point generator
    #--------------------------------------------------------------------------
    ctypedef struct yakfMerweSt:
        yakfSigmaSt base
        yakfFloat alpha
        yakfFloat beta
        yakfFloat kappa

    cdef const yakfSigmaMethodsSt yakf_merwe_spm

#------------------------------------------------------------------------------
cdef extern from "yakf_math.c":
    cdef void yakfm_set_u(yakfInt sz, yakfFloat *res, yakfFloat *u)

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
ctypedef union yakfPyEkfBaseUn:
    yakfBaseSt           base
    yakfAdaptiveSt       adaptive
    yakfRobustSt         robust
    yakfAdaptiveRobustSt ada_rob

#------------------------------------------------------------------------------
# Kalman filter C-structure with Python callback info
#------------------------------------------------------------------------------
ctypedef struct yakfPyBaseEkfSt:
    # Kalman filter base union
    yakfPyEkfBaseUn base

    # Python/Cython self
    void * py_self

#------------------------------------------------------------------------------
cdef int _U_sz(int dim_u):
    return max(1, (dim_u * (dim_u - 1))//2)

#------------------------------------------------------------------------------
#                             Basic Filter class
#------------------------------------------------------------------------------
cdef class yakfExtendedBase:
    # Kalman filter C-self
    cdef yakfPyBaseEkfSt c_self
    
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
        self.c_self.base.base.f = <yakfFuncP> yakf_py_ekf_fx
        self._fx = fx

        if not callable(jfx):
            raise ValueError('jfx must be callable!')
        self.c_self.base.base.jf = <yakfFuncP> yakf_py_ekf_jfx
        self._jfx = jfx

        if not callable(hx):
            raise ValueError('hx must be callable!')
        self.c_self.base.base.h = <yakfFuncP> yakf_py_ekf_hx
        self._hx = hx

        if not callable(jhx):
            raise ValueError('jhx must be callable!')
        self.c_self.base.base.jh = <yakfFuncP> yakf_py_ekf_jhx
        self._jhx = jhx

        if residual_z:   
            if not callable(residual_z):
                raise ValueError('residual_z must be callable!')
            self.c_self.base.base.zrf = <yakfResFuncP> yakf_py_ekf_zrf
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
        raise AttributeError('yakfExtendedBase does not support this!')
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
        raise NotImplementedError('yakfExtendedBase is the base class!')
    
    def update(self, z, **hx_args):

        self._z[:] = z           
        self._hx_args = hx_args
        self._update()
        
#==============================================================================
#                             Basic C-callbacks
#==============================================================================
# State transition function 
cdef void yakf_py_ekf_fx(yakfPyBaseEkfSt * self):
    
    py_self = <yakfExtendedBase>(self.py_self)
    if not isinstance(py_self, yakfExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfExtendedBase)!')
    
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
cdef void yakf_py_ekf_jfx(yakfPyBaseEkfSt * self):
    
    py_self = <yakfExtendedBase>(self.py_self)
    if not isinstance(py_self, yakfExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfExtendedBase)!')
       
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
cdef void yakf_py_ekf_hx(yakfPyBaseEkfSt * self):
    
    py_self = <yakfExtendedBase>(self.py_self)
    if not isinstance(py_self, yakfExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfExtendedBase)!')
    
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
cdef void yakf_py_ekf_jhx(yakfPyBaseEkfSt * self):
    
    py_self = <yakfExtendedBase>(self.py_self)
    if not isinstance(py_self, yakfExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfExtendedBase)!')
    
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
cdef void yakf_py_ekf_zrf(yakfPyBaseEkfSt * self, yakfFloat * zp):
    
    py_self = <yakfExtendedBase>(self.py_self)
    if not isinstance(py_self, yakfExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfExtendedBase)!')
    
    zrf = py_self._residual_z
    if not callable(zrf):
        raise ValueError('jhx must be callable!')
    
    #How about handling exceptions here???
    py_self._y[:] = zrf(py_self._z, py_self._y)
    
#==============================================================================
cdef class Bierman(yakfExtendedBase):
    def _update(self):
        yakf_bierman_update(&self.c_self.base.base, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class Joseph(yakfExtendedBase):
    def _update(self):
        yakf_joseph_update(&self.c_self.base.base, &self.v_z[0])

#==============================================================================
#                        Adaptive filter basic class
#==============================================================================
cdef class yakfAdaptiveBase(yakfExtendedBase):
    def __init__(self, int dim_x, int dim_z, yakfFloat dt, \
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
        self.c_self.base.adaptive.chi2 = <yakfFloat>value
#==============================================================================
cdef class AdaptiveBierman(yakfAdaptiveBase):
    def _update(self):
        yakf_adaptive_bierman_update(&self.c_self.base.adaptive, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class AdaptiveJoseph(yakfAdaptiveBase):
    def _update(self):
        yakf_adaptive_joseph_update(&self.c_self.base.adaptive, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class DoNotUseThisFilter(yakfAdaptiveBase):
    """
    WARNING!!!
    DO NOT USE THIS variant of Adaptive Joseph filter !!!
    It was implemented to show some flaws of the corresponding algorithm!
    """
    def _update(self):
        yakf_do_not_use_this_update(&self.c_self.base.adaptive, &self.v_z[0])

#==============================================================================
#                        Robust filter basic class
#==============================================================================
cdef class yakfRobustBase(yakfExtendedBase):
    
    cdef object _gz
    cdef object _gdotz
    
    def __init__(self, int dim_x, int dim_z, yakfFloat dt, \
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
            
            self.c_self.base.robust.g    = <yakfRobFuncP>yakf_py_ekf_rob_gz
            self.c_self.base.robust.gdot = <yakfRobFuncP>yakf_py_ekf_rob_gdotz
            
        else:
            self.c_self.base.robust.g    = <yakfRobFuncP>0
            self.c_self.base.robust.gdot = <yakfRobFuncP>0

#------------------------------------------------------------------------------
# State transition function 
cdef yakfFloat yakf_py_ekf_rob_gz(yakfPyBaseEkfSt * self, yakfFloat nu):

    py_self = <yakfRobustBase>(self.py_self)
    if not isinstance(py_self, yakfExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfExtendedBase)!')

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

    return <yakfFloat>ret

#------------------------------------------------------------------------------
# State transition function Jacobian
cdef yakfFloat yakf_py_ekf_rob_gdotz(yakfPyBaseEkfSt * self, yakfFloat nu):

    py_self = <yakfRobustBase>(self.py_self)
    if not isinstance(py_self, yakfExtendedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfExtendedBase)!')

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

    return <yakfFloat>ret

#==============================================================================
cdef class RobustBierman(yakfRobustBase):
    def _update(self):
        yakf_robust_bierman_update(&self.c_self.base.robust, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class RobustJoseph(yakfRobustBase):
    def _update(self):
        yakf_robust_joseph_update(&self.c_self.base.robust, &self.v_z[0])

#==============================================================================
#                        Robust filter basic class
#==============================================================================
cdef class yakfAdaptiveRobustBase(yakfRobustBase):

    def __init__(self, int dim_x, int dim_z, yakfFloat dt, \
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
        self.c_self.base.ada_rob.chi2 = <yakfFloat>value

#==============================================================================
cdef class AdaptiveRobustBierman(yakfAdaptiveRobustBase):
    def _update(self):
        yakf_adaptive_robust_bierman_update(&self.c_self.base.ada_rob, \
                                            &self.v_z[0])

#------------------------------------------------------------------------------
cdef class AdaptiveRobustJoseph(yakfAdaptiveRobustBase):
    def _update(self):
        yakf_adaptive_robust_joseph_update(&self.c_self.base.ada_rob, \
                                           &self.v_z[0])

#==============================================================================
#                          UD-factorized UKF API
#==============================================================================
#------------------------------------------------------------------------------
#                   Sigma points generator basic definitions
#------------------------------------------------------------------------------
ctypedef union yakfPySigmaBaseUn:
    yakfSigmaSt base
    yakfMerweSt merwe

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

    # Kalman filter memory views
    cdef yakfFloat [::1]    v_wc
    cdef yakfFloat [::1]    v_wm
    cdef yakfFloat [:, ::1] v_sigmas_x
    cdef yakfFloat [:, ::1] v_sigmas_z

    # Kalman filter numpy arrays
    cdef np.ndarray  _wc
    cdef np.ndarray  _wm
    cdef np.ndarray  _sigmas_x
    cdef np.ndarray  _sigmas_z

    # Callback info
    cdef object _addf

    #The object will be Extensible
    cdef dict __dict__

    def __init__(self, yakfInt dim_x, yakfInt dim_z, addf=None):
        #Setup callbacks
        self.c_self.py_self = <void *>self

        if addf:
            if not callable(addf):
                raise ValueError('addf must be callable!')
            self.c_self.base.base.addf = <yakfSigmaAddP>yakf_py_sigma_addf
            self._addf = addf
        else:
            self.c_self.base.base.addf = <yakfSigmaAddP>0

        pnum = self.get_np(dim_x)

        self.c_self.base.base.np = pnum

        self._wc  = np.zeros((pnum,), dtype=np.float64)
        self.v_wc = self._wc
        self.c_self.base.base.wc = &self.v_wc[0]

        self._wm  = np.zeros((pnum,), dtype=np.float64)
        self.v_wm = self._wm
        self.c_self.base.base.wm = &self.v_wm[0]

        self._sigmas_x  = np.zeros((pnum, dim_x), dtype=np.float64)
        self.v_sigmas_x = self._sigmas_x

        self._sigmas_z  = np.zeros((pnum, dim_z), dtype=np.float64)
        self.v_sigmas_z = self._sigmas_z

    cdef yakfInt get_np(self, int dim_x):
        raise NotImplementedError('yakfSigmaBase is the base class!')

    cdef const yakfSigmaMethodsSt * get_spm(self):
        raise NotImplementedError('yakfSigmaBase is the base class!')
        
    #==========================================================================
    #Decorators
    @property
    def sigmas_x(self):
        return self._sigmas_x

    @property
    def sigmas_z(self):
        return self._sigmas_z
    
    @property
    def wc(self):
        return self._wc
    
    @property
    def wm(self):
        return self._wm
#------------------------------------------------------------------------------
#                         UD-factorized UKF definitions
#------------------------------------------------------------------------------
ctypedef union yakfPyUnscentedBaseUn:
    yakfUnscentedSt base

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
    cdef yakfFloat [::1]    v_z

    cdef yakfFloat [::1]    v_x
    cdef yakfFloat [::1]    v_zp

    cdef yakfFloat [::1]    v_Up
    cdef yakfFloat [::1]    v_Dp

    cdef yakfFloat [::1]    v_Us
    cdef yakfFloat [::1]    v_Ds

    cdef yakfFloat [:, ::1] v_Pzx

    cdef yakfFloat [::1]    v_Uq
    cdef yakfFloat [::1]    v_Dq

    cdef yakfFloat [::1]    v_Ur
    cdef yakfFloat [::1]    v_Dr

    cdef yakfFloat [::1]    v_Sx
    cdef yakfFloat [::1]    v_Sz

    # Kalman filter numpy arrays
    cdef np.ndarray  _z

    cdef np.ndarray  _x
    cdef np.ndarray  _zp

    cdef np.ndarray  _Up
    cdef np.ndarray  _Dp

    cdef np.ndarray  _Us
    cdef np.ndarray  _Ds

    cdef np.ndarray  _Pzx

    cdef np.ndarray _Uq
    cdef np.ndarray _Dq

    cdef np.ndarray _Ur
    cdef np.ndarray _Dr

    cdef np.ndarray _Sx
    cdef np.ndarray _Sz

    # Callback info
    cdef object    _points

    cdef yakfFloat _dt
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

    def __init__(self, int dim_x, int dim_z, yakfFloat dt, hx, fx, points,\
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
        self.c_self.base.base.f = <yakfUnscentedFuncP> yakf_py_ukf_fx
        self._fx = fx

        if not callable(hx):
            raise ValueError('hx must be callable!')
        self.c_self.base.base.h = <yakfUnscentedFuncP> yakf_py_ukf_hx
        self._hx = hx

        if x_mean_fn:
            if not callable(x_mean_fn):
                raise ValueError('x_mean_fn must be callable!')
            self.c_self.base.base.xmf = <yakfUnscentedFuncP> yakf_py_ukf_xmf
            self._mean_x = x_mean_fn
        else:
            self.c_self.base.base.xmf = <yakfUnscentedFuncP> 0
            self._mean_x = None

        if residual_x:
            if not callable(residual_x):
                raise ValueError('residual_x must be callable!')
            self.c_self.base.base.xrf = <yakfUnscentedResFuncP> yakf_py_ukf_xrf
            self._residual_x = residual_x
        else:
            self.c_self.base.base.xrf = <yakfUnscentedResFuncP> 0
            self._residual_x = None

        if z_mean_fn:
            if not callable(z_mean_fn):
                raise ValueError('z_mean_fn must be callable!')
            self.c_self.base.base.zmf = <yakfUnscentedFuncP> yakf_py_ukf_zmf
            self._mean_z = z_mean_fn
        else:
            self.c_self.base.base.zmf = <yakfUnscentedFuncP> 0
            self._mean_z = None

        if residual_z:
            if not callable(residual_z):
                raise ValueError('residual_z must be callable!')
            self.c_self.base.base.zrf = <yakfUnscentedResFuncP> yakf_py_ukf_zrf
            self._residual_z = residual_z
        else:
            self.c_self.base.base.zrf = <yakfUnscentedResFuncP> 0
            self._residual_z = None

        # Allocate memories and setup the rest of c_self
        self._z  = np.zeros((dim_z,), dtype=np.float64)
        self.v_z = self._z

        self._x  = np.zeros((dim_x,), dtype=np.float64)
        self.v_x = self._x
        self.c_self.base.base.x = &self.v_x[0]

        self._zp  = np.zeros((dim_z,), dtype=np.float64)
        self.v_zp = self._zp
        self.c_self.base.base.zp = &self.v_zp[0]

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

        self._Sz  = np.zeros((dim_z,), dtype=np.float64)
        self.v_Sz = self._Sz
        self.c_self.base.base.Sz = &self.v_Sz[0]

        #Setup sigma points generator
        if not isinstance(points, yakfSigmaBase):
            raise ValueError('Invalid points type (must be subclass of yakfSigmaBase)!')

        _points = <yakfSigmaBase>points

        self._points = _points

        self.c_self.base.base.sigmas_x = &_points.v_sigmas_x[0, 0]
        self.c_self.base.base.sigmas_z = &_points.v_sigmas_z[0, 0]
        self.c_self.base.base.points   = &_points.c_self.base.base
        self.c_self.base.base.spm      = _points.get_spm()
        
        yakf_unscented_post_init(&self.c_self.base.base)

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
        return self._Sz

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
        yakf_unscented_predict(&self.c_self.base.base)

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
        raise NotImplementedError('yakfUnscentedBase is the base class!')

    def update(self, z, **hx_args):

        self._z[:] = z
        self._hx_args = hx_args
        self._update()

#==============================================================================
cdef void yakf_py_sigma_addf(yakfPyUnscentedSt * self, yakfFloat * delta, \
                             yakfFloat * pivot, yakfFloat mult):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

    _addf = py_self._points._addf
    if not callable(_addf):
        raise ValueError('_addf must be callable!')

    nx = self.base.base.Nx
    if nx <= 0:
        raise ValueError('nx must be > 0!')

    _delta = np.asarray(<yakfFloat[:nx]> delta) #np.float64_t
    _pivot = np.asarray(<yakfFloat[:nx]> pivot) #np.float64_t

    _delta[:] = _addf(_delta, _pivot, mult)

#==============================================================================
cdef void yakf_py_ukf_fx(yakfPyUnscentedSt * self, \
                         yakfFloat * new_x, yakfFloat * old_x):

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

    _new_x = np.asarray(<yakfFloat[:nx]> new_x) #np.float64_t
    _old_x = np.asarray(<yakfFloat[:nx]> old_x) #np.float64_t

    _new_x[:] = fx(_old_x, dt, **fx_args)

#------------------------------------------------------------------------------
cdef void yakf_py_ukf_xmf(yakfPyUnscentedSt * self, \
                          yakfFloat * res, yakfFloat * sigmas):

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

    _sigmas = np.asarray(<yakfFloat[:pnum, :nx]> sigmas) #np.float64_t
    _res    = np.asarray(<yakfFloat[:nx]> res)           #np.float64_t

    _res[:] = mean_x(_sigmas, _points._wm)

#------------------------------------------------------------------------------
cdef void yakf_py_ukf_xrf(yakfPyUnscentedSt * self, yakfFloat * res, \
                             yakfFloat * sigma, yakfFloat * pivot):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

    residual_x = py_self._residual_x
    if not callable(residual_x):
        raise ValueError('residual_x must be callable!')

    nx = self.base.base.Nx
    if nx <= 0:
        raise ValueError('nx must be > 0!')

    _res   = np.asarray(<yakfFloat[:nx]> res)   #np.float64_t
    _sigma = np.asarray(<yakfFloat[:nx]> sigma) #np.float64_t
    _pivot = np.asarray(<yakfFloat[:nx]> pivot) #np.float64_t

    _res[:] = residual_x(_sigma, _pivot)

#==============================================================================
cdef void yakf_py_ukf_hx(yakfPyUnscentedSt * self, \
                         yakfFloat * z, yakfFloat * x):

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

    _x = np.asarray(<yakfFloat[:nx]> x) #np.float64_t
    _z = np.asarray(<yakfFloat[:nz]> z) #np.float64_t
    #print(<np.int64_t>x - <np.int64_t>self.base.base.sigmas_x)
    #print(<np.int64_t>z - <np.int64_t>self.base.base.sigmas_z)

    _z[:] = hx(_x, **hx_args)

#------------------------------------------------------------------------------
cdef void yakf_py_ukf_zmf(yakfPyUnscentedSt * self, \
                          yakfFloat * res, yakfFloat * sigmas):

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

    _sigmas = np.asarray(<yakfFloat[:pnum, :nz]> sigmas) #np.float64_t
    _res    = np.asarray(<yakfFloat[:nz]> res)           #np.float64_t

    _res[:] = mean_z(_sigmas, _points._wm)

#------------------------------------------------------------------------------
cdef void yakf_py_ukf_zrf(yakfPyUnscentedSt * self, yakfFloat * res, \
                             yakfFloat * sigma, yakfFloat * pivot):

    py_self = <yakfUnscentedBase>(self.py_self)
    if not isinstance(py_self, yakfUnscentedBase):
        raise ValueError('Invalid py_self type (must be subclass of yakfUnscentedBase)!')

    residual_z = py_self._residual_z
    if not callable(residual_z):
        raise ValueError('residual_z must be callable!')

    nz = self.base.base.Nz
    if nz <= 0:
        raise ValueError('nx must be > 0!')

    _res   = np.asarray(<yakfFloat[:nz]> res)   #np.float64_t
    _sigma = np.asarray(<yakfFloat[:nz]> sigma) #np.float64_t
    _pivot = np.asarray(<yakfFloat[:nz]> pivot) #np.float64_t

    _res[:] = residual_z(_sigma, _pivot)

#==============================================================================
cdef class Unscented(yakfUnscentedBase):
    """
    UD-factorized UKF implementation
    """
    def _update(self):
        yakf_unscented_update(&self.c_self.base.base, &self.v_z[0])

#------------------------------------------------------------------------------
cdef class MerweSigmaPoints(yakfSigmaBase):
    """
    Van der Merwe sigma point generator implementation
    """
    def __init__(self, yakfInt dim_x, yakfInt dim_z, \
                 yakfFloat alpha, yakfFloat beta, yakfFloat kappa=0.0, **kwargs):
        super().__init__(dim_x, dim_z, **kwargs)
        self.c_self.base.merwe.alpha = alpha
        self.c_self.base.merwe.beta  = beta
        self.c_self.base.merwe.kappa = kappa
        
        
    cdef yakfInt get_np(self, int dim_x):
        return (2 * dim_x + 1)
        
    cdef const yakfSigmaMethodsSt * get_spm(self):
        return &yakf_merwe_spm

    #==========================================================================
    #Decorators
    @property
    def alpha(self):
        return self.c_self.base.merwe.alpha

    @property
    def beta(self):
        return self.c_self.base.merwe.beta
    
    @property
    def kappa(self):
        return self.c_self.base.merwe.kappa
