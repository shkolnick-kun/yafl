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
    and limitations under the License.
"""
#==============================================================================
#                                YAFL C API
#==============================================================================
#cython: language_level=3
#distutils: language=c
from libc cimport stdint
from libc.stdlib cimport malloc, free

#------------------------------------------------------------------------------
IF YAFLPY_USE_64_BIT:
    cdef extern from "yafl_config.h":
        ctypedef double         yaflFloat
        ctypedef stdint.int32_t yaflInt
ELSE:
    cdef extern from "yafl_config.h":
        ctypedef float          yaflFloat
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
        YAFL_ST_INV_ARG_12   = 0x1b0

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_set_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n)
    cdef yaflStatusEn yafl_math_add_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n)
    cdef yaflStatusEn yafl_math_sub_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n)

    cdef yaflStatusEn yafl_math_set_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n)
    cdef yaflStatusEn yafl_math_add_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n)
    cdef yaflStatusEn yafl_math_sub_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n)

    cdef yaflStatusEn yafl_math_set_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_vtv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n)
    cdef yaflStatusEn yafl_math_add_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n)
    cdef yaflStatusEn yafl_math_sub_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_set_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_set_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_add_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_sub_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b)

    cdef yaflStatusEn yafl_math_set_u(yaflInt sz, yaflFloat *res, yaflFloat *u)
    cdef yaflStatusEn yafl_math_add_u(yaflInt sz, yaflFloat *res, yaflFloat *u)
    cdef yaflStatusEn yafl_math_sub_u(yaflInt sz, yaflFloat *res, yaflFloat *u)

    #--------------------------------------------------------------------------
    cdef yaflFloat * _yafl_blk_ptr(yaflFloat * m, yaflInt nc, yaflInt r, yaflInt c)

    cdef yaflStatusEn yafl_math_bset_m(yaflInt nc, yaflFloat *res, yaflInt sr, yaflInt sc, yaflFloat *a)
    cdef yaflStatusEn yafl_math_badd_m(yaflInt nc, yaflFloat *res, yaflInt sr, yaflInt sc, yaflFloat *a)
    cdef yaflStatusEn yafl_math_bsub_m(yaflInt nc, yaflFloat *res, yaflInt sr, yaflInt sc, yaflFloat *a)

    cdef yaflStatusEn yafl_math_bset_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)
    cdef yaflStatusEn yafl_math_badd_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)
    cdef yaflStatusEn yafl_math_bsub_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)

    cdef yaflStatusEn yafl_math_bset_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)
    cdef yaflStatusEn yafl_math_badd_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)
    cdef yaflStatusEn yafl_math_bsub_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)

    cdef yaflStatusEn yafl_math_bset_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v)
    cdef yaflStatusEn yafl_math_badd_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v)
    cdef yaflStatusEn yafl_math_bsub_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_bset_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_badd_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_bsub_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_bset_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_badd_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_bsub_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_bset_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_badd_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b)
    cdef yaflStatusEn yafl_math_bsub_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_ruv(yaflInt sz, yaflFloat *res, yaflFloat *u)
    cdef yaflStatusEn yafl_math_rutv(yaflInt sz, yaflFloat *res, yaflFloat *u)
    cdef yaflStatusEn yafl_math_rum(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *u)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_mwgsu(yaflInt nr, yaflInt nc, yaflFloat *res_u, yaflFloat *res_d, yaflFloat *w, yaflFloat *d)

    #--------------------------------------------------------------------------
    cdef yaflStatusEn yafl_math_udu_up(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v)
    cdef yaflStatusEn yafl_math_udu_down(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v)

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

    ctypedef yaflStatusEn (* yaflKalmanUpdateCBP)(yaflKalmanBaseSt *)
    ctypedef yaflStatusEn (* yaflKalmanUpdateCBP2)(yaflKalmanBaseSt *, \
                                                   yaflFloat *)

    ctypedef yaflFloat (* yaflKalmanRobFuncP)(yaflKalmanBaseSt *, yaflFloat)

    ctypedef struct _yaflKalmanBaseSt:
        yaflKalmanFuncP f       #
        yaflKalmanFuncP h       #
        yaflKalmanResFuncP  zrf #
        yaflKalmanUpdateCBP rcb #

        yaflFloat * x    #
        yaflFloat * y    #

        yaflFloat * Up   #
        yaflFloat * Dp   #

        yaflFloat * Uq   #
        yaflFloat * Dq   #

        yaflFloat * Ur   #
        yaflFloat * Dr   #

        yaflFloat * l    #

        yaflFloat rff    #

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

        yaflFloat qff      #

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
    cdef yaflStatusEn \
    yafl_ekf_do_not_use_this_update_scalar(yaflKalmanBaseSt * self, yaflInt i)

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
        yaflFloat * Sz
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

    #==========================================================================
    #                  Julier sigma point generator
    #==========================================================================
    ctypedef struct yaflUKFJulierSt:
        yaflUKFSigmaSt base
        yaflFloat kappa

    cdef const yaflUKFSigmaMethodsSt yafl_ukf_julier_spm

    #==========================================================================
    #                  Interructing multiple model
    #==========================================================================
    ctypedef struct yaflFilterBankItemSt:
        yaflKalmanBaseSt   * filter
        yaflKalmanUpdateCBP  predict
        yaflKalmanUpdateCBP2 update
        yaflFloat          * Us
        yaflFloat          * Ds
        yaflFloat          * Xs

    ctypedef struct yaflIMMCBSt:
        yaflFilterBankItemSt * bank
        yaflFloat            * mu
        yaflFloat            * M
        yaflFloat            * Up
        yaflFloat            * Dp
        yaflFloat            * x
        yaflFloat            * cbar
        yaflFloat            * omega
        yaflFloat            * y
        yaflFloat            * W
        yaflFloat            * D
        yaflInt                Nb

    cdef yaflStatusEn yafl_imm_post_init(yaflIMMCBSt * self)

    cdef yaflStatusEn yafl_imm_predict(yaflIMMCBSt * self)

    cdef yaflStatusEn yafl_imm_update(yaflIMMCBSt * self, yaflFloat * z)

cdef extern from "immconst.c":
    # EKF function pointers
    cdef const yaflKalmanUpdateCBP  imm_yafl_ekf_base_predict

    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_bierman_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_joseph_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_adaptive_bierman_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_adaptive_joseph_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_robust_bierman_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_robust_joseph_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_adaptive_robust_bierman_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ekf_adaptive_robust_joseph_update

    #==========================================================================
    # EKF function pointers
    cdef const yaflKalmanUpdateCBP  imm_yafl_ukf_base_predict

    # Bierman variants
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ukf_bierman_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ukf_adaptive_bierman_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ukf_robust_bierman_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ukf_adaptive_robust_bierman_update
    # Full variants
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ukf_update
    cdef const yaflKalmanUpdateCBP2 imm_yafl_ukf_adaptive_update

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
IF YAFLPY_USE_64_BIT:
    NP_DTYPE = np.float64
ELSE:
    NP_DTYPE = np.float32

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
ST_INV_ARG_12 = YAFL_ST_INV_ARG_12

#==============================================================================
#                            Helper functions
#==============================================================================
cdef int _U_sz(int dim_u):
    assert 0 < dim_u
    return max(1, (dim_u * (dim_u - 1))//2)

#==============================================================================
#                            Test functions
#==============================================================================
def _set_vxn(v, yaflFloat n):
    assert isinstance(v, np.ndarray)
    assert NP_DTYPE == v.dtype
    assert 1 == len(v.shape)

    sz = v.shape[0]
    res = np.zeros((sz,), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_v = v
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_set_vxn(sz, &v_r[0], &v_v[0], n)
    return status, res

#------------------------------------------------------------------------------
def _add_vxn(a, b, yaflFloat n):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz

    res = a.copy()

    cdef yaflFloat [::1]    v_v = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_add_vxn(sz, &v_r[0], &v_v[0], n)
    return status, res

#------------------------------------------------------------------------------
def _sub_vxn(a, b, yaflFloat n):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz

    res = a.copy()

    cdef yaflFloat [::1]    v_v = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_sub_vxn(sz, &v_r[0], &v_v[0], n)
    return status, res

#------------------------------------------------------------------------------
def _set_vrn(v, yaflFloat n):
    assert isinstance(v, np.ndarray)
    assert NP_DTYPE == v.dtype
    assert 1 == len(v.shape)

    sz = v.shape[0]
    res = np.zeros((sz,), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_v = v
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_set_vrn(sz, &v_r[0], &v_v[0], n)
    return status, res

#------------------------------------------------------------------------------
def _add_vrn(a, b, yaflFloat n):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz

    res = a.copy()

    cdef yaflFloat [::1]    v_v = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_add_vrn(sz, &v_r[0], &v_v[0], n)
    return status, res

#------------------------------------------------------------------------------
def _sub_vrn(a, b, yaflFloat n):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz

    res = a.copy()

    cdef yaflFloat [::1]    v_v = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_sub_vrn(sz, &v_r[0], &v_v[0], n)
    return status, res

#------------------------------------------------------------------------------
def _set_vxv(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz

    res = np.zeros((sz,), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_set_vxv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _add_vxv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz
    assert res.shape[0] == sz

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_add_vxv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _sub_vxv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz
    assert res.shape[0] == sz

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_sub_vxv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _set_vrv(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz

    res = np.zeros((sz,), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_set_vrv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _add_vrv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz
    assert res.shape[0] == sz

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_add_vrv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _sub_vrv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz
    assert res.shape[0] == sz

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_sub_vrv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _vtv(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert b.shape[0] == sz

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat res = 0

    status = yafl_math_vtv(sz, &res, &v_a[0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _set_vvt(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    nr = a.shape[0]
    nc = b.shape[0]

    res = np.zeros((nr,nc), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_set_vvt(nr, nc, &v_r[0,0], &v_a[0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _add_vvt(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr = a.shape[0]
    nc = b.shape[0]
    assert res.shape[0] == nr
    assert res.shape[1] == nc

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_add_vvt(nr, nc, &v_r[0,0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _sub_vvt(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr = a.shape[0]
    nc = b.shape[0]
    assert res.shape[0] == nr
    assert res.shape[1] == nc

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_sub_vvt(nr, nc, &v_r[0,0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _set_vvtxn(a, b, yaflFloat n):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    nr = a.shape[0]
    nc = b.shape[0]

    res = np.zeros((nr,nc), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_set_vvtxn(nr, nc, &v_r[0,0], &v_a[0], &v_b[0], n)
    return status, res

#------------------------------------------------------------------------------
def _add_vvtxn(res, a, b, yaflFloat n):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr = a.shape[0]
    nc = b.shape[0]
    assert res.shape[0] == nr
    assert res.shape[1] == nc

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_add_vvtxn(nr, nc, &v_r[0,0], &v_a[0], &v_b[0], n)
    return status

#------------------------------------------------------------------------------
def _sub_vvtxn(res, a, b, yaflFloat n):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr = a.shape[0]
    nc = b.shape[0]
    assert res.shape[0] == nr
    assert res.shape[1] == nc

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_sub_vvtxn(nr, nc, &v_r[0,0], &v_a[0], &v_b[0], n)
    return status

#------------------------------------------------------------------------------
def _set_mv(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    nr = a.shape[0]
    nc = a.shape[1]

    assert b.shape[0] == nc

    res = np.zeros((nr,), dtype=NP_DTYPE)

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_set_mv(nr, nc, &v_r[0], &v_a[0,0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _add_mv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    nr = a.shape[0]
    nc = a.shape[1]

    assert b.shape[0]   == nc
    assert res.shape[0] == nr

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_add_mv(nr, nc, &v_r[0], &v_a[0,0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _sub_mv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    nr = a.shape[0]
    nc = a.shape[1]

    assert b.shape[0]   == nc
    assert res.shape[0] == nr

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_sub_mv(nr, nc, &v_r[0], &v_a[0,0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _set_vtm(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 2 == len(b.shape)

    nr = b.shape[0]
    nc = b.shape[1]

    assert a.shape[0] == nr

    res = np.zeros((nc,), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [:, ::1] v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_set_vtm(nr, nc, &v_r[0], &v_a[0], &v_b[0,0])
    return status, res

#------------------------------------------------------------------------------
def _add_vtm(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 2 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    nr = b.shape[0]
    nc = b.shape[1]

    assert a.shape[0]   == nr
    assert res.shape[0] == nc

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [:, ::1] v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_add_vtm(nr, nc, &v_r[0], &v_a[0], &v_b[0,0])
    return status

#------------------------------------------------------------------------------
def _sub_vtm(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 2 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    nr = b.shape[0]
    nc = b.shape[1]

    assert a.shape[0]   == nr
    assert res.shape[0] == nc

    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [:, ::1] v_b = b
    cdef yaflFloat [::1]    v_r = res

    status = yafl_math_sub_vtm(nr, nc, &v_r[0], &v_a[0], &v_b[0,0])
    return status

#------------------------------------------------------------------------------
def _set_mm(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 2 == len(b.shape)

    nr =  a.shape[0]
    ncr = a.shape[1]

    assert b.shape[0] == ncr
    nc =  b.shape[1]

    res = np.zeros((nr,nc), dtype=NP_DTYPE)

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [:, ::1] v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_set_mm(nr, ncr, nc, &v_r[0,0], &v_a[0,0], &v_b[0,0])
    return status, res

#------------------------------------------------------------------------------
def _add_mm(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 2 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr =  a.shape[0]
    ncr = a.shape[1]

    assert b.shape[0] == ncr
    nc =  b.shape[1]

    assert res.shape[0] == nr
    assert res.shape[1] == nc

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [:, ::1] v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_add_mm(nr, ncr, nc, &v_r[0,0], &v_a[0,0], &v_b[0,0])
    return status

#------------------------------------------------------------------------------
def _sub_mm(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 2 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr =  a.shape[0]
    ncr = a.shape[1]

    assert b.shape[0] == ncr
    nc =  b.shape[1]

    assert res.shape[0] == nr
    assert res.shape[1] == nc

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [:, ::1] v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_sub_mm(nr, ncr, nc, &v_r[0,0], &v_a[0,0], &v_b[0,0])
    return status

#------------------------------------------------------------------------------
def _set_u(u):
    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    n = u.shape[0]
    sz = int(round((1. + (1. + 8. * n)**0.5)/2.))

    assert _U_sz(sz) == n

    res = np.zeros((sz,sz), dtype=NP_DTYPE)

    cdef yaflFloat [::1]    v_u = u
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_set_u(sz, &v_r[0,0], &v_u[0])
    return status, res

#------------------------------------------------------------------------------
def _add_u(res, u):
    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    n = u.shape[0]
    sz = res.shape[0]
    assert res.shape[1] == sz
    assert _U_sz(sz) == n

    cdef yaflFloat [::1]    v_u = u
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_add_u(sz, &v_r[0,0], &v_u[0])
    return status

#------------------------------------------------------------------------------
def _sub_u(res, u):
    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    n = u.shape[0]
    sz = res.shape[0]
    assert res.shape[1] == sz
    assert _U_sz(sz) == n

    cdef yaflFloat [::1]    v_u = u
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_sub_u(sz, &v_r[0,0], &v_u[0])
    return status

#------------------------------------------------------------------------------
def _set_vtu(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = a.shape[0]
    assert _U_sz(sz) == b.shape[0]

    res = np.zeros((sz,), dtype=NP_DTYPE)

    cdef yaflFloat [::1] v_a = a
    cdef yaflFloat [::1] v_b = b
    cdef yaflFloat [::1] v_r = res

    status = yafl_math_set_vtu(sz, &v_r[0], &v_a[0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _add_vtu(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = a.shape[0]
    assert _U_sz(sz) == b.shape[0]
    assert res.shape[0] == sz

    cdef yaflFloat [::1] v_a = a
    cdef yaflFloat [::1] v_b = b
    cdef yaflFloat [::1] v_r = res

    status = yafl_math_add_vtu(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _sub_vtu(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = a.shape[0]
    assert _U_sz(sz) == b.shape[0]
    assert res.shape[0] == sz

    cdef yaflFloat [::1] v_a = a
    cdef yaflFloat [::1] v_b = b
    cdef yaflFloat [::1] v_r = res

    status = yafl_math_sub_vtu(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _set_uv(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    sz = b.shape[0]
    assert _U_sz(sz) == a.shape[0]

    res = np.zeros((sz,), dtype=NP_DTYPE)

    cdef yaflFloat [::1] v_a = a
    cdef yaflFloat [::1] v_b = b
    cdef yaflFloat [::1] v_r = res

    status = yafl_math_set_uv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _add_uv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = b.shape[0]
    assert _U_sz(sz) == a.shape[0]
    assert res.shape[0] == sz

    cdef yaflFloat [::1] v_a = a
    cdef yaflFloat [::1] v_b = b
    cdef yaflFloat [::1] v_r = res

    status = yafl_math_add_uv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _sub_uv(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 1 == len(res.shape)

    sz = b.shape[0]
    assert _U_sz(sz) == a.shape[0]
    assert res.shape[0] == sz

    cdef yaflFloat [::1] v_a = a
    cdef yaflFloat [::1] v_b = b
    cdef yaflFloat [::1] v_r = res

    status = yafl_math_sub_uv(sz, &v_r[0], &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _set_mu(a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    nr = a.shape[0]
    nc = a.shape[1]
    assert _U_sz(nc) == b.shape[0]

    res = np.zeros((nr,nc), dtype=NP_DTYPE)

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_set_mu(nr, nc, &v_r[0,0], &v_a[0,0], &v_b[0])
    return status, res

#------------------------------------------------------------------------------
def _add_mu(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr = a.shape[0]
    nc = a.shape[1]
    assert _U_sz(nc) == b.shape[0]
    assert res.shape == a.shape

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_add_mu(nr, nc, &v_r[0,0], &v_a[0,0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _sub_mu(res, a, b):
    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    nr = a.shape[0]
    nc = a.shape[1]
    assert _U_sz(nc) == b.shape[0]
    assert res.shape == a.shape

    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b
    cdef yaflFloat [:, ::1] v_r = res

    status = yafl_math_sub_mu(nr, nc, &v_r[0,0], &v_a[0,0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _bset_m(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    nr = a.shape[0]
    nc = a.shape[1]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a

    status = yafl_math_bset_m(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              nr, nc, &v_a[0,0])
    return status

#------------------------------------------------------------------------------
def _badd_m(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    nr = a.shape[0]
    nc = a.shape[1]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a

    status = yafl_math_badd_m(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              nr, nc, &v_a[0,0])
    return status

#------------------------------------------------------------------------------
def _bsub_m(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    nr = a.shape[0]
    nc = a.shape[1]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a

    status = yafl_math_bsub_m(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              nr, nc, &v_a[0,0])
    return status

#------------------------------------------------------------------------------
def _bset_u(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    n = a.shape[0]
    sz = int(round((1. + (1. + 8. * n)**0.5)/2.))

    assert _U_sz(sz) == n

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_bset_u(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _badd_u(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    n = a.shape[0]
    sz = int(round((1. + (1. + 8. * n)**0.5)/2.))

    assert _U_sz(sz) == n

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_badd_u(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _bsub_u(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    n = a.shape[0]
    sz = int(round((1. + (1. + 8. * n)**0.5)/2.))

    assert _U_sz(sz) == n

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_bsub_u(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _bset_ut(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    n = a.shape[0]
    sz = int(round((1. + (1. + 8. * n)**0.5)/2.))

    assert _U_sz(sz) == n

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_bset_ut(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _badd_ut(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    n = a.shape[0]
    sz = int(round((1. + (1. + 8. * n)**0.5)/2.))

    assert _U_sz(sz) == n

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_badd_ut(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _bsub_ut(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    n = a.shape[0]
    sz = int(round((1. + (1. + 8. * n)**0.5)/2.))

    assert _U_sz(sz) == n

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_bsub_ut(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _bset_v(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    sz = a.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + 1  <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_bset_v(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _badd_v(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    sz = a.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + 1  <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_badd_v(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _bsub_v(res, a, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    sz = a.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + 1  <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a

    status = yafl_math_bsub_v(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0])
    return status

#------------------------------------------------------------------------------
def _bset_vvt(res, a, b, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    sz = a.shape[0]
    assert b.shape[0] == sz

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_bset_vvt(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _badd_vvt(res, a, b, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    sz = a.shape[0]
    assert b.shape[0] == sz

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_badd_vvt(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _bsub_vvt(res, a, b, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 1 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    sz = a.shape[0]
    assert b.shape[0] == sz

    assert r >= 0
    assert c >= 0
    assert r + sz <= rnr
    assert c + sz <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_bsub_vvt(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              sz, &v_a[0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _bset_mu(res, a, b, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    nr = a.shape[0]
    nc = a.shape[1]

    assert _U_sz(nc) == b.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_bset_mu(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              nr, nc, &v_a[0,0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _badd_mu(res, a, b, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    nr = a.shape[0]
    nc = a.shape[1]

    assert _U_sz(nc) == b.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_badd_mu(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              nr, nc, &v_a[0,0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _bsub_mu(res, a, b, yaflInt r, yaflInt c):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    nr = a.shape[0]
    nc = a.shape[1]

    assert _U_sz(nc) == b.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_bsub_mu(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), \
                              nr, nc, &v_a[0,0], &v_b[0])
    return status

#------------------------------------------------------------------------------
def _bset_bu(res, a, b, yaflInt r, yaflInt c, yaflInt ar, yaflInt ac, yaflInt nr, yaflInt nc):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    anr = a.shape[0]
    anc = a.shape[1]

    assert _U_sz(nc) == b.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc
    assert ar + nr <= anr
    assert ac + nc <= anc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_bset_bu(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), nr, nc, \
                               anc, _yafl_blk_ptr(&v_a[0,0], anc, ar, ac), &v_b[0])
    return status

#------------------------------------------------------------------------------
def _badd_bu(res, a, b, yaflInt r, yaflInt c, yaflInt ar, yaflInt ac, yaflInt nr, yaflInt nc):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    anr = a.shape[0]
    anc = a.shape[1]

    assert _U_sz(nc) == b.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc
    assert ar + nr <= anr
    assert ac + nc <= anc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_badd_bu(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), nr, nc, \
                               anc, _yafl_blk_ptr(&v_a[0,0], anc, ar, ac), &v_b[0])
    return status

#------------------------------------------------------------------------------
def _bsub_bu(res, a, b, yaflInt r, yaflInt c, yaflInt ar, yaflInt ac, yaflInt nr, yaflInt nc):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(a, np.ndarray)
    assert NP_DTYPE == a.dtype
    assert 2 == len(a.shape)

    assert isinstance(b, np.ndarray)
    assert NP_DTYPE == b.dtype
    assert 1 == len(b.shape)

    rnr = res.shape[0]
    rnc = res.shape[1]

    anr = a.shape[0]
    anc = a.shape[1]

    assert _U_sz(nc) == b.shape[0]

    assert r >= 0
    assert c >= 0
    assert r + nr <= rnr
    assert c + nc <= rnc
    assert ar + nr <= anr
    assert ac + nc <= anc

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [:, ::1] v_a = a
    cdef yaflFloat [::1]    v_b = b

    status = yafl_math_bsub_bu(rnc, _yafl_blk_ptr(&v_r[0,0], rnc, r, c), nr, nc, \
                               anc, _yafl_blk_ptr(&v_a[0,0], anc, ar, ac), &v_b[0])
    return status

#------------------------------------------------------------------------------
def _ruv(u, v):
    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    assert isinstance(v, np.ndarray)
    assert NP_DTYPE == v.dtype
    assert 1 == len(v.shape)
    sz = v.shape[0]

    assert u.shape[0] == _U_sz(sz)

    v = v.copy()
    cdef yaflFloat [::1]    v_v = v
    cdef yaflFloat [::1]    v_u = u
    res = yafl_math_ruv(sz, &v_v[0], &v_u[0])
    return res, v

#------------------------------------------------------------------------------
def _rutv(u, v):
    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    assert isinstance(v, np.ndarray)
    assert NP_DTYPE == v.dtype
    assert 1 == len(v.shape)
    sz = v.shape[0]

    assert u.shape[0] == _U_sz(sz)

    v = v.copy()
    cdef yaflFloat [::1]    v_v = v
    cdef yaflFloat [::1]    v_u = u
    res = yafl_math_rutv(sz, &v_v[0], &v_u[0])
    return res, v

#------------------------------------------------------------------------------
def _rum(res, u):
    assert isinstance(res, np.ndarray)
    assert NP_DTYPE == res.dtype
    assert 2 == len(res.shape)

    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    nr = res.shape[0]
    nc = res.shape[1]
    assert _U_sz(nr) == u.shape[0]

    cdef yaflFloat [:, ::1] v_r = res
    cdef yaflFloat [::1]    v_u = u

    status = yafl_math_rum(nr, nc, &v_r[0,0], &v_u[0])
    return status

#------------------------------------------------------------------------------
def _mwgs(W, D):
    assert isinstance(W, np.ndarray)
    assert NP_DTYPE == W.dtype
    assert 2 == len(W.shape)
    nr = W.shape[0]
    nc = W.shape[1]

    assert isinstance(D, np.ndarray)
    assert NP_DTYPE == D.dtype
    assert 1 == len(D.shape)
    assert D.shape[0] == nc

    u = np.zeros((_U_sz(nr),), dtype=NP_DTYPE)
    d = np.ones((nr,), dtype=NP_DTYPE)

    cdef yaflFloat [:, ::1] v_W = W
    cdef yaflFloat [::1]    v_D = D
    cdef yaflFloat [::1]    v_u = u
    cdef yaflFloat [::1]    v_d = d

    res = yafl_math_mwgsu(nr, nc, &v_u[0], &v_d[0], &v_W[0,0], &v_D[0])

    return res, u, d

#------------------------------------------------------------------------------
def _udu_up(u,d, yaflFloat alpha, v):
    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    assert isinstance(d, np.ndarray)
    assert NP_DTYPE == d.dtype
    assert 1 == len(d.shape)

    assert isinstance(v, np.ndarray)
    assert NP_DTYPE == v.dtype
    assert 1 == len(v.shape)

    assert alpha > 0.0

    sz = d.shape[0]
    assert _U_sz(sz) == u.shape[0]
    assert v.shape[0] == sz

    cdef yaflFloat [::1]    v_u = u
    cdef yaflFloat [::1]    v_d = d
    cdef yaflFloat [::1]    v_v = v

    return yafl_math_udu_up(sz, &v_u[0], &v_d[0], alpha, &v_v[0])

#------------------------------------------------------------------------------
def _udu_down(u,d, yaflFloat alpha, v):
    assert isinstance(u, np.ndarray)
    assert NP_DTYPE == u.dtype
    assert 1 == len(u.shape)

    assert isinstance(d, np.ndarray)
    assert NP_DTYPE == d.dtype
    assert 1 == len(d.shape)

    assert isinstance(v, np.ndarray)
    assert NP_DTYPE == v.dtype
    assert 1 == len(v.shape)

    assert alpha > 0.0

    sz = d.shape[0]
    assert _U_sz(sz) == u.shape[0]
    assert v.shape[0] == sz

    cdef yaflFloat [::1]    v_u = u
    cdef yaflFloat [::1]    v_d = d
    cdef yaflFloat [::1]    v_v = v

    return yafl_math_udu_down(sz, &v_u[0], &v_d[0], alpha, &v_v[0])

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
#                             Basic Filter class
#------------------------------------------------------------------------------
cdef class yaflKalmanBase:
    # Kalman filter C-self
    cdef yaflPyKalmanBaseSt c_self

    cdef yaflFloat _l

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

    cdef object    _rcb

    cdef dict      _hx_args
    cdef object    _hx

    cdef object    _residual_z

    def __init__(self, int dim_x, int dim_z, yaflFloat dt, \
                 fx, hx, residual_z = None, rcb = None):

        #Store dimensions
        self.c_self.base.base.Nx = dim_x
        self.c_self.base.base.Nz = dim_z

        #Set likelihood storage
        self.c_self.base.base.l = &self._l

        #Set forgetting factors
        self.c_self.base.base.rff = 0.0

        #Setup callbacks
        self.c_self.py_self = <void *>self

        self._dt = dt
        self._fx_args = {}
        self._hx_args = {}

        if fx:
            if not callable(fx):
                raise ValueError('fx must be callable!')
            self.c_self.base.base.f = <yaflKalmanFuncP>yafl_py_kalman_fx
            self._fx = fx
        else:
            self.c_self.base.base.f = <yaflKalmanFuncP>0
            self._fx = None

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

        if rcb:
            if not callable(rcb):
                raise ValueError('rcb must be callable!')
            self.c_self.base.base.rcb = <yaflKalmanUpdateCBP>yafl_py_kalman_rcb
            self._rcb = rcb
        else:
            self.c_self.base.base.rcb = <yaflKalmanUpdateCBP>0
            self._rcb = None

        # Allocate memories and setup the rest of c_self
        self._z  = np.zeros((dim_z,), dtype=NP_DTYPE)
        self.v_z = self._z

        self._x  = np.zeros((dim_x,), dtype=NP_DTYPE)
        self.v_x = self._x
        self.c_self.base.base.x = &self.v_x[0]

        self._y  = np.zeros((dim_z,), dtype=NP_DTYPE)
        self.v_y = self._y
        self.c_self.base.base.y = &self.v_y[0]

        self._Up  = np.zeros((_U_sz(dim_x),), dtype=NP_DTYPE)
        self.v_Up = self._Up
        self.c_self.base.base.Up = &self.v_Up[0]

        self._Dp  = np.ones((dim_x,), dtype=NP_DTYPE)
        self.v_Dp = self._Dp
        self.c_self.base.base.Dp = &self.v_Dp[0]

        self._Uq  = np.zeros((_U_sz(dim_x),), dtype=NP_DTYPE)
        self.v_Uq = self._Uq
        self.c_self.base.base.Uq = &self.v_Uq[0]

        self._Dq  = np.ones((dim_x,), dtype=NP_DTYPE)
        self.v_Dq = self._Dq
        self.c_self.base.base.Dq = &self.v_Dq[0]

        self._Ur  = np.zeros((_U_sz(dim_z),), dtype=NP_DTYPE)
        self.v_Ur = self._Ur
        self.c_self.base.base.Ur = &self.v_Ur[0]

        self._Dr  = np.ones((dim_z,), dtype=NP_DTYPE)
        self.v_Dr = self._Dr
        self.c_self.base.base.Dr = &self.v_Dr[0]

    #==========================================================================
    cdef yaflPyKalmanBaseUn * cbase(self):
        return &self.c_self.base

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

    #--------------------------------------------------------------------------
    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, value):
        raise AttributeError('yaflKalmanBase does not support this!')

    #--------------------------------------------------------------------------
    @property
    def P(self):
        #
        if self.c_self.base.base.Nx > 1:
            _,u = _set_u(self._Up)
            return u.dot(np.diag(self._Dp).dot(u.T))
        #
        return self._Dp.copy().reshape((1,1))

    @P.setter
    def P(self, value):
        raise AttributeError('yaflKalmanBase does not support this!')

    #--------------------------------------------------------------------------
    @property
    def Q(self):
        #
        if self.c_self.base.base.Nx > 1:
            _,u = _set_u(self._Uq)
            return u.dot(np.diag(self._Dq).dot(u.T))
        #
        return self._Dq.copy().reshape((1,1))

    @Q.setter
    def Q(self, value):
        raise AttributeError('yaflKalmanBase does not support this!')

    #--------------------------------------------------------------------------
    @property
    def R(self):
        #
        if self.c_self.base.base.Nz > 1:
            _,u = _set_u(self._Ur)
            return u.dot(np.diag(self._Dr).dot(u.T))
        #
        return self._Dr.copy().reshape((1,1))

    @R.setter
    def R(self, value):
        raise AttributeError('yaflKalmanBase does not support this!')

    #--------------------------------------------------------------------------
    @property
    def rff(self):
        return self.c_self.base.base.rff

    @rff.setter
    def rff(self, yaflFloat value):
        if value < 0.0 or value >= 1.0:
            raise ValueError('rff value must be in [0.0, 1.0)')

        self.c_self.base.base.rff = value

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

        _new_x = np.asarray(<yaflFloat[:nx]> new_x) #
        _old_x = np.asarray(<yaflFloat[:nx]> old_x) #

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
        if nz <= 0:
            raise ValueError('nz must be > 0!')

        _x = np.asarray(<yaflFloat[:nx]> x) #
        _z = np.asarray(<yaflFloat[:nz]> z) #

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
            raise ValueError('nz must be > 0!')

        _res   = np.asarray(<yaflFloat[:nz]> res)   #
        _sigma = np.asarray(<yaflFloat[:nz]> sigma) #
        _pivot = np.asarray(<yaflFloat[:nz]> pivot) #

        _res[:] = residual_z(_sigma, _pivot)

        return YAFL_ST_OK

    except Exception as e:
        print(tb.format_exc())
        return YAFL_ST_INV_ARG_1

#------------------------------------------------------------------------------
cdef yaflStatusEn yafl_py_kalman_rcb(yaflPyKalmanBaseSt * self):
    try:
        if not isinstance(<object>(self.py_self), yaflKalmanBase):
            raise ValueError('Invalid py_self type (must be subclass of yaflKalmanBase)!')

        py_self = <yaflKalmanBase>(self.py_self)

        rcb = py_self._rcb
        if not callable(rcb):
            raise ValueError('residual_z must be callable!')

        hx_args = py_self._hx_args
        if not isinstance(hx_args, dict):
            raise ValueError('Invalid hx_args type (must be dict)!')

        rcb(**hx_args)

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
    cdef yaflFloat [::1]    v_W
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

        self.c_self.base.ekf.qff = 0.0

        if not callable(jfx):
            raise ValueError('jfx must be callable!')
        self.c_self.base.ekf.jf = <yaflKalmanFuncP>yafl_py_ekf_jfx
        self._jfx = jfx

        if not callable(jhx):
            raise ValueError('jhx must be callable!')
        self.c_self.base.ekf.jh = <yaflKalmanFuncP>yafl_py_ekf_jhx
        self._jhx = jhx


        # Allocate memories and setup the rest of c_self
        self._H  = np.zeros((dim_z, dim_x), dtype=NP_DTYPE)
        self.v_H = self._H
        self.c_self.base.ekf.H = &self.v_H[0, 0]

        # We may need some memory for P and R updates
        n = max(2 * dim_x * dim_x, (dim_x + dim_z) * dim_z)
        self._W  = np.zeros((n,), dtype=NP_DTYPE)
        self.v_W = self._W
        self.c_self.base.ekf.W = &self.v_W[0]

        n = max(2 * dim_x, dim_x + dim_z)
        self._D  = np.ones((n,), dtype=NP_DTYPE)
        self.v_D = self._D
        self.c_self.base.ekf.D = &self.v_D[0]

    #--------------------------------------------------------------------------
    @property
    def qff(self):
        return self.c_self.base.ekf.qff

    @qff.setter
    def qff(self, yaflFloat value):
        if value < 0.0 or value >= 1.0:
            raise ValueError('qff value must be in [0.0, 1.0)')

        self.c_self.base.ekf.qff = value

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

        # Read x size
        dim_x = self.base.base.Nx

        # Interpret x
        _x = np.asarray(<yaflFloat[:dim_x]> x)

        # Interpret w.shape as (dim_x, 2 * dim_x)
        _w = np.asarray(<yaflFloat[:dim_x, :2 * dim_x]> w)

        # Calculate Jacobian
        _w[:, :dim_x] = jfx(_x, dt, **fx_args)

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

        # Read x and z size
        dim_x = self.base.base.Nx
        dim_z = self.base.base.Nz

        # Interpret x
        _x = np.asarray(<yaflFloat[:dim_x]> x)

        # Interpret h.shape as (dim_x, 2 * dim_x)
        _h = np.asarray(<yaflFloat[:dim_z, :dim_x]> h)

        # Calculate Jacobian
        _h[:,:] = jhx(_x, **hx_args)

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
cdef class DoNotUseThisFilter(yaflAdaptiveBase):
    """
    WARNING!!!
    DO NOT USE THIS variant of Adaptive Joseph filter !!!
    It was implemented to show some flaws of the corresponding algorithm!
    """
    def _update(self):
        return yafl_ekf_base_update(&self.c_self.base.base, &self.v_z[0], \
                                    yafl_ekf_do_not_use_this_update_scalar)

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
    yaflUKFJulierSt julier

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
    cdef yaflFloat [::1]    v_Sz
    cdef yaflFloat [:, ::1] v_Pzx

    cdef yaflFloat [::1]    v_sigmas_x
    cdef yaflFloat [:, ::1] v_sigmas_z
    cdef yaflFloat [::1]    v_wm
    cdef yaflFloat [::1]    v_wc

    # Kalman filter numpy arrays
    cdef np.ndarray  _zp

    cdef np.ndarray  _Sx
    cdef np.ndarray  _Sz
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

        self._wm  = np.zeros((pnum,), dtype=NP_DTYPE)
        self.v_wm = self._wm
        self.c_self.base.ukf.wm = &self.v_wm[0]

        self._wc  = np.zeros((pnum,), dtype=NP_DTYPE)
        self.v_wc = self._wc
        self.c_self.base.ukf.wc = &self.v_wc[0]

        #Memory pool for sigmas_z and R updates
        n = max(pnum * dim_x, (dim_x + dim_z) * (dim_z + 1))
        self._sigmas_x  = np.zeros((n,), dtype=NP_DTYPE)
        self.v_sigmas_x = self._sigmas_x
        self.c_self.base.ukf.sigmas_x = &self.v_sigmas_x[0]

        self._sigmas_z  = np.zeros((pnum, dim_z), dtype=NP_DTYPE)
        self.v_sigmas_z = self._sigmas_z
        self.c_self.base.ukf.sigmas_z = &self.v_sigmas_z[0, 0]

        #Rest of rhe UKF
        self._zp  = np.zeros((dim_z,), dtype=NP_DTYPE)
        self.v_zp = self._zp
        self.c_self.base.ukf.zp = &self.v_zp[0]

        self._Pzx  = np.zeros((dim_z, dim_x), dtype=NP_DTYPE)
        self.v_Pzx = self._Pzx
        self.c_self.base.ukf.Pzx = &self.v_Pzx[0, 0]

        self._Sx  = np.zeros((dim_x,), dtype=NP_DTYPE)
        self.v_Sx = self._Sx
        self.c_self.base.ukf.Sx = &self.v_Sx[0]

        self._Sz  = np.zeros((dim_z,), dtype=NP_DTYPE)
        self.v_Sz = self._Sz
        self.c_self.base.ukf.Sz = &self.v_Sz[0]

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
        dim_x = self.c_self.base.base.Nx
        pnum  = self._points.pnum
        return self._sigmas_x[:pnum * dim_x].copy().reshape((pnum, dim_x))

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

        _delta = np.asarray(<yaflFloat[:nx]> delta) #
        _pivot = np.asarray(<yaflFloat[:nx]> pivot) #

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

        _sigmas = np.asarray(<yaflFloat[:pnum, :nx]> sigmas) #
        _res    = np.asarray(<yaflFloat[:nx]> res)           #

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

        _res   = np.asarray(<yaflFloat[:nx]> res)   #
        _sigma = np.asarray(<yaflFloat[:nx]> sigma) #
        _pivot = np.asarray(<yaflFloat[:nx]> pivot) #

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

        _sigmas = np.asarray(<yaflFloat[:pnum, :nz]> sigmas) #
        _res    = np.asarray(<yaflFloat[:nz]> res)           #

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

        self._Us  = np.zeros((_U_sz(dim_z),), dtype=NP_DTYPE)
        self.v_Us = self._Us
        self.c_self.base.ukf_full.Us = &self.v_Us[0]

        self._Ds  = np.ones((dim_z,), dtype=NP_DTYPE)
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

#==============================================================================
cdef class JulierSigmaPoints(yaflSigmaBase):
    """
    Van der Merwe sigma point generator implementation
    """
    def __init__(self, yaflInt dim_x, yaflFloat kappa=0.0, **kwargs):
        super().__init__(dim_x, **kwargs)
        self.c_self.base.julier.kappa = kappa

    cdef yaflInt get_np(self, int dim_x):
        return (2 * dim_x + 1)

    cdef const yaflUKFSigmaMethodsSt * get_spm(self):
        return &yafl_ukf_julier_spm

    #==========================================================================
    #Decorators
    @property
    def kappa(self):
        return self.c_self.base.merwe.kappa

    @kappa.setter
    def kappa(self, value):
        raise AttributeError('JulierSigmaPoints does not support this!')

#==============================================================================
cdef class IMMEstimator:
    #
    cdef yaflIMMCBSt c_self

    #Memoryviews
    cdef yaflFloat [::1]    v_mu
    cdef yaflFloat [:, ::1] v_M
    cdef yaflFloat [::1]    v_Up
    cdef yaflFloat [::1]    v_Dp
    cdef yaflFloat [::1]    v_x
    cdef yaflFloat [::1]    v_cbar
    cdef yaflFloat [:, ::1] v_omega
    cdef yaflFloat [::1]    v_y
    cdef yaflFloat [:, ::1] v_W
    cdef yaflFloat [::1]    v_D

    cdef np.ndarray  _mu
    cdef np.ndarray  _M
    cdef np.ndarray  _Up
    cdef np.ndarray  _Dp
    cdef np.ndarray  _x
    cdef np.ndarray  _cbar
    cdef np.ndarray  _omega
    cdef np.ndarray  _y
    cdef np.ndarray  _W
    cdef np.ndarray  _D

    cdef np.ndarray       _z
    cdef yaflFloat [::1] v_z

    cdef yaflFloat _dt
    cdef list      _filters

    #==========================================================================
    def __cinit__(self, filters, mu, M, yaflFloat dt):
        assert isinstance(filters, list)
        assert len(filters) > 1

        cdef yaflPyKalmanBaseUn * base

        self.c_self.bank = <yaflFilterBankItemSt *>malloc(sizeof(yaflFilterBankItemSt) * len(filters))

        for i,f in enumerate(filters):
            if isinstance(f, Bierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_bierman_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, Joseph):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_joseph_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, AdaptiveBierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_adaptive_bierman_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, AdaptiveJoseph):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_adaptive_joseph_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, RobustBierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_robust_bierman_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, RobustJoseph):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_robust_joseph_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, AdaptiveRobustBierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_adaptive_robust_bierman_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, AdaptiveRobustJoseph):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ekf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ekf_adaptive_robust_joseph_update
                self.c_self.bank[i].Us      = base.ekf.W
                self.c_self.bank[i].Ds      = base.ekf.D
                self.c_self.bank[i].Xs      = base.ekf.H

            elif isinstance(f, UnscentedBierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ukf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ukf_bierman_update
                self.c_self.bank[i].Us      = base.ukf.sigmas_x
                self.c_self.bank[i].Ds      = base.ukf.Sx
                self.c_self.bank[i].Xs      = base.ukf.Pzx

            elif isinstance(f, UnscentedAdaptiveBierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ukf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ukf_adaptive_bierman_update
                self.c_self.bank[i].Us      = base.ukf.sigmas_x
                self.c_self.bank[i].Ds      = base.ukf.Sx
                self.c_self.bank[i].Xs      = base.ukf.Pzx

            elif isinstance(f, UnscentedRobustBierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ukf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ukf_robust_bierman_update
                self.c_self.bank[i].Us      = base.ukf.sigmas_x
                self.c_self.bank[i].Ds      = base.ukf.Sx
                self.c_self.bank[i].Xs      = base.ukf.Pzx

            elif isinstance(f, UnscentedAdaptiveRobustBierman):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ukf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ukf_adaptive_robust_bierman_update
                self.c_self.bank[i].Us      = base.ukf.sigmas_x
                self.c_self.bank[i].Ds      = base.ukf.Sx
                self.c_self.bank[i].Xs      = base.ukf.Pzx

            elif isinstance(f, Unscented):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ukf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ukf_update
                self.c_self.bank[i].Us      = base.ukf.sigmas_x
                self.c_self.bank[i].Ds      = base.ukf.Sx
                self.c_self.bank[i].Xs      = base.ukf.Pzx

            elif isinstance(f, UnscentedAdaptive):
                base = (<yaflKalmanBase>f).cbase()
                self.c_self.bank[i].predict = imm_yafl_ukf_base_predict
                self.c_self.bank[i].update  = imm_yafl_ukf_adaptive_update
                self.c_self.bank[i].Us      = base.ukf.sigmas_x
                self.c_self.bank[i].Ds      = base.ukf.Sx
                self.c_self.bank[i].Xs      = base.ukf.Pzx

            else:
                raise ValueError('List filters must contain only yaflpy filters')

            self.c_self.bank[i].filter = &(<yaflKalmanBase>f).cbase().base

        self.c_self.Nb = len(filters)

    #==========================================================================
    def __init__(self, filters, mu, M, yaflFloat dt):
        
        cdef yaflPyKalmanBaseUn * base = (<yaflKalmanBase>filters[0]).cbase()

        nx = base.base.Nx
        nz = base.base.Nz
        hx = (<yaflKalmanBase>filters[0])._hx

        for f in filters:
            base = (<yaflKalmanBase>f).cbase()
            assert base.base.Nx == nx
            assert base.base.Nz == nz
            assert (<yaflKalmanBase>f)._hx == hx

        self._filters = filters

        assert isinstance(mu, np.ndarray)
        assert NP_DTYPE == mu.dtype
        assert 1 == len(mu.shape)

        assert isinstance(M, np.ndarray)
        assert NP_DTYPE == M.dtype
        assert 2 == len(M.shape)

        assert mu.shape[0] == self.c_self.Nb
        assert M.shape[0]  == self.c_self.Nb
        assert M.shape[1]  == self.c_self.Nb

        self._mu  = mu.copy()
        self.v_mu = self._mu
        self.c_self.mu = &self.v_mu[0]

        self._M  = M.copy()
        self.v_M = self._M
        self.c_self.M = &self.v_M[0, 0]

        self._Up  = np.zeros_like(filters[0].Up)
        self.v_Up = self._Up
        self.c_self.Up = &self.v_Up[0]

        self._Dp  = np.zeros_like(filters[0].Dp)
        self.v_Dp = self._Dp
        self.c_self.Dp = &self.v_Dp[0]

        self._x  = np.zeros_like(filters[0].x)
        self.v_x = self._x
        self.c_self.x = &self.v_x[0]

        self._cbar  = np.zeros_like(mu)
        self.v_cbar = self._cbar
        self.c_self.cbar = &self.v_cbar[0]

        self._omega  = np.zeros_like(M)
        self.v_omega = self._omega
        self.c_self.omega = &self.v_omega[0, 0]

        self._y  = np.zeros_like(filters[0].x)
        self.v_y = self._y
        self.c_self.y = &self.v_y[0]

        self._W  = np.zeros((2 * nx, nx), dtype=NP_DTYPE)
        self.v_W = self._W
        self.c_self.W = &self.v_W[0, 0]

        self._D  = np.zeros((2 * nx,), dtype=NP_DTYPE)
        self.v_D = self._D
        self.c_self.D = &self.v_D[0]

        assert YAFL_ST_OK == yafl_imm_post_init(&self.c_self)
        
        self._z  = np.zeros((nz,), dtype=NP_DTYPE)
        self.v_z = self._z

        self._dt = dt
        

    #==========================================================================
    def __dealloc__(self):
        free(self.c_self.bank)

    #==========================================================================
    def predict(self, dt=None, **fx_args):
        old_dt = self._dt

        if dt:
            if np.isnan(dt):
                raise ValueError('Invalid dt value (nan)!')
            self._dt = <yaflFloat>dt

        for f in self._filters:
            f._fx_args = fx_args
            f._dt = self._dt

        res = yafl_imm_predict(&self.c_self)
        if res > YAFL_ST_ERR_THR:
            raise ValueError('Bad return value on IMMEstimator.predict!')

        self._dt = old_dt
        return res

    #==========================================================================
    def update(self, z, **hx_args):

        self._z[:] = z

        for f in self._filters:
            f._hx_args = hx_args

        res = yafl_imm_update(&self.c_self, &self.v_z[0])
        if res > YAFL_ST_ERR_THR:
            raise ValueError('Bad return value on IMMEstimator.update!')

        return res
    
    #==========================================================================
    #Decorators
    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, value):
        raise AttributeError('IMMEstimator does not support this!')
    #--------------------------------------------------------------------------
    @property
    def M(self):
        return self._M

    @M.setter
    def M(self, value):
        raise AttributeError('IMMEstimator does not support this!')
    #--------------------------------------------------------------------------
    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, value):
        raise AttributeError('IMMEstimator does not support this!')
    #--------------------------------------------------------------------------
    @property
    def cbar(self):
        return self._cbar

    @cbar.setter
    def cbar(self, value):
        raise AttributeError('IMMEstimator does not support this!')
    #--------------------------------------------------------------------------
    @property
    def W(self):
        return self._W

    @W.setter
    def W(self, value):
        raise AttributeError('IMMEstimator does not support this!')
    #--------------------------------------------------------------------------
    @property
    def D(self):
        return self._D

    @D.setter
    def D(self, value):
        raise AttributeError('IMMEstimator does not support this!')
    #--------------------------------------------------------------------------
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x[:] = value
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
    def P(self):
        #
        if len(self._x) > 1:
            _,u = _set_u(self._Up)
            return u.dot(np.diag(self._Dp).dot(u.T))
        #
        return self._Dp.copy().reshape((1,1))

    @P.setter
    def P(self, value):
        raise AttributeError('IMMEstimator does not support this!')
