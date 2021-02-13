/*******************************************************************************
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
******************************************************************************/

#ifndef YAFL_H
#define YAFL_H

#include <yafl_config.h>
#include "yafl_math.h"

/*=============================================================================
               Basic UD-factorized Kalman filter Definitions
=============================================================================*/
typedef struct _yaflKalmanBaseSt yaflKalmanBaseSt;

typedef yaflStatusEn (* yaflKalmanFuncP)(yaflKalmanBaseSt *, yaflFloat *,    \
                                         yaflFloat *);

typedef yaflStatusEn (* yaflKalmanResFuncP)(yaflKalmanBaseSt *, yaflFloat *, \
                                            yaflFloat *, yaflFloat *);

typedef yaflStatusEn (* yaflKalmanScalarUpdateP)(yaflKalmanBaseSt *, yaflInt);

/*
Function pointer for robust versions of filters.

Based on:
1. West M., "Robust Sequential Approximate Bayesian Estimation",
   J. R. Statist. Soc. B (1981), 43, No. 2, pp. 157-166

2. Gaver, Donald Paul; Jacobs, Patricia A., "Robustifying the Kalman filter",
   Naval Postgraduate School technical report. 1987
   http://hdl.handle.net/10945/30147
*/
typedef yaflFloat (* yaflKalmanRobFuncP)(yaflKalmanBaseSt *, yaflFloat);

struct _yaflKalmanBaseSt {
    yaflKalmanFuncP f;       /*A state transition function*/
    yaflKalmanFuncP h;       /*A measurement function*/
    yaflKalmanResFuncP zrf;  /*Measurement residual function function*/

    yaflFloat * x;  /*State vector*/
    yaflFloat * y;  /*Innovation vector*/

    yaflFloat * Up; /*Upper triangular part of P*/
    yaflFloat * Dp; /*Diagonal part of P*/

    yaflFloat * Uq; /*Upper triangular part of Q*/
    yaflFloat * Dq; /*Diagonal part of Q*/

    yaflFloat * Ur; /*Upper triangular part of R*/
    yaflFloat * Dr; /*Diagonal part of R*/

    yaflInt   Nx;   /*State vector size*/
    yaflInt   Nz;   /*Measurement vector size*/
};

/*---------------------------------------------------------------------------*/
#define YAFL_KALMAN_BASE_MEMORY_MIXIN(nx, nz) \
    yaflFloat x[nx];                          \
    yaflFloat y[nz];                          \
                                              \
    yaflFloat Up[((nx - 1) * nx)/2];          \
    yaflFloat Dp[nx];                         \
                                              \
    yaflFloat Uq[((nx - 1) * nx)/2];          \
    yaflFloat Dq[nx];                         \
                                              \
    yaflFloat Ur[((nz - 1) * nz)/2];          \
    yaflFloat Dr[nz]

/*---------------------------------------------------------------------------*/
#define YAFL_KALMAN_BASE_INITIALIZER(_f, _h, _zrf, _nx, _nz, _mem) \
{                                                                  \
    .f   = (yaflKalmanFuncP)_f,                                    \
    .h   = (yaflKalmanFuncP)_h,                                    \
    .zrf = (yaflKalmanResFuncP)_zrf,                               \
                                                                   \
    .x   = _mem.x,                                                 \
    .y   = _mem.y,                                                 \
                                                                   \
    .Up  = _mem.Up,                                                \
    .Dp  = _mem.Dp,                                                \
                                                                   \
    .Uq  = _mem.Uq,                                                \
    .Dq  = _mem.Dq,                                                \
                                                                   \
    .Ur  = _mem.Ur,                                                \
    .Dr  = _mem.Dr,                                                \
                                                                   \
    .Nx  = _nx,                                                    \
    .Nz  = _nz                                                     \
}

/*=============================================================================
                    Basic UD-factorized EKF definitions
=============================================================================*/
typedef struct _yaflEKFBaseSt yaflEKFBaseSt;

struct _yaflEKFBaseSt {
    yaflKalmanBaseSt base; /*Base type*/

    yaflKalmanFuncP jf; /*Jacobian of a state transition function*/
    yaflKalmanFuncP jh; /*Jacobian of a measurement function*/

    yaflFloat * H;   /*Measurement Jacobian values*/
    yaflFloat * W;   /*Scratchpad memory block matrix*/
    yaflFloat * D;   /*Scratchpad memory diagonal matrix*/
};

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_BASE_MEMORY_MIXIN(nx, nz) \
    YAFL_KALMAN_BASE_MEMORY_MIXIN(nx, nz); \
                                           \
    yaflFloat H[nz * nx];                  \
    yaflFloat W[2 * nx * nx];              \
    yaflFloat D[2 * nx]

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem) \
{                                                                         \
    .base = YAFL_KALMAN_BASE_INITIALIZER(_f, _h, _zrf, _nx, _nz, _mem),   \
                                                                          \
    .jf  = (yaflKalmanFuncP)_jf,                                          \
    .jh  = (yaflKalmanFuncP)_jh,                                          \
                                                                          \
    .H   = _mem.H,                                                        \
    .W   = _mem.W,                                                        \
    .D   = _mem.D                                                         \
}

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ekf_base_predict(yaflKalmanBaseSt * self);

yaflStatusEn yafl_ekf_base_update(yaflKalmanBaseSt * self, yaflFloat * z, \
                                  yaflKalmanScalarUpdateP scalar_update);

/*---------------------------------------------------------------------------*/
#define YAFL_KALMAN_PREDICT_WRAPPER(sname, stype) \
static inline yaflStatusEn _yafl_##sname##_predict_wrapper(stype * self) \
{                                                                        \
    return yafl_ekf_base_predict((yaflKalmanBaseSt *)self);              \
}

/*---------------------------------------------------------------------------*/
#define YAFL_KALMAN_UPDATE_IMPL(bname, sname, stype)                                           \
yaflStatusEn yafl_##bname##_##sname##_scalar_update(yaflKalmanBaseSt * self, yaflInt i);       \
static inline yaflStatusEn yafl_##bname##_##sname##_update(stype * self, yaflFloat * z)        \
{                                                                                              \
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, yafl_##bname##_##sname##_scalar_update); \
}

/*---------------------------------------------------------------------------*/
//static inline yaflStatusEn _yafl_ekf_predict_wrapper(yaflEKFBaseSt * self)
//{
//    return yafl_ekf_base_predict((yaflKalmanBaseSt *)self);
//}
YAFL_KALMAN_PREDICT_WRAPPER(ekf, yaflEKFBaseSt)
/*-----------------------------------------------------------------------------
                               Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_BIERMAN_PREDICT _yafl_ekf_predict_wrapper

//yaflStatusEn yafl_ekf_bierman_scalar_update(yaflKalmanBaseSt * self, yaflInt i);
//
//static inline yaflStatusEn yafl_ekf_bierman_update(yaflEKFBaseSt * self, \
//                                                   yaflFloat * z)
//{
//    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
//                                yafl_ekf_bierman_scalar_update);
//}
YAFL_KALMAN_UPDATE_IMPL(ekf, bierman, yaflEKFBaseSt)
/*-----------------------------------------------------------------------------
                               Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_JOSEPH_PREDICT _yafl_ekf_predict_wrapper

yaflStatusEn yafl_ekf_joseph_scalar_update(yaflKalmanBaseSt * self, yaflInt i);

static inline yaflStatusEn yafl_ekf_joseph_update(yaflEKFBaseSt * self, \
                                                  yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
                                yafl_ekf_joseph_scalar_update);
}

/*=============================================================================
                    Adaptive UD-factorized EKF definitions
=============================================================================*/
typedef struct {
    yaflEKFBaseSt base;
    yaflFloat chi2; /*Divergence test threshold (chi-squared criteria)*/
} yaflEKFAdaptiveSt; /*Adaptive kKalman-Hinfinity filter structure*/

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem)  \
{                                                                              \
    .base = YAFL_EKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem), \
    .chi2 = 10.8275662                                                         \
}
/*
Default value for chi2 is:
  scipy.stats.chi2.ppf(0.999, 1)
*/

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _yafl_ada_ekf_predict_wrapper(yaflEKFAdaptiveSt * self)
{
    return yafl_ekf_base_predict((yaflKalmanBaseSt *)self);
}

/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_BIERAMN_PREDICT _yafl_ada_ekf_predict_wrapper

yaflStatusEn yafl_ekf_adaptive_bierman_scalar_update(yaflKalmanBaseSt * self, \
                                                     yaflInt i);

static inline yaflStatusEn \
    yafl_ekf_adaptive_bierman_update(yaflEKFAdaptiveSt * self, yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
                                yafl_ekf_adaptive_bierman_scalar_update);
}

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_JOSEPH_PREDICT _yafl_ada_ekf_predict_wrapper

yaflStatusEn yafl_ekf_adaptive_joseph_scalar_update(yaflKalmanBaseSt * self, \
                                                    yaflInt i);

static inline yaflStatusEn \
    yafl_ekf_adaptive_joseph_update(yaflEKFAdaptiveSt * self, yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
                                yafl_ekf_adaptive_joseph_scalar_update);
}

/*-----------------------------------------------------------------------------
                                 WARNING!!!

             DO NOT USE THIS variant of Adaptive Joseph filter !!!

     It was implemented to show some flaws of the corresponding algorithm!
-----------------------------------------------------------------------------*/
yaflStatusEn \
    yafl_ekf_do_not_use_this_scalar_update(yaflKalmanBaseSt * self, yaflInt i);


static inline yaflStatusEn \
    yafl_ekf_do_not_use_this_update(yaflEKFAdaptiveSt * self, yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
                                yafl_ekf_do_not_use_this_scalar_update);
}

/*=============================================================================
                    Robust UD-factorized EKF definitions
=============================================================================*/
/*
Based on:
1. West M., "Robust Sequential Approximate Bayesian Estimation",
   J. R. Statist. Soc. B (1981), 43, No. 2, pp. 157-166

2. Gaver, Donald Paul; Jacobs, Patricia A., "Robustifying the Kalman filter",
   Naval Postgraduate School technical report. 1987
   http://hdl.handle.net/10945/30147
*/
typedef struct {
    yaflEKFBaseSt base;
    yaflKalmanRobFuncP g;    /* g = -d(ln(pdf(y))) / dy */
    yaflKalmanRobFuncP gdot; /* gdot = G = d(g) / dy */
} yaflEKFRobustSt; /*Robust EKF*/

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf, _g, _gdot,  \
                                    _nx, _nz, _mem)                     \
{                                                                       \
    .base = YAFL_EKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, \
                                      _mem),                            \
    .g    = _g,                                                         \
    .gdot = _gdot                                                       \
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _yafl_rob_ekf_predict_wrapper(yaflEKFRobustSt * self)
{
    return yafl_ekf_base_predict((yaflKalmanBaseSt *)self);
}

/*-----------------------------------------------------------------------------
                           Robust Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ROBUST_BIERAMN_PREDICT _yafl_rob_ekf_predict_wrapper

yaflStatusEn \
    yafl_ekf_robust_bierman_scalar_update(yaflKalmanBaseSt * self, yaflInt i);

static inline yaflStatusEn \
    yafl_ekf_robust_bierman_update(yaflEKFRobustSt * self, yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z,
                                yafl_ekf_robust_bierman_scalar_update);
}

/*-----------------------------------------------------------------------------
                           Robust Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ROBUST_JOSEPH_PREDICT _yafl_rob_ekf_predict_wrapper

yaflStatusEn \
    yafl_ekf_robust_joseph_scalar_update(yaflKalmanBaseSt * self, yaflInt i);

static inline yaflStatusEn \
    yafl_ekf_robust_joseph_update(yaflEKFRobustSt * self, yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
                                yafl_ekf_robust_joseph_scalar_update);
}


/*=============================================================================
                 Adaptive robust UD-factorized EKF definitions
=============================================================================*/
typedef struct {
    yaflEKFRobustSt base;
    yaflFloat chi2; /*Divergence test threshold (chi-squared criteria)*/
} yaflEKFAdaptiveRobustSt; /*Robust EKF*/

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf,    \
                                             _g, _gdot, _nx, _nz, _mem) \
{                                                                       \
    .base = YAFL_EKF_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf,         \
                                        _g, _gdot, _nx, _nz, _mem),     \
    .chi2 = 8.8074684                                                   \
}
/*
Default value for chi2 is:
  scipy.stats.chi2.ppf(0.997, 1)
*/

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn \
    _yafl_ada_rob_predict_wrapper(yaflEKFAdaptiveRobustSt * self)
{
    return yafl_ekf_base_predict((yaflKalmanBaseSt *)self);
}

/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_ROBUST_BIERAMN_PREDICT _yafl_ada_rob_predict_wrapper

yaflStatusEn \
    yafl_ekf_adaptive_robust_bierman_scalar_update(yaflKalmanBaseSt * self, \
                                               yaflInt i);

static inline yaflStatusEn \
    yafl_ekf_adaptive_robust_bierman_update(yaflEKFAdaptiveRobustSt * self, \
                                            yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
                                yafl_ekf_adaptive_robust_bierman_scalar_update);
}

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_ROBUST_JOSEPH_PREDICT(self) _yafl_ada_rob_predict_wrapper

yaflStatusEn \
    yafl_ekf_adaptive_robust_joseph_scalar_update(yaflKalmanBaseSt * self, \
                                              yaflInt i);

static inline yaflStatusEn \
    yafl_ekf_adaptive_robust_joseph_update(yaflEKFAdaptiveRobustSt * self, \
                                            yaflFloat * z)
{
    return yafl_ekf_base_update((yaflKalmanBaseSt *)self, z, \
                                yafl_ekf_adaptive_robust_joseph_scalar_update);
}

/*=============================================================================
                    Basic UD-factorized UKF definitions
=============================================================================*/
typedef struct _yaflUKFBaseSt yaflUKFBaseSt;        /* The UKF base type */
/*
Used to add delta vectors to initial point in sigma point generation.
Parameters:
yaflFloat * x0
yaflFloat * delta_x
yaflFloat factor

Does:
delta_x = x0 + factor * delta_x
*/
typedef yaflStatusEn (* yaflUKFSigmaAddP)(yaflUKFBaseSt *, yaflFloat *, yaflFloat *, yaflFloat);

/*Sigma point generator info base type*/
typedef struct _yaflUKFSigmaSt {
    yaflInt         np; /* The number of sigma points          */
    yaflUKFSigmaAddP addf; /* Sigma point addition function       */
} yaflUKFSigmaSt;

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_SIGMA_BASE_INITIALIZER(_np, _addf, _mem) \
{                                                         \
    .np   = (yaflInt)_np,                                 \
    .addf = (yaflUKFSigmaAddP)_addf,                      \
}

/*---------------------------------------------------------------------------*/
/* Computes sigma points weights */
typedef yaflStatusEn (* yafkUKFSigmaGenWeigthsP)(yaflUKFBaseSt *);

/* Generates sigma points */
typedef yaflStatusEn (* yaflUKFSigmaGenSigmasP)(yaflUKFBaseSt *);

typedef struct _yaflUKFSigmaMethodsSt {
    yafkUKFSigmaGenWeigthsP   wf; /* Weight function                */
    yaflUKFSigmaGenSigmasP  spgf; /* Sigma point generator function */
} yaflUKFSigmaMethodsSt;

/*---------------------------------------------------------------------------*/
struct _yaflUKFBaseSt {

    yaflKalmanBaseSt base;

    /* A pointer to the sigma point generator structure */
    yaflUKFSigmaSt               * sp_info;
    /* A sigma point generator method table pointer     */
    const yaflUKFSigmaMethodsSt  * sp_meth;

    yaflKalmanFuncP    xmf; /* State mean function              */
    yaflKalmanResFuncP xrf; /* State residual function function */

    yaflKalmanFuncP    zmf; /* Measurement mean function function    */
    yaflFloat * zp; /* Predicted measurement vector */

    /*Scratchpad memory*/
    yaflFloat * Sx;  /* State       */
    yaflFloat * Pzx; /* Pzx cross covariance matrix */

    yaflFloat * sigmas_x; /* State sigma points       */
    yaflFloat * sigmas_z; /* Measurement sigma points */
    yaflFloat * wm;       /* Weights for mean calculations       */
    yaflFloat * wc;       /* Weights for covariance calculations */
};

/*---------------------------------------------------------------------------*/
/*
Warning: sigmas_x and _sigmas_z aren't defined in this mixin, see
         sigma points generators mixins!!!
*/
#define YAFL_UKF_BASE_MEMORY_MIXIN(nx, nz) \
    YAFL_KALMAN_BASE_MEMORY_MIXIN(nx, nz);  \
    yaflFloat zp[nz];                       \
                                            \
    yaflFloat Pzx[nz * nx];                 \
                                            \
    yaflFloat Sx[nx]

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_BASE_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf,          \
                                   _zrf, _nx, _nz, _mem)                      \
{                                                                             \
    .base = YAFL_KALMAN_BASE_INITIALIZER(_f, _h, _zrf, _nx, _nz, _mem),       \
                                                                              \
    .sp_info = _p,                                                            \
    .sp_meth = _pm,                                                           \
                                                                              \
    .xmf = (yaflUKFFuncP)_xmf,                                                \
    .xrf = (yaflUKFResFuncP)_xrf,                                             \
                                                                              \
    .zmf = (yaflUKFFuncP)_zmf,                                                \
                                                                              \
    .zp  = _mem.zp,                                                           \
                                                                              \
    .Sx  = _mem.Sx,                                                           \
    .Pzx = _mem.Pzx,                                                          \
                                                                              \
    .sigmas_x  = _mem.sigmas_x,                                               \
    .sigmas_z  = _mem.sigmas_z,                                               \
                                                                              \
    .wm   = _mem.wm,                                                          \
    .wc   = _mem.wc                                                           \
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn yafl_ukf_post_init(yaflUKFBaseSt * self)
{
    YAFL_CHECK(self,              YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->sp_meth,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->sp_meth->wf, YAFL_ST_INV_ARG_1);
    return self->sp_meth->wf(self); /*Need to compute weights before start*/
}

static inline yaflStatusEn yafl_ukf_gen_sigmas(yaflUKFBaseSt * self)
{
    YAFL_CHECK(self,                YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->sp_meth,       YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->sp_meth->spgf, YAFL_ST_INV_ARG_1);
    self->sp_meth->spgf(self);
    return self->sp_meth->wf(self);
}

yaflStatusEn yafl_ukf_base_predict(yaflUKFBaseSt * self);

yaflStatusEn yafl_ukf_base_update(yaflUKFBaseSt * self, yaflFloat * z, \
                                  yaflKalmanScalarUpdateP scalar_update);

/*=============================================================================
                               Bierman UKF
=============================================================================*/
#define YAFL_UKF_BIERMAN_PREDICT(self) \
    yafl_ukf_base_predict((yaflUKFBaseSt *)self)

yaflStatusEn yafl_ukf_bierman_update_scalar(yaflKalmanBaseSt * self, yaflInt i);

static inline yaflStatusEn \
    yafl_ukf_bierman_update(yaflUKFBaseSt * self, yaflFloat * z)
{
    return yafl_ukf_base_update(self, z, yafl_ukf_bierman_update_scalar);
}

/*=============================================================================
                          Adaptive Bierman UKF
=============================================================================*/
typedef struct {
    yaflUKFBaseSt base;
    yaflFloat chi2;
} yaflUKFAdaptivedSt;

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_ADAPTIVE_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h,       \
                                            _zmf, _zrf, _nx, _nz, _mem)  \
{                                                                        \
    .base = YAFL_UKF_BASE_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf, \
                                       _zrf, _nx, _nz, _mem) ,           \
    .chi2 = 10.8275662                                                   \
}

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_ADAPTIVE_BIERMAN_PREDICT(self) \
    yafl_ukf_base_predict((yaflUKFBaseSt *)self)

yaflStatusEn yafl_ukf_adaptive_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                     yaflInt i);

static inline yaflStatusEn \
    yafl_ukf_adaptive_bierman_update(yaflUKFAdaptivedSt * self, yaflFloat * z)
{
    return yafl_ukf_base_update((yaflUKFBaseSt *)self, z, \
                                yafl_ukf_adaptive_bierman_update_scalar);
}

/*=============================================================================
                           Robust Bierman UKF
=============================================================================*/
typedef yaflFloat (* yaflUKFRobFuncP)(yaflUKFBaseSt *, yaflFloat);

typedef struct {
    yaflUKFBaseSt base;
    yaflKalmanRobFuncP g;    /* g = -d(ln(pdf(y))) / dy */
    yaflKalmanRobFuncP gdot; /* gdot = G = d(g) / dy */
} yaflUKFRobustSt; /*Robust UKF*/

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_ROBUST_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf, _zrf, \
                                    _g, _gdot, _nx, _nz, _mem)               \
{                                                                            \
    .base = YAFL_UKF_BASE_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf,     \
                                       _zrf, _nx, _nz, _mem) ,               \
    .g    = _g,                                                              \
    .gdot = _gdot                                                            \
}

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ukf_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                     yaflInt i);

static inline yaflStatusEn \
    yafl_ukf_robust_bierman_update(yaflUKFAdaptivedSt * self, yaflFloat * z)
{
    return yafl_ukf_base_update((yaflUKFBaseSt *)self, z, \
                                yafl_ukf_robust_bierman_update_scalar);
}

/*=============================================================================
                        Adaptive robust Bierman UKF
=============================================================================*/
typedef struct {
    yaflUKFRobustSt base;
    yaflFloat chi2; /*Divergence test threshold (chi-squared criteria)*/
} yaflUKFAdaptiveRobustSt; /*Robust EKF*/

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_ADAPTIVE_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf,    \
                                             _g, _gdot, _nx, _nz, _mem) \
{                                                                       \
    .base = YAFL_UKF_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf,         \
                                        _g, _gdot, _nx, _nz, _mem),     \
    .chi2 = 8.8074684                                                   \
}

/*---------------------------------------------------------------------------*/
yaflStatusEn \
    yafl_ukf_adaptive_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                   yaflInt i);

static inline yaflStatusEn \
    yafl_ukf_adaptive_robust_bierman_update(yaflUKFAdaptivedSt * self, \
                                            yaflFloat * z)
{
    return yafl_ukf_base_update((yaflUKFBaseSt *)self, z, \
                                yafl_ukf_adaptive_robust_bierman_update_scalar);
}

/*=============================================================================
            Full UKF, not sequential square root version of UKF
=============================================================================*/
typedef struct {
    yaflUKFBaseSt base; /* Base type                   */
    yaflFloat *   Us;   /* Upper triangular part of S  */
    yaflFloat *   Ds;   /* Diagonal part of S          */
} yaflUKFSt;

/*---------------------------------------------------------------------------*/
/*
Warning: sigmas_x and _sigmas_z aren't defined in this mixin, see
         sigma points generators mixins!!!
*/
#define YAFL_UKF_MEMORY_MIXIN(nx, nz)   \
    YAFL_UKF_BASE_MEMORY_MIXIN(nx, nz); \
    yaflFloat Us[((nz - 1) * nz)/2];    \
    yaflFloat Ds[nz]

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf,               \
                                   _zrf, _nx, _nz, _mem)                      \
{                                                                             \
    .base = YAFL_UKF_BASE_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf,      \
                                      _zrf, _nx, _nz, _mem),                  \
    .Us  = _mem.Us,                                                           \
    .Ds  = _mem.Ds                                                            \
}

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_PREDICT(self) \
    yafl_ukf_base_predict((yaflUKFBaseSt *)self)

yaflStatusEn yafl_ukf_update(yaflUKFBaseSt * self, yaflFloat * z);

/*=============================================================================
                     Van der Merwe sigma point generator
=============================================================================*/
typedef struct _yaflUKFMerweSt {
    yaflUKFSigmaSt base;
    yaflFloat alpha;
    yaflFloat beta;
    yaflFloat kappa;
} yaflUKFMerweSt;

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_MERWE_MEMORY_MIXIN(nx, nz) \
    yaflFloat wm[2 * nx + 1];               \
    yaflFloat wc[2 * nx + 1];               \
    yaflFloat sigmas_x[(2 * nx + 1) * nx];  \
    yaflFloat sigmas_z[(2 * nx + 1) * nz]

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_MERWE_INITIALIZER(_nx, _addf, _alpha, _beta, _kappa, _mem) \
{                                                                           \
    .base  = YAFL_UKF_SIGMA_BASE_INITIALIZER((2 * _nx + 1), _addf, _mem),   \
    .alpha = _alpha,                                                        \
    .beta  = _beta,                                                         \
    .kappa = _kappa                                                         \
}

/*---------------------------------------------------------------------------*/
const yaflUKFSigmaMethodsSt yafl_ukf_merwe_spm;

#endif // YAFL_H
