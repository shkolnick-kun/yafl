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
                    Basic UD-factorized EKF definitions
=============================================================================*/
typedef struct _yaflEKFBaseSt yaflEKFBaseSt;

typedef void (* yaflEKFFuncP)(yaflEKFBaseSt *);
typedef void (* yaflEKFResFuncP)(yaflEKFBaseSt *, yaflFloat *);
typedef yaflStatusEn (* yaflEKFScalarUpdateP)(yaflEKFBaseSt *, yaflInt);

struct _yaflEKFBaseSt {
    yaflEKFFuncP f;    /*A state transition function*/
    yaflEKFFuncP jf;   /*Jacobian of a state transition function*/

    yaflEKFFuncP h;    /*A measurement function*/
    yaflEKFFuncP jh;   /*Jacobian of a measurement function*/

    yaflEKFResFuncP zrf;  /*Measurement residual function function*/

    yaflFloat * x;  /*State vector*/
    yaflFloat * y;  /*Innovation vector*/
    yaflFloat * H;  /*Measurement Jacobian values*/

    yaflFloat * Up; /*Upper triangular part of P*/
    yaflFloat * Dp; /*Diagonal part of P*/

    yaflFloat * Uq; /*Upper triangular part of Q*/
    yaflFloat * Dq; /*Diagonal part of Q*/

    yaflFloat * Ur; /*Upper triangular part of R*/
    yaflFloat * Dr; /*Diagonal part of R*/

    yaflFloat * W;  /*Scratchpad memory block matrix*/
    yaflFloat * D;  /*Scratchpad memory diagonal matrix*/

    yaflInt   Nx;   /*State vector size*/
    yaflInt   Nz;   /*Measurement vector size*/
};

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_BASE_MEMORY_MIXIN(nx, nz) \
    yaflFloat x[nx];                   \
    yaflFloat y[nz];                   \
    yaflFloat H[nz * nx];              \
                                       \
    yaflFloat Up[((nx - 1) * nx)/2];   \
    yaflFloat Dp[nx];                  \
                                       \
    yaflFloat Uq[((nx - 1) * nx)/2];   \
    yaflFloat Dq[nx];                  \
                                       \
    yaflFloat Ur[((nz - 1) * nz)/2];   \
    yaflFloat Dr[nz];                  \
                                       \
    yaflFloat W[2 * nx * nx];          \
    yaflFloat D[2 * nx];

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem)\
{                                                                    \
    .f   = (yaflEKFFuncP)_f,                                         \
    .jf  = (yaflEKFFuncP)_jf,                                        \
                                                                     \
    .h   = (yaflEKFFuncP)_h,                                         \
    .jh  = (yaflEKFFuncP)_jh,                                        \
                                                                     \
    .zrf = (yaflEKFResFuncP)_zrf,                                    \
                                                                     \
    .x   = _mem.x,                                                   \
    .y   = _mem.y,                                                   \
    .H   = _mem.H,                                                   \
                                                                     \
    .Up  = _mem.Up,                                                  \
    .Dp  = _mem.Dp,                                                  \
                                                                     \
    .Uq  = _mem.Uq,                                                  \
    .Dq  = _mem.Dq,                                                  \
                                                                     \
    .Ur  = _mem.Ur,                                                  \
    .Dr  = _mem.Dr,                                                  \
                                                                     \
    .W   = _mem.W,                                                   \
    .D   = _mem.D,                                                   \
                                                                     \
    .Nx  = _nx,                                                      \
    .Nz  = _nz                                                       \
}

/*---------------------------------------------------------------------------*/
void yafl_ekf_base_predict(yaflEKFBaseSt * self);
void yafl_ekf_base_update(yaflEKFBaseSt * self, yaflFloat * z, yaflEKFScalarUpdateP scalar_update);

/*-----------------------------------------------------------------------------
                               Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_BIERMAN_PREDICT yafl_ekf_base_predict
void yafl_ekf_bierman_update(yaflEKFBaseSt * self, yaflFloat * z);

/*-----------------------------------------------------------------------------
                               Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_JOSEPH_PREDICT yafl_ekf_base_predict
void yafl_ekf_joseph_update(yaflEKFBaseSt * self, yaflFloat * z);

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

/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_BIERAMN_PREDICT(self) yafl_ekf_base_predict((yaflEKFBaseSt *)self)
void yafl_ekf_adaptive_bierman_update(yaflEKFAdaptiveSt * self, yaflFloat * z);

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_JOSEPH_PREDICT(self) yafl_ekf_base_predict((yaflEKFBaseSt *)self)
void yafl_ekf_adaptive_joseph_update(yaflEKFAdaptiveSt * self, yaflFloat * z);

/*-----------------------------------------------------------------------------
                                 WARNING!!!

             DO NOT USE THIS variant of Adaptive Joseph filter !!!

     It was implemented to show some flaws of the corresponding algorithm!
-----------------------------------------------------------------------------*/
void yafl_ekf_do_not_use_this_update(yaflEKFAdaptiveSt * self, yaflFloat * z);

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
typedef yaflFloat (* yaflEKFRobFuncP)(yaflEKFBaseSt *, yaflFloat);

typedef struct {
    yaflEKFBaseSt base;
    yaflEKFRobFuncP g;    /* g = -d(ln(pdf(y))) / dy */
    yaflEKFRobFuncP gdot; /* gdot = G = d(g) / dy */
} yaflEKFRobustSt; /*Robust EKF*/

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem)    \
{                                                                              \
    .base = YAFL_EKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem), \
    .g    = (yaflEKFRobFuncP)0,                                                \
    .gdot = (yaflEKFRobFuncP)0                                                 \
}

/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ROBUST_BIERAMN_PREDICT(self) yafl_ekf_base_predict((yaflEKFBaseSt *)self)
void yafl_ekf_robust_bierman_update(yaflEKFRobustSt * self, yaflFloat * z);

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ROBUST_JOSEPH_PREDICT(self) yafl_ekf_base_predict((yaflEKFBaseSt *)self)
void yafl_ekf_robust_joseph_update(yaflEKFRobustSt * self, yaflFloat * z);

/*=============================================================================
                 Adaptive robust UD-factorized EKF definitions
=============================================================================*/
typedef struct {
    yaflEKFRobustSt base;
    yaflFloat chi2; /*Divergence test threshold (chi-squared criteria)*/
} yaflEKFAdaptiveRobustSt; /*Robust EKF*/

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem) \
{                                                                                    \
    .base = YAFL_EKF_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem),     \
    .chi2 = 8.8074684                                                                \
}
/*
Default value for chi2 is:
  scipy.stats.chi2.ppf(0.997, 1)
*/
/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_ROBUST_BIERAMN_PREDICT(self) \
    yafl_ekf_base_predict((yaflEKFBaseSt *)self)

void yafl_ekf_adaptive_robust_bierman_update(yaflEKFAdaptiveRobustSt * self, \
                                         yaflFloat * z);

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAFL_EKF_ADAPTIVE_ROBUST_JOSEPH_PREDICT(self) \
    yafl_ekf_base_predict((yaflEKFBaseSt *)self)

void yafl_ekf_adaptive_robust_joseph_update(yaflEKFAdaptiveRobustSt * self, \
                                        yaflFloat * z);

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
typedef void (* yaflUKFSigmaAddP)(yaflUKFBaseSt *, yaflFloat *, yaflFloat *, yaflFloat);

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
typedef void (* yafkUKFSigmaGenWeigthsP)(yaflUKFBaseSt *);

/* Generates sigma points */
typedef void (* yaflUKFSigmaGenSigmasP)(yaflUKFBaseSt *);

typedef struct _yaflUKFSigmaMethodsSt {
    yafkUKFSigmaGenWeigthsP   wf; /* Weight function                */
    yaflUKFSigmaGenSigmasP  spgf; /* Sigma point generator function */
} yaflUKFSigmaMethodsSt;

/*---------------------------------------------------------------------------*/
typedef yaflStatusEn (* yaflUKFScalarUpdateP)(yaflUKFBaseSt *, yaflInt);

typedef void (* yaflUKFFuncP)(yaflUKFBaseSt *, yaflFloat *, yaflFloat *);

typedef void (* yaflUKFResFuncP)(yaflUKFBaseSt *, yaflFloat *, yaflFloat *, \
                                 yaflFloat *);

struct _yaflUKFBaseSt {
    /* A pointer to the sigma point generator structure */
    yaflUKFSigmaSt               * sp_info;
    /* A sigma point generator method table pointer     */
    const yaflUKFSigmaMethodsSt  * sp_meth;

    yaflUKFFuncP      f; /* A state transition function      */
    yaflUKFFuncP    xmf; /* State mean function              */
    yaflUKFResFuncP xrf; /* State residual function function */

    yaflUKFFuncP      h; /* A measurement function                */
    yaflUKFFuncP    zmf; /* Measurement mean function function    */
    yaflUKFResFuncP zrf; /* Measurement residual function function*/

    yaflFloat * x;  /* State vector                 */
    yaflFloat * zp; /* Predicted measurement vector */
    yaflFloat * y;  /* Innovation */

    yaflFloat * Up;  /* Upper triangular part of P  */
    yaflFloat * Dp;  /* Diagonal part of P          */

    yaflFloat * Us;  /* Upper triangular part of S  */
    yaflFloat * Ds;  /* Diagonal part of S          */

    yaflFloat * Pzx; /* Pzx cross covariance matrix */

    yaflFloat * Uq;  /* Upper triangular part of Q   */
    yaflFloat * Dq;  /* Diagonal part of Q           */

    yaflFloat * Ur;  /* Upper triangular part of R   */
    yaflFloat * Dr;  /* Diagonal part of R           */

    yaflFloat * sigmas_x; /* State sigma points       */
    yaflFloat * sigmas_z; /* Measurement sigma points */
    yaflFloat * wm;       /* Weights for mean calculations       */
    yaflFloat * wc;       /* Weights for covariance calculations */

    /*Scratchpad memory*/
    yaflFloat * Sx;       /* State       */


    yaflInt   Nx; /* State vector size       */
    yaflInt   Nz; /* Measurement vector size */
};

/*---------------------------------------------------------------------------*/
/*
Warning: sigmas_x and _sigmas_z aren't defined in this mixin, see
         sigma points generators mixins!!!
*/
#define YAFL_UKF_MEMORY_MIXIN(nx, nz) \
    yaflFloat x[nx];                        \
    yaflFloat zp[nz];                       \
    yaflFloat y[nz];                        \
                                            \
    yaflFloat Up[((nx - 1) * nx)/2];        \
    yaflFloat Dp[nx];                       \
                                            \
    yaflFloat Us[((nz - 1) * nz)/2];        \
    yaflFloat Ds[nz];                       \
                                            \
    yaflFloat Pzx[nz * nx];                 \
                                            \
    yaflFloat Uq[((nx - 1) * nx)/2];        \
    yaflFloat Dq[nx];                       \
                                            \
    yaflFloat Ur[((nz - 1) * nz)/2];        \
    yaflFloat Dr[nz];                       \
                                            \
    yaflFloat Sx[nx];

/*---------------------------------------------------------------------------*/
#define YAFL_UKF_BASE_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf,          \
                                   _zrf, _nx, _nz, _mem)                      \
{                                                                             \
    .sp_info = _p,                                                            \
    .sp_meth = _pm,                                                           \
                                                                              \
    .f   = (yaflUKFFuncP)_f,                                                  \
    .xmf = (yaflUKFFuncP)_xmf,                                                \
    .xrf = (yaflUKFResFuncP)_xrf,                                             \
                                                                              \
    .h   = (yaflUKFFuncP)_h,                                                  \
    .xmf = (yaflUKFFuncP)_xmf,                                                \
    .xrf = (yaflUKFResFuncP)_xrf,                                             \
                                                                              \
    .x   = _mem.x,                                                            \
    .zp  = _mem.zp,                                                           \
    .y   = _mem.y,                                                            \
                                                                              \
    .Up  = _mem.Up,                                                           \
    .Dp  = _mem.Dp,                                                           \
                                                                              \
    .Us  = _mem.Us,                                                           \
    .Ds  = _mem.Ds,                                                           \
                                                                              \
    .Pzx = _mem.Pzx,                                                          \
                                                                              \
    .Uq  = _mem.Uq,                                                           \
    .Dq  = _mem.Dq,                                                           \
                                                                              \
    .Ur  = _mem.Ur,                                                           \
    .Dr  = _mem.Dr,                                                           \
                                                                              \
    .sigmas_x  = _mem.sigmas_x,                                               \
    .sigmas_z  = _mem.sigmas_z,                                               \
                                                                              \
    .wm   = _mem.wm,                                                          \
    .wc   = _mem.wc                                                           \
                                                                              \
    .Sx  = _mem.Sx,                                                           \
                                                                              \
    .Nx  = _nx,                                                               \
    .Nz  = _nz                                                                \
}

/*---------------------------------------------------------------------------*/
static inline void yafl_ukf_post_init(yaflUKFBaseSt * self)
{
    YAFL_ASSERT(self);
    YAFL_ASSERT(self->sp_meth);
    YAFL_ASSERT(self->sp_meth->wf);
    self->sp_meth->wf(self); /*Need to compute weights before start*/
}

static inline void yafl_ukf_gen_sigmas(yaflUKFBaseSt * self)
{
    YAFL_ASSERT(self);
    YAFL_ASSERT(self->sp_meth);
    YAFL_ASSERT(self->sp_meth->spgf);
    self->sp_meth->spgf(self);
}

void yafl_ukf_predict(yaflUKFBaseSt * self);
void yafl_ukf_update(yaflUKFBaseSt * self, yaflFloat * z);

/*===========================================================================*/
void yafl_ukf_base_update(yaflUKFBaseSt * self, yaflFloat * z, \
                               yaflUKFScalarUpdateP scalar_update);

/*===========================================================================*/
void yafl_ukf_bierman_update(yaflUKFBaseSt * self, yaflFloat * z);


/*===========================================================================*/
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

void yafl_ukf_adaptive_bierman_update(yaflUKFAdaptivedSt * self, yaflFloat * z);

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
    yaflFloat sigmas_z[(2 * nx + 1) * nz];

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
