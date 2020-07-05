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

#ifndef YAKF_H
#define YAKF_H

#include <yakf_config.h>
#include "yakf_math.h"

/*=============================================================================
                    Basic UD-factorized EKF definitions
=============================================================================*/
typedef struct _yakfBaseSt yakfBaseSt;

typedef void (* yakfFuncP)(yakfBaseSt *);
typedef void (* yakfResFuncP)(yakfBaseSt *, yakfFloat *);
typedef void (* yakfScalarUpdateP)(yakfBaseSt *, yakfInt);

struct _yakfBaseSt {
    yakfFuncP f;    /*A state transition function*/
    yakfFuncP jf;   /*Jacobian of a state transition function*/

    yakfFuncP h;    /*A measurement function*/
    yakfFuncP jh;   /*Jacobian of a measurement function*/

    yakfResFuncP zrf;  /*Measurement residual function function*/

    yakfFloat * x;  /*State vector*/
    yakfFloat * y;  /*Innovation vector*/
    yakfFloat * H;  /*Measurement Jacobian values*/

    yakfFloat * Up; /*Upper triangular part of P*/
    yakfFloat * Dp; /*Diagonal part of P*/

    yakfFloat * Uq; /*Upper triangular part of Q*/
    yakfFloat * Dq; /*Diagonal part of Q*/

    yakfFloat * Ur; /*Upper triangular part of R*/
    yakfFloat * Dr; /*Diagonal part of R*/

    yakfFloat * W;  /*Scratchpad memory block matrix*/
    yakfFloat * D;  /*Scratchpad memory diagonal matrix*/

    yakfInt   Nx;   /*State vector size*/
    yakfInt   Nz;   /*Measurement vector size*/
};

/*---------------------------------------------------------------------------*/
#define YAKF_BASE_MEMORY_MIXIN(nx, nz) \
    yakfFloat x[nx];                   \
    yakfFloat y[nz];                   \
    yakfFloat H[nz * nx];              \
                                       \
    yakfFloat Up[((nx - 1) * nx)/2];   \
    yakfFloat Dp[nx];                  \
                                       \
    yakfFloat Uq[((nx - 1) * nx)/2];   \
    yakfFloat Dq[nx];                  \
                                       \
    yakfFloat Ur[((nz - 1) * nz)/2];   \
    yakfFloat Dr[nz];                  \
                                       \
    yakfFloat W[2 * nx * nx];          \
    yakfFloat D[2 * nx];

/*---------------------------------------------------------------------------*/
#define YAKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem)\
{                                                                    \
    .f   = (yakfFuncP)_f,                                            \
    .jf  = (yakfFuncP)_jf,                                           \
                                                                     \
    .h   = (yakfFuncP)_h,                                            \
    .jh  = (yakfFuncP)_jh,                                           \
                                                                     \
    .zrf = (yakfResFuncP)_zrf,                                       \
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
void yakf_base_predict(yakfBaseSt * self);
void yakf_base_update(yakfBaseSt * self, yakfFloat * z, yakfScalarUpdateP scalar_update);

/*-----------------------------------------------------------------------------
                               Bierman filter
-----------------------------------------------------------------------------*/
#define YAKF_BIERMAN_PREDICT yakf_base_predict
void yakf_bierman_update(yakfBaseSt * self, yakfFloat * z);

/*-----------------------------------------------------------------------------
                               Joseph filter
-----------------------------------------------------------------------------*/
#define YAKF_JOSEPH_PREDICT yakf_base_predict
void yakf_joseph_update(yakfBaseSt * self, yakfFloat * z);

/*=============================================================================
                    Adaptive UD-factorized EKF definitions
=============================================================================*/
typedef struct {
    yakfBaseSt base;
    yakfFloat chi2; /*Divergence test threshold (chi-squared criteria)*/
} yakfAdaptiveSt; /*Adaptive kKalman-Hinfinity filter structure*/

/*---------------------------------------------------------------------------*/
#define YAKF_ADAPTIVE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem)  \
{                                                                          \
    .base = YAKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem), \
    .chi2 = 10.8275662                                                     \
}
/*
Default value for chi2 is:
  scipy.stats.chi2.ppf(0.999, 1)
*/

/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAKF_ADAPTIVE_BIERAMN_PREDICT(self) yakf_base_predict((yakfBaseSt *)self)
void yakf_adaptive_bierman_update(yakfAdaptiveSt * self, yakfFloat * z);

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAKF_ADAPTIVE_JOSEPH_PREDICT(self) yakf_base_predict((yakfBaseSt *)self)
void yakf_adaptive_joseph_update(yakfAdaptiveSt * self, yakfFloat * z);

/*-----------------------------------------------------------------------------
                                 WARNING!!!

             DO NOT USE THIS variant of Adaptive Joseph filter !!!

     It was implemented to show some flaws of the corresponding algorithm!
-----------------------------------------------------------------------------*/
void yakf_do_not_use_this_update(yakfAdaptiveSt * self, yakfFloat * z);

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
typedef yakfFloat (* yakfRobFuncP)(yakfBaseSt *, yakfFloat);

typedef struct {
    yakfBaseSt base;
    yakfRobFuncP g;    /* g = -d(ln(pdf(y))) / dy */
    yakfRobFuncP gdot; /* gdot = G = d(g) / dy */
} yakfRobustSt; /*Robust EKF*/

/*---------------------------------------------------------------------------*/
#define YAKF_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem)    \
{                                                                          \
    .base = YAKF_BASE_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem), \
    .g    = (yakfRobFuncP)0,                                               \
    .gdot = (yakfRobFuncP)0                                                \
}

/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAKF_ROBUST_BIERAMN_PREDICT(self) yakf_base_predict((yakfBaseSt *)self)
void yakf_robust_bierman_update(yakfRobustSt * self, yakfFloat * z);

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAKF_ROBUST_JOSEPH_PREDICT(self) yakf_base_predict((yakfBaseSt *)self)
void yakf_robust_joseph_update(yakfRobustSt * self, yakfFloat * z);

/*=============================================================================
                 Adaptive robust UD-factorized EKF definitions
=============================================================================*/
typedef struct {
    yakfRobustSt base;
    yakfFloat chi2; /*Divergence test threshold (chi-squared criteria)*/
} yakfAdaptiveRobustSt; /*Robust EKF*/

/*---------------------------------------------------------------------------*/
#define YAKF_ADAPTIVE_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem) \
{                                                                                \
    .base = YAKF_ROBUST_INITIALIZER(_f, _jf, _h, _jh, _zrf, _nx, _nz, _mem),     \
    .chi2 = 8.8074684                                                            \
}
/*
Default value for chi2 is:
  scipy.stats.chi2.ppf(0.997, 1)
*/
/*-----------------------------------------------------------------------------
                           Adaptive Bierman filter
-----------------------------------------------------------------------------*/
#define YAKF_ADAPTIVE_ROBUST_BIERAMN_PREDICT(self) \
    yakf_base_predict((yakfBaseSt *)self)

void yakf_adaptive_robust_bierman_update(yakfAdaptiveRobustSt * self, \
                                         yakfFloat * z);

/*-----------------------------------------------------------------------------
                           Adaptive Joseph filter
-----------------------------------------------------------------------------*/
#define YAKF_ADAPTIVE_ROBUST_JOSEPH_PREDICT(self) \
    yakf_base_predict((yakfBaseSt *)self)

void yakf_adaptive_robust_joseph_update(yakfAdaptiveRobustSt * self, \
                                        yakfFloat * z);

/*=============================================================================
                    Basic UD-factorized UKF definitions
=============================================================================*/
typedef struct _yakfUnscentedSt yakfUnscentedSt;        /* The UKF base type */
/*
Used to add delta vectors to initial point in sigma point generation.
Parameters:
yakfFloat * x0
yakfFloat * delta_x
yakfFloat factor

Does:
delta_x = x0 + factor * delta_x
*/
typedef void (* yakfSigmaAddP)(yakfUnscentedSt *, yakfFloat *, yakfFloat *, yakfFloat);

/*Sigma points base type*/
typedef struct _yakfSigmaSt {
    yakfInt     np;     /* The number of sigma points          */
    yakfFloat * wm;     /* Weights for mean calculations       */
    yakfFloat * wc;     /* Weights for covariance calculations */
    yakfSigmaAddP addf; /* Sigma point addition function       */
} yakfSigmaSt;

/*---------------------------------------------------------------------------*/
#define YAKF_SIGMA_BASE_INITIALIZER(_np, _addf, _mem) \
{                                                     \
    .np   = (yakfInt)_np,                             \
    .addf = (yakfSigmaAddP)_addf,                     \
    .wm   = _mem.wm,                                  \
    .wc   = _mem.wc                                   \
}

/*---------------------------------------------------------------------------*/
/* Computes sigma points weights */
typedef void (* yakfSigmaGenWeigthsP)(yakfUnscentedSt *);

/* Generates sigma points */
typedef void (* yakfSigmaGenSigmasP)(yakfUnscentedSt *);

typedef struct _yakfSigmaMethodsSt {
    yakfSigmaGenWeigthsP   wf; /* Weight function                */
    yakfSigmaGenSigmasP  spgf; /* Sigma point generator function */
} yakfSigmaMethodsSt;

/*---------------------------------------------------------------------------*/
typedef void (* yakfUnscentedFuncP)(yakfUnscentedSt *, yakfFloat *, yakfFloat *);

typedef void (* yakfUnscentedResFuncP)(yakfUnscentedSt *, yakfFloat *, \
                                       yakfFloat *, yakfFloat *);

struct _yakfUnscentedSt {
    /* A pointer to the sigma point generator structure */
    yakfSigmaSt           * points;
    /* A sigma point generator method table pointer     */
    const yakfSigmaMethodsSt * spm;

    yakfUnscentedFuncP      f; /* A state transition function      */
    yakfUnscentedFuncP    xmf; /* State mean function              */
    yakfUnscentedResFuncP xrf; /* State residual function function */

    yakfUnscentedFuncP      h; /* A measurement function                */
    yakfUnscentedFuncP    zmf; /* Measurement mean function function    */
    yakfUnscentedResFuncP zrf; /* Measurement residual function function*/

    yakfFloat * x;  /* State vector                 */
    yakfFloat * zp; /* Predicted measurement vector */

    yakfFloat * Up;  /* Upper triangular part of P  */
    yakfFloat * Dp;  /* Diagonal part of P          */

    yakfFloat * Us;  /* Upper triangular part of S  */
    yakfFloat * Ds;  /* Diagonal part of S          */

    yakfFloat * Pzx; /* Pzx cross covariance matrix */

    yakfFloat * Uq;  /* Upper triangular part of Q   */
    yakfFloat * Dq;  /* Diagonal part of Q           */

    yakfFloat * Ur;  /* Upper triangular part of R   */
    yakfFloat * Dr;  /* Diagonal part of R           */

    yakfFloat * sigmas_x; /* State sigma points       */
    yakfFloat * sigmas_z; /* Measurement sigma points */

    /*Scratchpad memory*/
    yakfFloat * Sx;       /* State       */
    yakfFloat * Sz;       /* Measurement */

    yakfInt   Nx; /* State vector size       */
    yakfInt   Nz; /* Measurement vector size */
};

/*---------------------------------------------------------------------------*/
/*
Warning: sigmas_x and _sigmas_z aren't defined in this mixin, see
         sigma points generators mixins!!!
*/
#define YAKF_UNSCENTED_MEMORY_MIXIN(nx, nz) \
    yakfFloat x[nx];                        \
    yakfFloat zp[nz];                       \
                                            \
    yakfFloat Up[((nx - 1) * nx)/2];        \
    yakfFloat Dp[nx];                       \
                                            \
    yakfFloat Us[((nz - 1) * nz)/2];        \
    yakfFloat Ds[nz];                       \
                                            \
    yakfFloat Pzx[nz * nx];                 \
                                            \
    yakfFloat Uq[((nx - 1) * nx)/2];        \
    yakfFloat Dq[nx];                       \
                                            \
    yakfFloat Ur[((nz - 1) * nz)/2];        \
    yakfFloat Dr[nz];                       \
                                            \
    yakfFloat Sx[nx];                       \
    yakfFloat Sz[nz];

/*---------------------------------------------------------------------------*/
#define YAKF_UNSCENTED_INITIALIZER(_p, _pm, _f, _xmf, _xrf, _h, _zmf,         \
                                   _zrf, _nx, _nz, _mem)                      \
{                                                                             \
    .f   = (yakfUnscentedFuncP)_f,                                            \
    .xmf = (yakfUnscentedFuncP)_xmf,                                          \
    .xrf = (yakfUnscentedResFuncP)_xrf,                                       \
                                                                              \
    .h   = (yakfUnscentedFuncP)_h,                                            \
    .xmf = (yakfUnscentedFuncP)_xmf,                                          \
    .xrf = (yakfUnscentedResFuncP)_xrf,                                       \
                                                                              \
    .x   = _mem.x,                                                            \
    .zp  = _mem.zp,                                                           \
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
    .Sx  = _mem.Sx,                                                           \
    .Sz  = _mem.Sz,                                                           \
                                                                              \
    .Nx  = _nx,                                                               \
    .Nz  = _nz                                                                \
}

/*---------------------------------------------------------------------------*/
static inline void yakf_unscented_post_init(yakfUnscentedSt * self)
{
    YAKF_ASSERT(self);
    YAKF_ASSERT(self->spm);
    YAKF_ASSERT(self->spm->wf);
    self->spm->wf(self); /*Need to compute weights before start*/
}

static inline void yakf_unscented_gen_sigmas(yakfUnscentedSt * self)
{
    YAKF_ASSERT(self);
    YAKF_ASSERT(self->spm);
    YAKF_ASSERT(self->spm->spgf);
    self->spm->spgf(self);
}

void yakf_unscented_predict(yakfUnscentedSt * self);
void yakf_unscented_update(yakfUnscentedSt * self, yakfFloat * z);

/*-----------------------------------------------------------------------------
                     Van der Merwe sigma point generator
-----------------------------------------------------------------------------*/
typedef struct _yakfMerweSt {
    yakfSigmaSt base;
    yakfFloat alpha;
    yakfFloat beta;
    yakfFloat kappa;
} yakfMerweSt;

/*---------------------------------------------------------------------------*/
#define YAKF_MERWE_MEMORY_MIXIN(nx, nz)    \
    yakfFloat wm[2 * nx + 1];              \
    yakfFloat wc[2 * nx + 1];              \
    yakfFloat sigmas_x[(2 * nx + 1) * nx]; \
    yakfFloat sigmas_z[(2 * nx + 1) * nz];

/*---------------------------------------------------------------------------*/
#define YAKF_MERWE_INITIALIZER(_nx, _addf, _alpha, _beta, _kappa, _mem) \
{                                                                       \
    .base  = YAKF_SIGMA_BASE_INITIALIZER((2 * _nx + 1), _addf, _mem),   \
    .alpha = _alpha,                                                    \
    .beta  = _beta,                                                     \
    .kappa = _kappa                                                     \
}

/*---------------------------------------------------------------------------*/
const yakfSigmaMethodsSt yakf_merwe_spm;

#endif // YAKF_H
