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

#endif // YAKF_H
