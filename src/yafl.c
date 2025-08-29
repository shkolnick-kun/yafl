/*******************************************************************************
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
******************************************************************************/
#include <string.h>

#include "yafl.h"

/*np.log(2. * np.pi)*/
#define YAFL_L2PI (1.8378770664093453)

#define _fx  (self->f)
#define _hx  (self->h)
#define _zrf (self->zrf)
#define _rcb (self->rcb)

#define _x   (self->x)
#define _y   (self->y)

#define _up  (self->Up)
#define _dp  (self->Dp)

#define _uq  (self->Uq)
#define _dq  (self->Dq)

#define _ur  (self->Ur)
#define _dr  (self->Dr)

#define _l   (self->l)

#define _nx  (self->Nx)
#define _nz  (self->Nz)

/*=============================================================================
                                  Base UDEKF
=============================================================================*/
#define _jfx (((yaflEKFBaseSt *)self)->jf)
#define _jhx (((yaflEKFBaseSt *)self)->jh)

#define _hy  (((yaflEKFBaseSt *)self)->H)
#define _w   (((yaflEKFBaseSt *)self)->W)
#define _d   (((yaflEKFBaseSt *)self)->D)

#define _qff (((yaflEKFBaseSt *)self)->qff)

yaflStatusEn yafl_ekf_base_predict(yaflKalmanBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt i;
    yaflInt nx2;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_up,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_dp,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_uq,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_dq,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_nx > 0, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_w,      YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_d,      YAFL_ST_INV_ARG_1);

    nx2 = _nx * 2;

    /*Default f(x) = x*/
    if (0 == _fx)
    {
        YAFL_CHECK(0 == _jfx, YAFL_ST_INV_ARG_1);

        for (i = 0; i < _nx; i++)
        {
            yaflInt j;
            yaflInt nci;

            nci = nx2 * i;
            for (j = 0; j < _nx; j++)
            {
                _w[nci + j] = (i != j) ? 0.0 : 1.0;
            }
        }
    }
    else
    {
        //yaflFloat * x;

        /*Must have some Jacobian function*/
        YAFL_CHECK(_jfx, YAFL_ST_INV_ARG_1);

        //x = self->x;
        YAFL_TRY(status,  _fx(self, _x, _x));  /* x = f(x_old, ...) */
        YAFL_TRY(status, _jfx(self, _w, _x));  /* Place F(x, ...)=df/dx to W  */
    }
    /* Now W = (F|***) */
    YAFL_TRY(status, \
             YAFL_MATH_BSET_BU(nx2, 0, _nx, _w, _nx, _nx, nx2, 0, 0, _w, _up));
    /* Now W = (F|FUp) */
    YAFL_TRY(status, yafl_math_bset_u(nx2, _w, _nx, _uq));
    /* Now W = (Uq|FUp) */

    /* D = concatenate([Dq, Dp]) */
    i = _nx*sizeof(yaflFloat);
    memcpy((void *)       _d,  (void *)_dq, i);
    memcpy((void *)(_d + _nx), (void *)_dp, i);

    /* Up, Dp = MWGSU(w, d)*/
    YAFL_TRY(status, yafl_math_mwgsu(_nx, nx2, _up, _dp, _w, _d));

    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn \
    _yafl_r_update(yaflInt nx, yaflInt nz, yaflFloat rff, yaflFloat * dp, \
                   yaflFloat * ur, yaflFloat * dr, yaflFloat * w, yaflFloat * d, \
                   yaflFloat * res_z)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt      nxz = nx + nz;

    YAFL_CHECK(nx > 0,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(nz > 0,     YAFL_ST_INV_ARG_2);

    YAFL_CHECK(rff >  0.0, YAFL_ST_INV_ARG_3);
    YAFL_CHECK(rff <  1.0, YAFL_ST_INV_ARG_3);

    YAFL_CHECK(dp,         YAFL_ST_INV_ARG_4);
    YAFL_CHECK(ur,         YAFL_ST_INV_ARG_5);
    YAFL_CHECK(dr,         YAFL_ST_INV_ARG_6);
    YAFL_CHECK(w,          YAFL_ST_INV_ARG_7);
    YAFL_CHECK(d,          YAFL_ST_INV_ARG_8);
    YAFL_CHECK(res_z,      YAFL_ST_INV_ARG_9);

    /*Will end R update, now W = (HUp|***) */
    YAFL_TRY(status, YAFL_MATH_BSET_U(nxz, 0, nx, w, nz, ur));
    /*Now W = (HUp|Ur) */

    /*Now D=(***)*/
    YAFL_TRY(status, yafl_math_set_vxn(nx,      d, dp, rff      ));
    /*Now D=(rff*Dp|***)*/
    YAFL_TRY(status, yafl_math_set_vxn(nz, d + nx, dr, 1.0 - rff));
    /* Now D=(rff * Dp|(1-rff) * dr) */

    /* Ur, Dr = MWGSU(w, d)*/
    YAFL_TRY(status, yafl_math_mwgsu(nz, nxz, ur, dr, w, d));
    /*Now R+ = (1.0 - rff) * R- + rff * H.dot(P.dot(H.T))*/
    YAFL_TRY(status, yafl_math_udu_up(nz, ur, dr, rff, res_z));
    /*Now R+ = (1.0 - rff) * R- + rff * (res_z.dot(res_z.T) + H.dot(P.dot(H.T)))*/

    return status;
}

/*---------------------------------------------------------------------------*/
static yaflStatusEn _yafl_ekf_compute_error(yaflKalmanBaseSt * self, yaflFloat * z)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt j;

    YAFL_CHECK(self,    YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_hx,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_x,      YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_y,      YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_nx > 0, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_nz > 0, YAFL_ST_INV_ARG_1);

    YAFL_TRY(status, _hx(self, _y,  _x)); /* self.y = h(x,...) */

    /*Compute measurement error*/
    if (0 == _zrf)
    {
        /*Default residual*/
        for (j = 0; j < _nz; j++)
        {
            _y[j] = z[j] - _y[j];
        }
    }
    else
    {
        /*zrf must be aware of self internal structure*/
        YAFL_TRY(status, _zrf(self, _y, z, _y)); /* self.y = zrf(z, h(x,...)) */
    }

    return status;
}

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ekf_base_update(yaflKalmanBaseSt * self, yaflFloat * z, yaflKalmanScalarUpdateP scalar_update)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt j;

    YAFL_CHECK(self,   YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_hx,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_x,      YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_y,      YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_ur,     YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_nx > 0, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_nz > 0, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_jhx,    YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_hy,     YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_w,      YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_d,      YAFL_ST_INV_ARG_1);

    YAFL_CHECK(z,      YAFL_ST_INV_ARG_2);
    YAFL_CHECK(scalar_update, YAFL_ST_INV_ARG_3);

    /*
    _y contains residual which we will need for possible R update.
    Will compute innovation later.
    */
    YAFL_TRY(status, _jhx(self, _hy, _x)); /* self.H = jh(x,...) */

    /*Update R if needed!*/
    if (self->rff > 0.0)
    {
        yaflStatusEn status_r = status; /*For quiet regularization*/
        /*Check input data*/
        YAFL_CHECK(_up, YAFL_ST_INV_ARG_1);
        YAFL_CHECK(_dp, YAFL_ST_INV_ARG_1);
        YAFL_CHECK(_dr, YAFL_ST_INV_ARG_1);
        /*Start R update*/
        YAFL_TRY(status_r, yafl_math_bset_mu(_nx + _nz, _w, _nz, _nx, _hy, _up));
        /*Now W = (HUp|***) */
        YAFL_TRY(status_r, \
                 _yafl_r_update(_nx, _nz, self->rff, _dp, _ur, _dr, _w, _d, _y));
        /*Updated R*/
    }
    /*R update call back*/
    if (_rcb)
    {
        _rcb(self);
    }
    /*Now _y can be used for innovation!*/

    /*Compute innovation*/
    YAFL_TRY(status, _yafl_ekf_compute_error(self, z));

    /* Decorrelate measurement noise */
    YAFL_TRY(status, yafl_math_ruv(_nz,      _y,  _ur));
    YAFL_TRY(status, yafl_math_rum(_nz, _nx, _hy, _ur));

    /*Start Q update if needed!*/
    if (_qff > 0.0)
    {
        yaflStatusEn status_q = status; /*For quiet regularization*/
        YAFL_CHECK(_dq, YAFL_ST_INV_ARG_1);
        YAFL_TRY(status_q, yafl_math_set_vxn(_nx, _dq, _dq, 1.0 - _qff));
    }

    /*Start log likelihood computation*/
    *self->l = _nz * YAFL_L2PI;

    /* Do scalar updates */
    for (j = 0; j < _nz; j++)
    {
        YAFL_TRY(status, scalar_update(self, j)); /*Q is updated here if needed!*/
    }

    /*Finalize log likelihood value*/
    *self->l *= -0.5;

    if (self->rff > 0.0)
    {
        /*Compute residual if needed*/
        YAFL_TRY(status, _yafl_ekf_compute_error(self, z));
    }
    return status;
}
/*=============================================================================
                                Bierman filter
=============================================================================*/
static inline yaflStatusEn \
    _bierman_update_body(yaflInt    nx, yaflFloat * x, yaflFloat * u, \
                        yaflFloat * d, yaflFloat * f, yaflFloat * v, \
                        yaflFloat   r, yaflFloat  nu, yaflFloat  ac, \
                        yaflFloat   gdot, yaflFloat * l)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt j;
    yaflInt k;
    yaflInt nxk;

    YAFL_CHECK(x, YAFL_ST_INV_ARG_2);
    YAFL_CHECK(u, YAFL_ST_INV_ARG_3);
    YAFL_CHECK(d, YAFL_ST_INV_ARG_4);
    YAFL_CHECK(f, YAFL_ST_INV_ARG_5);
    YAFL_CHECK(v, YAFL_ST_INV_ARG_6);

    for (k = 0, nxk = 0; k < nx; nxk += k++)
    {
        yaflFloat a;
        yaflFloat fk;
        yaflFloat vk;

        fk = gdot * f[k];
        /*Correct v in place*/
        vk = ac * v[k];
        v[k] = vk;
        a = r + fk * vk;
        /*Correct d in place*/
        d[k] *= ac * r / a;
#define p fk /*No need for separate p variable*/
        p = - fk / r;
        for (j = 0; j < k; j++)
        {
            yaflFloat ujk;
            yaflFloat vj;

            ujk = u[j + nxk];
            vj  = v[j];

            u[j + nxk] = ujk +   p * vj;
            v[j]       = vj  + ujk * vk;
        }
#undef  p /*Don't need p any more...*/
        r = a;
    }

    /*Compute log likelihood*/
    *l += nu * (nu / r) + YAFL_LOG(r);

    /*
    Now we must do:
    x += K * nu

    Since:
    r == a

    then we have:
    K == v / a == v / r

    and so:
    K * nu == (v / r) * nu == v / r * nu == v * (nu / r)

    Finally we get:
    x += v * (nu / r)
    */
    nu /= r;
    for (j=0; j < nx; j++)
    {
        v[j] *= nu;
        x[j] += v[j];
    }
    return status;
}

/*---------------------------------------------------------------------------*/
#define YAFL_Q_SCALAR_UPDATE(qff, knu)                                 \
do {                                                                   \
    if (qff > 0.0)                                                     \
    {                                                                  \
        yaflStatusEn status_q = status;                                \
        YAFL_CHECK(_uq, YAFL_ST_INV_ARG_1);                            \
        YAFL_CHECK(_dq, YAFL_ST_INV_ARG_1);                            \
        YAFL_TRY(status_q, yafl_math_udu_up(_nx, _uq, _dq, qff, knu)); \
    }                                                                  \
} while (0)

/*---------------------------------------------------------------------------*/
#define YAFL_SCALAR_UPDATE_ARGS_CHECKS()          \
do {                                              \
    YAFL_CHECK(self,         YAFL_ST_INV_ARG_1);  \
    YAFL_CHECK(self->Nz > i, YAFL_ST_INV_ARG_2);  \
} while (0)

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_BIERMAN_SELF_INTERNALS_CHECKS()  \
do {                                          \
    YAFL_CHECK(_nx > 0,   YAFL_ST_INV_ARG_1); \
    YAFL_CHECK(_up, YAFL_ST_INV_ARG_1);       \
    YAFL_CHECK(_dp, YAFL_ST_INV_ARG_1);       \
    YAFL_CHECK(_hy, YAFL_ST_INV_ARG_1);       \
    YAFL_CHECK(_y,  YAFL_ST_INV_ARG_1);       \
    YAFL_CHECK(_dr, YAFL_ST_INV_ARG_1);       \
    YAFL_CHECK(_d,  YAFL_ST_INV_ARG_1);       \
} while (0)

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ekf_bierman_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat * h;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_BIERMAN_SELF_INTERNALS_CHECKS();

    h = _hy + _nx * i;
    /* f = h.dot(Up) */
#   define f _d
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#define v h /*Don't need h any more, use it to store v*/
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));
    YAFL_TRY(status, \
             _bierman_update_body(_nx, _x, _up, _dp, f, v, _dr[i], _y[i], \
                                  1.0, 1.0, _l));
    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, v);
#   undef v /*Don't nee v any more*/
#   undef f

    return status;
}

/*=============================================================================
                                Joseph filter
=============================================================================*/
static inline yaflStatusEn \
    _joseph_update_body(yaflInt nx,    yaflFloat * x, yaflFloat * u, \
                        yaflFloat * d, yaflFloat * f, yaflFloat * v, \
                        yaflFloat * k, yaflFloat * w, yaflFloat  nu, \
                        yaflFloat  r, yaflFloat   s, yaflFloat  ac, \
                        yaflFloat  gdot, yaflFloat * l)
{
    yaflStatusEn status = YAFL_ST_OK;

    yaflInt nx1;

    YAFL_CHECK(x, YAFL_ST_INV_ARG_2);
    YAFL_CHECK(u, YAFL_ST_INV_ARG_3);
    YAFL_CHECK(d, YAFL_ST_INV_ARG_4);
    YAFL_CHECK(f, YAFL_ST_INV_ARG_5);
    YAFL_CHECK(v, YAFL_ST_INV_ARG_6);
    YAFL_CHECK(k, YAFL_ST_INV_ARG_7);
    YAFL_CHECK(w, YAFL_ST_INV_ARG_8);
    YAFL_CHECK(s > 0, YAFL_ST_INV_ARG_11);

    nx1 = nx + 1;

    /* k = Up.dot(v * ac / s) = Up.dot(v) * (ac / s) */
    /*May be used in place*/
    YAFL_TRY(status, yafl_math_set_vxn(nx, v, v, ac / s));
    YAFL_TRY(status, yafl_math_set_uv(nx, k, u, v));

#   define D v
    /*Set W and D*/
    /*May be used in place*/
    YAFL_TRY(status, yafl_math_set_vxn(nx, f, f, gdot));
    /*How about yafl_math_bset_vvtxn ?*/
    YAFL_TRY(status, yafl_math_bset_vvt(nx1, w, nx, k, f));
    YAFL_TRY(status, yafl_math_bsub_u(nx1, w, nx, u));

    /* Now w is (gdot*kf - Up|***) */
    YAFL_TRY(status, YAFL_MATH_BSET_V(nx1, 0, nx, w, nx, k));
    /* Now w is (gdot*kf - Up|k) */

    /* D = concatenate([ac * Dp, np.array([gdot * r])]) */
    YAFL_TRY(status, yafl_math_set_vxn(nx, D, d, ac));
    D[nx] = gdot * r;

    /* Up, Dp = MWGSU(W, D)*/
    YAFL_TRY(status, yafl_math_mwgsu(nx, nx1, u, d, w, D));

    /*Compute log likelihood*/
    *l += nu * (nu / s) + YAFL_LOG(s);

    /* x += k * nu */
#   define j nx1
    for (j=0; j < nx; j++)
    {
        k[j] *= nu;
        x[j] += k[j];
    }
#   undef j  /*Don't nee j any more*/
#   undef D  /*Don't nee D any more*/
    return status;
}

/*---------------------------------------------------------------------------*/
#define YAFL_EKF_JOSEPH_SELF_INTERNALS_CHECKS() \
do {                                            \
    YAFL_CHECK(_nx > 0, YAFL_ST_INV_ARG_1);     \
    YAFL_CHECK(_up,     YAFL_ST_INV_ARG_1);     \
    YAFL_CHECK(_dp,     YAFL_ST_INV_ARG_1);     \
    YAFL_CHECK(_hy,     YAFL_ST_INV_ARG_1);     \
    YAFL_CHECK(_y,      YAFL_ST_INV_ARG_1);     \
    YAFL_CHECK(_dr,     YAFL_ST_INV_ARG_1);     \
    YAFL_CHECK(_w,      YAFL_ST_INV_ARG_1);     \
    YAFL_CHECK(_d,      YAFL_ST_INV_ARG_1);     \
} while (0)

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ekf_joseph_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat s = 0.0;
    yaflFloat  * f;
    yaflFloat  * h;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_JOSEPH_SELF_INTERNALS_CHECKS();

#   define v _d
    f = v + _nx;
    h = _hy + _nx * i;

    /* f = h.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));

#   define r _dr[i]
    /* s = r + f.dot(v)*/
    YAFL_TRY(status, yafl_math_vtv(_nx, &s, f, v));
    s += r;

    /*K = Up.dot(v/s) = Up.dot(v)/s*/
#   define K h /*Don't need h any more, use it to store K*/
    YAFL_TRY(status, \
             _joseph_update_body(_nx, _x, _up, _dp, f, v, K, _w, _y[i], r, \
                                 s, 1.0, 1.0, _l));

    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, K);
#   undef K /*Don't nee K any more*/
#   undef r
#   undef v

    return status;
}

/*=============================================================================
                          Adaptive Bierman filter
=============================================================================*/
static inline yaflStatusEn \
    _adaptive_correction(yaflInt   nx, yaflFloat * res_ac, yaflFloat * res_s, \
                         yaflFloat * f, yaflFloat *      v, yaflFloat      r, \
                         yaflFloat  nu, yaflFloat     gdot, yaflFloat chi2)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat c  = 0.0;
    yaflFloat ac;
    yaflFloat s;

    YAFL_CHECK(res_ac, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(f,      YAFL_ST_INV_ARG_3);
    YAFL_CHECK(v,      YAFL_ST_INV_ARG_4);
    YAFL_CHECK(chi2 > 0,   YAFL_ST_INV_ARG_8);

    /* s = alpha**2 + gdot * f.dot(v)*/
    YAFL_TRY(status, yafl_math_vtv(nx, &c, f, v));
    c *= gdot;
    s = r + c;

    /* Divergence test */
    ac = gdot * (nu * (nu / chi2)) - s;
    if (ac > 0.0)
    {
        /*Adaptive correction factor*/
        ac = ac / c + 1.0;

        /*Corrected s*/
        s = ac * c + r;

        status |= YAFL_ST_MSK_ANOMALY; /*Anomaly detected!*/
    }
    else
    {
        ac = 1.0;
    }

    *res_ac = ac;

    if (res_s)
    {
        *res_s = s;
    }

    return status;
}

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ekf_adaptive_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                     yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat ac = 1.0;
    yaflFloat nu = 0.0;
    yaflFloat * h;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_BIERMAN_SELF_INTERNALS_CHECKS();

    nu = _y[i];
    h = _hy + _nx * i;

    /* f = h.dot(Up) */
#   define f _d
    //f = self->D;
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#   define v _hy /*Don't need h any more, use it to store v*/
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));

#   define r _dr[i]
    YAFL_TRY(status, \
             _adaptive_correction(_nx, &ac, 0, f, v, r, nu, 1.0, \
                                  ((yaflEKFAdaptiveSt *)self)->chi2));

    YAFL_TRY(status, \
             _bierman_update_body(_nx, _x, _up, _dp, f, v, r, nu, ac, 1.0, _l));

    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, v);
#   undef r
#   undef v /*Don't nee v any more*/
#   undef f

    return status;
}

/*=============================================================================
                           Adaptive Joseph filter
=============================================================================*/
yaflStatusEn \
    yafl_ekf_adaptive_joseph_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat s  = 1.0;
    yaflFloat ac = 1.0;
    yaflFloat nu = 0.0;
    yaflFloat * h;
    yaflFloat * f;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_JOSEPH_SELF_INTERNALS_CHECKS();

    h = _hy + _nx * i;

#   define v _d
    f = v + _nx;

    nu = _y[i];

    /* f = h.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));

#   define r _dr[i]
    YAFL_TRY(status, \
             _adaptive_correction(_nx, &ac, &s, f, v, r, nu, 1.0, \
                                  ((yaflEKFAdaptiveSt *)self)->chi2));

    /* K = Up.dot(v * ac / s) = Up.dot(v) * (ac / s) */
#   define K h /*Don't need h any more, use it to store K*/
    YAFL_TRY(status, \
             _joseph_update_body(_nx, _x, _up, _dp, f, v, K, _w, nu, r, s, \
                                 ac, 1.0, _l));

    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, K);
#   undef K /*Don't nee K any more*/
#   undef r
#   undef v

    return status;
}

/*=============================================================================
                                 WARNING!!!

             DO NOT USE THIS variant of Adaptive Joseph filter !!!

    It was implemented to show some flaws of the corresponding algorithm!
=============================================================================*/
yaflStatusEn \
    yafl_ekf_do_not_use_this_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat c = 0.0;
    yaflFloat s = 0.0;
    yaflInt nx;
    yaflInt nx1;
    yaflFloat * d;
    yaflFloat * u;
    yaflFloat * h;
    yaflFloat * f;
    yaflFloat * v;
    yaflFloat * w;
    yaflFloat nu;
    yaflFloat r;
    yaflFloat ac;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();

    nx = _nx;
    YAFL_EKF_JOSEPH_SELF_INTERNALS_CHECKS();

    nx1 = nx + 1;

    d = _dp;
    u = _up;

    h = _hy + nx * i;

    v = _d;
    f = v + nx;

    w = _w;

    nu = _y[i];
    r  = _dr[i];

    /* f = h.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(nx, f, h, u));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(nx, v, d, f));

    /* s = r + f.dot(v)*/
    YAFL_TRY(status, yafl_math_vtv(nx, &c, f, v));
    s = c + r;

    /* Divergence test */
    ac = (nu * (nu / (((yaflEKFAdaptiveSt *)self)->chi2))) - s;
    if (ac > 0.0)
    {
        /*Adaptive correction with no limitations approach*/
        YAFL_TRY(status, yafl_math_set_uv(nx, f, u, v));
        YAFL_TRY(status, yafl_math_udu_up(nx, u, d, (ac / c) / c, f));

        /*Recompute f,v,s*/
        /* f = h.dot(Up) */
        YAFL_TRY(status, yafl_math_set_vtu(nx, f, h, u));

        /* v = f.dot(Dp).T = Dp.dot(f.T).T */
        YAFL_TRY(status, YAFL_MATH_SET_DV(nx, v, d, f));

        /* s = r + f.dot(v)*/
        YAFL_TRY(status, yafl_math_vtv(nx, &c, f, v));
        s = c + r;
    }

    /* K = Up.dot(v * ac / s) = Up.dot(v) * (ac / s) */
#define K h /*Don't need h any more, use it to store K*/
#define D v
    YAFL_TRY(status, yafl_math_set_uv(nx, K, u, v));
    YAFL_TRY(status, yafl_math_set_vrn(nx, K, K, s)); /*May be used in place*/

    /*Set W and D*/
    YAFL_TRY(status, yafl_math_bset_vvt(nx1, w, nx, K, f));
    YAFL_TRY(status, yafl_math_bsub_u(nx1, w, nx, u));

    /* Now w is (Kf - Up|***) */
    YAFL_TRY(status, YAFL_MATH_BSET_V(nx1, 0, nx, w, nx, K));
    /* Now w is (Kf - Up|K) */

    /* D = concatenate([Dp, np.array([r])]) */
    memcpy((void *)D, (void *)d, nx * sizeof(yaflFloat));
    D[nx] = r;

    /* Up, Dp = MWGSU(W, D)*/
    YAFL_TRY(status, yafl_math_mwgsu(nx, nx1, u, d, w, D));

    /*Compute log likelihood*/
    *self->l += nu * (nu / s) + YAFL_LOG(s);

    /* self.x += K * nu */
    YAFL_TRY(status, yafl_math_add_vxn(nx, _x, K, nu));
#undef D /*Don't nee D any more*/
#undef K /*Don't nee K any more*/
    return status;
}

/*=============================================================================
                            Robust Bierman filter
=============================================================================*/
static inline yaflStatusEn \
    _scalar_robustify(yaflKalmanBaseSt * self, \
                      yaflKalmanRobFuncP g, yaflKalmanRobFuncP gdot, \
                      yaflFloat * gdot_res, yaflFloat * nu, yaflFloat r)
{
    yaflFloat tmp;

    YAFL_CHECK(self,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(gdot_res, YAFL_ST_INV_ARG_3);
    YAFL_CHECK(nu,       YAFL_ST_INV_ARG_4);

    if (g)
    {
        r = YAFL_SQRT(r);/*We need sqrt here*/

        tmp = *nu / r;
        *nu = r * g(self, tmp);

        YAFL_CHECK(gdot, YAFL_ST_INV_ARG_1);

        tmp = gdot(self, tmp);
    }
    else
    {
        tmp = 1.0;
    }

    *gdot_res = tmp;

    /*Detect glitches and return*/
    if (tmp < YAFL_EPS)
    {
        return YAFL_ST_MSK_GLITCH_LARGE;
    }

    if (tmp < (1.0 - 2.0*YAFL_EPS))
    {
        return YAFL_ST_MSK_GLITCH_SMALL;
    }

    return YAFL_ST_OK;
}

#define YAFL_EKF_SCALAR_ROBUSTIFY(_self, _gdot_res, _nu, _r) \
    _scalar_robustify(_self,                                 \
                      ((yaflEKFRobustSt *)_self)->g,         \
                      ((yaflEKFRobustSt *)_self)->gdot,      \
                      _gdot_res, _nu, _r)

/*---------------------------------------------------------------------------*/
yaflStatusEn \
    yafl_ekf_robust_bierman_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat gdot = 1.0;
    yaflFloat nu   = 0.0;
    yaflFloat r;
    yaflFloat * h;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_BIERMAN_SELF_INTERNALS_CHECKS();

    r  = _dr[i];
    nu = _y[i];

    YAFL_TRY(status, YAFL_EKF_SCALAR_ROBUSTIFY(self, &gdot, &nu, r));

    h = _hy + _nx * i;
    /* f = h.dot(Up) */
#   define f _d
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#   define v h /*Don't need h any more, use it to store v*/
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));

    YAFL_TRY(status, \
             _bierman_update_body(_nx, _x, _up, _dp, f, v, r, nu, \
                                  1.0, gdot, _l));

    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, v);
#   undef v  /*Don't nee v any more*/
#   undef f

    return status;
}

/*=============================================================================
                            Robust Joseph filter
=============================================================================*/
yaflStatusEn \
    yafl_ekf_robust_joseph_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat gdot = 0.0;
    yaflFloat s    = 0.0;
    yaflFloat nu   = 0.0;
    yaflFloat r;
    yaflFloat * h;
    yaflFloat * f;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_JOSEPH_SELF_INTERNALS_CHECKS();

    r   = _dr[i];
    nu  = _y[i];

    YAFL_TRY(status, YAFL_EKF_SCALAR_ROBUSTIFY(self, &gdot, &nu, r));

    h = _hy + _nx * i;
#   define v _d
    f = v + _nx;

    /* f = h.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));

    /* s = r + gdot * f.dot(v)*/
    YAFL_TRY(status, yafl_math_vtv(_nx, &s, f, v));
    s = r + gdot * s;

    /*K = Up.dot(v/s) = Up.dot(v)/s*/
#   define K h /*Don't need h any more, use it to store K*/
    YAFL_TRY(status, \
             _joseph_update_body(_nx, _x, _up, _dp, f, v, K, _w, nu, r, s, \
                                 1.0, gdot, _l));

    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, K);
#   undef K  /*Don't nee K any more*/
#   undef v  /*Don't nee v any more*/
    return status;
}

/*=============================================================================
                        Adaptive robust Bierman filter
=============================================================================*/
yaflStatusEn \
    yafl_ekf_adaptive_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                   yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat gdot = 1.0;
    yaflFloat ac   = 1.0;
    yaflFloat nu   = 0.0;
    yaflFloat r;
    yaflFloat * h;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_BIERMAN_SELF_INTERNALS_CHECKS();

    r  = _dr[i];
    nu = _y[i];

    YAFL_TRY(status, YAFL_EKF_SCALAR_ROBUSTIFY(self, &gdot, &nu, r));

    h = _hy + _nx * i;

    /* f = h.dot(Up) */
#   define f _d
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#   define v h /*Don't need h any more, use it to store v*/
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));

    YAFL_TRY(status, \
             _adaptive_correction(_nx, &ac, 0, f, v, r, nu, gdot, \
                                  ((yaflEKFAdaptiveRobustSt *)self)->chi2));

    YAFL_TRY(status, _bierman_update_body(_nx, _x, _up, _dp, f, v, r, nu, \
                                          ac, gdot, _l));

    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, v);
#   undef v  /*Don't nee v any more*/
#   undef f

    return status;
}

/*=============================================================================
                        Adaptive robust Joseph filter
=============================================================================*/
yaflStatusEn \
    yafl_ekf_adaptive_robust_joseph_update_scalar(yaflKalmanBaseSt * self, \
                                                  yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat gdot = 1.0;
    yaflFloat   ac = 1.0;
    yaflFloat    s = 0.0;
    yaflFloat * h;
    yaflFloat * f;
    yaflFloat r;
    yaflFloat nu;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    YAFL_EKF_JOSEPH_SELF_INTERNALS_CHECKS();

    r  = _dr[i]; /* alpha = r**0.5 is stored in Dr*/
    nu = _y[i];

    YAFL_TRY(status, YAFL_EKF_SCALAR_ROBUSTIFY(self, &gdot, &nu, r));

    h = _hy + _nx * i;
#   define v _d
    f = v + _nx;

    /* f = h.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(_nx, f, h, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(_nx, v, _dp, f));

    YAFL_TRY(status, \
             _adaptive_correction(_nx, &ac, &s, f, v, r, nu, gdot, \
                                  ((yaflEKFAdaptiveRobustSt *)self)->chi2));

    /* K = Up.dot(v * ac / s) = Up.dot(v) * (ac / s) */
#   define K h /*Don't need h any more, use it to store K*/
    YAFL_TRY(status, \
             _joseph_update_body(_nx, _x, _up, _dp, f, v, K, _w, nu, r, s, \
                                 ac, gdot, _l));

    /*Update Q if needed!*/
    YAFL_Q_SCALAR_UPDATE(_qff, K);
#   undef K  /*Don't nee K any more*/
#   undef v

    return status;
}

/*------------------------------------------------------------------------------
                                 Undef EKF stuff
------------------------------------------------------------------------------*/
#undef _jfx
#undef _jhx

#undef _hy
#undef _w
#undef _d
#undef _qff

/*=============================================================================
                    Basic UD-factorized UKF functions
=============================================================================*/
#define _kalman_self ((yaflKalmanBaseSt *)self)

#define _ufx  (_kalman_self->f)
#define _uhx  (_kalman_self->h)
#define _uzrf (_kalman_self->zrf)
#define _urcb (_kalman_self->rcb)

#define _ux   (_kalman_self->x)
#define _uy   (_kalman_self->y)

#define _uup  (_kalman_self->Up)
#define _udp  (_kalman_self->Dp)

#define _uuq  (_kalman_self->Uq)
#define _udq  (_kalman_self->Dq)

#define _uur  (_kalman_self->Ur)
#define _udr  (_kalman_self->Dr)

#define _ul   (_kalman_self->l)

#define _unx  (_kalman_self->Nx)
#define _unz  (_kalman_self->Nz)

/*UKF stuff*/
#define _xrf       (self->xrf)
#define _xmf       (self->xmf)

#define _zmf       (self->zmf)

#define _zp        (self->zp)

#define _sx        (self->Sx)
#define _sz        (self->Sz)

#define _pzx       (self->Pzx)

#define _sigmas_x  (self->sigmas_x)
#define _sigmas_z  (self->sigmas_z)

#define _wm        (self->wm)
#define _wc        (self->wc)

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _compute_res(yaflKalmanBaseSt * self, yaflInt sz,    \
                                        yaflKalmanResFuncP rf, yaflFloat * res, \
                                        yaflFloat * sigma, yaflFloat * pivot)
{
    yaflStatusEn status = YAFL_ST_OK;

    YAFL_CHECK(self,  YAFL_ST_INV_ARG_1);
    YAFL_CHECK(res,   YAFL_ST_INV_ARG_4);
    YAFL_CHECK(sigma, YAFL_ST_INV_ARG_5);
    YAFL_CHECK(pivot, YAFL_ST_INV_ARG_6);

    if (rf)
    {
        /*rf must be aware of sp and the current transform*/
        /* res = self.rf(sigma, pivot) */
        YAFL_TRY(status, rf(self, res, sigma, pivot));
    }
    else
    {
        yaflInt j;
        for (j = 0; j < sz; j++)
        {
            res[j] = sigma[j] - pivot[j];
        }
    }
    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _unscented_mean(yaflUKFBaseSt * self, \
                                           yaflInt    sz,        \
                                           yaflFloat * mean,     \
                                           yaflInt np,           \
                                           yaflFloat * sigmas,   \
                                           yaflKalmanFuncP    mf)
{
    yaflStatusEn status = YAFL_ST_OK;
    /*Args get checked later, or don't need to be checked*/
    if (mf)
    {
        /*mf must be aware of the current transform details...*/
        YAFL_CHECK(self,   YAFL_ST_INV_ARG_1);
        YAFL_CHECK(_wm,    YAFL_ST_INV_ARG_1);
        YAFL_CHECK(mean,   YAFL_ST_INV_ARG_3);
        YAFL_CHECK(sigmas, YAFL_ST_INV_ARG_5);
        YAFL_TRY(status, mf(_kalman_self, mean, sigmas));
    }
    else
    {
        YAFL_TRY(status, yafl_math_set_vtm(np, sz, mean, _wm, sigmas));
    }
    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _unscented_update(yaflUKFBaseSt * self,  \
                                             yaflInt sz,            \
                                             yaflFloat * u,         \
                                             yaflFloat * d,         \
                                             yaflFloat * pivot,     \
                                             yaflInt np,            \
                                             yaflFloat * sigmas,    \
                                             yaflFloat * sp,        \
                                             yaflKalmanResFuncP rf, \
                                             yaflFloat n)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt i;
    yaflFloat * wc;

    YAFL_CHECK(self,       YAFL_ST_INV_ARG_1);

    wc = _wc;
    YAFL_CHECK(wc,       YAFL_ST_INV_ARG_1);

    /*Other args get checked later, or don't need to be checked*/
    for (i = 0; i < np; i++)
    {
        yaflFloat wci;

        YAFL_TRY(status, _compute_res(_kalman_self, sz, rf, sp, sigmas + sz * i, pivot));
        /*Update res_u and res_d*/
        /*wc should be sorted in descending order*/
        wci = wc[i] * n;
        if (wci >= 0.0)
        {
            YAFL_TRY(status, yafl_math_udu_up(sz, u, d, wci, sp));
        }
        else
        {
            YAFL_TRY(status, yafl_math_udu_down(sz, u, d, -wci, sp));
        }
    }
    return status;
}

/*---------------------------------------------------------------------------*/
static yaflStatusEn _unscented_transform(yaflUKFBaseSt * self,   \
        yaflInt    res_sz,      \
        yaflFloat * res_v,      \
        yaflFloat * res_u,      \
        yaflFloat * res_d,      \
        yaflFloat * sp,         \
        yaflFloat * sigmas,     \
        yaflFloat * noise_u,    \
        yaflFloat * noise_d,    \
        yaflKalmanFuncP    mf,  \
        yaflKalmanResFuncP rf)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt np;
    yaflInt i;
    yaflUKFSigmaSt * sp_info;

    YAFL_CHECK(self,       YAFL_ST_INV_ARG_1);

    YAFL_CHECK(res_sz > 0, YAFL_ST_INV_ARG_2);
    YAFL_CHECK(res_u,      YAFL_ST_INV_ARG_4);
    YAFL_CHECK(res_d,      YAFL_ST_INV_ARG_5);
    YAFL_CHECK(sp,         YAFL_ST_INV_ARG_6);

    if (noise_u)
    {
        YAFL_CHECK(noise_d, YAFL_ST_INV_ARG_9);
    }

    YAFL_CHECK(_wm, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_wc, YAFL_ST_INV_ARG_1);

    sp_info = self->sp_info;
    YAFL_CHECK(sp_info, YAFL_ST_INV_ARG_1);

    np = sp_info->np;
    YAFL_CHECK(np > 1, YAFL_ST_INV_ARG_1);

    YAFL_TRY(status, _unscented_mean(self, res_sz, res_v, np, sigmas, mf));

    if (noise_u)
    {
        YAFL_CHECK(noise_d, YAFL_ST_INV_ARG_8);
        /*res_u, res_d = noise_u.copy(), noise_d.copy()*/
        memcpy((void *)res_u, (void *)noise_u, (res_sz * (res_sz - 1)) / 2 * sizeof(yaflFloat));
        memcpy((void *)res_d, (void *)noise_d, res_sz * sizeof(yaflFloat));
    }
    else
    {
        /* Zero res_u and res_d */
        for (i = res_sz - 1; i >= 0; i--)
        {
            res_d[i] = 0.0;
        }
        for (i = ((res_sz * (res_sz - 1)) / 2) - 1; i >= 0; i--)
        {
            res_u[i] = 0.0;
        }
    }

    YAFL_TRY(status, _unscented_update(self, res_sz, res_u, res_d, res_v, \
                                       np, sigmas, sp, rf, 1.0));
    return status;
}

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ukf_base_predict(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;

    yaflInt np;
    yaflInt nx;
    yaflInt i;
    yaflUKFSigmaSt * sp_info;

    /*Check some params and generate sigma points*/
    YAFL_TRY(status, yafl_ukf_gen_sigmas(self)); /*Self is checked here*/

    YAFL_CHECK(_unx, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_sigmas_x, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(self->sp_info, YAFL_ST_INV_ARG_1);
    sp_info = self->sp_info;

    np = sp_info->np;
    YAFL_CHECK(np > 1, YAFL_ST_INV_ARG_1);

    nx = _unx;
    YAFL_CHECK(nx > 0, YAFL_ST_INV_ARG_1);

    /*Compute process sigmas*/
    if (_ufx)
    {
        for (i = 0; i < np; i++)
        {
            yaflFloat * sigmai;
            sigmai = _sigmas_x + nx * i;
            YAFL_TRY(status, _ufx(_kalman_self, sigmai, sigmai));
        }
    }

    /*Predict x, Up, Dp*/
    YAFL_TRY(status, \
             _unscented_transform(self, nx, _ux, _uup, _udp, _sx, _sigmas_x, \
                                  _uuq, _udq, _xmf, _xrf));
    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _yafl_ukf_compute_sigmas_z_and_zp(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt i;
    yaflInt nx;
    yaflInt nz;
    yaflInt np;
    yaflUKFSigmaSt * sp_info; /*Sigma point generator info*/

    YAFL_CHECK(self,      YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_sigmas_x, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_sigmas_z, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_zp,       YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_wm,       YAFL_ST_INV_ARG_1);

    nx = _unx;
    YAFL_CHECK(nx, YAFL_ST_INV_ARG_1);

    nz = _unz;
    YAFL_CHECK(nz, YAFL_ST_INV_ARG_1);

    sp_info = self->sp_info;
    YAFL_CHECK(sp_info, YAFL_ST_INV_ARG_1);

    np = sp_info->np;
    YAFL_CHECK(np > 1, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_uhx, YAFL_ST_INV_ARG_1);

    for (i = 0; i < np; i++)
    {
        YAFL_TRY(status, _uhx(_kalman_self, _sigmas_z + nz * i, \
                              _sigmas_x + nx * i));
    }

    YAFL_TRY(status, _unscented_mean(self, nz, _zp, np, _sigmas_z, _zmf));

    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _yafl_ukf_compute_residual(yaflUKFBaseSt * self, \
                                                      yaflFloat * z)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat rff;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);

    rff = _kalman_self->rff;
    if (rff  > 0.0)
    {
        YAFL_CHECK(rff < 1.0, YAFL_ST_INV_ARG_1);
        /*Compute residual if needed (need to generate new sigmas)*/
        YAFL_TRY(status, yafl_ukf_gen_sigmas(self));
        YAFL_TRY(status, _yafl_ukf_compute_sigmas_z_and_zp(self));
        YAFL_TRY(status, _compute_res(_kalman_self, _unz, _uzrf, _uy, z, _zp));
    }
    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _yafl_ukf_compute_pzx(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;

    yaflUKFSigmaSt * sp_info;
    yaflKalmanResFuncP zrf;
    yaflKalmanResFuncP xrf;

    yaflFloat * sigmas_x;
    yaflFloat * sigmas_z;
    yaflFloat * pzx;
    yaflFloat * sx;
    yaflFloat * sz;
    yaflFloat * x;
    yaflFloat * zp;
    yaflFloat * wc;

    yaflInt np;
    yaflInt nz;
    yaflInt nx;
    yaflInt i;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);

    nx = _unx;
    YAFL_CHECK(nx, YAFL_ST_INV_ARG_1);

    nz = _unz;
    YAFL_CHECK(nz, YAFL_ST_INV_ARG_1);

    sp_info = self->sp_info;
    YAFL_CHECK(sp_info, YAFL_ST_INV_ARG_1);

    np = sp_info->np;
    YAFL_CHECK(np > 1, YAFL_ST_INV_ARG_1);

    wc = _wc;
    YAFL_CHECK(wc, YAFL_ST_INV_ARG_1);

    /*Will be checked by _compute_res*/
    zrf      = _uzrf;
    xrf      = _xrf;
    pzx      = _pzx;
    sx       = _sx;
    sz       = _sz;
    x        = _ux;
    zp       = _zp;

    sigmas_x = _sigmas_x;
    sigmas_z = _sigmas_z;

    /* Compute Pzx */
    YAFL_TRY(status, _compute_res(_kalman_self, nz, zrf, sz, sigmas_z, zp));
    YAFL_TRY(status, _compute_res(_kalman_self, nx, xrf, sx, sigmas_x, x));
    YAFL_TRY(status, yafl_math_set_vvtxn(nz, nx, pzx, sz, sx, wc[0]));

    for (i = 1; i < np; i++)
    {
        YAFL_TRY(status, _compute_res(_kalman_self, nz, zrf, sz, sigmas_z + nz * i, zp));
        YAFL_TRY(status, _compute_res(_kalman_self, nx, xrf, sx, sigmas_x + nx * i, x));
        YAFL_TRY(status, yafl_math_add_vvtxn(nz, nx, pzx, sz, sx, wc[i]));
    }
    return status;
}

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ukf_base_update(yaflUKFBaseSt * self, yaflFloat * z, \
                                  yaflKalmanScalarUpdateP scalar_update)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt np;
    yaflInt nx;
    yaflInt nz;
    yaflInt i;
    yaflUKFSigmaSt * sp_info; /*Sigma point generator info*/

    YAFL_CHECK(scalar_update, YAFL_ST_INV_ARG_3);

    YAFL_CHECK(self,     YAFL_ST_INV_ARG_1);

    nx = _unx;
    YAFL_CHECK(nx, YAFL_ST_INV_ARG_1);

    nz = _unz;
    YAFL_CHECK(nz, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_pzx, YAFL_ST_INV_ARG_1);

    sp_info = self->sp_info;
    YAFL_CHECK(sp_info, YAFL_ST_INV_ARG_1);

    np = sp_info->np;
    YAFL_CHECK(np > 1, YAFL_ST_INV_ARG_1);

    /* Compute measurement sigmas and zp*/
    YAFL_TRY(status, _yafl_ukf_compute_sigmas_z_and_zp(self));

    /* Compute Pzx = H.dot(Up.dot(Dp.dot(Up.T))*/
    YAFL_TRY(status, _yafl_ukf_compute_pzx(self));
    /*Now _sigmas_z and _sigmas_x may be spoiled*/

    /*Compute H.dot(Up)*/
    for (i = 0; i < nz; i++)
    {
        yaflFloat * h;
        h = _pzx + nx * i;
        YAFL_TRY(status, yafl_math_ruv    (nx, h, _uup));
        YAFL_TRY(status, YAFL_MATH_SET_RDV(nx, h, _udp, h));
    }

    /*Update R if needed!*/
    if (_kalman_self->rff > 0.0)
    {
        yaflStatusEn status_r = status; /*For quiet regularization*/
        yaflFloat * d;
        yaflFloat * w;
        /*
        WARNING:

        1. _sigmas_x and _sigmas_z must be parts of the larger
            memory pool which has minimum size of (nx + nz) * (nz + 1)
            or np * (nx + nz) where np is number of sigma points

        2. _sigmas_x must be at start of this pool
        */
        d = _sigmas_x;
        w = _sigmas_x + nx + nz;

        /*Check input data*/
        YAFL_CHECK(_udr, YAFL_ST_INV_ARG_1);
        /*Start R update*/
        YAFL_TRY(status_r, yafl_math_bset_m(nx + nz, w, nz, nx, _pzx));
        /*Now W = (HUp|***) */
        YAFL_TRY(status_r, \
                 _yafl_r_update(nx, nz, _kalman_self->rff, _udp, \
                                _uur, _udr, w, d, _uy));
        /*Updated R, now _uy may be spoiled*/
    }

    /*R update call back*/
    if (_urcb)
    {
        _urcb(_kalman_self);
    }

    /*Compute H*/
    for (i = 0; i < nz; i++)
    {
        YAFL_TRY(status, yafl_math_rutv(nx, _pzx + nx * i, _uup));
    }
    /*Now _pzx = H*/

    /*Compute innovation*/
    YAFL_TRY(status, _compute_res(_kalman_self, nz, _uzrf, _uy, z, _zp));

    /* Decorrelate measurements*/
    YAFL_TRY(status, yafl_math_ruv(nz,      _uy, _uur));
    YAFL_TRY(status, yafl_math_rum(nz, nx, _pzx, _uur));

    /*Start log likelihood computation*/
    *_kalman_self->l = nz * YAFL_L2PI;

    /*Now we can do scalar updates*/
    for (i = 0; i < nz; i++)
    {
        YAFL_TRY(status, scalar_update(_kalman_self, i));
    }

    /*Finalize log likelihood value*/
    *_kalman_self->l *= -0.5;

    YAFL_TRY(status, _yafl_ukf_compute_residual(self, z));
    return status;
}

/*=============================================================================
                                 Bierman UKF
=============================================================================*/
#define _ukf_self ((yaflUKFBaseSt *)self)

#define _upzx (_ukf_self->Pzx)
#define _usx  (_ukf_self->Sx)

#define YAFL_UKF_BIERMAN_SELF_INTERNALS_CHECKS() \
do {                                             \
    YAFL_CHECK(_nx > 0, YAFL_ST_INV_ARG_1);      \
    YAFL_CHECK(_up,     YAFL_ST_INV_ARG_1);      \
    YAFL_CHECK(_dp,     YAFL_ST_INV_ARG_1);      \
    YAFL_CHECK(_upzx,   YAFL_ST_INV_ARG_1);      \
    YAFL_CHECK(_usx,    YAFL_ST_INV_ARG_1);      \
    YAFL_CHECK(_y,      YAFL_ST_INV_ARG_1);      \
    YAFL_CHECK(_dr,     YAFL_ST_INV_ARG_1);      \
} while (0)

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ukf_bierman_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt nx;
    yaflFloat * v;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();

    nx = _nx;
    YAFL_UKF_BIERMAN_SELF_INTERNALS_CHECKS();

    v = _upzx + nx * i;/*h is stored in v*/

#   define f _usx
    /* f = h.T.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(nx, f, v, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(nx, v, _dp, f));

    YAFL_TRY(status, \
             _bierman_update_body(nx, _x, _up, _dp, f, v, _dr[i], _y[i], \
                                  1.0, 1.0, _ul));
#   undef f
    return status;
}

/*=============================================================================
                             Adaptive Bierman UKF
=============================================================================*/
yaflStatusEn yafl_ukf_adaptive_bierman_update_scalar(yaflKalmanBaseSt * self, yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt nx;
    yaflFloat * v;
    yaflFloat nu;
    yaflFloat r;
    yaflFloat ac;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();
    nx = _nx;
    YAFL_UKF_BIERMAN_SELF_INTERNALS_CHECKS();

    v = _upzx + nx * i;/*h is stored in v*/

#   define f _usx
    /* f = h.T.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(nx, f, v, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(nx, v, _dp, f));

    nu = self->y[i];
    r  = self->Dr[i];

    YAFL_TRY(status, \
             _adaptive_correction(nx, &ac, 0, f, v, r, nu, 1.0, \
                                  ((yaflUKFAdaptivedSt *)self)->chi2));

    YAFL_TRY(status, \
             _bierman_update_body(nx, _x, _up, _dp, f, v, r, _y[i], ac, 1.0, _ul));
#   undef f
    return status;
}

/*=============================================================================
                           Robust Bierman UKF
=============================================================================*/
#define YAFL_UKF_SCALAR_ROBUSTIFY(_self, _gdot_res, _nu, _r)   \
    _scalar_robustify(self,                            \
                      ((yaflUKFRobustSt *)self)->g,    \
                      ((yaflUKFRobustSt *)self)->gdot, \
                      _gdot_res, _nu, _r)

yaflStatusEn yafl_ukf_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                   yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt nx;
    yaflFloat gdot = 1.0;
    yaflFloat nu   = 0.0;
    yaflFloat r;
    yaflFloat * v;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();

    nx = _nx;
    YAFL_UKF_BIERMAN_SELF_INTERNALS_CHECKS();

    r  = _dr[i];
    nu = _y[i];

    YAFL_TRY(status, YAFL_UKF_SCALAR_ROBUSTIFY(self, &gdot, &nu, r));

    v = _upzx + nx * i;/*h is stored in v*/

#   define f _usx
    /* f = h.T.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(nx, f, v, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(nx, v, _dp, f));

    YAFL_TRY(status, \
             _bierman_update_body(nx, _x, _up, _dp, f, v, r, nu, \
                                  1.0, gdot, _ul));
#   undef f
    return status;
}

/*=============================================================================
                         Adaptive robust Bierman UKF
=============================================================================*/
yaflStatusEn \
    yafl_ukf_adaptive_robust_bierman_update_scalar(yaflKalmanBaseSt * self, \
                                                   yaflInt i)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt nx;
    yaflFloat gdot = 1.0;
    yaflFloat ac   = 1.0;
    yaflFloat nu   = 0.0;
    yaflFloat r;
    yaflFloat * v;

    YAFL_SCALAR_UPDATE_ARGS_CHECKS();

    nx = _nx;
    YAFL_UKF_BIERMAN_SELF_INTERNALS_CHECKS();

    r  = _dr[i];
    nu = _y[i];

    YAFL_TRY(status, YAFL_UKF_SCALAR_ROBUSTIFY(self, &gdot, &nu, r));

    v = _upzx + nx * i;/*h is stored in v*/

#   define f _usx
    /* f = h.T.dot(Up) */
    YAFL_TRY(status, yafl_math_set_vtu(nx, f, v, _up));

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAFL_TRY(status, YAFL_MATH_SET_DV(nx, v, _dp, f));


    YAFL_TRY(status, \
             _adaptive_correction(nx, &ac, 0, f, v, r, nu, gdot, \
                                  ((yaflUKFAdaptiveRobustSt *)self)->chi2));

    YAFL_TRY(status, \
             _bierman_update_body(nx, _x, _up, _dp, f, v, r, nu, ac, gdot, _ul));
#   undef f
    return status;
}

/*=============================================================================
            Full UKF, not a sequential square root version of UKF
=============================================================================*/
#define _uus (((yaflUKFSt *)self)->Us)
#define _uds (((yaflUKFSt *)self)->Ds)

static inline yaflStatusEn _yafl_ukf_compute_ms_zp_s(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;

    yaflUKFSigmaSt * sp_info;

    yaflFloat * sigmas_x;
    yaflFloat * sigmas_z;
    yaflFloat * dr;
    yaflFloat * ur;

    yaflFloat rff;

    yaflInt np;
    yaflInt nz;
    yaflInt nx;
    yaflInt i;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);

    nx = _unx;
    YAFL_CHECK(nx, YAFL_ST_INV_ARG_1);

    nz = _unz;
    YAFL_CHECK(nz, YAFL_ST_INV_ARG_1);

    sp_info = self->sp_info;
    YAFL_CHECK(sp_info, YAFL_ST_INV_ARG_1);

    np = sp_info->np;
    YAFL_CHECK(np > 1, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_uhx, YAFL_ST_INV_ARG_1);

    sigmas_x = _sigmas_x;
    YAFL_CHECK(sigmas_x, YAFL_ST_INV_ARG_1);

    sigmas_z = _sigmas_z;
    YAFL_CHECK(sigmas_z, YAFL_ST_INV_ARG_1);


    /* Compute measurement sigmas */
    for (i = 0; i < np; i++)
    {
        YAFL_TRY(status, _uhx(_kalman_self, sigmas_z + nz * i, sigmas_x + nx * i));
    }

    ur  = _uur;
    dr  = _udr;

    rff = _kalman_self->rff;
    if (rff > 0.0)
    {
        /* Compute zp, update Ur, Dr, compute Us, Ds */
        yaflFloat tmp;
        yaflStatusEn r_status = YAFL_ST_OK; /*For silent R regularization*/

        YAFL_CHECK(rff < 1.0,  YAFL_ST_INV_ARG_1);
        YAFL_CHECK(dr,         YAFL_ST_INV_ARG_1);

        YAFL_TRY(status, _unscented_mean(self, nz, _zp, np, sigmas_z, _zmf));

        /*Update Ur, Dr*/
        /* R *= 1.0- rff */
        tmp = 1.0 - rff;
        for (i=0; i < nz; i++)
        {
            dr[i] *= tmp;
        }

        /* R += rff * y.dot(y.T) */
        YAFL_TRY(r_status, yafl_math_udu_up(nz, ur, dr, rff, _uy));

        /* R += rff * Pzz*/
        YAFL_TRY(r_status, _unscented_update(self, nz, ur, dr, _zp, \
                                             np, sigmas_z, _sz, _uzrf, rff));

        /*R update call back*/
        if (_urcb)
        {
            _urcb(_kalman_self);
        }

        /*Compute Us, Ds*/
        memcpy((void *)_uus, (void *)ur, (nz * (nz - 1)) / 2 * sizeof(yaflFloat));
        memcpy((void *)_uds, (void *)dr, nz * sizeof(yaflFloat));

        YAFL_TRY(status, _unscented_update(self, nz, _uus, _uds, _zp, \
                                       np, sigmas_z, _sz, _uzrf, 1.0));
    }
    else
    {
        /* Compute zp, Us, Ds */
        YAFL_TRY(status, \
                 _unscented_transform(self, nz, _zp, _uus, _uds, \
                                      _sz, sigmas_z, ur, dr, _zmf, _uzrf));
    }
    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn yafl_ukf_update_epilogue(yaflUKFBaseSt * self, yaflFloat * z)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt nz;
    yaflInt nx;
    yaflInt i;

    yaflFloat * y;
    yaflFloat * ds;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);

    y = _uy;
    YAFL_CHECK(y,  YAFL_ST_INV_ARG_1);

    nx = _unx;
    YAFL_CHECK(nx, YAFL_ST_INV_ARG_1);

    nz = _unz;
    YAFL_CHECK(nz, YAFL_ST_INV_ARG_1);

    ds = _uds;
    YAFL_CHECK(ds, YAFL_ST_INV_ARG_1);

    /* Compute Pzx */
    YAFL_TRY(status, _yafl_ukf_compute_pzx(self));

    /* Decorrelate measurements part 2*/
    YAFL_TRY(status, yafl_math_rum(nz, nx, _pzx, _uus));

    /*Start log likelihood computation*/
    *_kalman_self->l = nz * YAFL_L2PI;

    /*Now we can do scalar updates*/
    for (i = 0; i < nz; i++)
    {
        yaflFloat * pzxi;
        pzxi = _pzx + nx * i;
        /*
        self.x += K * y[i]

        K * y[i] = Pzx[i].T / ds[i] * y[i] = Pzx[i].T * (y[i] / ds[i])

        self.x += Pzx[i].T * (y[i] / ds[i])
        */
        YAFL_TRY(status, yafl_math_add_vxn(nx, _ux, pzxi, y[i] / ds[i]));

        /*
        P -= K.dot(S.dot(K.T))
        K.dot(S.dot(K.T)) = (Pzx[i].T / ds[i] * ds[i]).outer(Pzx[i] / ds[i]))
        K.dot(S.dot(K.T)) = (Pzx[i].T).outer(Pzx[i]) / ds[i]
        P -= (Pzx[i].T).outer(Pzx[i]) / ds[i]
        Up, Dp = udu(P)
        */
        YAFL_TRY(status, yafl_math_udu_down(nx, _uup, _udp, 1.0 / ds[i], pzxi));

        /*Compute log likelihood*/
        *_kalman_self->l += y[i] * (y[i] / ds[i]) + YAFL_LOG(ds[i]);
    }

    /*Finalize log likelihood value*/
    *_kalman_self->l *= -0.5;

    YAFL_TRY(status, _yafl_ukf_compute_residual(self, z));
    return status;
}

/*---------------------------------------------------------------------------*/
yaflStatusEn yafl_ukf_update(yaflUKFBaseSt * self, yaflFloat * z)
{
    yaflStatusEn status = YAFL_ST_OK;

    /* Compute measurement sigmas, zp, Us, Ds*/
    YAFL_TRY(status, _yafl_ukf_compute_ms_zp_s(self));

    /*Compute innovation*/
    YAFL_TRY(status, _compute_res(_kalman_self, _unz, _uzrf, _uy, z, _zp));

    /* Decorrelate measurements part 1*/
    YAFL_TRY(status, yafl_math_ruv(_unz, _uy, _uus));

    /*
    Decorrelate measurements part 2, compute Pzx,
    update P and x, compute residual if needed
    */
    YAFL_TRY(status, yafl_ukf_update_epilogue(self, z));

    return status;
}

/*=============================================================================
       Full adaptive UKF, not a sequential square root version of UKF
=============================================================================*/
static inline yaflStatusEn _ukf_compute_md(yaflUKFBaseSt * self, yaflFloat * z, yaflFloat * md)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt i;
    yaflInt nz;
    yaflFloat dist;
    yaflFloat * y;
    yaflFloat * ds;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_uus, YAFL_ST_INV_ARG_1);
    ds = _uds;
    YAFL_CHECK(ds,   YAFL_ST_INV_ARG_1);

    nz = _unz;
    YAFL_CHECK(nz, YAFL_ST_INV_ARG_1);

    y = _uy;
    YAFL_CHECK(y,    YAFL_ST_INV_ARG_1);
    YAFL_CHECK(z,    YAFL_ST_INV_ARG_2);

    YAFL_CHECK(md,   YAFL_ST_INV_ARG_3);

    YAFL_TRY(status, _compute_res(_kalman_self, nz, _uzrf, y, z, _zp));
    YAFL_TRY(status, yafl_math_ruv(nz, y, _uus));

    dist = 0;
    for (i = 0; i < nz; i++)
    {
        /* Since ds is positive definite we won't check values. */
        dist += y[i] * y[i] / ds[i];
    }

    *md = dist;

    return status;
}

yaflStatusEn yafl_ukf_adaptive_update(yaflUKFBaseSt * self, yaflFloat * z)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFloat delta;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);

    /* Compute measurement sigmas, zp, Us, Ds*/
    YAFL_TRY(status, _yafl_ukf_compute_ms_zp_s(self));

    /* Compute md, innovation, decorrelate measurements part 1*/
    YAFL_TRY(status, _ukf_compute_md(self, z, &delta));

#   define _chi2 (((yaflUKFFullAdapiveSt *)self)->chi2)
    if (delta > _chi2)
    {
        yaflInt nz;

        /* Adaptive correction */
        yaflFloat ac;

        nz = _unz;
        YAFL_CHECK(nz, YAFL_ST_INV_ARG_1);

        /* Compute correction factor, we don't need old _zp and _sigmas_z now*/
        YAFL_TRY(status, \
                 _unscented_transform(self, nz, _zp, _uus, _uds, _sz, _sigmas_z, \
                                      0, 0, _zmf, _uzrf));

        YAFL_TRY(status, _ukf_compute_md(self, z, &ac));
        ac *= 1.0 / _chi2 - 1.0 / delta;

        /* Correct _udp*/
        YAFL_TRY(status, yafl_math_set_vxn(_unx, _udp, _udp, 1.0 + ac));

        /* Generate new sigmas */
        YAFL_TRY(status, yafl_ukf_gen_sigmas(self));

        /* Now restart update with new sigmas */
        YAFL_TRY(status, _yafl_ukf_compute_ms_zp_s(self));

        /*  Recompute innovation*/
        YAFL_TRY(status, _compute_res(_kalman_self, nz, _uzrf, _uy, z, _zp));

        /*  Decorrelate measurements part 1*/
        YAFL_TRY(status, yafl_math_ruv(nz, _uy, _uus));
    }
#   undef _chi2

    /*
    Decorrelate measurements part 2, compute Pzx,
    update P and x, compute residual if needed
    */
    YAFL_TRY(status, yafl_ukf_update_epilogue(self, z));

    return status;
}

/*=============================================================================
                    Van der Merwe sigma points generator
=============================================================================*/
static yaflStatusEn _merwe_compute_weights(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt np;
    yaflInt nx;
    yaflInt i;
    yaflFloat * wc;
    yaflFloat * wm;
    yaflUKFSigmaSt * sp_info;
    yaflFloat lambda;
    yaflFloat alpha;
    yaflFloat c;
    yaflFloat d;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);
    nx = _unx;
    YAFL_CHECK(_unx, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_wm, YAFL_ST_INV_ARG_1);
    wm = _wm;

    YAFL_CHECK(_wc, YAFL_ST_INV_ARG_1);
    wc = _wc;

    YAFL_CHECK(self->sp_info, YAFL_ST_INV_ARG_1);
    sp_info = self->sp_info;

    YAFL_CHECK(sp_info->np, YAFL_ST_INV_ARG_1);
    np = sp_info->np - 1;    /*Achtung!!!*/

    alpha = ((yaflUKFMerweSt *)sp_info)->alpha;
    alpha *= alpha;
#define _alpha2 alpha

    lambda = _alpha2 * (nx + ((yaflUKFMerweSt *)sp_info)->kappa) - nx;

    d = lambda / (nx + lambda);
    wm[np] = d;
    wc[np] = d + (1.0 - _alpha2 + ((yaflUKFMerweSt *)sp_info)->beta);

    c = 0.5 / (nx + lambda);
    for (i = np - 1; i >= 0; i--)
    {
        wm[i] = c;
        wc[i] = c;
    }
#undef _alpha2
    return status;
}

/*---------------------------------------------------------------------------*/
static inline yaflStatusEn _add_delta(yaflUKFBaseSt * self,                  \
                                      yaflUKFSigmaAddP addf, yaflInt sz,     \
                                      yaflFloat * delta, yaflFloat * pivot,  \
                                      yaflFloat mult)
{
    yaflStatusEn status = YAFL_ST_OK;
    if (addf)
    {
        /* addf must be aware of self internal structure*/
        /* delta = self.addf(delta, pivot, mult) */
        YAFL_TRY(status, addf(self, delta, pivot, mult));
    }
    else
    {
        yaflInt j;
        for (j = 0; j < sz; j++)
        {
            delta[j] = pivot[j] + mult * delta[j];
        }
    }
    return status;
}

/*---------------------------------------------------------------------------*/
static yaflStatusEn _merwe_generate_points(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt nx;
    yaflInt i;
    yaflFloat * x;
    yaflFloat * sigmas_x;
    yaflFloat * dp;
    yaflUKFSigmaSt * sp_info;
    yaflUKFSigmaAddP addf;
    yaflFloat lambda_p_n;
    yaflFloat alpha;

    YAFL_CHECK(self,     YAFL_ST_INV_ARG_1);
    nx = _unx;
    YAFL_CHECK(nx, YAFL_ST_INV_ARG_1);

    x = _ux;
    YAFL_CHECK(x, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_uup, YAFL_ST_INV_ARG_1);

    dp = _udp;
    YAFL_CHECK(dp, YAFL_ST_INV_ARG_1);

    sigmas_x = _sigmas_x;
    YAFL_CHECK(sigmas_x, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(self->sp_info, YAFL_ST_INV_ARG_1);
    sp_info = self->sp_info;

    alpha = ((yaflUKFMerweSt *)sp_info)->alpha;
    lambda_p_n = alpha * alpha * (nx + ((yaflUKFMerweSt *)sp_info)->kappa);

    YAFL_TRY(status, yafl_math_bset_ut(nx, sigmas_x, nx, _uup));
    memcpy((void *)(sigmas_x + nx * nx), (void *)sigmas_x, \
           nx * nx * sizeof(yaflFloat));

    addf = sp_info->addf;
    for (i = 0; i < nx; i++)
    {
        yaflFloat mult;
        mult = YAFL_SQRT(dp[i] * lambda_p_n);
        YAFL_TRY(status, \
                 _add_delta(self, addf, nx, sigmas_x + nx * i,        x,   mult));
        YAFL_TRY(status, \
                 _add_delta(self, addf, nx, sigmas_x + nx * (nx + i), x, - mult));
    }
    memcpy((void *)(sigmas_x + 2 * nx * nx), (void *)x, nx * sizeof(yaflFloat));
    return status;
}

/*---------------------------------------------------------------------------*/
const yaflUKFSigmaMethodsSt yafl_ukf_merwe_spm =
{
    .wf   = _merwe_compute_weights,
    .spgf = _merwe_generate_points
};

/*=============================================================================
                    Julier sigma points generator
=============================================================================*/
static yaflStatusEn _julier_compute_weights(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt np;
    yaflInt nx;
    yaflInt i;
    yaflFloat * wc;
    yaflFloat * wm;
    yaflUKFSigmaSt * sp_info;
    yaflFloat kappa;
    yaflFloat c;

    YAFL_CHECK(self, YAFL_ST_INV_ARG_1);
    nx = _unx;
    YAFL_CHECK(_unx, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_wm, YAFL_ST_INV_ARG_1);
    wm = _wm;

    YAFL_CHECK(_wc, YAFL_ST_INV_ARG_1);
    wc = _wc;

    YAFL_CHECK(self->sp_info, YAFL_ST_INV_ARG_1);
    sp_info = self->sp_info;

    YAFL_CHECK(sp_info->np, YAFL_ST_INV_ARG_1);
    np = sp_info->np - 1;    /*Achtung!!!*/

    kappa = ((yaflUKFJulierSt *)sp_info)->kappa;

    wm[np] = kappa / (nx + kappa);
    wc[np] = wm[np];

    c = 0.5 / (nx + kappa);
    for (i = np - 1; i >= 0; i--)
    {
        wm[i] = c;
        wc[i] = c;
    }
    return status;
}

/*---------------------------------------------------------------------------*/
static yaflStatusEn _julier_generate_points(yaflUKFBaseSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflInt nx;
    yaflInt i;
    yaflFloat * x;
    yaflFloat * sigmas_x;
    yaflFloat * dp;
    yaflUKFSigmaSt * sp_info;
    yaflUKFSigmaAddP addf;
    yaflFloat kappa_p_n;

    YAFL_CHECK(self,     YAFL_ST_INV_ARG_1);
    nx = _unx;
    YAFL_CHECK(nx, YAFL_ST_INV_ARG_1);

    x = _ux;
    YAFL_CHECK(x, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(_uup, YAFL_ST_INV_ARG_1);

    dp = _udp;
    YAFL_CHECK(dp, YAFL_ST_INV_ARG_1);

    sigmas_x = _sigmas_x;
    YAFL_CHECK(sigmas_x, YAFL_ST_INV_ARG_1);

    YAFL_CHECK(self->sp_info, YAFL_ST_INV_ARG_1);
    sp_info = self->sp_info;

    kappa_p_n = ((yaflUKFJulierSt *)sp_info)->kappa + nx;

    YAFL_TRY(status, yafl_math_bset_ut(nx, sigmas_x, nx, _uup));
    memcpy((void *)(sigmas_x + nx * nx), (void *)sigmas_x, \
           nx * nx * sizeof(yaflFloat));

    addf = sp_info->addf;
    for (i = 0; i < nx; i++)
    {
        yaflFloat mult;
        mult = YAFL_SQRT(dp[i] * kappa_p_n);
        YAFL_TRY(status, \
                 _add_delta(self, addf, nx, sigmas_x + nx * i,        x,   mult));
        YAFL_TRY(status, \
                 _add_delta(self, addf, nx, sigmas_x + nx * (nx + i), x, - mult));
    }
    memcpy((void *)(sigmas_x + 2 * nx * nx), (void *)x, nx * sizeof(yaflFloat));
    return status;
}

/*---------------------------------------------------------------------------*/
const yaflUKFSigmaMethodsSt yafl_ukf_julier_spm =
{
    .wf   = _julier_compute_weights,
    .spgf = _julier_generate_points
};

/*=============================================================================
                          Undef UKF stuff
=============================================================================*/
#undef _kalman_self


#undef _ufx
#undef _uhx
#undef _uzrf
#undef _urcb

#undef _ux
#undef _uy

#undef _uup
#undef _udp

#undef _uuq
#undef _udq

#undef _uur
#undef _udr

#undef _ul

#undef _unx
#undef _unz

#undef _xrf
#undef _xmf

#undef _zmf

#undef _zp

#undef _sx
#undef _sz

#undef _pzx

#undef _sigmas_x
#undef _sigmas_z

#undef _wm
#undef _wc

/*----------------------------------------------------------------------------*/
#undef _ukf_self

#undef _upzx
#undef _usx

/*----------------------------------------------------------------------------*/
#undef _uus
#undef _uds
/*=============================================================================
                          Undef Kalman filter stuff
=============================================================================*/
#undef _fx
#undef _hx
#undef _zrf
#undef _rcb

#undef _x
#undef _y

#undef _up
#undef _dp

#undef _uq
#undef _dq

#undef _ur
#undef _dr

#undef _l

#undef _nx
#undef _nz

yaflStatusEn yafl_imm_post_init(yaflIMMCBSt * self)
{
    yaflInt i;
    yaflInt nb;
    yaflInt nx = 0;
    yaflInt nz = 0;

    YAFL_CHECK(self,        YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->bank,  YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->mu,    YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->M,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->Up,    YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->Dp,    YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->x,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->cbar,  YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->omega, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->y,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->W,     YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->D,     YAFL_ST_INV_ARG_1);

    nb = self->Nb;
    YAFL_CHECK(nb > 1,     YAFL_ST_INV_ARG_1);

    for (i=0; i<nb; i++)
    {
        /*Bank item*/
        yaflFilterBankItemSt * bi = self->bank + i;

        YAFL_CHECK(bi,              YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter,      YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter->x,   YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter->Up,  YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter->Dp,  YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->predict,     YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->update,      YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->Us,          YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->Ds,          YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->Xs,          YAFL_ST_INV_ARG_1);

        if (YAFL_UNLIKELY(0 == i))
        {
            nx = bi->filter->Nx;
            nz = bi->filter->Nz;
        }

        YAFL_CHECK((bi->filter->Nx == nx), YAFL_ST_INV_ARG_1);
        YAFL_CHECK((bi->filter->Nz == nz), YAFL_ST_INV_ARG_1);
    }

    return YAFL_ST_OK;
}

yaflStatusEn yafl_imm_predict(yaflIMMCBSt * self)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFilterBankItemSt * bi;
    yaflFilterBankItemSt * bj;
    yaflInt i;
    yaflInt j;
    yaflInt k;
    yaflInt nb;
    yaflInt nbj;
    yaflInt nu;
    yaflInt nx;
    yaflInt nx2;
    yaflInt szu;
    yaflInt szx;

    YAFL_CHECK(self,               YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->bank,         YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->bank->filter, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->omega,        YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->y,            YAFL_ST_INV_ARG_1);

    nb = self->Nb;
    YAFL_CHECK(nb > 1,     YAFL_ST_INV_ARG_1);

    nx = self->bank->filter->Nx;
    nx2 = nx * 2;
    nu = (((nx - 1) * nx) / 2);
    szx = nx * sizeof(yaflFloat);
    szu = nu * sizeof(yaflFloat);

    YAFL_TRY(status, yafl_math_set_vtm(nb, nb, self->cbar, self->mu, self->M));
    for (j = 0; j < nb; j++)
    {
        yaflInt nbj = nb * j;
        bj = self->bank + j;

        YAFL_CHECK(bj,             YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bj->filter,     YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bj->filter->Up, YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bj->filter->Dp, YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bj->filter->x,  YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bj->Us,         YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bj->Ds,         YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bj->Xs,         YAFL_ST_INV_ARG_1);

        for (i = 0; i < nb; i++)
        {
            self->omega[nbj + i] = self->M[nbj + i] * self->mu[j] / self->cbar[i];
        }
    }

    /*Calculate mixed states and covariances*/
    for (i = 0; i < nb; i++)
    {
        bi = self->bank + i;

        /*Mixed state*/
        for (j = 0; j < nb; j++)
        {
            if (YAFL_UNLIKELY(0 == j))
            {
                /*XS[i] = omega[0, i] * f[0].x*/
                YAFL_TRY(status, yafl_math_set_vxn(nx, bi->Xs, self->bank[0].filter->x, self->omega[i]));
                continue;
            }
            /*XS[i] += omega[j, i] * f[j].x*/
            YAFL_TRY(status, yafl_math_add_vxn(nx, bi->Xs, self->bank[j].filter->x, self->omega[nb * j + i]));
        }

        /*Mixed covariance*/
        for (j = 0; j < nb; j++)
        {
            nbj = nb * j;
            bj = self->bank + j;

            if (YAFL_UNLIKELY(0 == j))
            {
                /*PS[i] = omega[0, i] * f[j].P*/
                memcpy((void *)(bi->Us), (void *)(bj->filter->Up), szu);
                YAFL_TRY(status, yafl_math_set_vxn(nx, bi->Ds, bj->filter->Dp, self->omega[i]));
            }
            else
            {
                /*Must use MWGS update here*/
                YAFL_TRY(status, YAFL_MATH_BSET_U(nx2, 0,  0, self->W, nx, bi->Us        ));
                /* Now W = (Us[i]|***) */
                YAFL_TRY(status, YAFL_MATH_BSET_U(nx2, 0, nx, self->W, nx, bj->filter->Up));
                /* Now W = (Us[i]|f[j].Up) */

                /* D = concatenate([Ds[i], omega[j, i] * f[j].Dp]) */
                memcpy((void *)(self->D),      (void *)(bi->Ds), szx);
                YAFL_TRY(status, yafl_math_set_vxn(nx, self->D + nx, bj->filter->Dp, self->omega[nbj + i]));

                /*PS[i] += omega[j, i] * f[j].P*/
                YAFL_TRY(status, yafl_math_mwgsu(nx, nx2, bi->Us, bi->Ds, self->W, self->D));
            }

            /*y = f[j].x - XS[i]*/
            for (k = 0; k < nx; k++)
            {
                self->y[k] = bj->filter->x[k] - bi->Xs[k];
            }

            /*PS[i] += omega[j, i] * outer(y.T, y)*/
            YAFL_TRY(status, yafl_math_udu_up(nx, bi->Us, bi->Ds, self->omega[nbj + i], self->y));
        }
    }

    for (i = 0; i < nb; i++)
    {
        bi = self->bank + i;

        /*Copy states from scratch memory before predict*/
        memcpy((void *)(bi->filter->Up), (void *)(bi->Us), szu);
        memcpy((void *)(bi->filter->Dp), (void *)(bi->Ds), szx);
        memcpy((void *)(bi->filter->x),  (void *)(bi->Xs), szx);

        /*Now we can predict*/
        YAFL_TRY(status, bi->predict(bi->filter));
    }

    return status;
}

yaflStatusEn yafl_imm_update(yaflIMMCBSt * self, yaflFloat * z)
{
    yaflStatusEn status = YAFL_ST_OK;
    yaflFilterBankItemSt * bi;

    yaflFloat s;

    yaflInt i;
    yaflInt j;
    yaflInt nb;
    yaflInt nu;
    yaflInt nx;
    yaflInt nx2;
    yaflInt szu;
    yaflInt szx;

    YAFL_CHECK(self,               YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->bank,         YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->bank->filter, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->Up,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->Dp,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->x,            YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->y,            YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->cbar,         YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->mu,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(self->D,            YAFL_ST_INV_ARG_1);

    nb = self->Nb;
    YAFL_CHECK(nb > 1,     YAFL_ST_INV_ARG_1);

    YAFL_CHECK(z,          YAFL_ST_INV_ARG_2);

    nx = self->bank->filter->Nx;
    nx2 = nx * 2;
    nu = (((nx - 1) * nx) / 2);
    szx = nx * sizeof(yaflFloat);
    szu = nu * sizeof(yaflFloat);

    s = 0;
    for (i = 0; i < nb; i++)
    {
        bi = self->bank + i;
        /*Check filters*/
        YAFL_CHECK(bi,             YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->update,     YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter,     YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter->Up, YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter->Dp, YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter->x,  YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->filter->l,  YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->Us,         YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->Ds,         YAFL_ST_INV_ARG_1);
        YAFL_CHECK(bi->Xs,         YAFL_ST_INV_ARG_1);

        /*Update filters*/
        YAFL_TRY(status, bi->update(bi->filter, z));

        /*Start mu update using filters log likelihoods*/
        yaflFloat li = YAFL_EXP(*self->bank[i].filter->l);
        if (li < YAFL_EPS)
        {
            li = YAFL_EPS;
            status |= YAFL_ST_MSK_REGULARIZED;
        }

        self->mu[i] = self->cbar[i] * li;
        s += self->mu[i];
    }
    /*Finalize mode probabilities*/
    for (i = 0; i < nb; i++)
    {
        self->mu[i] /= s;
    }

    /*Calculate updated state*/
    for (i = 0; i < nb; i++)
    {
        if (YAFL_UNLIKELY(0 == i))
        {
            /*x = f[0].x * mu[0]*/
            YAFL_TRY(status, yafl_math_set_vxn(nx, self->x, self->bank[0].filter->x, self->mu[0]));
            continue;
        }
        /*x += f[i].x * mu[i]*/
        YAFL_TRY(status, yafl_math_add_vxn(nx, self->x, self->bank[i].filter->x, self->mu[i]));
    }

    /*Calculate updated covariance*/
    for (i = 0; i < nb; i++)
    {
        bi = self->bank + i;

        if (YAFL_UNLIKELY(0 == i))
        {
            /*P = mu[0] * f[0].P*/
            memcpy((void *)(self->Up), (void *)(bi->filter->Up), szu);
            YAFL_TRY(status, yafl_math_set_vxn(nx, self->Dp, bi->filter->Dp, self->mu[0]));
        }
        else
        {
            /*Must use MWGS update here*/
            YAFL_TRY(status, YAFL_MATH_BSET_U(nx2, 0,  0, self->W, nx, self->Up      ));
            /* Now W = (Up|***) */
            YAFL_TRY(status, YAFL_MATH_BSET_U(nx2, 0, nx, self->W, nx, bi->filter->Up));
            /* Now W = (Up|f[i].Up) */

            /* D = concatenate([Dp, mu[i] * f[i].Dp]) */
            memcpy((void *)(self->D), (void *)(self->Dp), szx);
            YAFL_TRY(status, yafl_math_set_vxn(nx, self->D + nx, bi->filter->Dp, self->mu[i]));

            /*P += mu[i] * f[i].P*/
            YAFL_TRY(status, yafl_math_mwgsu(nx, nx2, self->Up, self->Dp, self->W, self->D));
        }

        /*y = f[i].x - x*/
        for (j = 0; j < nx; j++)
        {
            self->y[j] = bi->filter->x[j] - self->x[j];
        }

        /*P += mu[i] * outer(y.T, y)*/
        YAFL_TRY(status, yafl_math_udu_up(nx, self->Up, self->Dp, self->mu[i], self->y));
    }

    return status;
}

