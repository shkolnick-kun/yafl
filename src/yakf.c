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

#include "yakf.h"
#include <string.h>

/*=============================================================================
                                  Base UDEKF
=============================================================================*/
void yakf_base_predict(yakfBaseSt * self)
{
    yakfInt i;
    yakfInt nx;
    yakfInt nx2;
    yakfFloat * w;
    yakfFloat * d;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->Uq);
    YAKF_ASSERT(self->Dq);
    YAKF_ASSERT(self->W);
    YAKF_ASSERT(self->D);

    nx2 = nx * 2;
    w  = self->W;

    /*Default f(x) = x*/
    if (0 == self->f)
    {
        YAKF_ASSERT(0 == self->jf);

        for (i = 0; i < nx; i++)
        {
            yakfInt j;
            yakfInt nci;

            nci = nx2 * i;
            for (j = 0; j < nx; j++)
            {
                w[nci + j] = (i != j) ? 0.0 : 1.0;
            }
        }
    }
    else
    {
        /*Must have some Jacobian function*/
        YAKF_ASSERT(self->jf);

        self->f(self);  /* x = f(x_old, ...) */
        self->jf(self); /* Place F(x, ...)=df/dx to W  */
    }
    /* Now W = (F|***) */
    YAKFM_BSET_BU(nx2, 0, nx, w, nx, nx, nx2, 0, 0, w, self->Up);
    /* Now W = (F|FUp) */
    yakfm_bset_u(nx2, w, nx, self->Uq);
    /* Now W = (Uq|FUp) */

    /* D = concatenate([Dq, Dp]) */
    d = self->D;
    i = nx*sizeof(yakfFloat);
    memcpy((void *)       d, (void *)self->Dq, i);
    memcpy((void *)(d + nx), (void *)self->Dp, i);

    /* Up, Dp = MWGSU(w, d)*/
    yakfm_mwgsu(nx, nx2, self->Up, self->Dp, w, d);
}

void yakf_base_update(yakfBaseSt * self, yakfFloat * z, yakfScalarUpdateP scalar_update)
{
    yakfInt j;
    yakfInt nx;
    yakfInt nz;
    yakfFloat * y;
    yakfFloat * ur;
    yakfFloat * h;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Ur);

    YAKF_ASSERT(self->h);
    YAKF_ASSERT(self->jh);
    YAKF_ASSERT(self->H);

    nz = self->Nz;
    YAKF_ASSERT(nz > 0);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(z);

    y  = self->y;

    self->h(self);  /* self.y =  h(x,...) */
    self->jh(self); /* self.H = jh(x,...) */

    if (0 == self->zrf)
    {
        /*Default residual*/
        for (j = 0; j < nz; j++)
        {
            y[j] = z[j] - y[j];
        }
    }
    else
    {
        self->zrf(self, z); /* self.y = zrf(z, h(x,...)) */
    }

    /* Decorrelate measurement noise */
    ur = self->Ur;
    h  = self->H;
    yakfm_ruv(nz,     y, ur);
    yakfm_rum(nz, nx, h, ur);

    /* Do scalar updates */
    for (j = 0; j < nz; j++)
    {
        scalar_update(self, j);
    }
}
/*=============================================================================
                                Bierman filter
=============================================================================*/
static void _bierman_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt j;
    yakfInt k;
    yakfInt nx;
    yakfInt nxk;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat r;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->D);

    u = self->Up;
    d = self->Dp;

    h = self->H + nx * i;

    /* f = h.dot(Up) */
    f = self->D;
    yakfm_set_vtu(nx, f, h, u);

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#define v h /*Don't need h any more, use it to store v*/
    YAKFM_SET_DV(nx, v, d, f);

    r = self->Dr[i];
    for (k = 0, nxk = 0; k < nx; nxk += k++)
    {
        yakfFloat a;
        yakfFloat fk;
        yakfFloat vk;

        fk = f[k];
        vk = v[k];
        a = r + fk * vk;
        d[k] *= r / a;
#define p fk /*No need for separate p variable*/
        p = - fk / r;
        for (j = 0; j < k; j++)
        {
            yakfFloat ujk;
            yakfFloat vj;

            ujk = u[j + nxk];
            vj  = v[j];

            u[j + nxk] = ujk +   p * vj;
            v[j]       = vj  + ujk * vk;
        }
#undef  p /*Don't need p any more...*/
        r = a;
    }
    /*
    Now we must do:
    self.x += K * y[i]

    Since:
    r == a

    then we have:
    K == v / a == v / r

    and so:
    K * y[i] == (v / r) * y[i] == v / r * y[i] == v * (y[i] / r)

    Finally we get:
    self.x += v * (y[i] / r)
    */
    yakfm_add_vxn(nx, self->x, v, self->y[i] / r);
#undef v /*Don't nee v any more*/
}

void yakf_bierman_update(yakfBaseSt * self, yakfFloat * z)
{
    yakf_base_update(self, z, _bierman_scalar_update);
}

/*=============================================================================
                                Joseph filter
=============================================================================*/
static void _joseph_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt nx;
    yakfInt nx1;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat * v;
    yakfFloat * w;
    yakfFloat r;
    yakfFloat s;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->W);
    YAKF_ASSERT(self->D);

    nx1 = nx + 1;

    d = self->Dp;
    u = self->Up;

    h = self->H + nx * i;

    v = self->D;
    f = v + nx;

    w = self->W;

    /* f = h.dot(Up) */
    yakfm_set_vtu(nx, f, h, u);

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAKFM_SET_DV(nx, v, d, f);

    r = self->Dr[i];

    /* s = r + f.dot(v)*/
    s = r + yakfm_vtv(nx, f, v);

    /*K = Up.dot(v/s) = Up.dot(v)/s*/
#define K h /*Don't need h any more, use it to store K*/
    yakfm_set_uv(nx, K, u, v);
    yakfm_set_vrn(nx, K, K, s); /*May be used in place*/

    /*Set W and D*/
    yakfm_bset_vvt(nx1, w, nx, K, f);
    yakfm_bsub_u(nx1, w, nx, u);

    /* Now w is (Kf - Up|***) */
    YAKFM_BSET_V(nx1, 0, nx, w, nx, K);
    /* Now w is (Kf - Up|K) */
#define D v
    /* D = concatenate([Dp, np.array([r])]) */
    memcpy((void *)D, (void *)d, nx * sizeof(yakfFloat));
    D[nx] = r;

    /* Up, Dp = MWGSU(W, D)*/
    yakfm_mwgsu(nx, nx1, u, d, w, D);
#undef D
    /* self.x += K * y[i] */
    yakfm_add_vxn(nx, self->x, K, self->y[i]);
#undef K /*Don't nee K any more*/
}

void yakf_joseph_update(yakfBaseSt * self, yakfFloat * z)
{
    yakf_base_update(self, z, _joseph_scalar_update);
}

/*=============================================================================
                          Adaptive Bierman filter
=============================================================================*/
static void _adaptive_bierman_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt j;
    yakfInt k;
    yakfInt nx;
    yakfInt nxk;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat nu;
    yakfFloat r;
    yakfFloat c;
    yakfFloat s;
    yakfFloat ac;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->D);

    u = self->Up;
    d = self->Dp;

    h = self->H + nx * i;

    nu = self->y[i];
    r  = self->Dr[i];

    /* f = h.dot(Up) */
    f = self->D;
    yakfm_set_vtu(nx, f, h, u);

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#define v h /*Don't need h any more, use it to store v*/
    YAKFM_SET_DV(nx, v, d, f);

    /* s = r + f.dot(v)*/
    c = yakfm_vtv(nx, f, v);
    s = c + r;

    /* Divergence test */
    ac = (nu * (nu / (((yakfAdaptiveSt *)self)->chi2))) - s;
    if (ac > 0.0)
    {
        /*Adaptive correction factor*/
        ac = ac / c + 1.0;
    }
    else
    {
        ac = 1.0;
    }

    for (k = 0, nxk = 0; k < nx; nxk += k++)
    {
        yakfFloat a;
        yakfFloat fk;
        yakfFloat vk;

        fk = f[k];
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
            yakfFloat ujk;
            yakfFloat vj;

            ujk = u[j + nxk];
            vj  = v[j];

            u[j + nxk] = ujk +   p * vj;
            v[j]       = vj  + ujk * vk;
        }
#undef  p /*Don't need p any more...*/
        r = a;
    }
    /*
    Now we must do:
    self.x += K * nu

    Since:
    r == a

    then we have:
    K == v / a == v / r

    and so:
    K * nu == (v / r) * nu == v / r * nu == v * (nu / r)

    Finally we get:
    self.x += v * (nu / r)
    */
    yakfm_add_vxn(nx, self->x, v, nu / r);
#undef v /*Don't nee v any more*/
}

void yakf_adaptive_bierman_update(yakfAdaptiveSt * self, yakfFloat * z)
{
    yakf_base_update((yakfBaseSt *)self, z, _adaptive_bierman_scalar_update);
}

/*=============================================================================
                           Adaptive Joseph filter
=============================================================================*/
static void _adaptive_joseph_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt nx;
    yakfInt nx1;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat * v;
    yakfFloat * w;
    yakfFloat nu;
    yakfFloat r;
    yakfFloat c;
    yakfFloat s;
    yakfFloat ac;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->W);
    YAKF_ASSERT(self->D);

    nx1 = nx + 1;

    d = self->Dp;
    u = self->Up;

    h = self->H + nx * i;

    v = self->D;
    f = v + nx;

    w = self->W;

    nu = self->y[i];
    r  = self->Dr[i];

    /* f = h.dot(Up) */
    yakfm_set_vtu(nx, f, h, u);

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAKFM_SET_DV(nx, v, d, f);

    /* s = r + f.dot(v)*/
    c = yakfm_vtv(nx, f, v);
    s = c + r;

    /* Divergence test */
    ac = (nu * (nu / (((yakfAdaptiveSt *)self)->chi2))) - s;
    if (ac > 0.0)
    {
        /*Adaptive correction factor*/
        ac = ac / c + 1.0;

        /*Corrected s*/
        s  = ac * c + r;
    }
    else
    {
        ac = 1.0;
    }

    /* K = Up.dot(v * ac / s) = Up.dot(v) * (ac / s) */
#define K h /*Don't need h any more, use it to store K*/
#define D v
    yakfm_set_vxn(nx, v, v, ac / s); /*May be used in place*/
    yakfm_set_uv(nx, K, u, v);

    /*Set W and D*/
    yakfm_bset_vvt(nx1, w, nx, K, f);
    yakfm_bsub_u(nx1, w, nx, u);

    /* Now w is (Kf - Up|***) */
    YAKFM_BSET_V(nx1, 0, nx, w, nx, K);
    /* Now w is (Kf - Up|K) */

    /* D = concatenate([ac * Dp, np.array([r])]) */
    yakfm_set_vxn(nx, D, d, ac);
    D[nx] = r;

    /* Up, Dp = MWGSU(W, D)*/
    yakfm_mwgsu(nx, nx1, u, d, w, D);

    /* self.x += K * nu */
    yakfm_add_vxn(nx, self->x, K, nu);
#undef D /*Don't nee D any more*/
#undef K /*Don't nee K any more*/
}

void yakf_adaptive_joseph_update(yakfAdaptiveSt * self, yakfFloat * z)
{
    yakf_base_update((yakfBaseSt *)self, z, _adaptive_joseph_scalar_update);
}

/*=============================================================================
                                 WARNING!!!

             DO NOT USE THIS variant of Adaptive Joseph filter !!!

    It was implemented to show some flaws of the corresponding algorithm!
=============================================================================*/
static void _do_not_use_this_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt nx;
    yakfInt nx1;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat * v;
    yakfFloat * w;
    yakfFloat nu;
    yakfFloat r;
    yakfFloat c;
    yakfFloat s;
    yakfFloat ac;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->W);
    YAKF_ASSERT(self->D);

    nx1 = nx + 1;

    d = self->Dp;
    u = self->Up;

    h = self->H + nx * i;

    v = self->D;
    f = v + nx;

    w = self->W;

    nu = self->y[i];
    r  = self->Dr[i];

    /* f = h.dot(Up) */
    yakfm_set_vtu(nx, f, h, u);

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
    YAKFM_SET_DV(nx, v, d, f);

    /* s = r + f.dot(v)*/
    c = yakfm_vtv(nx, f, v);
    s = c + r;

    /* Divergence test */
    ac = (nu * (nu / (((yakfAdaptiveSt *)self)->chi2))) - s;
    if (ac > 0.0)
    {
        /*Adaptive correction with no limitations approach*/
        yakfm_set_uv(nx, f, u, v);
        yakfm_udu_up(nx, u, d, (ac / c) / c, f);

        /*Recompute f,v,s*/
        /* f = h.dot(Up) */
        yakfm_set_vtu(nx, f, h, u);

        /* v = f.dot(Dp).T = Dp.dot(f.T).T */
        YAKFM_SET_DV(nx, v, d, f);

        /* s = r + f.dot(v)*/
        s  = r + yakfm_vtv(nx, f, v);;
    }

    /* K = Up.dot(v * ac / s) = Up.dot(v) * (ac / s) */
#define K h /*Don't need h any more, use it to store K*/
#define D v
    yakfm_set_uv(nx, K, u, v);
    yakfm_set_vrn(nx, K, K, s); /*May be used in place*/

    /*Set W and D*/
    yakfm_bset_vvt(nx1, w, nx, K, f);
    yakfm_bsub_u(nx1, w, nx, u);

    /* Now w is (Kf - Up|***) */
    YAKFM_BSET_V(nx1, 0, nx, w, nx, K);
    /* Now w is (Kf - Up|K) */

    /* D = concatenate([Dp, np.array([r])]) */
    memcpy((void *)D, (void *)d, nx * sizeof(yakfFloat));
    D[nx] = r;

    /* Up, Dp = MWGSU(W, D)*/
    yakfm_mwgsu(nx, nx1, u, d, w, D);

    /* self.x += K * nu */
    yakfm_add_vxn(nx, self->x, K, nu);
#undef D /*Don't nee D any more*/
#undef K /*Don't nee K any more*/
}

void yakf_do_not_use_this_update(yakfAdaptiveSt * self, yakfFloat * z)
{
    yakf_base_update((yakfBaseSt *)self, z, _do_not_use_this_update);
}
