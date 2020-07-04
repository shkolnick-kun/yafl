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
        /*zrf must be aware of self internal structure*/
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

    /* self.x += K * y[i] */
    yakfm_add_vxn(nx, self->x, K, self->y[i]);
#undef D /*Don't nee D any more*/
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

/*=============================================================================
                            Robust Bierman filter
=============================================================================*/
static void _robust_bierman_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt j;
    yakfInt k;
    yakfInt nx;
    yakfInt nxk;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat r05;
    yakfFloat gdot;
    yakfFloat y;
    yakfRobFuncP g;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->D);

    r05 = self->Dr[i]; /* alpha = r**0.5 is stored in Dr*/
    y   = self->y[i];
    g   = ((yakfRobustSt *)self)->g;
    if (g)
    {
        gdot = y / r05; /*Use gdot as temp variable*/
        y = r05 * g(self, gdot);

        g = ((yakfRobustSt *)self)->gdot;
        YAKF_ASSERT(g);

        gdot = g(self, gdot);
    }
    else
    {
        gdot = 1.0;
    }

    r05 *= r05;
#define A2 r05 /*Now it is r = alpha**2 */

    u = self->Up;
    d = self->Dp;

    h = self->H + nx * i;

    /* f = h.dot(Up) */
    f = self->D;
    yakfm_set_vtu(nx, f, h, u);

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#define v h /*Don't need h any more, use it to store v*/
    YAKFM_SET_DV(nx, v, d, f);

    for (k = 0, nxk = 0; k < nx; nxk += k++)
    {
        yakfFloat a;
        yakfFloat fk;
        yakfFloat vk;

        fk = gdot * f[k];
        vk = v[k];
        a = A2 + fk * vk;
        d[k] *= A2 / a;
#define p fk /*No need for separate p variable*/
        p = - fk / A2;
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
        A2 = a;
    }
    /*
    Now we must do:
    y = alpha * g(y[i] / alpha)
    self.x += K * y

    Since:
    A2 == a

    then we have:
    K == v / a == v / A2

    and so:
    K * nu == (v / A2) * y == v / A2 * y == v * (y / A2)

    Finally we get:
    self.x += v * (y / A2)
    */
    yakfm_add_vxn(nx, self->x, v, y / A2);
#undef v  /*Don't nee v any more*/
#undef A2 /*Don't nee A2 any more*/
}

void yakf_robust_bierman_update(yakfRobustSt * self, yakfFloat * z)
{
    yakf_base_update((yakfBaseSt *)self, z, _robust_bierman_scalar_update);
}

/*=============================================================================
                            Robust Joseph filter
=============================================================================*/
static void _robust_joseph_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt nx;
    yakfInt nx1;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat * v;
    yakfFloat * w;
    yakfFloat r05;
    yakfFloat gdot;
    yakfFloat s;
    yakfFloat y;
    yakfRobFuncP g;

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

    r05 = self->Dr[i]; /* alpha = r**0.5 is stored in Dr*/
    y   = self->y[i];
    g   = ((yakfRobustSt *)self)->g;
    if (g)
    {
        s = y / r05;
        y = r05 * g(self, s);

        g = ((yakfRobustSt *)self)->gdot;
        YAKF_ASSERT(g);

        gdot = g(self, s);
    }
    else
    {
        gdot = 1.0;
    }

    r05 *= r05;
#define A2 r05 /*Now it is r = alpha**2 */

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

    /* s = alpha**2 + gdot * f.dot(v)*/
    s = A2 + gdot * yakfm_vtv(nx, f, v);

    /*K = Up.dot(v/s) = Up.dot(v)/s*/
#define K h /*Don't need h any more, use it to store K*/
#define D v
    yakfm_set_uv(nx, K, u, v);
    yakfm_set_vrn(nx, K, K, s); /*May be used in place*/

    /*Set W and D*/
    yakfm_set_vxn(nx, f, f, gdot); /*May be used in place*/
    yakfm_bset_vvt(nx1, w, nx, K, f); /*How about yakfm_bset_vvtxn ?*/
    yakfm_bsub_u(nx1, w, nx, u);

    /* Now w is (gdot*Kf - Up|***) */
    YAKFM_BSET_V(nx1, 0, nx, w, nx, K);
    /* Now w is (gdot*Kf - Up|K) */

    /* D = concatenate([Dp, np.array([gdot * alpha**2])]) */
    memcpy((void *)D, (void *)d, nx * sizeof(yakfFloat));
    D[nx] = gdot * A2;

    /* Up, Dp = MWGSU(W, D)*/
    yakfm_mwgsu(nx, nx1, u, d, w, D);

    /* self.x += K * alpha * g(y[i] / alpha) */
    yakfm_add_vxn(nx, self->x, K, y);
#undef D  /*Don't nee D any more*/
#undef K  /*Don't nee K any more*/
#undef A2 /*Don't nee A2 any more*/
}

void yakf_robust_joseph_update(yakfRobustSt * self, yakfFloat * z)
{
    yakf_base_update((yakfBaseSt *)self, z, _robust_joseph_scalar_update);
}

/*=============================================================================
                        Adaptive robust Bierman filter
=============================================================================*/
static void _ada_rob_bierman_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt j;
    yakfInt k;
    yakfInt nx;
    yakfInt nxk;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat r05;
    yakfFloat gdot;
    yakfFloat nu;
    yakfFloat c;
    yakfFloat s;
    yakfFloat ac;
    yakfRobFuncP g;

    YAKF_ASSERT(self);

    nx = self->Nx;
    YAKF_ASSERT(nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->y);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->D);

    r05 = self->Dr[i]; /* alpha = r**0.5 is stored in Dr*/
    nu   = self->y[i];
    g   = ((yakfRobustSt *)self)->g;
    if (g)
    {
        gdot = nu / r05; /*Use gdot as temp variable*/
        nu = r05 * g(self, gdot);

        g = ((yakfRobustSt *)self)->gdot;
        YAKF_ASSERT(g);

        gdot = g(self, gdot);
    }
    else
    {
        gdot = 1.0;
    }

    r05 *= r05;
#define A2 r05 /*Now it is r = alpha**2 */

    u = self->Up;
    d = self->Dp;

    h = self->H + nx * i;

    /* f = h.dot(Up) */
    f = self->D;
    yakfm_set_vtu(nx, f, h, u);

    /* v = f.dot(Dp).T = Dp.dot(f.T).T */
#define v h /*Don't need h any more, use it to store v*/
    YAKFM_SET_DV(nx, v, d, f);

    /* s = alpha**2 + gdot * f.dot(v)*/
    c = gdot * yakfm_vtv(nx, f, v);
    s = A2 + c;

    /* Divergence test */
    ac = gdot * (nu * (nu / (((yakfAdaptiveRobustSt *)self)->chi2))) - s;
    if (ac > 0.0)
    {
        /*Adaptive correction factor*/
        ac = ac / c + 1.0;

        /*Corrected s*/
        s  = A2 + ac * c;
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

        fk = gdot * f[k];
        /*Correct v in place*/
        vk = ac * v[k];
        v[k] = vk;
        a = A2 + fk * vk;
        /*Correct d in place*/
        d[k] *= ac * A2 / a;
#define p fk /*No need for separate p variable*/
        p = - fk / A2;
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
        A2 = a;
    }
    /*
    Now we must do:
    nu = alpha * g(y[i] / alpha)
    self.x += K * nu

    Since:
    A2 == a

    then we have:
    K == v / a == v / A2

    and so:
    K * nu == (v / A2) * nu == v / A2 * nu == v * (nu / A2)

    Finally we get:
    self.x += v * (nu / A2)
    */
    yakfm_add_vxn(nx, self->x, v, nu / A2);
#undef v  /*Don't nee v any more*/
#undef A2 /*Don't nee A2 any more*/
}

void yakf_adaptive_robust_bierman_update(yakfAdaptiveRobustSt * self, yakfFloat * z)
{
    yakf_base_update((yakfBaseSt *)self, z, _ada_rob_bierman_scalar_update);
}

/*=============================================================================
                        Adaptive robust Joseph filter
=============================================================================*/
static void _ada_rob_joseph_scalar_update(yakfBaseSt * self, yakfInt i)
{
    yakfInt nx;
    yakfInt nx1;
    yakfFloat * d;
    yakfFloat * u;
    yakfFloat * h;
    yakfFloat * f;
    yakfFloat * v;
    yakfFloat * w;
    yakfFloat r05;
    yakfFloat gdot;
    yakfFloat ac;
    yakfFloat c;
    yakfFloat s;
    yakfFloat nu;
    yakfRobFuncP g;

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

    r05 = self->Dr[i]; /* alpha = r**0.5 is stored in Dr*/
    nu  = self->y[i];
    g   = ((yakfRobustSt *)self)->g;
    if (g)
    {
        s  = nu / r05;
        nu = r05 * g(self, s); /*nu = alpha * g(y[i] / alpha)*/

        g = ((yakfRobustSt *)self)->gdot;
        YAKF_ASSERT(g);

        gdot = g(self, s);
    }
    else
    {
        gdot = 1.0;
    }

    r05 *= r05;
#define A2 r05 /*Now it is r = alpha**2 */

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

    /* s = alpha**2 + gdot * f.dot(v)*/
    c = gdot * yakfm_vtv(nx, f, v);
    s = A2 + c;

    /* Divergence test */
    ac = gdot * (nu * (nu / (((yakfAdaptiveRobustSt *)self)->chi2))) - s;
    if (ac > 0.0)
    {
        /*Adaptive correction factor*/
        ac = ac / c + 1.0;

        /*Corrected s*/
        s  = A2 + ac * c;
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
    yakfm_set_vxn(nx, f, f, gdot); /*May be used in place*/
    yakfm_bset_vvt(nx1, w, nx, K, f); /*How about yakfm_bset_vvtxn ?*/
    yakfm_bsub_u(nx1, w, nx, u);

    /* Now w is (gdot*Kf - Up|***) */
    YAKFM_BSET_V(nx1, 0, nx, w, nx, K);
    /* Now w is (gdot*Kf - Up|K) */

    /* D = concatenate([ac * Dp, np.array([gdot * alpha**2])]) */
    yakfm_set_vxn(nx, D, d, ac);
    D[nx] = gdot * A2;

    /* Up, Dp = MWGSU(W, D)*/
    yakfm_mwgsu(nx, nx1, u, d, w, D);

    /* self.x += K * nu */
    yakfm_add_vxn(nx, self->x, K, nu);
#undef D  /*Don't nee D any more*/
#undef K  /*Don't nee K any more*/
#undef A2 /*Don't nee A2 any more*/
}

void yakf_adaptive_robust_joseph_update(yakfAdaptiveRobustSt * self, yakfFloat * z)
{
    yakf_base_update((yakfBaseSt *)self, z, _ada_rob_joseph_scalar_update);
}

/*=============================================================================
                    Basic UD-factorized UKF functions
=============================================================================*/
static inline void _compute_res(yakfUnscentedSt * self, yakfInt sz,          \
                                yakfUnscentedResFuncP rf, yakfFloat * sigma, \
                                yakfFloat * pivot, yakfFloat * res)
{
    if (rf)
    {
        /*rf must be aware of sp and the current transform*/
        rf(self, res, sigma, pivot); /* sp = self.rf(sigmas[i], res_v) */
    }
    else
    {
        yakfInt j;
        for (j = 0; j < sz; j++)
        {
            res[j] = sigma[j] - pivot[j];
        }
    }
}

static void _unscented_transform(yakfUnscentedSt * self, \
                                 yakfInt    res_sz,      \
                                 yakfFloat * res_v,      \
                                 yakfFloat * res_u,      \
                                 yakfFloat * res_d,      \
                                 yakfFloat * sp,         \
                                 yakfFloat * sigmas,     \
                                 yakfFloat * noise_u,    \
                                 yakfFloat * noise_d,    \
                                 yakfUnscentedFuncP mf,  \
                                 yakfUnscentedResFuncP rf)
{
    yakfInt np;
    yakfInt i;
    yakfSigmaSt * points;
    yakfFloat * wc;

    YAKF_ASSERT(self);

    YAKF_ASSERT(res_sz > 0);
    YAKF_ASSERT(res_v);
    YAKF_ASSERT(res_u);
    YAKF_ASSERT(res_d);
    YAKF_ASSERT(sp);

    if (noise_u)
    {
        YAKF_ASSERT(noise_d);
    }

    YAKF_ASSERT(self->points);
    points = self->points;

    YAKF_ASSERT(points->np > 1);
    np = points->np;

    YAKF_ASSERT(points->wm);
    YAKF_ASSERT(points->wc);
    wc = points->wc;

    if (mf)
    {
        mf(self, res_v, sigmas); /*mf must be aware of the current transform details...*/
    }
    else
    {
        yakfm_set_vtm(np, res_sz, res_v, points->wm, sigmas);
    }

    if (noise_u)
    {
        /*res_u, res_d = noise_u.copy(), noise_d.copy()*/
        memcpy((void *)res_u, (void *)noise_u, (res_sz * (res_sz - 1)) / 2 * sizeof(yakfFloat));
        memcpy((void *)res_d, (void *)noise_d, res_sz * sizeof(yakfFloat));

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

    for (i = 0; i < np; i++)
    {
        _compute_res(self, res_sz, rf, sigmas + res_sz * i, res_v, sp);
        /*Update res_u and res_d*/
        /*wc should be sorted in descending order*/
        if (wc[i] >= 0.0)
        {
            yakfm_udu_up(res_sz, res_u, res_d, wc[i], sp);
        }
        else
        {
            yakfm_udu_down(res_sz, res_u, res_d, -wc[i], sp);
        }
    }
}

void yakf_unscented_predict(yakfUnscentedSt * self)
{
    yakfInt np;
    yakfInt nx;
    yakfInt i;
    yakfSigmaSt * points;
    yakfUnscentedFuncP fx; /*State transition function*/
    yakfFloat * sigmas_x;

    /*Check some params and generate sigma points*/
    yakf_unscented_gen_sigmas(self); /*Self is checked here*/

    YAKF_ASSERT(self->Nx);
    nx = self->Nx;

    YAKF_ASSERT(self->sigmas_x);
    sigmas_x = self->sigmas_x;

    YAKF_ASSERT(self->points);
    points = self->points;

    YAKF_ASSERT(points->np > 1);
    np = points->np;

    /*Compute process sigmas*/
    YAKF_ASSERT(self->f);
    fx = self->f;
    if (fx)
    {
        for (i = 0; i < np; i++)
        {
            yakfFloat * sigmai;
            sigmai = sigmas_x + nx * i;
            fx(self, sigmai, sigmai);
        }
    }

    /*Predict x, Up, Dp*/
    _unscented_transform(self, nx, self->x, self->Up, self->Dp, self->Sx, \
                         sigmas_x, self->Uq, self->Dq, self->xmf, self->xrf);
}

void yakf_unscented_update(yakfUnscentedSt * self, yakfFloat * z)
{
    yakfInt np;
    yakfInt nx;
    yakfInt nz;
    yakfInt i;
    yakfSigmaSt * points;      /*Sigma points generator*/
    yakfUnscentedFuncP hx;     /*State transition function*/
    yakfUnscentedResFuncP xrf; /*State residual function*/
    yakfUnscentedResFuncP zrf; /*Measurement residual function*/
    yakfFloat * sigmas_x;      /*State sigma points*/
    yakfFloat * sigmas_z;      /*Measurement sigma points*/
    /*Scratchpad memory*/
    yakfFloat * spx;           /*For state residuals*/
    yakfFloat * spz;           /*For measurement residuals*/

    yakfFloat * x;             /*State*/
    yakfFloat * zp;            /*Predicted measurement*/
    yakfFloat * wc;            /*Covariance computation weights*/
    yakfFloat * pzx;           /*Pzx cross covariance matrix*/

    /*Pzz covariance*/
    yakfFloat * us;            /*Unit upper triangular factor*/
    yakfFloat * ds;            /*Diagonal factor*/

    /*P covariance*/
    yakfFloat * up;            /*Unit upper triangular factor*/
    yakfFloat * dp;            /*Diagonal factor*/

    YAKF_ASSERT(self);
    YAKF_ASSERT(self->Nx);
    nx = self->Nx;

    YAKF_ASSERT(self->sigmas_x);
    sigmas_x = self->sigmas_x;

    YAKF_ASSERT(self->Nz);
    nz = self->Nz;

    YAKF_ASSERT(self->sigmas_z);
    sigmas_z = self->sigmas_z;

    YAKF_ASSERT(self->points);
    points = self->points;

    YAKF_ASSERT(points->np > 1);
    np = points->np;

    YAKF_ASSERT(points->wc);
    wc = points->wc;

    YAKF_ASSERT(self->h);
    hx = self->h;

    YAKF_ASSERT(self->Pzx);
    pzx = self->Pzx;

    YAKF_ASSERT(self->Sx);
    spx = self->Sx;

    YAKF_ASSERT(self->Sz);
    spz = self->Sz;

    YAKF_ASSERT(self->zp);
    zp = self->zp;

    YAKF_ASSERT(self->x);
    x = self->x;

    YAKF_ASSERT(self->Us);
    us = self->Us;

    YAKF_ASSERT(self->Ds);
    ds = self->Ds;

    YAKF_ASSERT(self->Up);
    up = self->Up;

    YAKF_ASSERT(self->Dp);
    dp = self->Dp;

    /* Compute measurement sigmas */
    for (i = 0; i < np; i++)
    {
        hx(self, sigmas_z + nz * i, sigmas_x + nx * i);
    }

    /* Compute zp, Us, Ds */
    zrf = self->zrf;
    _unscented_transform(self, nz, zp, us, ds, spz, \
                         sigmas_z, self->Ur, self->Dr, self->zmf, zrf);

    /* Compute Pzx */
    xrf = self->xrf;
    _compute_res(self, nz, zrf, sigmas_z, zp, spz);
    _compute_res(self, nx, xrf, sigmas_x,  x, spx);
    yakfm_set_vvtxn(nz, nx, pzx, spz, spx, wc[0]);

    for (i = 1; i < np; i++)
    {
        _compute_res(self, nz, zrf, sigmas_z + nz * i, zp, spz);
        _compute_res(self, nx, xrf, sigmas_x + nx * i,  x, spx);
        yakfm_add_vvtxn(nz, nx, pzx, spz, spx, wc[i]);
    }

    /*Compute innovation*/
#define Y spx
    _compute_res(self, nz, zrf, z, zp, Y);

    /* Decorrelate measurements*/
    yakfm_ruv(nz,       Y, us);
    yakfm_rum(nz, nx, pzx, us);

    /*Now we can do scalar updates*/
    for (i = 0; i < nz; i++)
    {
        yakfFloat * pzxi;
        pzxi = pzx + nx * i;
        /*
        self.x += K * Y[i]

        K * Y[i] = Pzx[i].T / ds[i] * Y[i] = Pzx[i].T * (Y[i] / ds[i])

        self.x += Pzx[i].T * (Y[i] / ds[i])
        */
        yakfm_add_vxn(nx, x, pzxi, Y[i] / ds[i]);

        /*
        P -= K.dot(S.dot(K.T))
        K.dot(S.dot(K.T)) = (Pzx[i].T / ds[i] * ds[i]).outer(Pzx[i] / ds[i]))
        K.dot(S.dot(K.T)) = (Pzx[i].T).outer(Pzx[i]) / ds[i]
        P -= (Pzx[i].T).outer(Pzx[i]) / ds[i]
        Up, Dp = udu(P)
        */
        yakfm_udu_down(nx, up, dp, 1.0 / ds[i], pzxi);
    }
#undef Y
}