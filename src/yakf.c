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
    YAKFM_BSET_U(nx2, 0, 0, w, nx, self->Uq);
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
    YAKF_ASSERT(self->Nx > 1);
    YAKF_ASSERT(self->Up);
    YAKF_ASSERT(self->Dp);
    YAKF_ASSERT(self->H);
    YAKF_ASSERT(self->Dr);
    YAKF_ASSERT(self->D);

    nx = self->Nx;

    u = self->Up;
    d = self->Dp;

    h = self->H + nx * i;

    /* f = h.dot(Up) */
    f = self->D;
    yakfm_set_vtu(nx, f, h, self->Up);

    /* v = Dp.dot(f) */
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
    yakfm_add_nv(nx, self->x, v, self->y[i] / r);
#undef v /*Don't nee v any more*/
}

void yakf_bierman_update(yakfBaseSt * self, yakfFloat * z)
{
    yakf_base_update(self, z, _bierman_scalar_update);
}
