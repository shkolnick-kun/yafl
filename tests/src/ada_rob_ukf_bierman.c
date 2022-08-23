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
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <hdf5/serial/hdf5.h>
#include <hdf5utils.h>
#include <yafl.h>

/*-----------------------------------------------------------------------------
                            Kalman filter things
-----------------------------------------------------------------------------*/
#define NX 4
#define NZ 2

yaflStatusEn fx(yaflKalmanBaseSt * self, yaflFloat * x, yaflFloat * xz)
{
    (void)xz;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(4 == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_2);

    x[0] += 0.1 * x[2];
    x[1] += 0.1 * x[3];
    return YAFL_ST_OK;
}

yaflStatusEn hx(yaflKalmanBaseSt * self, yaflFloat * y, yaflFloat * x)
{
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(2 == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(y,             YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_3);

    y[0] = x[0];
    y[1] = x[1];
    return YAFL_ST_OK;
}

/*---------------------------------------------------------------------------*/
yaflFloat g(yaflKalmanBaseSt * self, yaflFloat x)
{
    (void)self;

    if (3.0 >= fabs(x))
    {
        return x;
    }

    if (6.0 >= fabs(x))
    {
        return x/3.0;
    }

    return (x > 0.0) ? 1.0 : -1.0;
}

yaflFloat gdot(yaflKalmanBaseSt * self, yaflFloat x)
{
    (void)self;

    if (3.0 >= fabs(x))
    {
        return 1.0;
    }

    if (6.0 >= fabs(x))
    {
        return 1.0/3.0;
    }

    return 0.0;
}

/*---------------------------------------------------------------------------*/
typedef struct
{
    YAFL_UKF_BASE_MEMORY_MIXIN(NX, NZ);
    YAFL_UKF_JULIER_MEMORY_MIXIN(NX, NZ);
    yaflFloat dummy[30];
} kfMemorySt;

#define DP (0.1)
#define DX (1.0e-6)
#define DZ (400)
kfMemorySt kf_memory =
{
    .x = {
        [0] = 50.0,
        [2] = 10.0
    },

    .Up = {
        0,
        0,0,
        0,0,0
    },
    .Dp = {DP, DP, DP, DP},

    .Uq = {
        0,
        0,0,
        0,0,0
    },
    .Dq = {DX, DX, DX, DX},

    .Ur = {0},
    .Dr = {DZ, DZ}
};

yaflUKFJulierSt         sp = YAFL_UKF_JULIER_INITIALIZER(NX, 0, 0.0, kf_memory);
yaflUKFAdaptiveRobustSt kf = YAFL_UKF_ADAPTIVE_ROBUST_INITIALIZER(&sp.base, &yafl_ukf_julier_spm, fx, 0, 0, hx, 0, 0, g, gdot, NX, NZ, 0.0, 8.807468393511947, kf_memory);

/*-----------------------------------------------------------------------------
                                  Test data
-----------------------------------------------------------------------------*/
#define IN_FILE  "../data/input.h5"
#define IN_DS    "noisy"

#define OUT_FILE "../data/output.h5"
#define OUT_DS   "kf_out"

/*---------------------------------------------------------------------------*/
int main (void)
{
    hid_t  file;
    herr_t status;
    int    i;

    file = H5Fopen(IN_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5UtilsMatSt mat = hdf5_utils_read_array(file, IN_DS);
    status = H5Fclose(file);
    printf("File %s closed with status: %d\n", IN_FILE, status);

    assert(mat.shape.dim.y == NZ);

    yafl_ukf_post_init(&kf.base.base);

    for (i=0; i<mat.shape.dim.x; i++)
    {
        yafl_ukf_adaptive_robust_bierman_predict(&kf);
        yafl_ukf_adaptive_robust_bierman_update(&kf, mat.data + NZ*i);
        mat.data[NZ*i + 0] = kf.base.base.base.x[0];
        mat.data[NZ*i + 1] = kf.base.base.base.x[1];
    }

    file = H5Fcreate(OUT_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hdf5_utils_write_array(file, OUT_DS, &mat);
    status = H5Fclose(file);
    printf("File %s closed with status: %d\n", OUT_FILE, status);

    free(mat.data);
    return 0;
}
