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

yaflStatusEn jfx(yaflKalmanBaseSt * self, yaflFloat * w, yaflFloat * x)
{
    yaflInt i;
    yaflInt nx;
    yaflInt nx2;

    (void)x;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(4 == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(w,             YAFL_ST_INV_ARG_2);

    nx  = self->Nx;
    nx2 = nx * 2;

    for (i = 0; i < nx; i++)
    {
        yaflInt j;
        yaflInt nci;

        nci = nx2 * i;
        for (j = 0; j < nx; j++)
        {
            w[nci + j] = (i != j) ? 0.0 : 1.0;
        }
    }

    w[nx2*0 + 2] = 0.1;
    w[nx2*1 + 3] = 0.1;
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

yaflStatusEn jhx(yaflKalmanBaseSt * self, yaflFloat * h, yaflFloat * x)
{
    yaflInt i;
    yaflInt nx;
    yaflInt nz;

    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(4 == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(2 == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(h,             YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_3);

    nx = self->Nx;
    nz = self->Nz;

    for (i = 0; i < nz; i++)
    {
        yaflInt j;
        yaflInt nci;

        nci = nx * i;
        for (j = 0; j < nx; j++)
        {
            h[nci + j] = 0.0;
        }
    }

    h[nx*0 + 0] = 1.0;
    h[nx*1 + 1] = 1.0;
    return YAFL_ST_OK;
}
/*---------------------------------------------------------------------------*/
typedef struct
{
    YAFL_EKF_BASE_MEMORY_MIXIN(NX, NZ);
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

yaflEKFBaseSt kf = YAFL_EKF_BASE_INITIALIZER(fx, jfx, hx, jhx, 0, NX, NZ, 0.0, 0.0, kf_memory);

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
    for (i=0; i<mat.shape.dim.x; i++)
    {
        YAFL_EKF_JOSEPH_PREDICT(&kf);
        yafl_ekf_joseph_update(&kf, mat.data + NZ*i);
        mat.data[NZ*i + 0] = kf.base.x[0];
        mat.data[NZ*i + 1] = kf.base.x[1];
    }

    file = H5Fcreate(OUT_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hdf5_utils_write_array(file, OUT_DS, &mat);
    status = H5Fclose(file);
    printf("File %s closed with status: %d\n", OUT_FILE, status);

    free(mat.data);
    return 0;
}
