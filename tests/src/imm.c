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
#define NX 3
#define NZ 1
#define DT (0.1)

yaflStatusEn cv(yaflKalmanBaseSt * self, yaflFloat * x, yaflFloat * xz)
{
    (void)xz;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(3 == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_2);

    x[0] += x[1] * DT;
    return YAFL_ST_OK;
}

yaflStatusEn jcv(yaflKalmanBaseSt * self, yaflFloat * w, yaflFloat * x)
{
    yaflInt i;
    yaflInt nx;
    yaflInt nx2;

    (void)x;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(3 == self->Nx, YAFL_ST_INV_ARG_1);
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

    w[nx2*0 + 1] = DT;
    return YAFL_ST_OK;
}

yaflStatusEn ca(yaflKalmanBaseSt * self, yaflFloat * x, yaflFloat * xz)
{
    (void)xz;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(3 == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_2);

    x[0] += (x[1] + 0.5 * x[2] * DT) * DT;
    x[1] += x[2] * DT;
    return YAFL_ST_OK;
}

yaflStatusEn jca(yaflKalmanBaseSt * self, yaflFloat * w, yaflFloat * x)
{
    yaflInt i;
    yaflInt nx;
    yaflInt nx2;

    (void)x;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(3 == self->Nx, YAFL_ST_INV_ARG_1);
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

    w[nx2*0 + 1] = DT;
    w[nx2*0 + 2] = DT * DT * 0.5;
    w[nx2*1 + 2] = DT;
    return YAFL_ST_OK;
}

yaflStatusEn hx(yaflKalmanBaseSt * self, yaflFloat * y, yaflFloat * x)
{
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(1 == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(y,             YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_3);

    y[0] = x[0];
    return YAFL_ST_OK;
}

yaflStatusEn jhx(yaflKalmanBaseSt * self, yaflFloat * h, yaflFloat * x)
{
    yaflInt i;
    yaflInt nx;
    yaflInt nz;

    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(3 == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(1 == self->Nz, YAFL_ST_INV_ARG_1);
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
    return YAFL_ST_OK;
}
/*---------------------------------------------------------------------------*/
typedef struct
{
    YAFL_EKF_BASE_MEMORY_MIXIN(NX, NZ);
} ekfMemorySt;

#define DP (0.1)
#define DX (1.0e-6)
#define DZ (1)

ekfMemorySt cv_memory =
{
    .x = {
        [0] = 0.0,
        [1] = 0.0,
        [2] = 0.0
    },

    .Up = {
        0,
        0,0,
    },
    .Dp = {DP, DP, DP},

    .Uq = {
        0,
        0,0,
    },
    .Dq = {DX, DX, DX},

    //.Ur = {0},
    .Dr = {DZ}
};

yaflEKFBaseSt cv_kf = YAFL_EKF_BASE_INITIALIZER(cv, jcv, hx, jhx, 0, NX, NZ, 0.0, 0.0, cv_memory);

/*---------------------------------------------------------------------------*/
typedef struct
{
    YAFL_UKF_MEMORY_MIXIN(NX, NZ);
    YAFL_UKF_JULIER_MEMORY_MIXIN(NX, NZ);
} ukfMemorySt;

ekfMemorySt ca_memory =
//ukfMemorySt ca_memory =
{
    .x = {
        [0] = 0.0,
        [1] = 0.0,
        [2] = 0.0
    },

    .Up = {
        0,
        0,0,
    },
    .Dp = {DP, DP, DP},

    .Uq = {
        0,
        0,0,
    },
    .Dq = {DX, DX, DX},

    //.Ur = {0},
    .Dr = {DZ}
};

yaflEKFBaseSt ca_kf = YAFL_EKF_BASE_INITIALIZER(ca, jca, hx, jhx, 0, NX, NZ, 0.0, 0.0, ca_memory);
//yaflUKFJulierSt sp    = YAFL_UKF_JULIER_INITIALIZER(NX, 0, 0.0, kf_memory);
//yaflUKFSt       ca_kf = YAFL_UKF_INITIALIZER(&sp.base, &yafl_ukf_julier_spm, ca, 0, 0, hx, 0, 0, NX, NZ, 0.0, ca_memory);

/*---------------------------------------------------------------------------*/
typedef struct {
    YAFL_IMM_MEMORY_MIXIN(2,3);
}immMemorySt;

yaflFilterBankItemSt imm_bank[2] = {
    [0] = YAFL_IMM_EKF_ITEM_INITIALIZER(&cv_kf, cv_memory, yafl_ekf_base_predict, yafl_ekf_bierman_update),
    [1] = YAFL_IMM_EKF_ITEM_INITIALIZER(&ca_kf, ca_memory, yafl_ekf_base_predict, yafl_ekf_bierman_update)
    //[1] = YAFL_IMM_UKF_ITEM_INITIALIZER(&ca_kf, ca_memory, yafl_ukf_base_predict, yafl_ukf_update)
};

immMemorySt imm_memory = {
    .mu = {0.5, 0.5},
    .M  = {
        0.95, 0.05,
        0.05, 0.95
    },
    .Up = {
        0,
        0,0,
    },
    .Dp = {DP, DP, DP},
    .x = {
        [0] = 0.0,
        [1] = 0.0,
        [2] = 0.0
    },
    .cbar  = {0.0, 0.0},
    .omega = {
        0.0, 0.0,
        0.0, 0.0
    },
    .y = {
        [0] = 0.0,
        [1] = 0.0,
        [2] = 0.0
    },
    .D = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .W = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        },
};

yaflIMMCBSt imm = YAFL_IMM_INITIALIZER(imm_bank, 2, imm_memory);

/*-----------------------------------------------------------------------------
                                  Test data
-----------------------------------------------------------------------------*/
#define IN_FILE  "../../data/input.h5"
#define IN_DS    "noisy"

#define OUT_FILE "../../data/output.h5"
#define OUT_DS   "imm_out"

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
    assert(YAFL_ST_OK == yafl_imm_post_init(&imm));

    for (i=0; i<mat.shape.dim.x; i++)
    {
        assert(YAFL_ST_ERR_THR > yafl_imm_predict(&imm));
        assert(YAFL_ST_ERR_THR > yafl_imm_update(&imm, mat.data + NZ*i));
        mat.data[NZ*i + 0] = imm.x[0];
    }

    file = H5Fcreate(OUT_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hdf5_utils_write_array(file, OUT_DS, &mat);
    status = H5Fclose(file);
    printf("File %s closed with status: %d\n", OUT_FILE, status);

    free(mat.data);
    return 0;
}
