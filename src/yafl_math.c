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

#include "yafl_math.h"

#define _DO_VXN(name, op)                                                \
yaflStatusEn name(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n) \
{                                                                        \
    yaflInt k;                                                           \
                                                                         \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                  \
    YAFL_CHECK(v,   YAFL_ST_INV_ARG_3);                                  \
                                                                         \
    for (k = 0; k < sz; k++)                                             \
    {                                                                    \
        res[k] op v[k] * n;                                              \
    }                                                                    \
    return YAFL_ST_OK;                                                   \
}

_DO_VXN(yafl_math_set_vxn,  =)
_DO_VXN(yafl_math_add_vxn, +=)
_DO_VXN(yafl_math_sub_vxn, -=)

#define _DO_VRN(name, op)                                                \
yaflStatusEn name(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n) \
{                                                                        \
    yaflInt k;                                                           \
                                                                         \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                  \
    YAFL_CHECK(v,   YAFL_ST_INV_ARG_3);                                  \
    YAFL_CHECK(n,   YAFL_ST_INV_ARG_4);                                  \
                                                                         \
    for (k = 0; k < sz; k++)                                             \
    {                                                                    \
        res[k] op v[k] / n;                                              \
    }                                                                    \
    return YAFL_ST_OK;                                                   \
}

_DO_VRN(yafl_math_set_vrn,  =)
_DO_VRN(yafl_math_add_vrn, +=)
_DO_VRN(yafl_math_sub_vrn, -=)

#define _DO_VXV(name, op)                                                \
yaflStatusEn name(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)\
{                                                                        \
    yaflInt k;                                                           \
                                                                         \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                  \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_3);                                  \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_4);                                  \
                                                                         \
    for (k = 0; k < sz; k++)                                             \
    {                                                                    \
        res[k] op a[k] * b[k];                                           \
    }                                                                    \
    return YAFL_ST_OK;                                                   \
}

_DO_VXV(yafl_math_set_vxv,  =)
_DO_VXV(yafl_math_add_vxv, +=)
_DO_VXV(yafl_math_sub_vxv, -=)

#define _DO_VRV(name, op)                                                \
yaflStatusEn name(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b)\
{                                                                        \
    yaflInt k;                                                           \
                                                                         \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                  \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_3);                                  \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_4);                                  \
                                                                         \
    for (k = 0; k < sz; k++)                                             \
    {                                                                    \
        res[k] op a[k] / b[k];                                           \
    }                                                                    \
    return YAFL_ST_OK;                                                   \
}

_DO_VRV(yafl_math_set_vrv,  =)
_DO_VRV(yafl_math_add_vrv, +=)
_DO_VRV(yafl_math_sub_vrv, -=)

yaflFloat yafl_math_vtv(yaflInt sz, yaflFloat *a, yaflFloat *b)
{
    yaflInt k;
    yaflFloat res;

    YAFL_CHECK(a, 0.0);
    YAFL_CHECK(b, 0.0);

    res = a[0] * b[0];
    for (k = 1; k < sz; k++)
    {
        res += a[k] * b[k];
    }
    return res;
}

#define _DO_VVT(name, op)                                                             \
yaflStatusEn name(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b) \
{                                                                                     \
    yaflInt j;                                                                        \
                                                                                      \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_3);                                               \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_4);                                               \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_5);                                               \
                                                                                      \
    for (j = 0; j < nr; j++)                                                          \
    {                                                                                 \
        yaflInt k;                                                                    \
        yaflInt ncj;                                                                  \
        yaflFloat aj;                                                                 \
                                                                                      \
        ncj = nc * j;                                                                 \
        aj  = a[j];                                                                   \
                                                                                      \
        for (k = 0; k < nc; k++)                                                      \
        {                                                                             \
            res[ncj + k] op aj * b[k];                                                \
        }                                                                             \
    }                                                                                 \
    return YAFL_ST_OK;                                                                \
}

_DO_VVT(yafl_math_set_vvt,  =)
_DO_VVT(yafl_math_add_vvt, +=)
_DO_VVT(yafl_math_sub_vvt, -=)

#define _DO_VVTXN(name, op)                                                                        \
yaflStatusEn name(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n) \
{                                                                                                  \
    yaflInt j;                                                                                     \
                                                                                                   \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_3);                                                            \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_4);                                                            \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_5);                                                            \
                                                                                                   \
    for (j = 0; j < nr; j++)                                                                       \
    {                                                                                              \
        yaflInt k;                                                                                 \
        yaflInt ncj;                                                                               \
        yaflFloat aj;                                                                              \
                                                                                                   \
        ncj = nc * j;                                                                              \
        aj  = a[j] * n;                                                                            \
                                                                                                   \
        for (k = 0; k < nc; k++)                                                                   \
        {                                                                                          \
            res[ncj + k] op aj * b[k];                                                             \
        }                                                                                          \
    }                                                                                              \
    return YAFL_ST_OK;                                                                             \
}

_DO_VVTXN(yafl_math_set_vvtxn,  =)
_DO_VVTXN(yafl_math_add_vvtxn, +=)
_DO_VVTXN(yafl_math_sub_vvtxn, -=)

#define _DO_MV(name, op1, op2)                                                        \
yaflStatusEn name(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b) \
{                                                                                     \
    yaflInt j;                                                                        \
                                                                                      \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_3);                                               \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_4);                                               \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_5);                                               \
                                                                                      \
    for (j = 0; j < nr; j++)                                                          \
    {                                                                                 \
        yaflInt k;                                                                    \
        yaflInt ncj;                                                                  \
        yaflFloat resj;                                                               \
                                                                                      \
        ncj  =   nc * j;                                                              \
        resj =   res[j];                                                              \
        resj op1 a[ncj] * b[0];                                                       \
                                                                                      \
        for (k = 1; k < nc; k++)                                                      \
        {                                                                             \
            resj op2 a[ncj + k] * b[k];                                               \
        }                                                                             \
        res[j] = resj;                                                                \
    }                                                                                 \
    return YAFL_ST_OK;                                                                \
}

_DO_MV(yafl_math_set_mv,  =, +=)
_DO_MV(yafl_math_add_mv, +=, +=)
_DO_MV(yafl_math_sub_mv, -=, -=)

#define _DO_VTM(name, op1, op2)                                                       \
yaflStatusEn name(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b) \
{                                                                                     \
    yaflInt j;                                                                        \
                                                                                      \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_3);                                               \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_4);                                               \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_5);                                               \
                                                                                      \
    for (j = 0; j < nc; j++)                                                          \
    {                                                                                 \
        res[j] op1 a[0] * b[j];                                                       \
    }                                                                                 \
                                                                                      \
    for (j = 1; j < nr; j++)                                                          \
    {                                                                                 \
        yaflInt k;                                                                    \
        yaflInt ncj;                                                                  \
        yaflFloat aj;                                                                 \
                                                                                      \
        ncj = nc * j;                                                                 \
        aj = a[j];                                                                    \
                                                                                      \
        for (k = 0; k < nc; k++)                                                      \
        {                                                                             \
            res[k] op2 aj * b[ncj + k];                                               \
        }                                                                             \
    }                                                                                 \
    return YAFL_ST_OK;                                                                \
}

_DO_VTM(yafl_math_set_vtm,  =, +=)
_DO_VTM(yafl_math_add_vtm, +=, +=)
_DO_VTM(yafl_math_sub_vtm, -=, -=)

/*This is right as it is OMP friendly style*/
#define _DO_MM(name, op1, op2)                                                                      \
yaflStatusEn name(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b) \
{                                                                                                   \
    yaflInt i;                                                                                      \
                                                                                                    \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_4);                                                             \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_5);                                                             \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_6);                                                             \
                                                                                                    \
    for (i = 0; i < nr; i++)                                                                        \
    {                                                                                               \
        yaflInt j;                                                                                  \
        yaflInt k;                                                                                  \
        yaflInt nci;                                                                                \
        yaflInt ncri;                                                                               \
                                                                                                    \
        nci = nc * i;                                                                               \
        ncri = ncr * i;                                                                             \
                                                                                                    \
        for (k = 0; k < nc; k++)                                                                    \
        {                                                                                           \
            res[nci + k] op1 a[ncri] * b[k];                                                        \
        }                                                                                           \
                                                                                                    \
        for (j = 1; j < ncr; j++)                                                                   \
        {                                                                                           \
            yaflInt ncj;                                                                            \
            yaflFloat aij;                                                                          \
                                                                                                    \
            ncj = nc * j;                                                                           \
            aij = a[ncri + j];                                                                      \
                                                                                                    \
            for (k = 0; k < nc; k++)                                                                \
            {                                                                                       \
                res[nci + k] op2 aij * b[ncj + k];                                                  \
            }                                                                                       \
        }                                                                                           \
    }                                                                                               \
    return YAFL_ST_OK;                                                                              \
}

_DO_MM(yafl_math_set_mm,  =, +=)
_DO_MM(yafl_math_add_mm, +=, +=)
_DO_MM(yafl_math_sub_mm, -=, -=)

#define _DO_VTU(name, op1, op2)                                           \
yaflStatusEn name(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b) \
{                                                                         \
    yaflInt j;                                                            \
    yaflInt szj;                                                          \
                                                                          \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                   \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_3);                                   \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_4);                                   \
                                                                          \
    for (j = 0; j < sz; j++)                                              \
    {                                                                     \
        res[j] op1 a[j];                                                  \
    }                                                                     \
                                                                          \
    for (szj = 0, j = 1; j < sz; szj += j++)                              \
    {                                                                     \
        yaflInt k;                                                        \
        yaflFloat resj;                                                   \
                                                                          \
        resj = res[j];                                                    \
        for (k = 0; k < j; k++)                                           \
        {                                                                 \
            resj op2 a[k] * b[k + szj];                                   \
        }                                                                 \
        res[j] = resj;                                                    \
    }                                                                     \
    return YAFL_ST_OK;                                                    \
}

_DO_VTU(yafl_math_set_vtu,  =, +=)
_DO_VTU(yafl_math_add_vtu, +=, +=)
_DO_VTU(yafl_math_sub_vtu, -=, -=)

#define _DO_UV(name, op1, op2)                                            \
yaflStatusEn name(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b) \
{                                                                         \
    yaflInt j;                                                            \
    yaflInt szj;                                                          \
                                                                          \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                   \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_3);                                   \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_4);                                   \
                                                                          \
    for (j = 0; j < sz; j++)                                              \
    {                                                                     \
        res[j] op1 b[j];                                                  \
    }                                                                     \
                                                                          \
    for (szj = 0, j = 1; j < sz; szj += j++)                              \
    {                                                                     \
        yaflInt k;                                                        \
        yaflFloat bj;                                                     \
                                                                          \
        bj = b[j];                                                        \
        for (k = 0; k < j; k++)                                           \
        {                                                                 \
            res[k] op2 a[k + szj] * bj;                                   \
        }                                                                 \
    }                                                                     \
    return YAFL_ST_OK;                                                    \
}

_DO_UV(yafl_math_set_uv,  =, +=)
_DO_UV(yafl_math_add_uv, +=, +=)
_DO_UV(yafl_math_sub_uv, -=, -=)

/*This is right as it is OMP friendly style*/
#define _DO_MU(name, op1, op2)                                                         \
yaflStatusEn name(yaflInt nr,  yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b) \
{                                                                                      \
    yaflInt i;                                                                         \
                                                                                       \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_3);                                                \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_4);                                                \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_5);                                                \
                                                                                       \
    for (i = 0; i < nr; i++)                                                           \
    {                                                                                  \
        yaflInt j;                                                                     \
        yaflInt nci;                                                                   \
        yaflInt ncj;                                                                   \
                                                                                       \
        nci = nc * i;                                                                  \
                                                                                       \
        for (j = 0; j < nc; j++)                                                       \
        {                                                                              \
            res[nci + j] op1 a[nci + j];                                               \
        }                                                                              \
                                                                                       \
        for (ncj = 0, j = 1; j < nc; ncj += j++)                                       \
        {                                                                              \
            yaflInt k;                                                                 \
            yaflFloat resij;                                                           \
                                                                                       \
            resij = res[nci + j];                                                      \
            for (k = 0; k < j; k++)                                                    \
            {                                                                          \
                resij op2 a[nci + k] * b[k + ncj];                                     \
            }                                                                          \
            res[nci + j] = resij;                                                      \
        }                                                                              \
    }                                                                                  \
    return YAFL_ST_OK;                                                                 \
}

_DO_MU(yafl_math_set_mu,  =, +=)
_DO_MU(yafl_math_add_mu, +=, +=)
_DO_MU(yafl_math_sub_mu, -=, -=)

yaflStatusEn yafl_math_set_u(yaflInt sz, yaflFloat *res, yaflFloat *u)
{
    yaflInt i;

    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_3);

    for (i = 0; i < sz; i++)
    {
        yaflInt szi;
        yaflInt j;

        szi = sz * i;
        for (j = 0; j < i; j++)
        {
            res[szi + j] = 0.0;
        }

        res[szi + j] = 1.0;

        for (j++; j < sz; j++)
        {
            res[szi + j] = u[i + ((j - 1) * j) / 2];
        }
    }
    return YAFL_ST_OK;
}

#define _DO_U(name, op)                                     \
yaflStatusEn name(yaflInt sz, yaflFloat *res, yaflFloat *u) \
{                                                           \
    yaflInt i;                                              \
                                                            \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                     \
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_3);                     \
                                                            \
    for (i = 0; i < sz; i++)                                \
    {                                                       \
        yaflInt szi;                                        \
        yaflInt j;                                          \
                                                            \
        res[(sz + 1) * i] op 1.0;                           \
                                                            \
        szi = sz * i;                                       \
        for (j = i + 1; j < sz; j++)                        \
        {                                                   \
            res[szi + j] op u[i + ((j - 1) * j) / 2];       \
        }                                                   \
    }                                                       \
    return YAFL_ST_OK;                                      \
}

_DO_U(yafl_math_add_u, +=)
_DO_U(yafl_math_sub_u, -=)

yaflStatusEn yafl_math_bset_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)
{
    yaflInt i;

    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_4);

    for (i = 0; i < sz; i++)
    {
        yaflInt nci;
        yaflInt j;

        nci = nc * i;
        for (j = 0; j < i; j++)
        {
            res[nci + j] = 0.0;
        }

        res[nci + j] = 1.0;

        for (j++; j < sz; j++)
        {
            res[nci + j] = u[i + ((j - 1) * j) / 2];
        }
    }
    return YAFL_ST_OK;
}

#define _DO_BU(name, op)                                                \
yaflStatusEn name(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u) \
{                                                                       \
    yaflInt i;                                                          \
                                                                        \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                 \
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_4);                                 \
                                                                        \
    for (i = 0; i < sz; i++)                                            \
    {                                                                   \
        yaflInt nci;                                                    \
        yaflInt j;                                                      \
                                                                        \
        res[(nc + 1) * i] op 1.0;                                       \
                                                                        \
        nci = nc * i;                                                   \
        for (j = i + 1; j < sz; j++)                                    \
        {                                                               \
            res[nci + j] op u[i + ((j - 1) * j) / 2];                   \
        }                                                               \
    }                                                                   \
    return YAFL_ST_OK;                                                  \
}

_DO_BU(yafl_math_badd_u, +=)
_DO_BU(yafl_math_bsub_u, -=)

yaflStatusEn yafl_math_bset_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u)
{
    yaflInt i;
    yaflInt szi;

    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_4);

    for (szi = 0, i = 0; i < sz; szi += i++)
    {
        yaflInt nci;
        yaflInt j;

        nci = nc * i;
        for (j = 0; j < i; j++)
        {
            res[nci + j] = u[j + szi];
        }

        res[nci + j] = 1.0;

        for (j++; j < sz; j++)
        {
            res[nci + j] = 0.0;
        }
    }
    return YAFL_ST_OK;
}

#define _DO_BUT(name, op)                                               \
yaflStatusEn name(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u) \
{                                                                       \
    yaflInt i;                                                          \
    yaflInt szi;                                                        \
                                                                        \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                 \
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_4);                                 \
                                                                        \
    for (szi = 0, i = 1; i < sz; szi += i++)                            \
    {                                                                   \
        yaflInt nci;                                                    \
        yaflInt j;                                                      \
                                                                        \
        nci = nc * i;                                                   \
        for (j = 0; j < i; j++)                                         \
        {                                                               \
            res[nci + j] op u[j + szi];                                 \
        }                                                               \
                                                                        \
        res[nci + j] op 1.0;                                            \
    }                                                                   \
    return YAFL_ST_OK;                                                  \
}

_DO_BUT(yafl_math_badd_ut, +=)
_DO_BUT(yafl_math_bsub_ut, -=)

#define _DO_BV(name, op)                                                \
yaflStatusEn name(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v) \
{                                                                       \
    yaflInt i;                                                          \
                                                                        \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                 \
    YAFL_CHECK(v,   YAFL_ST_INV_ARG_4);                                 \
                                                                        \
    for (i = 0; i < sz; i++)                                            \
    {                                                                   \
        res[nc * i] op v[i];                                            \
    }                                                                   \
    return YAFL_ST_OK;                                                  \
}

_DO_BV(yafl_math_bset_v,  =)
_DO_BV(yafl_math_badd_v, +=)
_DO_BV(yafl_math_bsub_v, -=)

#define _DO_BVVT(name, op)                                                            \
yaflStatusEn name(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b) \
{                                                                                     \
    yaflInt j;                                                                        \
                                                                                      \
    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);                                               \
    YAFL_CHECK(a,   YAFL_ST_INV_ARG_4);                                               \
    YAFL_CHECK(b,   YAFL_ST_INV_ARG_5);                                               \
                                                                                      \
    for (j = 0; j < sz; j++)                                                          \
    {                                                                                 \
        yaflInt k;                                                                    \
        yaflInt ncj;                                                                  \
        yaflFloat aj;                                                                 \
                                                                                      \
        ncj = nc * j;                                                                 \
        aj = a[j];                                                                    \
                                                                                      \
        for (k = 0; k < sz; k++)                                                      \
        {                                                                             \
            res[ncj + k] op aj * b[k];                                                \
        }                                                                             \
    }                                                                                 \
    return YAFL_ST_OK;                                                                \
}

_DO_BVVT(yafl_math_bset_vvt,  =)
_DO_BVVT(yafl_math_badd_vvt, +=)
_DO_BVVT(yafl_math_bsub_vvt, -=)

#define _DO_BMU(name, op1, op2)                                                                    \
yaflStatusEn name(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b) \
{                                                                                                  \
    yaflInt i;                                                                                     \
                                                                                                   \
    YAFL_CHECK(res,      YAFL_ST_INV_ARG_2);                                                       \
    YAFL_CHECK(a,        YAFL_ST_INV_ARG_5);                                                       \
    YAFL_CHECK(b,        YAFL_ST_INV_ARG_6);                                                       \
    YAFL_CHECK(rnc > nc, YAFL_ST_INV_ARG_1);                                                       \
                                                                                                   \
    for (i = 0; i < nr; i++)                                                                       \
    {                                                                                              \
        yaflInt j;                                                                                 \
        yaflInt nci;                                                                               \
        yaflInt rnci;                                                                              \
        yaflInt ncj;                                                                               \
                                                                                                   \
        nci  = nc * i;                                                                             \
        rnci = rnc * i;                                                                            \
                                                                                                   \
        for (j = 0; j < nc; j++)                                                                   \
        {                                                                                          \
            res[rnci + j] op1 a[nci + j];                                                          \
        }                                                                                          \
                                                                                                   \
        for (ncj = 0, j = 1; j < nc; ncj += j++)                                                   \
        {                                                                                          \
            yaflInt k;                                                                             \
            yaflFloat resij;                                                                       \
                                                                                                   \
            resij = res[rnci + j];                                                                 \
            for (k = 0; k < j; k++)                                                                \
            {                                                                                      \
                resij op2 a[nci + k] * b[k + ncj];                                                 \
            }                                                                                      \
            res[rnci + j] = resij;                                                                 \
        }                                                                                          \
    }                                                                                              \
    return YAFL_ST_OK;                                                                             \
}

_DO_BMU(yafl_math_bset_mu,  =, +=)
_DO_BMU(yafl_math_badd_mu, +=, +=)
_DO_BMU(yafl_math_bsub_mu, -=, -=)

#define _DO_BBU(name, op1, op2)                                                                                 \
yaflStatusEn name(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b) \
{                                                                                                               \
    yaflInt i;                                                                                                  \
                                                                                                                \
    YAFL_CHECK(res,      YAFL_ST_INV_ARG_2);                                                                    \
    YAFL_CHECK(a,        YAFL_ST_INV_ARG_6);                                                                    \
    YAFL_CHECK(b,        YAFL_ST_INV_ARG_7);                                                                    \
    YAFL_CHECK(rnc > nc, YAFL_ST_INV_ARG_1);                                                                    \
    YAFL_CHECK(anc > nc, YAFL_ST_INV_ARG_5);                                                                    \
                                                                                                                \
    for (i = 0; i < nr; i++)                                                                                    \
    {                                                                                                           \
        yaflInt j;                                                                                              \
        yaflInt anci;                                                                                           \
        yaflInt rnci;                                                                                           \
        yaflInt ncj;                                                                                            \
                                                                                                                \
        anci = anc * i;                                                                                         \
        rnci = rnc * i;                                                                                         \
                                                                                                                \
        for (j = 0; j < nc; j++)                                                                                \
        {                                                                                                       \
            res[rnci + j] op1 a[anci + j];                                                                      \
        }                                                                                                       \
                                                                                                                \
        for (ncj = 0, j = 1; j < nc; ncj += j++)                                                                \
        {                                                                                                       \
            yaflInt k;                                                                                          \
            yaflFloat resij;                                                                                    \
                                                                                                                \
            resij = res[rnci + j];                                                                              \
            for (k = 0; k < j; k++)                                                                             \
            {                                                                                                   \
                resij op2 a[anci + k] * b[k + ncj];                                                             \
            }                                                                                                   \
            res[rnci + j] = resij;                                                                              \
        }                                                                                                       \
    }                                                                                                           \
    return YAFL_ST_OK;                                                                                          \
}

_DO_BBU(yafl_math_bset_bu,  =, +=)
_DO_BBU(yafl_math_badd_bu, +=, +=)
_DO_BBU(yafl_math_bsub_bu, -=, -=)

yaflStatusEn yafl_math_ruv(yaflInt sz, yaflFloat *res, yaflFloat *u)
{
    yaflInt j;
    yaflInt szj;

    YAFL_CHECK(res, YAFL_ST_INV_ARG_2);
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_3);

    for (j = sz - 1, szj = ((j - 1) * j) / 2; j > 0; szj -= --j)
    /*for (j = sz - 1; j > 0; j--)*/
    {
        yaflInt i;
        yaflFloat resj;
        /*yaflInt szj;

        szj = ((j - 1)*j)/2;*/
        resj = res[j];

        for (i = j - 1; i >= 0; i--)
        {
            res[i] -= u[i + szj] * resj;
        }
    }
    return YAFL_ST_OK;
}

yaflStatusEn yafl_math_rum(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *u)
{
    yaflInt j;
    yaflInt nrj;

    YAFL_CHECK(res, YAFL_ST_INV_ARG_3);
    YAFL_CHECK(u,   YAFL_ST_INV_ARG_4);

    for (j = nr - 1, nrj = ((j - 1) * j) / 2; j > 0; nrj -= --j)
    {
        yaflInt ncj;
        yaflInt i;

        ncj = nc * j;

        for (i = j - 1; i >= 0; i--)
        {
            yaflInt k;
            yaflInt nci;

            yaflFloat uij;

            nci = nc * i;
            uij = u[i + nrj];

            for (k = nc - 1; k >= 0; k--)
            {
                res[nci + k] -= uij * res[ncj + k];
            }
        }
    }
    return YAFL_ST_OK;
}

yaflStatusEn yafl_math_mwgsu(yaflInt nr, yaflInt nc, yaflFloat *res_u, yaflFloat *res_d, yaflFloat *w, yaflFloat *d)
{
    yaflInt j;
    yaflInt nrj;

    YAFL_CHECK(res_u, YAFL_ST_INV_ARG_3);
    YAFL_CHECK(res_d, YAFL_ST_INV_ARG_4);
    YAFL_CHECK(w,     YAFL_ST_INV_ARG_5);
    YAFL_CHECK(d,     YAFL_ST_INV_ARG_6);

    for (j = nr - 1, nrj = ((j - 1) * j) / 2; j >= 0; nrj -= --j)
    {
        yaflInt   ncj;
        yaflInt   k;
        yaflFloat res_dj;
        yaflFloat wjk;

        ncj = nc * j;

        /*res_d[j] = w[j].dot(d * w[j])*/
        wjk     = w[ncj + nc - 1];
        wjk    *= wjk;
        res_dj  = wjk * d[nc - 1];
        for (k = nc - 2; k >= 0; k--)
        {
            wjk     = w[ncj + k];
            wjk    *= wjk;
            res_dj += wjk * d[k];
        }

        /*Bad Eigenvalue workaround*/
        if (res_dj < YAFL_UDU_EPS)
        {
            res_d[j] = YAFL_UDU_EPS;

            for (k = j - 1; k >= 0; k--)
            {
                res_u[k + nrj] = 0;
            }
            continue;
        }

        /*Good Eigenvalue*/
        res_d[j] = res_dj;

        for (k = j - 1; k >= 0; k--)
        {
            yaflInt   nck;
            yaflInt   i;
            yaflFloat res_ukj;

            nck = nc * k;

            /* res_u[k,j] = w[k].dot(d * w[j])/res_d[j] */
            res_ukj = w[nck] * d[0] * w[ncj] / res_dj;
            for (i = 1; i < nc; i++)
            {
                res_ukj += w[nck + i] * d[i] * w[ncj + i] / res_dj;
            }
            res_u[k + nrj] = res_ukj;

            /* w[k] -= res_u[k,j] * w[j] */
            w[nck + nc - 1] -= res_ukj * w[ncj + nc - 1];
            for (i = nc - 2; i >= 0; i--)
            {
                w[nck + i] -= res_ukj * w[ncj + i];
            }
        }
    }
    return YAFL_ST_OK;
}

yaflStatusEn yafl_math_udu_up(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v)
{
    yaflInt j;
    yaflInt szj;

    YAFL_CHECK(res_u,      YAFL_ST_INV_ARG_2);
    YAFL_CHECK(res_d,      YAFL_ST_INV_ARG_3);
    YAFL_CHECK(v,          YAFL_ST_INV_ARG_5);
    YAFL_CHECK(alpha >= 0, YAFL_ST_INV_ARG_4);

    for (j = sz - 1, szj = ((j - 1) * j) / 2; j >= 0; szj -= --j)
    {
        yaflInt k;
        yaflFloat pj;
        yaflFloat dj;
        yaflFloat res_dj;
        yaflFloat betaj;

        dj = res_d[j];

        pj = v[j];
        res_dj = dj + alpha * pj * pj;
        if (res_dj < YAFL_UDU_EPS)
        {
            res_dj = YAFL_UDU_EPS;
        }

        betaj = alpha * pj / res_dj;

        /*Update res_d[j]*/
        res_d[j] = res_dj;

        /*Update res_u[:j,j]*/
        for (k = j - 1; k >= 0; k--)
        {
            /*Buffer vars*/
            yaflFloat ukj;
            yaflFloat vk;

            ukj = res_u[k + szj];

            vk  = v[k] - pj * ukj;
            v[k] = vk;

            res_u[k + szj] = ukj + betaj * vk;
        }
        alpha *= dj / res_dj;
    }
    return YAFL_ST_OK;
}

yaflStatusEn yafl_math_udu_down(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v)
{
    yaflInt j;
    yaflInt szj;
    yaflFloat pj;
    yaflFloat dj;

    YAFL_CHECK(res_u,      YAFL_ST_INV_ARG_2);
    YAFL_CHECK(res_d,      YAFL_ST_INV_ARG_3);
    YAFL_CHECK(v,          YAFL_ST_INV_ARG_5);
    YAFL_CHECK(alpha >= 0, YAFL_ST_INV_ARG_4);

    /*Solve U*p = v*/
    yafl_math_ruv(sz, v, res_u);

    /*Compute alpha0*/
    pj = v[0];
    dj = pj * pj / res_d[0];
    for (j = 1; j < sz; j++)
    {
        if (res_d[j] < YAFL_UDU_EPS)
        {
            res_d[j] = YAFL_UDU_EPS;
        }

        pj  = v[j];
        dj += pj * pj / res_d[j];
    }

    dj = 1 - alpha * dj;

    if (dj < YAFL_UDU_EPS)
    {
        dj = YAFL_UDU_EPS;
    }

    alpha /= dj;

    /*Conmpute update*/
    for (szj = 0, j = 0; j < sz; szj += j++)
    {
        yaflInt k;
        yaflFloat res_dj;
        yaflFloat betaj;

        dj = res_d[j];
        pj = v[j]; /*We don't need to do v[j] = pj; as it was done in #L774*/

        res_dj = dj * dj / (dj + alpha * pj * pj);
        /*No need to check dj or res_d[j] as it was set big enough in #L783*/
        betaj = alpha * pj / dj; /*Sign "-" will be used in #L825*/

        /*Update res_d[j]*/
        res_d[j] = res_dj;

        for (k = 0; k < j; k++)
        {
            /*Buffer vars*/
            yaflFloat ukj;
            yaflFloat vk;

            ukj = res_u[k + szj];
            vk  = v[k];

            res_u[k + szj] = ukj - betaj * vk;
            v[k]           = vk  + pj    * ukj;
        }
        alpha *= res_dj / dj;
    }
    return YAFL_ST_OK;
}
