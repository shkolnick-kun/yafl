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

#include "yakf_math.h"

#define _DO_VXN(name, op)                                         \
void name(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n) \
{                                                                \
    yakfInt k;                                                   \
                                                                 \
    YAKF_ASSERT(res);                                            \
    YAKF_ASSERT(v);                                              \
                                                                 \
    for (k = 0; k < sz; k++)                                     \
    {                                                            \
        res[k] op v[k] * n;                                      \
    }                                                            \
}

_DO_VXN(yakfm_set_vxn,  =)
_DO_VXN(yakfm_add_vxn, +=)
_DO_VXN(yakfm_sub_vxn, -=)

#define _DO_VRN(name, op)                                        \
void name(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n) \
{                                                                \
    yakfInt k;                                                   \
                                                                 \
    YAKF_ASSERT(res);                                            \
    YAKF_ASSERT(v);                                              \
    YAKF_ASSERT(n);                                              \
                                                                 \
    for (k = 0; k < sz; k++)                                     \
    {                                                            \
        res[k] op v[k] / n;                                      \
    }                                                            \
}

_DO_VRN(yakfm_set_vrn,  =)
_DO_VRN(yakfm_add_vrn, +=)
_DO_VRN(yakfm_sub_vrn, -=)

#define _DO_VXV(name, op)                                        \
void name(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b)\
{                                                                \
    yakfInt k;                                                   \
                                                                 \
    YAKF_ASSERT(res);                                            \
    YAKF_ASSERT(a);                                              \
    YAKF_ASSERT(b);                                              \
                                                                 \
    for (k = 0; k < sz; k++)                                     \
    {                                                            \
        res[k] op a[k] * b[k];                                   \
    }                                                            \
}

_DO_VXV(yakfm_set_vxv,  =)
_DO_VXV(yakfm_add_vxv, +=)
_DO_VXV(yakfm_sub_vxv, -=)

#define _DO_VRV(name, op)                                        \
void name(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b)\
{                                                                \
    yakfInt k;                                                   \
                                                                 \
    YAKF_ASSERT(res);                                            \
    YAKF_ASSERT(a);                                              \
    YAKF_ASSERT(b);                                              \
                                                                 \
    for (k = 0; k < sz; k++)                                     \
    {                                                            \
        res[k] op a[k] / b[k];                                   \
    }                                                            \
}

_DO_VRV(yakfm_set_vrv,  =)
_DO_VRV(yakfm_add_vrv, +=)
_DO_VRV(yakfm_sub_vrv, -=)

yakfFloat yakfm_vtv(yakfInt sz, yakfFloat *a, yakfFloat *b)
{
    yakfInt k;
    yakfFloat res;

    YAKF_ASSERT(a);
    YAKF_ASSERT(b);

    res = a[0] * b[0];
    for (k = 1; k < sz; k++)
    {
        res += a[k] * b[k];
    }
    return res;
}

#define _DO_VVT(name, op)                                                     \
void name(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b) \
{                                                                             \
    yakfInt j;                                                                \
                                                                              \
    YAKF_ASSERT(res);                                                         \
    YAKF_ASSERT(a);                                                           \
    YAKF_ASSERT(b);                                                           \
                                                                              \
    for (j = 0; j < nr; j++)                                                  \
    {                                                                         \
        yakfInt k;                                                            \
        yakfInt ncj;                                                          \
        yakfFloat aj;                                                         \
                                                                              \
        ncj = nc * j;                                                         \
        aj  = a[j];                                                           \
                                                                              \
        for (k = 0; k < nc; k++)                                              \
        {                                                                     \
            res[ncj + k] op aj * b[k];                                        \
        }                                                                     \
    }                                                                         \
}

_DO_VVT(yakfm_set_vvt,  =)
_DO_VVT(yakfm_add_vvt, +=)
_DO_VVT(yakfm_sub_vvt, -=)

#define _DO_VVTXN(name, op)                                                                \
void name(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b, yakfFloat n) \
{                                                                                          \
    yakfInt j;                                                                             \
                                                                                           \
    YAKF_ASSERT(res);                                                                      \
    YAKF_ASSERT(a);                                                                        \
    YAKF_ASSERT(b);                                                                        \
                                                                                           \
    for (j = 0; j < nr; j++)                                                               \
    {                                                                                      \
        yakfInt k;                                                                         \
        yakfInt ncj;                                                                       \
        yakfFloat aj;                                                                      \
                                                                                           \
        ncj = nc * j;                                                                      \
        aj  = a[j] * n;                                                                    \
                                                                                           \
        for (k = 0; k < nc; k++)                                                           \
        {                                                                                  \
            res[ncj + k] op aj * b[k];                                                     \
        }                                                                                  \
    }                                                                                      \
}

_DO_VVTXN(yakfm_set_vvtxn,  =)
_DO_VVTXN(yakfm_add_vvtxn, +=)
_DO_VVTXN(yakfm_sub_vvtxn, -=)

#define _DO_MV(name, op1, op2)                                               \
void name(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b)\
{                                                                            \
    yakfInt j;                                                               \
                                                                             \
    YAKF_ASSERT(res);                                                        \
    YAKF_ASSERT(a);                                                          \
    YAKF_ASSERT(b);                                                          \
                                                                             \
    for (j = 0; j < nr; j++)                                                 \
    {                                                                        \
        yakfInt k;                                                           \
        yakfInt ncj;                                                         \
        yakfFloat resj;                                                      \
                                                                             \
        ncj  =   nc * j;                                                     \
        resj =   res[j];                                                     \
        resj op1 a[ncj] * b[0];                                              \
                                                                             \
        for (k = 1; k < nc; k++)                                             \
        {                                                                    \
            resj op2 a[ncj + k] * b[k];                                      \
        }                                                                    \
        res[j] = resj;                                                       \
    }                                                                        \
}

_DO_MV(yakfm_set_mv,  =, +=)
_DO_MV(yakfm_add_mv, +=, +=)
_DO_MV(yakfm_sub_mv, -=, -=)

#define _DO_VTM(name, op1, op2)                                              \
void name(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b)\
{                                                                            \
    yakfInt j;                                                               \
                                                                             \
    YAKF_ASSERT(res);                                                        \
    YAKF_ASSERT(a);                                                          \
    YAKF_ASSERT(b);                                                          \
                                                                             \
    for (j = 0; j < nc; j++)                                                 \
    {                                                                        \
        res[j] op1 a[0] * b[j];                                              \
    }                                                                        \
                                                                             \
    for (j = 1; j < nr; j++)                                                 \
    {                                                                        \
        yakfInt k;                                                           \
        yakfInt ncj;                                                         \
        yakfFloat aj;                                                        \
                                                                             \
        ncj = nc * j;                                                        \
        aj = a[j];                                                           \
                                                                             \
        for (k = 0; k < nc; k++)                                             \
        {                                                                    \
            res[k] op2 aj * b[ncj + k];                                      \
        }                                                                    \
    }                                                                        \
}

_DO_VTM(yakfm_set_vtm,  =, +=)
_DO_VTM(yakfm_add_vtm, +=, +=)
_DO_VTM(yakfm_sub_vtm, -=, -=)

/*This is right as it is OMP friendly style*/
#define _DO_MM(name, op1, op2)                                                             \
void name(yakfInt nr,  yakfInt ncr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b)\
{                                                                                          \
    yakfInt i;                                                                             \
                                                                                           \
    YAKF_ASSERT(res);                                                                      \
    YAKF_ASSERT(a);                                                                        \
    YAKF_ASSERT(b);                                                                        \
                                                                                           \
    for (i = 0; i < nr; i++)                                                               \
    {                                                                                      \
        yakfInt j;                                                                         \
        yakfInt k;                                                                         \
        yakfInt nci;                                                                       \
        yakfInt ncri;                                                                      \
                                                                                           \
        nci = nc * i;                                                                      \
        ncri = ncr * i;                                                                    \
                                                                                           \
        for (k = 0; k < nc; k++)                                                           \
        {                                                                                  \
            res[nci + k] op1 a[ncri] * b[k];                                               \
        }                                                                                  \
                                                                                           \
        for (j = 1; j < ncr; j++)                                                          \
        {                                                                                  \
            yakfInt ncj;                                                                   \
            yakfFloat aij;                                                                 \
                                                                                           \
            ncj = nc * j;                                                                  \
            aij = a[ncri + j];                                                             \
                                                                                           \
            for (k = 0; k < nc; k++)                                                       \
            {                                                                              \
                res[nci + k] op2 aij * b[ncj + k];                                         \
            }                                                                              \
        }                                                                                  \
    }                                                                                      \
}

_DO_MM(yakfm_set_mm,  =, +=)
_DO_MM(yakfm_add_mm, +=, +=)
_DO_MM(yakfm_sub_mm, -=, -=)

#define _DO_VTU(name, op1, op2)                                  \
void name(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b)\
{                                                                \
    yakfInt j;                                                   \
    yakfInt szj;                                                 \
                                                                 \
    YAKF_ASSERT(res);                                            \
    YAKF_ASSERT(a);                                              \
    YAKF_ASSERT(b);                                              \
                                                                 \
    for (j = 0; j < sz; j++)                                     \
    {                                                            \
        res[j] op1 a[j];                                         \
    }                                                            \
                                                                 \
    for (szj = 0, j = 1; j < sz; szj += j++)                     \
    {                                                            \
        yakfInt k;                                               \
        yakfFloat resj;                                          \
                                                                 \
        resj = res[j];                                           \
        for (k = 0; k < j; k++)                                  \
        {                                                        \
            resj op2 a[k] * b[k + szj];                          \
        }                                                        \
        res[j] = resj;                                           \
    }                                                            \
}

_DO_VTU(yakfm_set_vtu,  =, +=)
_DO_VTU(yakfm_add_vtu, +=, +=)
_DO_VTU(yakfm_sub_vtu, -=, -=)

#define _DO_UV(name, op1, op2)                                   \
void name(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b)\
{                                                                \
    yakfInt j;                                                   \
    yakfInt szj;                                                 \
                                                                 \
    YAKF_ASSERT(res);                                            \
    YAKF_ASSERT(a);                                              \
    YAKF_ASSERT(b);                                              \
                                                                 \
    for (j = 0; j < sz; j++)                                     \
    {                                                            \
        res[j] op1 b[j];                                         \
    }                                                            \
                                                                 \
    for (szj = 0, j = 1; j < sz; szj += j++)                     \
    {                                                            \
        yakfInt k;                                               \
        yakfFloat bj;                                            \
                                                                 \
        bj = b[j];                                               \
        for (k = 0; k < j; k++)                                  \
        {                                                        \
            res[k] op2 a[k + szj] * bj;                          \
        }                                                        \
    }                                                            \
}

_DO_UV(yakfm_set_uv,  =, +=)
_DO_UV(yakfm_add_uv, +=, +=)
_DO_UV(yakfm_sub_uv, -=, -=)

/*This is right as it is OMP friendly style*/
#define _DO_MU(name, op1, op2)                                                 \
void name(yakfInt nr,  yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b) \
{                                                                              \
    yakfInt i;                                                                 \
                                                                               \
    YAKF_ASSERT(res);                                                          \
    YAKF_ASSERT(a);                                                            \
    YAKF_ASSERT(b);                                                            \
                                                                               \
    for (i = 0; i < nr; i++)                                                   \
    {                                                                          \
        yakfInt j;                                                             \
        yakfInt nci;                                                           \
        yakfInt ncj;                                                           \
                                                                               \
        nci = nc * i;                                                          \
                                                                               \
        for (j = 0; j < nc; j++)                                               \
        {                                                                      \
            res[nci + j] op1 a[nci + j];                                       \
        }                                                                      \
                                                                               \
        for (ncj = 0, j = 1; j < nc; ncj += j++)                               \
        {                                                                      \
            yakfInt k;                                                         \
            yakfFloat resij;                                                   \
                                                                               \
            resij = res[nci + j];                                              \
            for (k = 0; k < j; k++)                                            \
            {                                                                  \
                resij op2 a[nci + k] * b[k + ncj];                             \
            }                                                                  \
            res[nci + j] = resij;                                              \
        }                                                                      \
    }                                                                          \
}

_DO_MU(yakfm_set_mu,  =, +=)
_DO_MU(yakfm_add_mu, +=, +=)
_DO_MU(yakfm_sub_mu, -=, -=)

void yakfm_set_u(yakfInt sz, yakfFloat *res, yakfFloat *u)
{
    yakfInt i;

    YAKF_ASSERT(res);
    YAKF_ASSERT(u);

    for (i = 0; i < sz; i++)
    {
        yakfInt szi;
        yakfInt j;

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
}

#define _DO_U(name, op)                               \
void name(yakfInt sz, yakfFloat *res, yakfFloat *u)   \
{                                                     \
    yakfInt i;                                        \
                                                      \
    YAKF_ASSERT(res);                                 \
    YAKF_ASSERT(u);                                   \
                                                      \
    for (i = 0; i < sz; i++)                          \
    {                                                 \
        yakfInt szi;                                  \
        yakfInt j;                                    \
                                                      \
        res[(sz + 1) * i] op 1.0;                     \
                                                      \
        szi = sz * i;                                 \
        for (j = i + 1; j < sz; j++)                  \
        {                                             \
            res[szi + j] op u[i + ((j - 1) * j) / 2]; \
        }                                             \
    }                                                 \
}

_DO_U(yakfm_add_u, +=)
_DO_U(yakfm_sub_u, -=)

void yakfm_bset_u(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *u)
{
    yakfInt i;

    YAKF_ASSERT(res);
    YAKF_ASSERT(u);

    for (i = 0; i < sz; i++)
    {
        yakfInt nci;
        yakfInt j;

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
}

#define _DO_BU(name, op)                                        \
void name(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *u) \
{                                                               \
    yakfInt i;                                                  \
                                                                \
    YAKF_ASSERT(res);                                           \
    YAKF_ASSERT(u);                                             \
                                                                \
    for (i = 0; i < sz; i++)                                    \
    {                                                           \
        yakfInt nci;                                            \
        yakfInt j;                                              \
                                                                \
        res[(nc + 1) * i] op 1.0;                               \
                                                                \
        nci = nc * i;                                           \
        for (j = i + 1; j < sz; j++)                            \
        {                                                       \
            res[nci + j] op u[i + ((j - 1) * j) / 2];           \
        }                                                       \
    }                                                           \
}

_DO_BU(yakfm_badd_u, +=)
_DO_BU(yakfm_bsub_u, -=)

#define _DO_BV(name, op)                                        \
void name(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *v) \
{                                                               \
    yakfInt i;                                                  \
                                                                \
    YAKF_ASSERT(res);                                           \
    YAKF_ASSERT(v);                                             \
                                                                \
    for (i = 0; i < sz; i++)                                    \
    {                                                           \
        res[nc * i] op v[i];                                    \
    }                                                           \
}

_DO_BV(yakfm_bset_v,  =)
_DO_BV(yakfm_badd_v, +=)
_DO_BV(yakfm_bsub_v, -=)

#define _DO_BVVT(name, op)                                                    \
void name(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *a, yakfFloat *b) \
{                                                                             \
    yakfInt j;                                                                \
                                                                              \
    YAKF_ASSERT(res);                                                         \
    YAKF_ASSERT(a);                                                           \
    YAKF_ASSERT(b);                                                           \
                                                                              \
    for (j = 0; j < sz; j++)                                                  \
    {                                                                         \
        yakfInt k;                                                            \
        yakfInt ncj;                                                          \
        yakfFloat aj;                                                         \
                                                                              \
        ncj = nc * j;                                                         \
        aj = a[j];                                                            \
                                                                              \
        for (k = 0; k < sz; k++)                                              \
        {                                                                     \
            res[ncj + k] op aj * b[k];                                        \
        }                                                                     \
    }                                                                         \
}

_DO_BVVT(yakfm_bset_vvt,  =)
_DO_BVVT(yakfm_badd_vvt, +=)
_DO_BVVT(yakfm_bsub_vvt, -=)

#define _DO_BMU(name, op1, op2)                                                            \
void name(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfFloat *a, yakfFloat *b) \
{                                                                                          \
    yakfInt i;                                                                             \
                                                                                           \
    YAKF_ASSERT(res);                                                                      \
    YAKF_ASSERT(a);                                                                        \
    YAKF_ASSERT(b);                                                                        \
    YAKF_ASSERT(rnc > nc);                                                                 \
                                                                                           \
    for (i = 0; i < nr; i++)                                                               \
    {                                                                                      \
        yakfInt j;                                                                         \
        yakfInt nci;                                                                       \
        yakfInt rnci;                                                                      \
        yakfInt ncj;                                                                       \
                                                                                           \
        nci  = nc * i;                                                                     \
        rnci = rnc * i;                                                                    \
                                                                                           \
        for (j = 0; j < nc; j++)                                                           \
        {                                                                                  \
            res[rnci + j] op1 a[nci + j];                                                  \
        }                                                                                  \
                                                                                           \
        for (ncj = 0, j = 1; j < nc; ncj += j++)                                           \
        {                                                                                  \
            yakfInt k;                                                                     \
            yakfFloat resij;                                                               \
                                                                                           \
            resij = res[rnci + j];                                                         \
            for (k = 0; k < j; k++)                                                        \
            {                                                                              \
                resij op2 a[nci + k] * b[k + ncj];                                         \
            }                                                                              \
            res[rnci + j] = resij;                                                         \
        }                                                                                  \
    }                                                                                      \
}

_DO_BMU(yakfm_bset_mu,  =, +=)
_DO_BMU(yakfm_badd_mu, +=, +=)
_DO_BMU(yakfm_bsub_mu, -=, -=)

#define _DO_BBU(name, op1, op2)                                                                         \
void name(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfInt anc, yakfFloat *a, yakfFloat *b) \
{                                                                                                       \
    yakfInt i;                                                                                          \
                                                                                                        \
    YAKF_ASSERT(res);                                                                                   \
    YAKF_ASSERT(a);                                                                                     \
    YAKF_ASSERT(b);                                                                                     \
    YAKF_ASSERT(anc > nc);                                                                              \
    YAKF_ASSERT(rnc > nc);                                                                              \
                                                                                                        \
    for (i = 0; i < nr; i++)                                                                            \
    {                                                                                                   \
        yakfInt j;                                                                                      \
        yakfInt anci;                                                                                   \
        yakfInt rnci;                                                                                   \
        yakfInt ncj;                                                                                    \
                                                                                                        \
        anci = anc * i;                                                                                 \
        rnci = rnc * i;                                                                                 \
                                                                                                        \
        for (j = 0; j < nc; j++)                                                                        \
        {                                                                                               \
            res[rnci + j] op1 a[anci + j];                                                              \
        }                                                                                               \
                                                                                                        \
        for (ncj = 0, j = 1; j < nc; ncj += j++)                                                        \
        {                                                                                               \
            yakfInt k;                                                                                  \
            yakfFloat resij;                                                                            \
                                                                                                        \
            resij = res[rnci + j];                                                                      \
            for (k = 0; k < j; k++)                                                                     \
            {                                                                                           \
                resij op2 a[anci + k] * b[k + ncj];                                                     \
            }                                                                                           \
            res[rnci + j] = resij;                                                                      \
        }                                                                                               \
    }                                                                                                   \
}

_DO_BBU(yakfm_bset_bu,  =, +=)
_DO_BBU(yakfm_badd_bu, +=, +=)
_DO_BBU(yakfm_bsub_bu, -=, -=)

void yakfm_ruv(yakfInt sz, yakfFloat *res, yakfFloat *u)
{
    yakfInt j;
    yakfInt szj;

    YAKF_ASSERT(res);
    YAKF_ASSERT(u);

    for (j = sz - 1, szj = ((j - 1) * j) / 2; j > 0; szj -= --j)
    /*for (j = sz - 1; j > 0; j--)*/
    {
        yakfInt i;
        yakfFloat resj;
        /*yakfInt szj;

        szj = ((j - 1)*j)/2;*/
        resj = res[j];

        for (i = j - 1; i >= 0; i--)
        {
            res[i] -= u[i + szj] * resj;
        }
    }
}

void yakfm_rum(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *u)
{
    yakfInt j;
    yakfInt nrj;

    YAKF_ASSERT(res);
    YAKF_ASSERT(u);

    for (j = nr - 1, nrj = ((j - 1) * j) / 2; j > 0; nrj -= --j)
    {
        yakfInt ncj;
        yakfInt i;

        ncj = nc * j;

        for (i = j - 1; i >= 0; i--)
        {
            yakfInt k;
            yakfInt nci;

            yakfFloat uij;

            nci = nc * i;
            uij = u[i + nrj];

            for (k = nc - 1; k >= 0; k--)
            {
                res[nci + k] -= uij * res[ncj + k];
            }
        }
    }
}

void yakfm_mwgsu(yakfInt nr, yakfInt nc, yakfFloat *res_u, yakfFloat *res_d, yakfFloat *w, yakfFloat *d)
{
    yakfInt j;
    yakfInt nrj;

    YAKF_ASSERT(res_u);
    YAKF_ASSERT(res_d);
    YAKF_ASSERT(w);
    YAKF_ASSERT(d);

    for (j = nr - 1, nrj = ((j - 1) * j) / 2; j >= 0; nrj -= --j)
    {
        yakfInt   ncj;
        yakfInt   k;
        yakfFloat res_dj;
        yakfFloat wjk;

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
        if (res_dj < YAKF_UDU_EPS)
        {
            res_d[j] = YAKF_UDU_EPS;

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
            yakfInt   nck;
            yakfInt   i;
            yakfFloat res_ukj;

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
}

void yakfm_udu_up(yakfInt sz, yakfFloat *res_u, yakfFloat *res_d, yakfFloat alpha, yakfFloat *v)
{
    yakfInt j;
    yakfInt szj;

    YAKF_ASSERT(res_u);
    YAKF_ASSERT(res_d);
    YAKF_ASSERT(v);
    YAKF_ASSERT(alpha >= 0);

    for (j = sz - 1, szj = ((j - 1) * j) / 2; j >= 0; szj -= --j)
    {
        yakfInt k;
        yakfFloat pj;
        yakfFloat dj;
        yakfFloat res_dj;
        yakfFloat betaj;

        dj = res_d[j];

        pj = v[j];
        res_dj = dj + alpha * pj * pj;
        if (res_dj < YAKF_UDU_EPS)
        {
            res_dj = YAKF_UDU_EPS;
        }

        betaj = alpha * pj / res_dj;

        /*Update res_d[j]*/
        res_d[j] = res_dj;

        /*Update res_u[:j,j]*/
        for (k = j - 1; k >= 0; k--)
        {
            /*Buffer vars*/
            yakfFloat ukj;
            yakfFloat vk;

            ukj = res_u[k + szj];

            vk  = v[k] - pj * ukj;
            v[k] = vk;

            res_u[k + szj] = ukj + betaj * vk;
        }
        alpha *= dj / res_dj;
    }
}

void yakfm_udu_down(yakfInt sz, yakfFloat *res_u, yakfFloat *res_d, yakfFloat alpha, yakfFloat *v)
{
    yakfInt j;
    yakfInt szj;
    yakfFloat pj;
    yakfFloat dj;

    YAKF_ASSERT(res_u);
    YAKF_ASSERT(res_d);
    YAKF_ASSERT(v);
    YAKF_ASSERT(alpha >= 0);

    /*Solve U*p = v*/
    yakfm_ruv(sz, v, res_u);

    /*Compute alpha0*/
    pj = v[0];
    dj = pj * pj / res_d[0];
    for (j = 1; j < sz; j++)
    {
        if (res_d[j] < YAKF_UDU_EPS)
        {
            res_d[j] = YAKF_UDU_EPS;
        }

        pj  = v[j];
        dj += pj * pj / res_d[j];
    }

    dj = 1 - alpha * dj;

    if (dj < YAKF_UDU_EPS)
    {
        dj = YAKF_UDU_EPS;
    }

    alpha /= dj;

    /*Conmpute update*/
    for (szj = 0, j = 0; j < sz; szj += j++)
    {
        yakfInt k;
        yakfFloat res_dj;
        yakfFloat betaj;

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
            yakfFloat ukj;
            yakfFloat vk;

            ukj = res_u[k + szj];
            vk  = v[k];

            res_u[k + szj] = ukj - betaj * vk;
            v[k]           = vk  + pj    * ukj;
        }
        alpha *= res_dj / dj;
    }
}
