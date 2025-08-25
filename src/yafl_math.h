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

#ifndef YAFL_MATH_H
#define YAFL_MATH_H

/*Config file must be included at first place*/
#include <yafl_config.h>

#define _YAFL_CHECK(cond, err, file, func, line)                           \
do {                                                                       \
    if (YAFL_UNLIKELY(!(cond)))                                            \
    {                                                                      \
        YAFL_DBG("YAFL:The expression (%s) is false in \n function: %s",   \
                 #cond, func);                                             \
        YAFL_DBG("\n file: %s\n line: %d\n will return: %s\n",             \
                 file, line, #err);                                        \
        return err;                                                        \
    }                                                                      \
} while (0)

#define YAFL_CHECK(cond, err) _YAFL_CHECK(cond, err, __FILE__, __func__, __LINE__)

typedef enum {
    /*Warning flag masks*/
    YAFL_ST_MSK_REGULARIZED  = 0x01,
    YAFL_ST_MSK_GLITCH_SMALL = 0x02,
    YAFL_ST_MSK_GLITCH_LARGE = 0x04,
    YAFL_ST_MSK_ANOMALY      = 0x08,
    /*Everything is OK*/
    YAFL_ST_OK           = 0x00,
    /*Warning statuses (see Warning flag masks)*/
    YAFL_ST_R            = 0x01,
    YAFL_ST_S            = 0x02,
    YAFL_ST_SR           = 0x03,
    YAFL_ST_L            = 0x04,
    YAFL_ST_LR           = 0x05,
    YAFL_ST_SL           = 0x06,
    YAFL_ST_SLR          = 0x07,
    YAFL_ST_A            = 0x08,
    YAFL_ST_AR           = 0x09,
    YAFL_ST_SA           = 0x0a,
    YAFL_ST_SAR          = 0x0b,
    YAFL_ST_LA           = 0x0c,
    YAFL_ST_LAR          = 0x0d,
    YAFL_ST_SLA          = 0x0e,
    YAFL_ST_SLAR         = 0x0f,
    /*Error threshold value (greater values are errors)*/
    YAFL_ST_ERR_THR      = 0x010,
    /*Invalid argument numer*/
    YAFL_ST_INV_ARG_1    = 0x100,
    YAFL_ST_INV_ARG_2    = 0x110,
    YAFL_ST_INV_ARG_3    = 0x120,
    YAFL_ST_INV_ARG_4    = 0x130,
    YAFL_ST_INV_ARG_5    = 0x140,
    YAFL_ST_INV_ARG_6    = 0x150,
    YAFL_ST_INV_ARG_7    = 0x160,
    YAFL_ST_INV_ARG_8    = 0x170,
    YAFL_ST_INV_ARG_9    = 0x180,
    YAFL_ST_INV_ARG_10   = 0x190,
    YAFL_ST_INV_ARG_11   = 0x1a0,
    YAFL_ST_INV_ARG_12   = 0x1b0,
    /*Other errors*/
    YAFL_ST_ERR_OTHER

} yaflStatusEn;

/* We need some way of printing human readable statuses */
char * yafl_fail_dsc(yaflStatusEn status);

#define _YAFL_TRY(status, exp, file, func, line)                              \
do {                                                                          \
    (status) |= (exp);                                                        \
    if (YAFL_UNLIKELY(YAFL_ST_ERR_THR <= (status)))                           \
    {                                                                         \
        YAFL_DBG("YAFL:The expression (%s) gave an error in \n function: %s", \
                 #exp, func);                                                 \
        YAFL_DBG("\n file: %s\n line: %d\n will return: %s\n",                \
                 file, line, yafl_fail_dsc(status));                          \
        return status;                                                        \
    }                                                                         \
} while (0)

#define YAFL_TRY(status, exp) _YAFL_TRY(status, exp, __FILE__, __func__, __LINE__)

/*=======================================================================================
                                    Basic operations
=======================================================================================*/
/*
Notation:
[nN] - number
[vV] - vector
[mM] - matrix
[bB] - matrix block
[dD] - diagonal matrix
[uU] - upper triangular unit matrix
[lL] - lower triangular unit matrix

res - result of operation

<*>[tT] - transposed *, e.g. vt is transposed vector
<*>[rR] - reversed   *, e.d. rn is 1/n
x - multiplication

sz - size of a vector or a square/triangular matrix

WARNING: Only a few functions/macros may be used for in place processing!!!
--------------------------------------------------------------------------------------------------------
         Function/Macro (May be used for in place processing)                              NumPy expr
------------------------------------------------------------------------------------------------------*/
yaflStatusEn yafl_math_set_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res  = n*v */
yaflStatusEn yafl_math_add_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res += n*v */
yaflStatusEn yafl_math_sub_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res -= n*v */

yaflStatusEn yafl_math_set_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res  = v/n */
yaflStatusEn yafl_math_add_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res += v/n */
yaflStatusEn yafl_math_sub_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res -= v/n */

yaflStatusEn yafl_math_set_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a*b */
yaflStatusEn yafl_math_add_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a*b */
yaflStatusEn yafl_math_sub_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a*b */

yaflStatusEn yafl_math_set_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a/b */
yaflStatusEn yafl_math_add_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a/b */
yaflStatusEn yafl_math_sub_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a/b */

#define YAFL_MATH_SET_DV yafl_math_set_vxv                                             /* res  = d.dot(v) */
#define YAFL_MATH_ADD_DV yafl_math_add_vxv                                             /* res += d.dot(v) */
#define YAFL_MATH_SUB_DV yafl_math_sub_vxv                                             /* res -= d.dot(v) */

#define YAFL_MATH_SET_RDV(sz, res, a, b) yafl_math_set_vrv(sz, res, b, a)              /* res  = (1/d).dot(v)*/
#define YAFL_MATH_ADD_RDV(sz, res, a, b) yafl_math_add_vrv(sz, res, b, a)              /* res += (1/d).dot(v)*/
#define YAFL_MATH_SUB_RDV(sz, res, a, b) yafl_math_sub_vrv(sz, res, b, a)              /* res += (1/d).dot(v)*/

#define YAFL_MATH_SET_VTD yafl_math_set_vxv                                            /* res  = v.T.dot(d) */
#define YAFL_MATH_ADD_VTD yafl_math_add_vxv                                            /* res += v.T.dot(d) */
#define YAFL_MATH_SUB_VTD yafl_math_sub_vxv                                            /* res -= v.T.dot(d) */

#define YAFL_MATH_SET_VTRD(sz, res, a, b) yafl_math_set_vrv(sz, res, b, a)             /* res  = v.T.dot(1/d)*/
#define YAFL_MATH_ADD_VTRD(sz, res, a, b) yafl_math_add_vrv(sz, res, b, a)             /* res += v.T.dot(1/d)*/
#define YAFL_MATH_SUB_VTRD(sz, res, a, b) yafl_math_sub_vrv(sz, res, b, a)             /* res += v.T.dot(1/d)*/

/*
m - rectangular matrix
v - vector

nc - number of columns in a matrix
nr - number of rows in a matrix
-----------------------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                                          NumPy expr
---------------------------------------------------------------------------------------------------------------------------------------------*/
yaflStatusEn yafl_math_vtv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);                                /* a.T.dot(b)             */

yaflStatusEn yafl_math_set_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);                /* res  = outer(a, b)     */
yaflStatusEn yafl_math_add_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);                /* res += outer(a, b)     */
yaflStatusEn yafl_math_sub_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);                /* res -= outer(a, b)     */

yaflStatusEn yafl_math_set_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n); /* res  = outer(a, b) * n */
yaflStatusEn yafl_math_add_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n); /* res += outer(a, b) * n */
yaflStatusEn yafl_math_sub_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n); /* res -= outer(a, b) * n */

/*
m - rectangular matrix
v - vector

nc - number of columns in a matrix
nr - number of rows in a matrix
----------------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                                      NumPy expr
--------------------------------------------------------------------------------------------------------------------------------------*/

yaflStatusEn yafl_math_set_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);               /* res  = a.dot(b)   */
yaflStatusEn yafl_math_add_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);               /* res += a.dot(b)   */
yaflStatusEn yafl_math_sub_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);               /* res -= a.dot(b)   */

yaflStatusEn yafl_math_set_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);              /* res  = a.T.dot(b) */
yaflStatusEn yafl_math_add_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);              /* res += a.T.dot(b) */
yaflStatusEn yafl_math_sub_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);              /* res -= a.T.dot(b) */

/*ncr - number of columns of in @a and of rows in @b*/
yaflStatusEn yafl_math_set_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a.dot(b)   */
yaflStatusEn yafl_math_add_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a.dot(b)   */
yaflStatusEn yafl_math_sub_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a.dot(b)   */

/*
m - rectangular matrix
u - upper triangular unit matrix
v - vector
--------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                         NumPy expr
------------------------------------------------------------------------------------------------------------------------*/
yaflStatusEn yafl_math_set_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);            /* res  = a.T.dot(b) */
yaflStatusEn yafl_math_add_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);            /* res += a.T.dot(b) */
yaflStatusEn yafl_math_sub_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);            /* res -= a.T.dot(b) */

yaflStatusEn yafl_math_set_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);             /* res  = a.dot(b)   */
yaflStatusEn yafl_math_add_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);             /* res += a.dot(b)   */
yaflStatusEn yafl_math_sub_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);             /* res -= a.dot(b)   */

yaflStatusEn yafl_math_set_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a.dot(b)   */
yaflStatusEn yafl_math_add_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a.dot(b)   */
yaflStatusEn yafl_math_sub_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a.dot(b)   */

yaflStatusEn yafl_math_set_u(yaflInt sz, yaflFloat *res, yaflFloat *u);                            /* res  = u          */
yaflStatusEn yafl_math_add_u(yaflInt sz, yaflFloat *res, yaflFloat *u);                            /* res += u          */
yaflStatusEn yafl_math_sub_u(yaflInt sz, yaflFloat *res, yaflFloat *u);                            /* res -= u          */

/*Matrix row size and block start address*/
#define YAFL_BLK_PTR(m,nc,r,c) ((yaflFloat *)m + nc * r + c)
#define YAFL_BLK(m,nc,r,c) nc, YAFL_BLK_PTR(m,nc,r,c)

/*For test purposes...*/
static inline yaflFloat * _yafl_blk_ptr(yaflFloat * m, yaflInt nc, yaflInt r, yaflInt c)
{
    return YAFL_BLK_PTR(m,nc,r,c);
}

/*
Block operations:

m - rectangular matrix
u - upper triangular unit matrix
v - vector
------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                   NumPy expr
----------------------------------------------------------------------------------------------------------------------*/
yaflStatusEn yafl_math_bset_m(yaflInt nc, yaflFloat *res, yaflInt sr, yaflInt sc, yaflFloat *a); /* res[:sr, :sc]  = a       */
yaflStatusEn yafl_math_badd_m(yaflInt nc, yaflFloat *res, yaflInt sr, yaflInt sc, yaflFloat *a); /* res[:sr, :sc] += a       */
yaflStatusEn yafl_math_bsub_m(yaflInt nc, yaflFloat *res, yaflInt sr, yaflInt sc, yaflFloat *a); /* res[:sr, :sc] -= a       */

#define YAFL_MATH_BSET_M(nc, r, c, m, sr, sc, a) yafl_math_bset_u(YAFL_BLK(m,nc,r,c), sr, sc, a) /* m[r:r+sz, c:c+sz]  = a   */
#define YAFL_MATH_BADD_M(nc, r, c, m, sr, sc, a) yafl_math_badd_u(YAFL_BLK(m,nc,r,c), sr, sc, a) /* m[r:r+sz, c:c+sz] += a   */
#define YAFL_MATH_BSUB_M(nc, r, c, m, sr, sc, a) yafl_math_bsub_u(YAFL_BLK(m,nc,r,c), sr, sc, a) /* m[r:r+sz, c:c+sz] -= a   */

yaflStatusEn yafl_math_bset_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);      /* res[:sz, :sz]  = u       */
yaflStatusEn yafl_math_badd_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);      /* res[:sz, :sz] += u       */
yaflStatusEn yafl_math_bsub_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);      /* res[:sz, :sz] -= u       */

#define YAFL_MATH_BSET_U(nc, r, c, m, sz, u) yafl_math_bset_u(YAFL_BLK(m,nc,r,c), sz, u)  /* m[r:r+sz, c:c+sz]  = u   */
#define YAFL_MATH_BADD_U(nc, r, c, m, sz, u) yafl_math_badd_u(YAFL_BLK(m,nc,r,c), sz, u)  /* m[r:r+sz, c:c+sz] += u   */
#define YAFL_MATH_BSUB_U(nc, r, c, m, sz, u) yafl_math_bsub_u(YAFL_BLK(m,nc,r,c), sz, u)  /* m[r:r+sz, c:c+sz] -= u   */

yaflStatusEn yafl_math_bset_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);     /* res[:sz, :sz]  = u.T     */
yaflStatusEn yafl_math_badd_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);     /* res[:sz, :sz] += u.T     */
yaflStatusEn yafl_math_bsub_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);     /* res[:sz, :sz] -= u.T     */

#define YAFL_MATH_BSET_UT(nc, r, c, m, sz, u) yafl_math_bset_u(YAFL_BLK(m,nc,r,c), sz, u) /* m[r:r+sz, c:c+sz]  = u.T */
#define YAFL_MATH_BADD_UT(nc, r, c, m, sz, u) yafl_math_badd_u(YAFL_BLK(m,nc,r,c), sz, u) /* m[r:r+sz, c:c+sz] += u.T */
#define YAFL_MATH_BSUB_UT(nc, r, c, m, sz, u) yafl_math_bsub_u(YAFL_BLK(m,nc,r,c), sz, u) /* m[r:r+sz, c:c+sz] -= u.T */

yaflStatusEn yafl_math_bset_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v);      /* res[:sz, 0]  = v         */
yaflStatusEn yafl_math_badd_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v);      /* res[:sz, 0] += v         */
yaflStatusEn yafl_math_bsub_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v);      /* res[:sz, 0] -= v         */

#define YAFL_MATH_BSET_V(nc, r, c, m, sz, v) yafl_math_bset_v(YAFL_BLK(m,nc,r,c), sz, v)  /* m[r:r+sz, c]  = v        */
#define YAFL_MATH_BADD_V(nc, r, c, m, sz, v) yafl_math_badd_v(YAFL_BLK(m,nc,r,c), sz, v)  /* m[r:r+sz, c] += v        */
#define YAFL_MATH_BSUB_V(nc, r, c, m, sz, v) yafl_math_bsub_v(YAFL_BLK(m,nc,r,c), sz, v)  /* m[r:r+sz, c] -= v        */

/*
Block operations:

m - rectangular matrix
v - vector
-------------------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                                  NumPy expr
-----------------------------------------------------------------------------------------------------------------------------------------*/
yaflStatusEn yafl_math_bset_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b); /* res[:sz, :sz]  = outer(a, b)     */
yaflStatusEn yafl_math_badd_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b); /* res[:sz, :sz] += outer(a, b)     */
yaflStatusEn yafl_math_bsub_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b); /* res[:sz, :sz] -= outer(a, b)     */

#define YAFL_MATH_BSET_VVT(nc, r, c, m, sz, a, b) yafl_math_bset_vvt(YAFL_BLK(m,nc,r,c), sz, a, b)   /* m[r:r+sz, c:c+sz]  = outer(a, b) */
#define YAFL_MATH_BADD_VVT(nc, r, c, m, sz, a, b) yafl_math_badd_vvt(YAFL_BLK(m,nc,r,c), sz, a, b)   /* m[r:r+sz, c:c+sz] += outer(a, b) */
#define YAFL_MATH_BSUB_VVT(nc, r, c, m, sz, a, b) yafl_math_bsub_vvt(YAFL_BLK(m,nc,r,c), sz, a, b)   /* m[r:r+sz, c:c+sz] -= outer(a, b) */

/*
Block operations:

m - rectangular matrix
u - upper triangular unit matrix

rnc - number of columns in result matrices
----------------------------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                                           NumPy expr
--------------------------------------------------------------------------------------------------------------------------------------------------*/
yaflStatusEn yafl_math_bset_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b); /* res[:nr, :nc]  = a.dot(b)     */
yaflStatusEn yafl_math_badd_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b); /* res[:nr, :nc] += a.dot(b)     */
yaflStatusEn yafl_math_bsub_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b); /* res[:nr, :nc] -= a.dot(b)     */

#define YAFL_MATH_BSET_MU(rnc, r, c, m, nr, nc, a, b) yafl_math_bset_mu(YAFL_BLK(m,rnc,r,c), nr, nc, a, b)       /* m[r:r+nr, c:c+nc]  = a.dot(b) */
#define YAFL_MATH_BADD_MU(rnc, r, c, m, nr, nc, a, b) yafl_math_badd_mu(YAFL_BLK(m,rnc,r,c), nr, nc, a, b)       /* m[r:r+nr, c:c+nc] += a.dot(b) */
#define YAFL_MATH_BSUB_MU(rnc, r, c, m, nr, nc, a, b) yafl_math_bsub_mu(YAFL_BLK(m,rnc,r,c), nr, nc, a, b)       /* m[r:r+nr, c:c+nc] -= a.dot(b) */

/*
Block operations:

u - upper triangular unit matrix

rnc - number of columns in result matrices
anc - number of columns in a matrices blocks

Warning:
Matrix blocks a and res must not overlap!
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                 Function/Macro                                                                                                   NumPy expr
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
yaflStatusEn yafl_math_bset_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b);               /* res[:nr, :nc]  = a[:nr, :nc].dot(b)               */
yaflStatusEn yafl_math_badd_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b);               /* res[:nr, :nc] += a[:nr, :nc].dot(b)               */
yaflStatusEn yafl_math_bsub_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b);               /* res[:nr, :nc] -= a[:nr, :nc].dot(b)               */

#define YAFL_MATH_BSET_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yafl_math_bset_bu(YAFL_BLK(m,rnc,r,c), nr, nc, YAFL_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc]  = a[ar:ar+nr, ac:ac+nc].dot(b) */
#define YAFL_MATH_BADD_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yafl_math_badd_bu(YAFL_BLK(m,rnc,r,c), nr, nc, YAFL_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc] += a[ar:ar+nr, ac:ac+nc].dot(b) */
#define YAFL_MATH_BSUB_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yafl_math_bsub_bu(YAFL_BLK(m,rnc,r,c), nr, nc, YAFL_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc] -= a[ar:ar+nr, ac:ac+nc].dot(b) */

/*
Back substitution:

m - rectangular matrix
u - upper triangular unit matrix
v - vector
------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                           NumPy expr
----------------------------------------------------------------------------------------------------------------*/
/*res is vector*/
yaflStatusEn yafl_math_ruv(yaflInt sz, yaflFloat *res, yaflFloat *u);             /*res = linalg.inv(u).dot(res)*/
yaflStatusEn yafl_math_rutv(yaflInt sz, yaflFloat *res, yaflFloat *u);            /*res = linalg.inv(u.T).dot(res)*/
/*res is matrix*/
yaflStatusEn yafl_math_rum(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *u); /*res = linalg.inv(u).dot(res)*/

/*
Modified Weighted Gram-Schmidt update (for UDU' decomposition)

Does in place:
res_u, res_d = MWGS(w,d)

Warning:
Matrix w is not valid after call.
*/
yaflStatusEn yafl_math_mwgsu(yaflInt nr, yaflInt nc, yaflFloat *res_u, yaflFloat *res_d, yaflFloat *w, yaflFloat *d);

/*
Rank 1 UDU' update.

Based on:
1. Gill, Murray, and Wright, "Practical Optimization", p42.
2. Bierman, "Factorization Methods for Discrete Sequential Estimation", p44.

Does in place:
p = u.dot(d.dot(u.T))
p += alpha * outer(v,v.T)
u,d = udu(p)

Warning:
Vector v is not valid after call.

TODO: add doc with derivation!
*/
yaflStatusEn yafl_math_udu_up(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v);

/*
Rank 1 UDU' downdate.
Based on:
1. Gill, Murray, and Wright, "Practical Optimization", p43.
2. Bierman, "Factorization Methods for Discrete Sequential Estimation", p44.

Does in place:
p = u.dot(d.dot(u.T))
p -= alpha * outer(v,v.T)
u,d = udu(p)

Warning:
Vector v is not valid after call.

TODO: add doc with derivation!
*/
yaflStatusEn yafl_math_udu_down(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v);
#endif // YAFL_MATH_H

