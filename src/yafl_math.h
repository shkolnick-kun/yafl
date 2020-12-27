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

#ifndef YAFL_MATH_H_
#define YAFL_MATH_H_

/*Config file must be included at first place*/
#include <yafl_config.h>

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
--------------------------------------------------------------------------------------------
         Function/Macro (May be used for in place processing)                  NumPy expr
------------------------------------------------------------------------------------------*/
void yaflm_set_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res  = n*v */
void yaflm_add_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res += n*v */
void yaflm_sub_vxn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res -= n*v */

void yaflm_set_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res  = v/n */
void yaflm_add_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res += v/n */
void yaflm_sub_vrn(yaflInt sz, yaflFloat *res, yaflFloat *v, yaflFloat n);  /* res -= v/n */

void yaflm_set_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a*b */
void yaflm_add_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a*b */
void yaflm_sub_vxv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a*b */

void yaflm_set_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a/b */
void yaflm_add_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a/b */
void yaflm_sub_vrv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a/b */

#define YAFLM_SET_DV yaflm_set_vxv                                          /* res  = d.dot(v) */
#define YAFLM_ADD_DV yaflm_add_vxv                                          /* res += d.dot(v) */
#define YAFLM_SUB_DV yaflm_sub_vxv                                          /* res -= d.dot(v) */

#define YAFLM_SET_RDV(sz, res, a, b) yaflm_set_vrv(sz, res, b, a)           /* res  = (1/d).dot(v)*/
#define YAFLM_ADD_RDV(sz, res, a, b) yaflm_add_vrv(sz, res, b, a)           /* res += (1/d).dot(v)*/
#define YAFLM_SUB_RDV(sz, res, a, b) yaflm_sub_vrv(sz, res, b, a)           /* res += (1/d).dot(v)*/

#define YAFLM_SET_VTD yaflm_set_vxv                                         /* res  = v.T.dot(d) */
#define YAFLM_ADD_VTD yaflm_add_vxv                                         /* res += v.T.dot(d) */
#define YAFLM_SUB_VTD yaflm_sub_vxv                                         /* res -= v.T.dot(d) */

#define YAFLM_SET_VTRD(sz, res, a, b) yaflm_set_vrv(sz, res, b, a)          /* res  = v.T.dot(1/d)*/
#define YAFLM_ADD_VTRD(sz, res, a, b) yaflm_add_vrv(sz, res, b, a)          /* res += v.T.dot(1/d)*/
#define YAFLM_SUB_VTRD(sz, res, a, b) yaflm_sub_vrv(sz, res, b, a)          /* res += v.T.dot(1/d)*/

/*
m - rectangular matrix
v - vector

nc - number of columns in a matrix
nr - number of rows in a matrix
-----------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                               NumPy expr
---------------------------------------------------------------------------------------------------------------------------------*/
yaflFloat yaflm_vtv(yaflInt sz, yaflFloat *a, yaflFloat *b);                                           /* a.T.dot(b)             */

void yaflm_set_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);                /* res  = outer(a, b)     */
void yaflm_add_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);                /* res += outer(a, b)     */
void yaflm_sub_vvt(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b);                /* res -= outer(a, b)     */

void yaflm_set_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n); /* res  = outer(a, b) * n */
void yaflm_add_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n); /* res += outer(a, b) * n */
void yaflm_sub_vvtxn(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b, yaflFloat n); /* res -= outer(a, b) * n */

/*
m - rectangular matrix
v - vector

nc - number of columns in a matrix
nr - number of rows in a matrix
------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                         NumPy expr
----------------------------------------------------------------------------------------------------------*/

void yaflm_set_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a.dot(b) */
void yaflm_add_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a.dot(b) */
void yaflm_sub_mv(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a.dot(b) */

void yaflm_set_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a.T.dot(b) */
void yaflm_add_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a.T.dot(b) */
void yaflm_sub_vtm(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a.T.dot(b) */

/*ncr - number of columns of in @a and of rows in @b*/
void yaflm_set_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a.dot(b) */
void yaflm_add_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a.dot(b) */
void yaflm_sub_mm(yaflInt nr,  yaflInt ncr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a.dot(b) */

/*
m - rectangular matrix
u - upper triangular unit matrix
v - vector
--------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                           NumPy expr
------------------------------------------------------------------------------------------------------------*/
void yaflm_set_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);            /* res  = a.T.dot(b) */
void yaflm_add_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);            /* res += a.T.dot(b) */
void yaflm_sub_vtu(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);            /* res -= a.T.dot(b) */

void yaflm_set_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);             /* res  = a.dot(b)   */
void yaflm_add_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);             /* res += a.dot(b)   */
void yaflm_sub_uv(yaflInt sz, yaflFloat *res, yaflFloat *a, yaflFloat *b);             /* res -= a.dot(b)   */

void yaflm_set_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res  = a.dot(b)   */
void yaflm_add_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res += a.dot(b)   */
void yaflm_sub_mu(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *a, yaflFloat *b); /* res -= a.dot(b)   */

void yaflm_set_u(yaflInt sz, yaflFloat *res, yaflFloat *u);                            /* res  = u          */
void yaflm_add_u(yaflInt sz, yaflFloat *res, yaflFloat *u);                            /* res += u          */
void yaflm_sub_u(yaflInt sz, yaflFloat *res, yaflFloat *u);                            /* res -= u          */

/*Matrix row size and block start address*/
#define YAFL_BLK(m,nc,r,c) nc, ((yaflFloat *)m + nc * r + c)
/*
Block operations:

m - rectangular matrix
u - upper triangular unit matrix
v - vector
----------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                           NumPy expr
--------------------------------------------------------------------------------------------------------------*/
void yaflm_bset_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);            /* res[:sz, :sz]  = u     */
void yaflm_badd_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);            /* res[:sz, :sz] += u     */
void yaflm_bsub_u(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);            /* res[:sz, :sz] -= u     */

#define YAFLM_BSET_U(nc, r, c, m, sz, u) yaflm_bset_u(YAFL_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz]  = u */
#define YAFLM_BADD_U(nc, r, c, m, sz, u) yaflm_badd_u(YAFL_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz] += u */
#define YAFLM_BSUB_U(nc, r, c, m, sz, u) yaflm_bsub_u(YAFL_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz] -= u */

void yaflm_bset_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);            /* res[:sz, :sz]  = u.T   */
void yaflm_badd_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);            /* res[:sz, :sz] += u.T   */
void yaflm_bsub_ut(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *u);            /* res[:sz, :sz] -= u.T   */

#define YAFLM_BSET_UT(nc, r, c, m, sz, u) yaflm_bset_u(YAFL_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz]  = u.T */
#define YAFLM_BADD_UT(nc, r, c, m, sz, u) yaflm_badd_u(YAFL_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz] += u.T */
#define YAFLM_BSUB_UT(nc, r, c, m, sz, u) yaflm_bsub_u(YAFL_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz] -= u.T */

void yaflm_bset_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v);            /* res[:sz, 0]  = v    */
void yaflm_badd_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v);            /* res[:sz, 0] += v    */
void yaflm_bsub_v(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *v);            /* res[:sz, 0] -= v    */

#define YAFLM_BSET_V(nc, r, c, m, sz, v) yaflm_bset_v(YAFL_BLK(m,nc,r,c), sz, v)    /* m[r:r+sz, c]  = v   */
#define YAFLM_BADD_V(nc, r, c, m, sz, v) yaflm_badd_v(YAFL_BLK(m,nc,r,c), sz, v)    /* m[r:r+sz, c] += v   */
#define YAFLM_BSUB_V(nc, r, c, m, sz, v) yaflm_bsub_v(YAFL_BLK(m,nc,r,c), sz, v)    /* m[r:r+sz, c] -= v   */

/*
Block operations:

m - rectangular matrix
v - vector
---------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                       NumPy expr
-------------------------------------------------------------------------------------------------------------------------------*/
void yaflm_bset_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b);   /* res[:sz, :sz]  = outer(a, b)     */
void yaflm_badd_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b);   /* res[:sz, :sz] += outer(a, b)     */
void yaflm_bsub_vvt(yaflInt nc, yaflFloat *res, yaflInt sz, yaflFloat *a, yaflFloat *b);   /* res[:sz, :sz] -= outer(a, b)     */

#define YAFLM_BSET_VVT(nc, r, c, m, sz, a, b) yaflm_bset_vvt(YAFL_BLK(m,nc,r,c), sz, a, b) /* m[r:r+sz, c:c+sz]  = outer(a, b) */
#define YAFLM_BADD_VVT(nc, r, c, m, sz, a, b) yaflm_badd_vvt(YAFL_BLK(m,nc,r,c), sz, a, b) /* m[r:r+sz, c:c+sz] += outer(a, b) */
#define YAFLM_BSUB_VVT(nc, r, c, m, sz, a, b) yaflm_bsub_vvt(YAFL_BLK(m,nc,r,c), sz, a, b) /* m[r:r+sz, c:c+sz] -= outer(a, b) */

/*
Block operations:

m - rectangular matrix
u - upper triangular unit matrix

rnc - number of columns in result matrices
----------------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                                NumPy expr
--------------------------------------------------------------------------------------------------------------------------------------*/
void yaflm_bset_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b); /* res[:nr, :nc]  = a.dot(b)     */
void yaflm_badd_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b); /* res[:nr, :nc] += a.dot(b)     */
void yaflm_bsub_mu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflFloat *a, yaflFloat *b); /* res[:nr, :nc] -= a.dot(b)     */

#define YAFLM_BSET_MU(rnc, r, c, m, nr, nc, a, b) yaflm_bset_mu(YAFL_BLK(m,rnc,r,c), nr, nc, a, b)   /* m[r:r+nr, c:c+nc]  = a.dot(b) */
#define YAFLM_BADD_MU(rnc, r, c, m, nr, nc, a, b) yaflm_badd_mu(YAFL_BLK(m,rnc,r,c), nr, nc, a, b)   /* m[r:r+nr, c:c+nc] += a.dot(b) */
#define YAFLM_BSUB_MU(rnc, r, c, m, nr, nc, a, b) yaflm_bsub_mu(YAFL_BLK(m,rnc,r,c), nr, nc, a, b)   /* m[r:r+nr, c:c+nc] -= a.dot(b) */

/*
Block operations:

u - upper triangular unit matrix

rnc - number of columns in result matrices
anc - number of columns in a matrices blocks

Warning:
Matrix blocks a and res must not overlap!
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                 Function/Macro                                                                                           NumPy expr
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
void yaflm_bset_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b);                   /* res[:nr, :nc]  = a[:nr, :nc].dot(b)               */
void yaflm_badd_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b);                   /* res[:nr, :nc] += a[:nr, :nc].dot(b)               */
void yaflm_bsub_bu(yaflInt rnc, yaflFloat *res, yaflInt nr, yaflInt nc, yaflInt anc, yaflFloat *a, yaflFloat *b);                   /* res[:nr, :nc] -= a[:nr, :nc].dot(b)               */

#define YAFLM_BSET_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yaflm_bset_bu(YAFL_BLK(m,rnc,r,c), nr, nc, YAFL_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc]  = a[ar:ar+nr, ac:ac+nc].dot(b) */
#define YAFLM_BADD_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yaflm_badd_bu(YAFL_BLK(m,rnc,r,c), nr, nc, YAFL_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc] += a[ar:ar+nr, ac:ac+nc].dot(b) */
#define YAFLM_BSUB_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yaflm_bsub_bu(YAFL_BLK(m,rnc,r,c), nr, nc, YAFL_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc] -= a[ar:ar+nr, ac:ac+nc].dot(b) */

/*
Back substitution:

m - rectangular matrix
u - upper triangular unit matrix
v - vector
------------------------------------------------------------------------------------------------------
                                   Function/Macro                                 NumPy expr
----------------------------------------------------------------------------------------------------*/
/*res is vector*/
void yaflm_ruv(yaflInt sz, yaflFloat *res, yaflFloat *u);             /*res = linalg.inv(u).dot(res)*/
/*res is matrix*/
void yaflm_rum(yaflInt nr, yaflInt nc, yaflFloat *res, yaflFloat *u); /*res = linalg.inv(u).dot(res)*/

/*
Modified Weighted Gram-Schmidt update (for UDU' decomposition)

Does in place:
res_u, res_d = MWGS(w,d)

Warning:
Matrix w is not valid after call.
*/
void yaflm_mwgsu(yaflInt nr, yaflInt nc, yaflFloat *res_u, yaflFloat *res_d, yaflFloat *w, yaflFloat *d);

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
void yaflm_udu_up(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v);

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
void yaflm_udu_down(yaflInt sz, yaflFloat *res_u, yaflFloat *res_d, yaflFloat alpha, yaflFloat *v);
#endif // YAFL_MATH_H_
