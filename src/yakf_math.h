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

#ifndef YAKF_MATH_H_
#define YAKF_MATH_H_

/*Config file must be included at first place*/
#include <yakf_config.h>

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
--------------------------------------------------------------------------------------------
                        Function/Macro                                         NumPy expr
------------------------------------------------------------------------------------------*/
void yakfm_set_nv(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n);   /* res  = n*v */
void yakfm_add_nv(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n);   /* res += n*v */
void yakfm_sub_nv(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n);   /* res -= n*v */

void yakfm_set_vrn(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n);  /* res  = v/n */
void yakfm_add_vrn(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n);  /* res += v/n */
void yakfm_sub_vrn(yakfInt sz, yakfFloat *res, yakfFloat *v, yakfFloat n);  /* res -= v/n */

void yakfm_set_vxv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res  = a*b */
void yakfm_add_vxv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res += a*b */
void yakfm_sub_vxv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res -= a*b */

void yakfm_set_vrv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res  = a/b */
void yakfm_add_vrv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res += a/b */
void yakfm_sub_vrv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res -= a/b */

#define YAKFM_SET_DV yakfm_set_vxv                                          /* res  = d.dot(v) */
#define YAKFM_ADD_DV yakfm_add_vxv                                          /* res += d.dot(v) */
#define YAKFM_SUB_DV yakfm_sub_vxv                                          /* res -= d.dot(v) */

#define YAKFM_SET_RDV(sz, res, a, b) yakfm_set_vrv(sz, res, b, a)           /* res  = (1/d).dot(v)*/
#define YAKFM_ADD_RDV(sz, res, a, b) yakfm_add_vrv(sz, res, b, a)           /* res += (1/d).dot(v)*/
#define YAKFM_SUB_RDV(sz, res, a, b) yakfm_sub_vrv(sz, res, b, a)           /* res += (1/d).dot(v)*/

#define YAKFM_SET_VTD yakfm_set_vxv                                         /* res  = v.T.dot(d) */
#define YAKFM_ADD_VTD yakfm_add_vxv                                         /* res += v.T.dot(d) */
#define YAKFM_SUB_VTD yakfm_sub_vxv                                         /* res -= v.T.dot(d) */

#define YAKFM_SET_VTRD(sz, res, a, b) yakfm_set_vrv(sz, res, b, a)          /* res  = v.T.dot(1/d)*/
#define YAKFM_ADD_VTRD(sz, res, a, b) yakfm_add_vrv(sz, res, b, a)          /* res += v.T.dot(1/d)*/
#define YAKFM_SUB_VTRD(sz, res, a, b) yakfm_sub_vrv(sz, res, b, a)          /* res += v.T.dot(1/d)*/

yakfFloat yakfm_vtv(yakfInt sz, yakfFloat *a, yakfFloat *b);                /* a.T.dot(b) */

void yakfm_set_vvt(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res  = outer(a, b) */
void yakfm_add_vvt(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res += outer(a, b) */
void yakfm_sub_vvt(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res -= outer(a, b) */
/*
m - rectangular matrix
v - vector

nc - number of columns in a matrix
nr - number of rows in a matrix
------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                         NumPy expr
----------------------------------------------------------------------------------------------------------*/
void yakfm_set_mv(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res  = a.dot(b) */
void yakfm_add_mv(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res += a.dot(b) */
void yakfm_sub_mv(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res -= a.dot(b) */

void yakfm_set_vtm(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res  = a.T.dot(b) */
void yakfm_add_vtm(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res += a.T.dot(b) */
void yakfm_sub_vtm(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res -= a.T.dot(b) */

/*ncr - number of columns of in @a and of rows in @b*/
void yakfm_set_mm(yakfInt nr,  yakfInt ncr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res  = a.dot(b) */
void yakfm_add_mm(yakfInt nr,  yakfInt ncr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res += a.dot(b) */
void yakfm_sub_mm(yakfInt nr,  yakfInt ncr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res -= a.dot(b) */

/*
m - rectangular matrix
u - upper triangular unit matrix
v - vector
--------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                           NumPy expr
------------------------------------------------------------------------------------------------------------*/
void yakfm_set_vtu(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b);            /* res  = a.T.dot(b) */
void yakfm_add_vtu(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b);            /* res += a.T.dot(b) */
void yakfm_sub_vtu(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b);            /* res -= a.T.dot(b) */

void yakfm_set_uv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b);             /* res  = a.dot(b)   */
void yakfm_add_uv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b);             /* res += a.dot(b)   */
void yakfm_sub_uv(yakfInt sz, yakfFloat *res, yakfFloat *a, yakfFloat *b);             /* res -= a.dot(b)   */

void yakfm_set_mu(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res  = a.dot(b)   */
void yakfm_add_mu(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res += a.dot(b)   */
void yakfm_sub_mu(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *a, yakfFloat *b); /* res -= a.dot(b)   */

void yakfm_set_u(yakfInt sz, yakfFloat *res, yakfFloat *u);                            /* res  = u          */
void yakfm_add_u(yakfInt sz, yakfFloat *res, yakfFloat *u);                            /* res += u          */
void yakfm_sub_u(yakfInt sz, yakfFloat *res, yakfFloat *u);                            /* res -= u          */

/*Matrix row size and block start address*/
#define YAKF_BLK(m,nc,r,c) nc, ((yakfFloat *)m + nc * r + c)
/*
Block operations:

m - rectangular matrix
u - upper triangular unit matrix
v - vector
----------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                           NumPy expr
--------------------------------------------------------------------------------------------------------------*/
void yakfm_bset_u(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *u);            /* res[:sz, :sz]  = u     */
void yakfm_badd_u(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *u);            /* res[:sz, :sz] += u     */
void yakfm_bsub_u(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *u);            /* res[:sz, :sz] -= u     */

#define YAKFM_BSET_U(nc, r, c, m, sz, u) yakfm_bset_u(YAKF_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz]  = u */
#define YAKFM_BADD_U(nc, r, c, m, sz, u) yakfm_badd_u(YAKF_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz] += u */
#define YAKFM_BSUB_U(nc, r, c, m, sz, u) yakfm_bsub_u(YAKF_BLK(m,nc,r,c), sz, u)    /* m[r:r+sz, c:c+sz] -= u */

void yakfm_bset_v(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *v);            /* res[:sz, 0]  = v    */
void yakfm_badd_v(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *v);            /* res[:sz, 0] += v    */
void yakfm_bsub_v(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *v);            /* res[:sz, 0] -= v    */

#define YAKFM_BSET_V(nc, r, c, m, sz, v) yakfm_bset_v(YAKF_BLK(m,nc,r,c), sz, v)    /* m[r:r+sz, c]  = v   */
#define YAKFM_BADD_V(nc, r, c, m, sz, v) yakfm_badd_v(YAKF_BLK(m,nc,r,c), sz, v)    /* m[r:r+sz, c] += v   */
#define YAKFM_BSUB_V(nc, r, c, m, sz, v) yakfm_bsub_v(YAKF_BLK(m,nc,r,c), sz, v)    /* m[r:r+sz, c] -= v   */

/*
Block operations:

m - rectangular matrix
v - vector
---------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                       NumPy expr
-------------------------------------------------------------------------------------------------------------------------------*/
void yakfm_bset_vvt(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *a, yakfFloat *b);   /* res[:sz, :sz]  = outer(a, b)     */
void yakfm_badd_vvt(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *a, yakfFloat *b);   /* res[:sz, :sz] += outer(a, b)     */
void yakfm_bsub_vvt(yakfInt nc, yakfFloat *res, yakfInt sz, yakfFloat *a, yakfFloat *b);   /* res[:sz, :sz] -= outer(a, b)     */

#define YAKFM_BSET_VVT(nc, r, c, m, sz, a, b) yakfm_bset_vvt(YAKF_BLK(m,nc,r,c), sz, a, b) /* m[r:r+sz, c:c+sz]  = outer(a, b) */
#define YAKFM_BADD_VVT(nc, r, c, m, sz, a, b) yakfm_badd_vvt(YAKF_BLK(m,nc,r,c), sz, a, b) /* m[r:r+sz, c:c+sz] += outer(a, b) */
#define YAKFM_BSUB_VVT(nc, r, c, m, sz, a, b) yakfm_bsub_vvt(YAKF_BLK(m,nc,r,c), sz, a, b) /* m[r:r+sz, c:c+sz] -= outer(a, b) */

/*
Block operations:

m - rectangular matrix
u - upper triangular unit matrix

rnc - number of columns in result matrices
----------------------------------------------------------------------------------------------------------------------------------------
                                   Function/Macro                                                                NumPy expr
--------------------------------------------------------------------------------------------------------------------------------------*/
void yakfm_bset_mu(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfFloat *a, yakfFloat *b); /* res[:nr, :nc]  = a.dot(b)     */
void yakfm_badd_mu(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfFloat *a, yakfFloat *b); /* res[:nr, :nc] += a.dot(b)     */
void yakfm_bsub_mu(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfFloat *a, yakfFloat *b); /* res[:nr, :nc] -= a.dot(b)     */

#define YAKFM_BSET_MU(rnc, r, c, m, nr, nc, a, b) yakfm_bset_mu(YAKF_BLK(m,rnc,r,c), nr, nc, a, b)   /* m[r:r+nr, c:c+nc]  = a.dot(b) */
#define YAKFM_BADD_MU(rnc, r, c, m, nr, nc, a, b) yakfm_badd_mu(YAKF_BLK(m,rnc,r,c), nr, nc, a, b)   /* m[r:r+nr, c:c+nc] += a.dot(b) */
#define YAKFM_BSUB_MU(rnc, r, c, m, nr, nc, a, b) yakfm_bsub_mu(YAKF_BLK(m,rnc,r,c), nr, nc, a, b)   /* m[r:r+nr, c:c+nc] -= a.dot(b) */

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
void yakfm_bset_bu(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfInt anc, yakfFloat *a, yakfFloat *b);                   /* res[:nr, :nc]  = a[:nr, :nc].dot(b)               */
void yakfm_badd_bu(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfInt anc, yakfFloat *a, yakfFloat *b);                   /* res[:nr, :nc] += a[:nr, :nc].dot(b)               */
void yakfm_bsub_bu(yakfInt rnc, yakfFloat *res, yakfInt nr, yakfInt nc, yakfInt anc, yakfFloat *a, yakfFloat *b);                   /* res[:nr, :nc] -= a[:nr, :nc].dot(b)               */

#define YAKFM_BSET_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yakfm_bset_bu(YAKF_BLK(m,rnc,r,c), nr, nc, YAKF_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc]  = a[ar:ar+nr, ac:ac+nc].dot(b) */
#define YAKFM_BADD_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yakfm_badd_bu(YAKF_BLK(m,rnc,r,c), nr, nc, YAKF_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc] += a[ar:ar+nr, ac:ac+nc].dot(b) */
#define YAKFM_BSUB_BU(rnc, r, c, m, nr, nc, anc, ar, ac, a, b) yakfm_bsub_bu(YAKF_BLK(m,rnc,r,c), nr, nc, YAKF_BLK(a,anc,ar,ac), b) /* m[r:r+nr, c:c+nc] -= a[ar:ar+nr, ac:ac+nc].dot(b) */

/*
Back substitution:

m - rectangular matrix
u - upper triangular unit matrix
v - vector
------------------------------------------------------------------------------------------------------
                                   Function/Macro                                 NumPy expr
----------------------------------------------------------------------------------------------------*/
/*res is vector*/
void yakfm_ruv(yakfInt sz, yakfFloat *res, yakfFloat *u);             /*res = linalg.inv(u).dot(res)*/
/*res is matrix*/
void yakfm_rum(yakfInt nr, yakfInt nc, yakfFloat *res, yakfFloat *u); /*res = linalg.inv(u).dot(res)*/

/*
Modified Weighted Gram-Schmidt update (for UDU' decomposition)

Does in place:
res_u, res_d = MWGS(w,d)

Warning:
Matrix w is not valid after call.
*/
void yakfm_mwgsu(yakfInt nr, yakfInt nc, yakfFloat *res_u, yakfFloat *res_d, yakfFloat *w, yakfFloat *d);

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
void yakfm_udu_up(yakfInt sz, yakfFloat *res_u, yakfFloat *res_d, yakfFloat alpha, yakfFloat *v);

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
void yakfm_udu_down(yakfInt sz, yakfFloat *res_u, yakfFloat *res_d, yakfFloat alpha, yakfFloat *v);
#endif // YAKF_MATH_H_
