# What we are doing here: the idea of adaptive correction
We have a system which can be cheracterised with the following set of equatins:

<a name="1"></a>$ \mathbf{x}_k = f(\mathbf{x}_{k-1}, \mathbf{u}_k, \mathbf{w}_k) \quad (1),$

<a name="2"></a>$\mathbf{z}_k = h(\mathbf{x}_k) + \mathbf{v}_k \quad (2),$
    
Where [(1)](#1) is state transition equation and [(2)](#2) is measurement equation, $\mathbf w_k$ and $\mathbf v_k$ are process and measurement noise respectively, $\mathbf u_k$ is control vector.

We use the Extended Kalman Filter (EKF) as our main estimation algorithm and the $H_\infty$ filter as backup algorithm to deal with Kalman filter divergenсe.

# The EKF algorithm
## Prediction step
<a name="3"></a>$\mathbf{x}_{k|k-1} = f(\mathbf{x}_{k-1|k-1}, \mathbf{u}_k, 0) \quad (3),$

<a name="4"></a>$\mathbf{P}_{k|k-1} =  \mathbf{F}_{k} \mathbf{P}_{k-1|k-1} \mathbf{F}_{k}^{\top} + \mathbf{B}_{k} \mathbf{Q}_{k} \mathbf{B}_{k}^{\top} \quad (4),$

<a name="5"></a>$\mathbf{y}_{k} = \mathbf{z}_{k} - h(\mathbf{x}_{k|k-1}) \quad (5),$

## Correction step
The convenient form of correction step is:

<a name="6"></a>$\mathbf{S}_{k} = \mathbf{H}_{k} \mathbf{P}_{k|k-1} \mathbf{H}_{k}^{\top} + \mathbf{R}_{k} \quad (6),$

<a name="7"></a>$\mathbf{K}_{k} = \mathbf{P}_{k|k-1} \mathbf{H}_{k}^{\top} \mathbf{S}_{k}^{-1} \quad (7),$

<a name="8"></a>$\mathbf{x}_{k|k} = \mathbf{x}_{k|k-1} + \mathbf{K}_{k} \mathbf{y}_{k} \quad (8),$

<a name="9"></a>$\mathbf{P}_{k|k} = \left( \mathbf{I} - \mathbf{K}_{k} \mathbf{H}_{k} \right) \mathbf{P}_{k|k-1} \left( \mathbf{I} - \mathbf{K}_{k} \mathbf{H}_{k} \right)^{\top} + \mathbf{K}_{k} \mathbf{R}_{k} \mathbf{K}_{k}^{\top} \quad (9),$

To compare the EKF with the $H_\infty$ filter we need to rewrite the correction step in the following form:

<a name="10"></a>$\mathbf{P}_{k|k} = \left( \mathbf{P}_{k|k-1}^{-1} + \mathbf{H}_{k}^{\top} \mathbf{R}_{k}^{-1} \mathbf{H}_{k} \right)^{-1} \quad (10),$

<a name="11"></a>$\mathbf{K}_{k} = \mathbf{P}_{k|k} \mathbf{H}_{k}^{\top} \mathbf{R}_{k}^{-1} \quad (11),$

<a name="12"></a>$\mathbf{x}_{k|k} = \mathbf{x}_{k|k-1} + \mathbf{K}_{k} \mathbf{y}_{k} \quad (12),$

## Where:
$ \mathbf{F}_{k} = \left . \frac{\partial f}{\partial \mathbf{x} } \right \vert _{\mathbf{x}_{k-1|k-1}, \mathbf{u}_{k}}, \quad $ 
$ \mathbf{B}_{k} = \left . \frac{\partial f}{\partial \mathbf{v} } \right \vert _{\mathbf{x}_{k-1|k-1}, \mathbf{u}_{k}}, \quad $
$ \mathbf{H}_{k} = \left . \frac{\partial h}{\partial \mathbf{x} } \right \vert _{\mathbf{x}_{k|k-1}}, \quad $
$ \mathbf{Q}_{k} \ge 0 $ is the process noise covariance
$ , \quad \mathbf{R}_{k} \ge 0 $ is the measurement noise covariance


# The $H_\infty$ filter algorithm
The $H_\infty$ filter algorithm is given in [**[Banavar1992]**](#banavar) this algorithm can be written as follows:

## Prediction step
<a name="13"></a>$\mathbf{x}_{k|k-1} = f(\mathbf{x}_{k-1|k-1}, \mathbf{u}_k, 0) \quad (13),$

<a name="14"></a>$\mathbf{P}_{k|k-1} =  \mathbf{F}_{k} \mathbf{P}_{k-1|k-1} \mathbf{F}_{k}^{\top} + \mathbf{B}_{k} \mathbf{Q}_{k} \mathbf{B}_{k}^{\top} \quad (14),$

<a name="15"></a>$\mathbf{y}_{k} = \mathbf{z}_{k} - h(\mathbf{x}_{k|k-1}) \quad (15),$

## Correction step
<a name="16"></a>$\mathbf{\tilde{S}}_{k} =  \mathbf{L}_{b}^{\top} \mathbf{\bar{S}}_{k} \mathbf{L}_{k} \quad (16),$

<a name="17"></a>$\mathbf{P}_{k|k} = \left( \mathbf{P}_{k|k-1}^{-1} - \theta \cdot \mathbf{\tilde{S}}_{k} + \mathbf{H}_{k}^{\top} \mathbf{R}_{k}^{-1} \mathbf{H}_{k} \right)^{-1} \quad (17),$

<a name="18"></a>$\mathbf{K}_{k} = \mathbf{P}_{k|k} \mathbf{H}_{k}^{\top} \mathbf{R}_{k}^{-1} \quad (18),$

<a name="19"></a>$\mathbf{x}_{k|k} = \mathbf{x}_{k|k-1} + \mathbf{K}_{k} \mathbf{y}_{k} \quad (19),$

## Where:
$ \mathbf{\bar{S}}_{k} \gt 0$ is user defined matrix, 
$ \mathbf{L}_{k} $ is full rank user defined state weigthing matrix,
$ \theta $ is user defined performance bounary value.

# Filtering algorithms comparizon
If we compare equations [(10)](#10) and [(17)](#17) then we cen see that the least equation can be computed in two steps:

<a name="20"></a>$\mathbf{P}_{b} = \left( \mathbf{P}_{b-1}^{-1} - \theta \cdot \mathbf{\tilde{S}}_{b} \right)^{-1} = \mathbf{P}_{b-1} + \Delta{\mathbf{Q}_{b}} \quad (20),$

<a name="21"></a>$\mathbf{P}_{d} = \left( \mathbf{P}_{d-1}^{-1} + \mathbf{H}_{d}^{\top} \mathbf{R}_{d}^{-1} \mathbf{H}_{d} \right)^{-1} \quad (21),$

The order of these steps may vary, but the result should be ste same. 

The step [(21)](#21) is exactly the same as [(10)](#10) which is used in the Kalman filter state covariance update.

The step [(20)](#20) can be viewed as an additional process noise. Note that $\Delta \mathbf{Q}_{b}$ must be positive defnite for the existence of the filter.

# Adaptive correction derivation
## Divergenсe criterion
Kalman filter divergence is detected when the folowing criterion is met:

<a name="22"></a> $ \mathbf{y}_{b}^{\top} \mathbf{S}_{b}^{\mathsf{-1}} \mathbf{y}_{b} \gt {\beta}_{n} \quad (22),$

where ${\beta}_{n}$ is threshold ${\chi}_{\alpha,n}^{2}$ value.

## Correction: Initial derivation

After adaptive correction step case we get:

<a name="23"></a> $ \mathbf{y}_{b}^{\top} \mathbf{S}_{corr}^{\mathsf{-1}} \mathbf{y}_{b} = {\beta}_{n} \quad (23),$

By multimlication of the equation [(23)](#23) by $\mathbf{y}_{b}$ on the left and by $\mathbf{y}_{b} ^ {\top}$ on the rigth sides we get:

<a name="24"></a> $ \mathbf{y}_{b} \mathbf{y}_{b}^{\top} \mathbf{S}_{corr}^{\mathsf{-1}} \mathbf{y}_{b} \mathbf{y}_{b}^{\top} = {\beta}_{n} \cdot \mathbf{y}_{b} \mathbf{y}_{b}^{\top} \quad (24),$

By compating left and igth sides of the equation [(24)](#24) we get:

$ \mathbf{y}_{b} {y}_{b}^{\top} \mathbf{S}_{corr}^{\mathsf{-1}} = {\beta}_{n} \cdot \mathbf{I} ,$

or

$ \mathbf{y}_{b} \mathbf{y}_{b}^{\top} = {\beta}_{n} \cdot \mathbf{S}_{corr} ,$

or

<a name="25"></a> $ \mathbf{S}_{corr} = \cfrac {\mathbf{y}_{b} \mathbf{y}_{b}^{\top}} {{\beta}_{n}} \quad (25),$

## Correction: No limits approach

By substitution of [(20)](#20) rigth side to [(16)](#16) we get:

<a name="26"></a> $\mathbf{S}_{corr} = {S}_{b} + \mathbf{H}_{b} {\Delta \mathbf{Q}_{b}} \mathbf{H}_{b}^{\top} \quad (26),$

and so substitution [(26)](#26) to [(25)](#25) gives us:

<a name="27"></a> $ {S}_{b} + \mathbf{H}_{b} {\Delta \mathbf{Q}_{b}} \mathbf{H}_{b}^{\top} = \cfrac {\mathbf{y}_{b} \mathbf{y}_{b}^{\top}} {{\beta}_{n}} \quad (27),$

or 

<a name="28"></a> $ {A}_{b} = \mathbf{H}_{b} {\Delta \mathbf{Q}_{b}} \mathbf{H}_{b}^{\top} = \cfrac {\mathbf{y}_{b} \mathbf{y}_{b}^{\top}} {{\beta}_{n}} - {S}_{b} \quad (28),$

on the other hand [(20)](#20) gives us:

$ \mathbf{P}_{b} = \left( \mathbf{P}_{b-1}^{-1} - \theta \cdot \mathbf{\tilde{S}}_{b} \right)^{-1} = \mathbf{P}_{b-1} + \Delta \mathbf{Q}_{b} $

or 

$ \mathbf{P}_{b} = \left( \mathbf{P}_{b-1}^{-1} - \theta \cdot \mathbf{L}_{b}^{\top} \mathbf{\bar{S}}_{b} \mathbf{L}_{b} \right)^{-1} = \left( \mathbf{P}_{b-1}^{-1} - \mathbf{L}_{b}^{\top} \mathbf{M}_{b}^{-1} \mathbf{L}_{b} \right)^{-1} \quad ,$ 
where $\mathbf{M}_{b}^{-1} =  \theta \cdot \mathbf{\bar{S}}_{b} $

or 

$ \mathbf{P}_{b} = \mathbf{P}_{b-1} + \mathbf{P}_{b-1}  \mathbf{L}_{b}^{\top} \left( \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right)^{-1} \mathbf{L}_{b} \mathbf{P}_{b-1} $

which means that:

<a name="29"></a> $ {\Delta \mathbf{Q}_{b}} = \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \left( \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right)^{-1} \mathbf{L}_{b} \mathbf{P}_{b-1} \quad (29),$

from [(28)](#28) we get:

$ \mathbf{A}_{b} = \mathbf{H}_{b} {\Delta \mathbf{Q}_{b}} \mathbf{H}_{b}^{\top} = \mathbf{H}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \left( \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right)^{-1} \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{H}_{b}^{\top} = \mathbf{C}_{b} \left( \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right)^{-1} \mathbf{C}_{b}^{\top} $

or 

$ \mathbf{A}_{b} = \mathbf{C}_{b} \left( \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right)^{-1} \mathbf{C}_{b}^{\top} $

or 

$ \mathbf{C}_{b}^{+} \mathbf{A}_{b} = \left( \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right)^{-1} \mathbf{C}_{b}^{\top} $

or

$ \left( \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right) \mathbf{C}_{b}^{+} = \mathbf{C}_{b}^{\top} \mathbf{A}_{b}^{-1}$

or

$ \mathbf{M}_{b} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} = \mathbf{C}_{b}^{\top} \mathbf{A}_{b}^{-1} \mathbf{C}_{b}$

or

<a name="30"></a> $ \mathbf{M}_{b} = \mathbf{C}_{b}^{\top} \mathbf{A}_{b}^{-1} \mathbf{C}_{b} + \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \quad (30),$


where:

<a name="31"></a> $ \mathbf{C}_{b} = \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{H}_{b}^{\top} \quad (31),$

so for $\Delta \mathbf{Q}_{b}$ we have:

$ {\Delta \mathbf{Q}_{b}} = \mathbf{P}_{b-1}  \mathbf{L}_{b}^{\top} \left( \mathbf{C}_{b}^{\top} \mathbf{A}_{b}^{-1} \mathbf{C}_{b} + \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} - \mathbf{L}_{b} \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \right)^{-1} \mathbf{L}_{b} \mathbf{P}_{b-1}  = \mathbf{P}_{b-1}  \mathbf{L}_{b}^{\top} \left( \mathbf{C}_{b}^{\top} \mathbf{A}_{b}^{-1} \mathbf{C}_{b} \right)^{-1} \mathbf{L}_{b} \mathbf{P}_{b-1} $

or 

<a name="32"></a> $ {\Delta \mathbf{Q}_{b}} = \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \mathbf{C}_{b}^{+} \mathbf{A}_{b} \mathbf{C}_{b}^{+ \top} \mathbf{L}_{b} \mathbf{P}_{b-1}  \quad (32),$

The Kalman filter with $H_\infty$ correction [(32)](#32) exists when:

$ {\Delta \mathbf{Q}_{b}} = \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \mathbf{C}_{b}^{+} \mathbf{A}_{b} \mathbf{C}_{b}^{+ \top} \mathbf{L}_{b} \mathbf{P}_{b-1}  \gt 0 ,$

According to [**[Horn et al]**](#horn_et_al) the above citeria is met when:
 * Maxtix $ \mathbf{Z}_{b} = \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \mathbf{C}_{b}^{+} $ has full rank and
 * $ \mathbf{A}_{b} \gt 0 $

Let us consider $ \mathbf{Z}_{b} = \mathbf{P}_{b-1} \mathbf{L}_{b}^{\top} \mathbf{C}_{b}^{+} $ :
 * as $ \mathbf{P}_{b-1} \gt 0 $ and $ \mathbf{L}_{b} $ has full rank then $\mathbf{P}_{b-1} \mathbf{L}_{b}$ also has full rank;
 
 * as both $\mathbf{H}_{b}$ and $\mathbf{L}_{b}$ have full rank, аnd $ \mathbf{P}_{b-1} \gt 0$, then $\mathbf{C}_{b-1}$ also has full rank, consequently: the
product $ \mathbf{C}_{b}^{\top} \mathbf{C}_{b} \gt 0 $ (see [**[Horn et al]**](#horn_et_al)), and so, $ \mathbf{C}_{b}^{+} = \left( \mathbf{C}_{b}^{\top} \mathbf{C}_{b} \right)^{-1} \mathbf{C}_{b}^{\top} $ also have full rank;

 * as both product $\mathbf{P}_{b-1} \mathbf{L}_{b}$ and the matrix $ \mathbf{C}_{b}^{+} $ have full rank, then $ \mathbf{Z}_{b} $ also has full rank.
 
Let us consider the matrix A b :
 * $\mathbf{S}_{b} \ge \mathbf{H}_{b} {\Delta \mathbf{P}_{b-1}} \mathbf{H}_{b}^{\top} \gt 0$ 
 
 * now we must consider the product $ \mathbf{x}^{\top} \left( \mathbf{y}_{b} \mathbf{y}_{b}^{\top} \right) \mathbf{x} $ for $ \mathbf{x} \ne 0 $:
  * if the dimensionality of x is greater than one, then $ \mathbf{x}^{\top} \left( \mathbf{y}_{b} \mathbf{y}_{b}^{\top} \right) \mathbf{x} = \left( \mathbf{x}^{\top} \mathbf{y}_{b} \right) \left( \mathbf{y}_{b}^{\top} \mathbf{x} \right) \ge 0 $, consequently, regarding to positive definiteness of S b , we can not guarante the existense of the filter;
  * if the dimensionality of x is one (when we implement sequential filter), then $ \mathbf{x}^{\top} \left( \mathbf{y}_{b} \mathbf{y}_{b}^{\top} \right) \mathbf{x} = {y}_{b}^{2} {x}^(2) \gt 0 $ ,
    * in such case the equation [(28)](#28) reduces to ${a}_{b} = \cfrac{{y}_{b}^{2}}{{\beta}_{n}} - {s}_{b}$
    * the criterion [(22)](#22) is met when $ \cfrac{{y}_{b}^{2}}{{\beta}_{1}} \gt {s}_{b} $ or $\cfrac{{y}_{b}^{2}}{{\beta}_{n}} - {s}_{b} \gt 0$
    * so when thecriterion [(22)](#22) is met we have $ \mathbf{A}_{b} = {a}_{b} \gt 0 $, 
    
    thus trere is at least one sequential filter implementation.

So we can select  $ \mathbf{L}_{b} $, $ \mathbf{\bar{S}}_{b} $, $ \theta $ parameters of $H_\infty$ filter so that we can use operation [(20)](#20) to correct the
Kalman filter divergence. Our selection is reduced to definition of $ \mathbf{L}_{b} $ and {y}_{b} to be used in certain
filtering algorithm.


# References

<a name="banavar"></a> **[Banavar1992]** R. Banavar, “A game theoretic approach to linear dynamic estimation”, Doctoral Dissertation,
University of Texas at Austin, May 1992.

<a name="horn_et_al"></a> **[Horn et al]** R.A. Horn, C.R Johnson, Charles R. “Matrix Analysis (2nd ed.).”, Cambridge University Press,
2013, ISBN 978-0-521-38632-6.
