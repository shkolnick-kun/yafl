# Yet Another Filtering library

**YAFL** means Yet Another Filtering Library. Our library is in **aplha** stage. So, if you need some mature lib then you should consider the solutions listed below.

There sre several libraries which implement Kalman filters for, e.g.:
* [TinyEKF](https://github.com/simondlevy/TinyEKF) which is intended for usage on FPU enabled platforms;
* [libfixkalman](https://github.com/sunsided/libfixkalman) which can be used without FPU.

There are also libraries for python:
* [FilterPy](https://github.com/rlabbe/filterpy);
* [pykalman](https://github.com/pykalman/pykalman).

## The library
Technically speaking all filters in YAFL are adaptive since all of them have at least a measurement noice covariance adjustment.
The term **Adaptive** is used in our docs for **Kalman filter** variants with H-infinity divergence correction.

In YAFL you can find these **Kalman filter** variants:

| Algorithm family | Basic        | Adaptive     | Robust       | Adaptive robust |
| :--------------- | :----------: | :----------: | :----------: | --------------: |
| **SUD EKF**      | ✓            | ✓            | ✓            | ✓               |
| **SUD UKF**      | ✓            | ✓            | ✓            | ✓               |
| **UD UKF**       | ✓            | ✓            |              |                 |

where:
* **SUD** means Sequential UD-factorized
* **UD**  means UD-factorized
* **EKF** means Extended Kalman Filter
* **UKF** means Unscented Kalman Filter
* **Basic** means basic algorithm
* **Adaptive** means a **Kalman filter** with [**adaptive divergence correction**](./doc/AdaptiveCorrection.pdf). We use H-infinity filter to correct the divergence
* **Robust**   means **Robustified Kalman filter**, see [**West1981**](#west_1981)

For all **EKF** variants we have **Bierman** and **Joseph** updates.
For sequential UD-factorized **UKF** only **Bierman** updates have been implemented.

And yes, we [**can actually**](./doc/UsingEKFTricksWithSPKF.pdf) use **EKF** tricks with **UKF**!

## Notes on process and measurement noice covariance adjustments
We used [this paper](https://arxiv.org/pdf/1702.00884.pdf) to implement optional Q and R adaptive adjustments.
Here are som notes on our implementation:
* All filters in this lib have the optional measurement noice covariance adjustment which can be enabled by setting `rff` to a small positive number e.g. `1e-4`.
* All EKF filters in this lib have the optional process noice covariance adjustment which can be enabled by setting `qff` to a small positive number e.g. `1e-4`.
* None of UKF filters have the optional process noice covariance adjustment as it leads to filters instability.

## Notes on implementation
The library is written in C and is intended for embedded systems usage:
* We use static memory allocation
* We use cache-friendly algorithms when available.
* Regularization techniques are used when necessary. The code is numerically stable.
* Depends only on C standard library.

## Using with C
To use the library you need to:
* go to our [Releases page](https://github.com/shkolnick-kun/yafl/releases),
* download and unpack the latest release source code archive,
* add the folowing files to you C/C++ project:
  * [src/yafl_math.c](./src/yafl_math.c)
  * [src/yafl_math.h](./src/yafl_math.h)
  * [src/yafl.c](./src/yafl.c)
  * [src/yafl.h](./src/yafl.h)
* write yafl_config.h file and add it to you project. For Cortex-M4F or similar the file may look ike this:

```C
/*yafl_config.h*/

/*yafl_config.h*/

#ifndef YAFL_CONFIG_H
#define YAFL_CONFIG_H

#include <math.h>
#include <stdint.h>

#ifdef DEBUG
    /*
    In this example we will use standard output.
    You can actually use any printf implementation you want.
    */
#   include <stdio.h>
#   define YAFL_DBG(...) fprintf(stderr, __VA_ARGS__)

    /*
    Using branch speculation may save some clocks...
    */
#   ifdef __GNUC__
#       define YAFL_UNLIKELY(x) __builtin_expect((x), 0)
#   else /*__GNUC__*/
#       define YAFL_UNLIKELY(x) (x)
#   endif/*__GNUC__*/
#else /*DEBUG*/
#   define YAFL_DBG(...) /*Do nothing here*/
    /*
    Here we have "Never" actually, but you can use some of above definitions if you want.
    */
#   define YAFL_UNLIKELY(x) (0)
#endif/*DEBUG*/

#define YAFL_EPS  (1.0e-6)
#define YAFL_SQRT sqrtf
#define YAFL_ABS  fabsf
#define YAFL_EXP  expf
#define YAFL_LOG  logf


typedef float   yaflFloat;
typedef int32_t   yaflInt;

#endif/*YAFL_CONFIG_H*/

```
* read the [C-Manual](./doc/C-Manual.md) for usage details,
* write some usefull code which use our library in you project.


## Using with Python
We also have a Python extension for prototyping purposes. Python 3.5+ with 64bit is supproted.

To use the extension you need to:
* go to [Releases](https://github.com/shkolnick-kun/yafl/releases),
* download latest yaflpy-\<latest version\>.tar.gz,
* install it:
```bash
# Cython, numpy, scipy, setuptools, wheel
# are needed at this point
pip install path_to/yaflpy-\<latest version\>.tar.gz
```
* read the [Python-Manual](./doc/Python-Manual.md) for usage details.
* import the extension:
```Python
import yaflpy
```
* write some code which use the extension.

## References
<a name="west_1981"> **\[West1981\]** M. West, "Robust Sequential Approximate Bayesian Estimation", J. R. Statist. Soc. B (1981), 43, No. 2, pp. 157-166 </a>
