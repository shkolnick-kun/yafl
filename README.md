# Yet Another Filtering library

**YAFL** means Yet Another Filtering Library.

## The library
Here you can find some **Kalman filter** variants:

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

For all **EKF** variants we have **Bierman** and **Joseph** updates.
For sequential UD-factorized **UKF** only **Bierman** updates are available.

The library is written in C and is intended for embedded systems usage:
 * We use static memory allocation
 * We use cache-friendly algorithms when available.
 * Regularization techniques used if necessary. The code is numerically stable.
 * Depends only on C standart library.

To use the libabry you need to:
 * download the source code and
 * add the folowing files to you C/C++ project:
  * [src/yafl_math.c](./src/yafl_math.c)
  * [src/yafl_math.h](./src/yafl_math.h)
  * [src/yafl.c](./src/yafl.c)
  * [src/yafl.h](./src/yafl.h)
 * write yafl_config.h file and add it to you project. For Cortex-M4F or similar the file may look ike this:

```C
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
#   define YAFL_LOG(...) fprintf(stderr, __VA_ARGS__)

    /*
    Using branch speculation may save some clocks...
    */
#   ifdef __GNUC__
#       define YAFL_UNLIKELY(x) __builtin_expect((x), 0)
#   else /*__GNUC__*/
#       define YAFL_UNLIKELY(x) (x)
#   endif/*__GNUC__*/

#else /*DEBUG*/

    /*
    Here we have "Never" actually, but you can use previous definition if you want
    */
#   define YAFL_UNLIKELY(x) (0)

#endif/*DEBUG*/

#define YAFL_EPS  (1.0e-7)

#define YAFL_SQRT sqrtf
#define YAFL_ABS  fabs

typedef float   yaflFloat;
typedef int32_t   yaflInt;

/* WARNING!!!
Fast UKF SSR updates may give dramatically incorrect results in case of adaptive Bierman filter
*/
//#define YAFL_USE_FAST_UKF

#endif // YAFL_CONFIG_H

```

## Using with Python
We also have a Python extension for prototyping purposes. Python 3.5+ with 64bit is supproted.
