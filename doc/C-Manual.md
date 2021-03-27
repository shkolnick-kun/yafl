# YAFL C Manual
## First steps
Well, they are:
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
    Here we have "Never" actually, but you can use some of above definitions if you want.
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
* read this manul for usage details,
* write some usefull code which use our library in you project.

## SO WUT NOW?!!
Now you should:
```C
#include <yafl.h>

/*OK! We can do some filtering...*/

#define NX 4 /* State vector dimension       */
#define NZ 2 /* Observation vectio dimention */

/* This is our state transition function */
yaflStatusEn fx(yaflKalmanBaseSt * self, yaflFloat * x, yaflFloat * xz)
{
    (void)xz;
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(4 == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_2);

    /*We have a linear uniform motion here:*/
    x[0] += 0.1 * x[2];
    x[1] += 0.1 * x[3];
    return YAFL_ST_OK;
}

/* This is state transition Jacobian function */
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

    /* W has: nx rows and 2*nx columns */
    /* We must put our Jacobian to left side of W */
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

    /*
    Now W is

    (JFX|*)

    where JFX is:

    | 1   0   .1  0  |
    | 0   1   0   .1 |
    | 0   0   1   0  |
    | 0   0   0   1  |
    */
    return YAFL_ST_OK;
}

/*This is our measurement function*/
yaflStatusEn hx(yaflKalmanBaseSt * self, yaflFloat * y, yaflFloat * x)
{
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(2 == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(y,             YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_3);

    /*We see only coordinates, not the velocities*/
    y[0] = x[0];
    y[1] = x[1];
    return YAFL_ST_OK;
}

/* This is hx Jacobian function */
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

    /* Here we will use a full H matrix which has nz rows and nx columns */
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

    /*
    Now h is:
    | 1   0   0   0 |
    | 0   1   0   0 |

    */
    return YAFL_ST_OK;
}
/*---------------------------------------------------------------------------*/
/*This is our filter memory structure*/
typedef struct
{
    YAFL_EKF_BASE_MEMORY_MIXIN(NX, NZ); /*This mixin is actually used to declare EKF data storage*/
    yaflFloat dummy[30]; /*We may need some additional memory...*/
} kfMemorySt;

/*Now we will initiate our memory*/
/*Some constants*/
#define DP (0.1)
#define DX (1.0e-6)
#define DZ (400)

kfMemorySt kf_memory =
{
    /*Initial state vector*/
    .x = {
        [0] = 50.0,
        [2] = 10.0
        /*Other values are zeros*/
    },

    /*State covariance components*/
    .Up = {
        /*
        Here we have a unit upper triangular matrix.
        We don't need to store ones, so, only upper parts of three columns are stored
        */
        0,    /*1st column*/
        0,0,  /*2nd column*/
        0,0,0 /*3rd column*/
    },
    .Dp = {DP, DP, DP, DP}, /*Diagonal matrix is stored in a vector*/

    /*State noise covariance components*/
    .Uq = {
        0,
        0,0,
        0,0,0
    },
    .Dq = {DX, DX, DX, DX},

    /*Measurement noise covariance components*/
    .Ur = {0},
    .Dr = {DZ, DZ}
};

/*This is EKF structure definition*/
/*                                           fx  jfx  hx  jhx zrf nx  nz   memory */
yaflEKFBaseSt kf = YAFL_EKF_BASE_INITIALIZER(fx, jfx, hx, jhx, 0, NX, NZ, kf_memory);

/*
Arguments of initializer macro:
fx     - state transition function pointer
jfx    - state transition Jacobian function pointer
hx     - measurement function pointer
jhx    - measurement Jacobian function pointer
zrf    - measurement Residual function pointer (needed to calculate the distance between forecast and measurement vectors)
memory - the name of the memeory structure.
*/
```

So we've just declared the **EKF** memory structure, initiated it, we also declared and initiated **EKF** control block.
Now we may want to filter some measurements:

```C
extern yaflStatusEn get_some_measurement(yaflFloat * measurement_vector);

yaflFloat zp[NZ]; /*Memory for measurement vectors*/

yaflStatusEn status;

while (some_condition)
{
    /* This is Joseph filter predict step */
    status = YAFL_EKF_JOSEPH_PREDICT(&kf);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }

    /*Now get one measurement vector*/
    status = get_some_measurement(&zp[0]);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }

    /*OK! We have a correct measurement at this point. Let's handle it...*/
    status = yafl_ekf_joseph_update(&kf, &zp[0]);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }
}
```
Yeah! Thats it! We've jut used Joseph filter to filter our measurements.

## OK! What else do we have in [yafl.h](./src/yafl.h)?
