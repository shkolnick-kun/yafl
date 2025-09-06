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
/*                                           fx  jfx  hx  jhx zrf nx  nz  rff  qff  memory */
yaflEKFBaseSt kf = YAFL_EKF_BASE_INITIALIZER(fx, jfx, hx, jhx, 0, NX, NZ, 0.0, 0.0, kf_memory);

/*
Arguments of initializer macro:
fx     - state transition function pointer
jfx    - state transition Jacobian function pointer
hx     - measurement function pointer
jhx    - measurement Jacobian function pointer
zrf    - measurement Residual function pointer (needed to calculate the distance between forecast and measurement vectors)
nx     - the dimension of state vector
nz     - the dimension of measurement vector
rff    - measurement noice covatiance forgetting factor
qff    - process noice covatiance forgetting factor
memory - the name of a memory structure.
*/
```

So we've just declared the **EKF** memory structure, initiated it, we also declared and initiated **EKF** control block.
Now we may want to filter some measurements:

```C
extern volatile int some_condition;
extern yaflStatusEn get_some_measurement(yaflFloat * measurement_vector);

yaflFloat z[NZ]; /*Memory for measurement vectors*/

yaflStatusEn status;

while (some_condition)
{
    /* This is Bierman filter predict step */
    status = YAFL_EKF_BIERMAN_PREDICT(&kf);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }

    /*Now get one measurement vector*/
    status = get_some_measurement(&z[0]);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }

    /*OK! We have a correct measurement at this point. Let's process it...*/
    status = yafl_ekf_bierman_update(&kf, &z[0]);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }
}
```
Yeah! Thats it! We've used Sequential Bierman filter to process our measurements.

## What else do we have in YAFL

Filter algorithms and data structures are declared in [src/yafl.h](./src/yafl.h).
To power the filtering code we need some math functions which are declared in [src/yafl_math](./src/yafl_math.h).

Functions declared in [src/yafl.h](./src/yafl.h) and [src/yafl_math](./src/yafl_math.h) return result code, defined in `yaflStatusEn`.
If everything is OK then `YAFL_ST_OK = 0` is returned. If something went wrong the result code is actually a bit mask.
Here is the bit field list:
* Warning bit fields:
  * `YAFL_ST_MSK_REGULARIZED`  - we had to do regularization to maintain the filter integrity
  * `YAFL_ST_MSK_GLITCH_SMALL` - a **small** error signal glitch detected: the error signal is large and its weigth is reduced. This field is used n **Robust** filters.
  * `YAFL_ST_MSK_GLITCH_LARGE` - a **large** error signal glitch detected, the error signal is too large and only its sign should be used for filter update. This field is used n **Robust** filters.
  * `YAFL_ST_MSK_ANOMALY` - an error signal anomaly detected, adaptive correction was done to handle it. This field is used in **Adaptive** filters.

* Error bit field values:
  * `YAFL_ST_INV_ARG_1` - invalid argument 1
  * `YAFL_ST_INV_ARG_2` - invalid argument 2
  * `YAFL_ST_INV_ARG_3` - invalid argument 3
  * `...`
  * `YAFL_ST_INV_ARG_11` - invalid argument 11

The execution of the YAFL code is stoped on the first error, so only one error may appear during YAFL function call.
We think that errors are actually unlikely to happen so every error check in our code is wrapped in `YAFL_UNLIKELY` macro.
If this macro doesn't mean **never** and if some **printf-like** function oe macro is used in `YAFL_DBG` macro, then you can get something like **stacktrace** in your console og log file if error happens.

On the other hand a single call to some YAFL function may yield many warning bits.
To disthiguish between warinings and errors `YAFL_ST_ERR_THR` is used, e.g.:

```C
    if (YAFL_ST_OK != status)
    {
        if (YAFL_ST_ERR_THR >= status)
        {
            /*Handle warnings*/
            (void)status;
        }
        else
        {
            /*Handle errors*/
            (void)status;
        }
    }
```

We prefer using static memory allocation so filter structure initializers look like this:
```C
/*yaflKalmanBaseSt is for internal usage only, so our code axemple is about EKF*/
yaflEKFBaseSt kf = YAFL_EKF_BASE_INITIALIZER(fx, jfx, hx, jhx, zrf, nx, nz, rff, qff, memory);
```
where:
* `fx`     - state transition function pointer
* `jfx`    - state transition Jacobian function pointer
* `hx`     - measurement function pointer
* `jhx`    - measurement Jacobian function pointer
* `zrf`    - measurement Residual function pointer (needed to calculate the distance between forecast and measurement vectors)
* `nx`     - the dimension of state vector
* `nz`     - the dimension of measurement vector
* `rff`    - measurement noice covatiance forgetting factor
* `qff`    - process noice covatiance forgetting factor
* `memory` - the name of the memory structure.

The filter control block does not store an data it has only pointers to data storage, so we must have som separate data storage called `memory`.

Actually **all** filter control block structures in YAFL have only **immutable** data so we can declare them as **constants** if needed.

Functions `fx`, `jfx`, `hx`, `jhx`, must have prototypes which must be similar to:

`yaflStatusEn fn(yaflKalmanBaseSt * self, yaflFloat *a, yaflFloat *b)`.

The function `zrf` must have the prototype which must be similar to:

`yaflStatusEn fn(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *measurement, yaflFloat *forecast)`

it is utilized to compute the error vectors as there are some cases when we can't use the Euclid distance.

As you can see all functions have `yaflKalmanBaseSt * self` as their first parameter this is done to enable passing additional data to these functions.
So all these functions are **"methods"** of some **"class"** which has `yaflKalmanBaseSt` as its **"basic class"**. As you can see we are using **OOP** approach here.

Now let's talk about memory. We use **mixins** to declare filter memory pools so users can put more fields at their filter memory pools.
This also enables us to do nested **mixin** macros calls while keeping memory pool strustures flat.

### Notes on Q and R adjustments
We used [this paper](https://arxiv.org/pdf/1702.00884.pdf) to implement optional Q and R adaptive adjustments.
Here are som notes on our implementation:
* All filters in this lib have the optional measurement noice covariance adjustment which can be enabled by setting `rff` to a small positive number e.g. `1e-4`.
* All EKF filters in this lib have the optional process noice covariance adjustment which can be enabled by setting `qff` to a small positive number e.g. `1e-4`.
* None of UKF filters have the optional process noice covariance adjustment as it leads to filters instability.

### EKF stuff
The basic type for all **EKF** control blocks is `yaflEKFBaseSt`. It is initialized with `YAFL_EKF_BASE_INITIALIZER` macro.
In the example above we gave the initial explanation of its work. Now let's go deeper.

We will use numpy code to explain stuff so you should learn something about numpy.

#### State transition function
Here is the basic definition:
```C
yaflStatusEn fx(yaflKalmanBaseSt * self, yaflFloat * new_x, yaflFloat * old_x)
{
    (void)xz;
    YAFL_CHECK(self,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NX == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(new_x,          YAFL_ST_INV_ARG_2);
    YAFL_CHECK(old_x,          YAFL_ST_INV_ARG_3);

    /*
    Here we must update state, numpy code is:

    new_x = self.fx(old_x)
    */

    return YAFL_ST_OK;
}
```
This definition applies to both **EKF** and **UKF**.

In **EKF** case `new_x == old_x` so you can just do `self.x = self.fx(self.x)`.

In **UKF** case `new_x != old_x`.

#### State transition function Jacobian
Tbe basic definition is:
```C
yaflStatusEn jfx(yaflKalmanBaseSt * self, yaflFloat * w, yaflFloat * x)
{
    yaflInt i;
    yaflInt nx;
    yaflInt nx2;

    YAFL_CHECK(self,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NX == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(w,              YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,              YAFL_ST_INV_ARG_3);

    nx  = self->Nx;
    nx2 = nx * 2;

    /* Here we must compute Jacobian of state transition function.

    w has: nx rows and 2*nx columns so we will do:

    w[:,0:nx] = self.jfx(self.x)
    */

    return YAFL_ST_OK;
}
```
This definition is **EKF** specific.

#### Measurement function
Tbe basic definition is:
```C
yaflStatusEn hx(yaflKalmanBaseSt * self, yaflFloat * y, yaflFloat * x)
{
    YAFL_CHECK(self,          YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NZ == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(y,             YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,             YAFL_ST_INV_ARG_3);

    /*
    Here we must do:

    y = self.hx(x)
    */
    return YAFL_ST_OK;
}
```
This definition applies to both **EKF** and **UKF**.

In **EKF** case so you actually do `self.y = self.hx(self.x)`.

For **UKF** `x` and `y` are state and measurement sigma points respectively.

#### Jacobian of measurement function
The basic defeinition is:
```C
yaflStatusEn jhx(yaflKalmanBaseSt * self, yaflFloat * h, yaflFloat * x)
{
    yaflInt i;
    yaflInt nx;
    yaflInt nz;

    YAFL_CHECK(self,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NX == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NZ == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(h,              YAFL_ST_INV_ARG_2);
    YAFL_CHECK(x,              YAFL_ST_INV_ARG_3);

    nx = self->Nx;
    nz = self->Nz;

    /* Here we will use a full H matrix which has nz rows and nx columns
    We must do:

    self.h = self.jhx(self.x)
    */
    return YAFL_ST_OK;
}
```
This definition is **EKF** specific.

#### Measurement residual function:
The basic definition is:

```C
yaflStatusEn zrf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *a, yaflFloat *b)
{
    yaflInt i;
    yaflInt nz;

    YAFL_CHECK(self,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NZ == self->Nz, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(result,         YAFL_ST_INV_ARG_2);
    YAFL_CHECK(a,              YAFL_ST_INV_ARG_2);
    YAFL_CHECK(b,              YAFL_ST_INV_ARG_3);

    /*
    Here we must do:

    result = self.zrf(a, b)

    for linear case:

    result = a - b
    */

    return YAFL_ST_OK;
}
```

This definition applies to both **EKF** and **UKF**.

In **EKF** case `result == b == self.y` so you actually do `self.y = self.zrf(measurement, self.y)`.

In **UKF** case `result == self.y`, `a == self.zp` and `b` is some of measurement sigma points respectively.

#### Memory
Memory pools are declared using mixins:
```C
typedef struct {
    YAFL_EKF_BASE_MEMORY_MIXIN(NX, NZ); /* This is for EKF storage. */
    /* You can declare som more storage here */
} someMemorySt;
```
Where:
* `NX` - A state vector size.
* `NZ` - A measurement vector size.

The most convenient way to initiate a filter memory is somthing like this:
```C
someMemorySt memory =
{
    /*Initial state vector*/
    .x = {
        /*Initiate x[0]...x[NX-1]*/
    },

    /*State covariance components*/
    .Up = {
        /*
        Here we have a unit upper triangular matrix.
        We don't need to store ones, so, only upper parts of three columns are stored
        */
        0,               /*1st column*/
        0,0,             /*2nd column*/
        0,0,0,           /*3rd column*/
        /*...*/
        0,0,0, /*...*/ 0 /*(NX - 1)-th column*/
    },
    /*     0    1   2       (NX - 1) elements*/
    .Dp = {DP, DP, DP, /*...*/ DP}, /*Diagonal matrix is stored in a vector*/

    /*State noise covariance components*/
    .Uq = {
        0,               /*1st column*/
        0,0,             /*2nd column*/
        0,0,0,           /*3rd column*/
        /*...*/
        0,0,0, /*...*/ 0 /*(NX - 1)-th column*/
    },
    /*     0    1   2       (NX - 1) elements*/
    .Dq = {DX, DX, DX, /*...*/ DX},

    /*Measurement noise covariance components*/
    .Ur = {
        0,               /*1st column*/
        0,0,             /*2nd column*/
        0,0,0,           /*3rd column*/
        /*...*/
        0,0,0, /*...*/ 0 /*(NZ - 1)-th column*/
    },
    /*     0    1   2       (NZ - 1) elements*/
    .Dr = {DX, DX, DX, /*...*/ DX},
};
```

#### Basic EKF variants
Basic **EKF** control block can be initialized with:
```C
yaflEKFBaseSt kf = YAFL_EKF_BASE_INITIALIZER(fx, jfx, hx, jhx, zrf, nx, nz, rff, qff, memory);
```

Predict macros are:
* `YAFL_EKF_BIERMAN_PREDICT(self)` for Bierman filter
* `YAFL_EKF_JOSEPH_PREDICT(self)` for Joseph sequential UD-factorized filter

Where `self` is a pointer to a filter.

These macros call `yaflStatusEn _yafl_ekf_predict_wrapper(yaflEKFBaseSt * self);` which is not supposed to be called directly by the user.

The example of filter predict call is:
```C
    status = YAFL_EKF_BIERMAN_PREDICT(&kf);
```

Update functions are:
* `yaflStatusEn yafl_ekf_bierman_update(yaflEKFBaseSt * self, yaflFloat *z);` for Bierman filter
* `yaflStatusEn yafl_ekf_joseph_update(yaflEKFBaseSt * self, yaflFloat *z);` for Joseph sequential UD-factorized filter

Where:
* `self` is a pointer to a filter,
* `z`    is a pointer to a measurement vector.

The example of filter predict call is:
```C
    status = yafl_ekf_bierman_update(&kf, &z[0]);
```

#### Adaptive EKF variants
In these EKF variants H-infinity filter algorithm is used to mitigate the Kalman filter divergence.
The Chi-square test is used to detect the filter divergence.

The filter control block type is:
```C
typedef struct {
    yaflEKFBaseSt base;
    yaflFloat chi2;
} yaflEKFAdaptiveSt;
```
The memory mixin used is `YAFL_EKF_BASE_MEMORY_MIXIN`.

The initilizer macro is:
```C
YAFL_EKF_ADAPTIVE_INITIALIZER(fx, jfx, hx, jhx, zrf, nx, nz, rff, qff, chi2, memory)
```

As you can see it takes the same parameter list as `YAFL_EKF_BASE_INITIALIZER`. the `chi2` field is set to `scipy.stats.chi2.ppf(0.999, 1)`

Predict macros are:
* `YAFL_EKF_ADAPTIVE_BIERAMN_PREDICT(self)` for Bierman filter
* `YAFL_EKF_ADAPTIVE_JOSEPH_PREDICT(self)` for Joseph sequential UD-factorized filter

Where `self` is a pointer to a filter.

These macros call `yaflStatusEn _yafl_ada_ekf_predict_wrapper(yaflEKFAdaptiveSt * self);` which is not supposed to be called directly by the user.

The example of filter predict call is:
```C
    status = YAFL_EKF_ADAPTIVE_BIERAMN_PREDICT(&kf);
```

Update functions are:
* `yaflStatusEn yafl_ekf_adaptive_bierman_update(yaflEKFAdaptiveSt * self, yaflFloat *z);` for Bierman filter
* `yaflStatusEn yafl_ekf_adaptive_joseph_update(yaflEKFAdaptiveSt * self, yaflFloat *z);` for Joseph sequential UD-factorized filter

Where:
* `self` is a pointer to a filter,
* `z`    is a pointer to a measurement vector.

The example of filter predict call is:
```C
    status = yafl_ekf_adaptive_bierman_update(&kf, &z[0]);
```

#### Robust EKF variants
In these EKF variants generalized linear model is used to deal with measurement glithces and non-Gaussian noise.

The influence function and it's first derivative are used to mitigate glitches.

The functions which we have used for test are:
```C

/*
The noise model is Normal with Poisson outliers.

See:

A.V. Chernodarov, G.I. Djandjgava and A.P. Rogalev (1999).
Monitoring and adaptive robust protection of the integrity of air data
inertial satellite navigation systems for maneuverable aircraft.

In: Proceedings of the RTO SCI International Conference on Integrated Navigation Systems,

held at "Electropribor", St. Petersburg, pp. 2111-10, RTO-MP-43 , NeuiIly-sur-Seine Cedex, France.
*/
yaflFloat psi(yaflKalmanBaseSt * self, yaflFloat normalized_error)
{
    (void)self;
    yaflFloat abs_ne = fabs(normalized_error);
    if (3.0 >= abs_ne)
    {
        return normalized_error;
    }

    if (6.0 >= abs_ne)
    {
        /*Small glitch*/
        return normalized_error / 3.0;
    }

    /*Large glitch*/
    return (normalized_error >= 0.0) ? 1.0 : -1.0;
}

yaflFloat psi_dot(yaflKalmanBaseSt * self, yaflFloat normalized_error)
{
    (void)self;
    yaflFloat abs_ne = fabs(normalized_error);
    if (3.0 >= abs_ne)
    {
        return 1.0;
    }

    if (6.0 >= abs_ne)
    {
        /*Small glitch*/
        return 1.0 / 3.0;
    }

    /*Large glitch*/
    return 0.0;
}

```

Where normalized error is `error[i]/Dr[i]`.

**Note that `Dr` vector now holds square roots of corresponding dispersions**

The filter control block type is:
```C
typedef struct {
    yaflEKFBaseSt base;
    yaflKalmanRobFuncP g;    /* g = -d(ln(pdf(y))) / dy */
    yaflKalmanRobFuncP gdot; /* gdot = G = d(g) / dy */
} yaflEKFRobustSt;
```
The memory mixin used is `YAFL_EKF_BASE_MEMORY_MIXIN`.

The initilizer macro is:
```C
YAFL_EKF_ROBUST_INITIALIZER(fx, jfx, hx, jhx, zrf, g, gdot, nx, nz, rff, qff, memory)
```

Note that we have two additional arguments:
* `g` is influence function
* `gdot` is influence function first derivative

The example for `YAFL_EKF_ROBUST_INITIALIZER` call is:
```C
yaflEKFRobustSt kf = YAFL_EKF_ROBUST_INITIALIZER(fx, jfx, hx, jhx, zrf, psi, psi_dot, nx, nz, rff, qff, memory);
```

Predict macros are:
* `YAFL_EKF_ROBUST_BIERAMN_PREDICT(self)` for Bierman filter
* `YAFL_EKF_ROBUST_JOSEPH_PREDICT(self)` for Joseph sequential UD-factorized filter

Where `self` is a pointer to a filter.

These macros call `yaflStatusEn _yafl_rob_ekf_predict_wrapper(yaflEKFRobustSt * self);` which is not supposed to be called directly by the user.

The example of filter predict call is:
```C
    status = YAFL_EKF_ROBUST_BIERAMN_PREDICT(&kf);
```

Update functions are:
* `yaflStatusEn yafl_ekf_robust_bierman_update(yaflEKFRobustSt * self, yaflFloat *z);` for Bierman filter
* `yaflStatusEn yafl_ekf_robust_joseph_update(yaflEKFRobustSt * self, yaflFloat *z);` for Joseph sequential UD-factorized filter

Where:
* `self` is a pointer to a filter,
* `z`    is a pointer to a measurement vector.

The example of filter predict call is:
```C
    status = yafl_ekf_robust_bierman_update(&kf, &z[0]);
```

#### Adaptive robust EKF variants
These EKF variants use both generalized linear model and adaptive correction.

The filter control block type is:
```C
typedef struct {
    yaflEKFRobustSt base;
    yaflFloat chi2;
} yaflEKFAdaptiveRobustSt;
```
The memory mixin used is `YAFL_EKF_BASE_MEMORY_MIXIN`.

The initilizer macro is:
```C
YAFL_EKF_ADAPTIVE_ROBUST_INITIALIZER(fx, jfx, hx, jhx, zrf, g, gdot, nx, nz, rff, qff, chi2, memory)
```

As you can see it takes the same parameter list as `YAFL_EKF_ROBUST_INITIALIZER`. the `chi2` field is set to `scipy.stats.chi2.ppf(0.997, 1)`

Predict macros are:
* `YAFL_EKF_ADAPTIVE_ROBUST_BIERAMN_PREDICT(self)` for Bierman filter
* `YAFL_EKF_ADAPTIVE_ROBUST_JOSEPH_PREDICT(self)` for Joseph sequential UD-factorized filter

Where `self` is a pointer to a filter.

These macros call `yaflStatusEn _yafl_ada_rob_predict_wrapper(yaflEKFAdaptiveRobustSt * self);` which is not supposed to be called directly by the user.

The example of filter predict call is:
```C
    status = YAFL_EKF_ADAPTIVE_ROBUST_BIERAMN_PREDICT(&kf);
```

Update functions are:
* `yaflStatusEn yafl_ekf_adaptive_robust_bierman_update(yaflEKFAdaptiveRobustSt * self, yaflFloat *z);` for Bierman filter
* `yaflStatusEn yafl_ekf_adaptive_robust_joseph_update(yaflEKFAdaptiveRobustSt * self, yaflFloat *z);` for Joseph sequential UD-factorized filter

Where:
* `self` is a pointer to a filter,
* `z`    is a pointer to a measurement vector.

The example of filter predict call is:
```C
    status = yafl_ekf_adaptive_robust_bierman_update(&kf, &z[0]);
```

### UKF stuff
The basic type for all **UKF** control blocks is `yaflUKFBaseSt`.

The memory mixin used is `YAFL_UKF_BASE_MEMORY_MIXIN` which is similar to other memory mixins in YAFL.

```C
YAFL_UKF_BASE_INITIALIZER(points, points_methods, fx, xmf, xrf, hx, zmf, zrf, nx, nz, rff, memory)
```
It is quite similar to **EKF** base initializer but there are no Jacobian function pointers and there are som new parameters:
* `points`         - a pointer to sigma point generator object
* `points_methods` - a pointer to sigma point generator vtable
* `xmf`            - a pointer to a state mean function
* `xrf`            - a pointer to a state residual function
* `zmf`            - a ponter to a measurement mean function

**Note that we don't have `qff` parameter since `Q` adjustment leads to numerical instability of Unscented filter variants!**

#### Siama points generators
Sigma point generator is runtime extension to `yaflUKFBaseSt` data type.
We desided to add such runtime extension to **UKF** types because user may
need to customize sigma point generation processes.

In `YAFL_UKF_BASE_INITIALIZER` call we must pass a pointer to Sigma point
generator object and a pointer to corresponding vtable.

Vtable is an object of `yaflUKFSigmaMethodsSt` type and it has two methods pointers to functions like:
* `yaflStatusEn wf(yaflUKFBaseSt * self);` which computes sigma points wegths
* `yaflStatusEn spgf(yaflUKFBaseSt * self);` which generates sigma points

As you can see the `self` parameter has `yaflUKFBaseSt` type, so these methods are actually **UKF** virtual methods.

On the other hand all sigma point generators have `yaflUKFSigmaSt` base type which is:
```C
typedef struct _yaflUKFSigmaSt {
    yaflInt         np;    /* The number of sigma points          */
    yaflUKFSigmaAddP addf; /* Sigma point addition function       */
} yaflUKFSigmaSt;
```

The actual `addf` function is something like this:
```C
yaflStatusEn my_cool_addf(yaflUKFBaseSt * self, yaflFloat * delta, yaflFloat * pivot, yaflFloat mult)
{
    YAFL_CHECK(self,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NX == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(delta,          YAFL_ST_INV_ARG_2);
    YAFL_CHECK(pivot,          YAFL_ST_INV_ARG_3);

    /*
    Here we something like this:

    delta = pivot + mult * delta;
    */
    return YAFL_ST_OK;
}
```

It is used by sigma point generator to add deltas to initial state vector to get state sigma points
in case when simple addition is not possible, e.g. in constrained **UKF** variants.

**Currently we suppport only Van der Merwe sigma points generator, more planned.**

To store sigma points we use corresponding **mixins**, e.g. `YAFL_UKF_MERWE_MEMORY_MIXIN`, see later.

#### Van der Merwe sigma points generator
Van der Merwe sigma points generator have `yaflUKFMerweSt` type.

The `yaflUKFMerweSt` object can be defined as:
```C
yaflUKFMerweSt my_sp_object = YAFL_UKF_MERWE_INITIALIZER(nx, addf, alpha, beta, kappa, memory);
```

Where:
* `nx` - is state vector size
* `alpha`, `beta`, `kappa` - are Van der Merwe sigma points parameters
* `memory` - is a name of **UKF** memory pool object

In case uf **UKF** memory pool structure declaration should look like this:
```C
typedef struct {
    YAFL_UKF_BASE_MEMORY_MIXIN(NX, NZ);
    YAFL_UKF_MERWE_MEMORY_MIXIN(NX, NZ);
    /*Other fields*/
} myUKFMemorySt;
```

Vam der Merwe sigma point generator methods are in
```C
const yaflUKFSigmaMethodsSt yafl_ukf_merwe_spm;
```

#### Julier sigma points generator
Julier sigma points generator have `yaflUKFJulierSt` type.

The `yaflUKFJulierSt` object can be defined as:
```C
yaflUKFJulierSt my_sp_object = YAFL_UKF_JULIER_INITIALIZER(nx, addf, kappa, memory);
```

Where:
* `nx` - is state vector size
* `kappa` - is Julier scaling parameter
* `memory` - is a name of **UKF** memory pool object

In case uf **UKF** memory pool structure declaration should look like this:
```C
typedef struct {
    YAFL_UKF_BASE_MEMORY_MIXIN(NX, NZ);
    YAFL_UKF_JULIER_MEMORY_MIXIN(NX, NZ);
    /*Other fields*/
} myUKFMemorySt;
```

Julier sigma points generator methods are in
```C
const yaflUKFSigmaMethodsSt yafl_ukf_julier_spm;
```

#### Mean and residual functions
Mean and residual functions are used by **Unscented transform** in cases when simple weighted mean and Euclidian distance are unusable.

State residual functions should look like this:
```C
yaflStatusEn xrf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *a, yaflFloat *b)
{
    yaflInt i;
    yaflInt nz;

    YAFL_CHECK(self,           YAFL_ST_INV_ARG_1);
    YAFL_CHECK(NX == self->Nx, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(result,         YAFL_ST_INV_ARG_2);
    YAFL_CHECK(a,              YAFL_ST_INV_ARG_3);
    YAFL_CHECK(b,              YAFL_ST_INV_ARG_4);

    /*
    Here we must do:

    result = self.xrf(a, b)

    for linear case:

    result = a - b
    */

    return YAFL_ST_OK;
}
```

State mean functions should look like this:
```C
#define _UKF_SELF ((yaflUKFBaseSt *)self) /*Or you can pass yaflUKFBaseSt * self as first argument*/
yaflStatusEn xmf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *sigmas_x)
{
    yaflInt nx;
    yaflInt np;

    YAFL_CHECK(self,                YAFL_ST_INV_ARG_1);
    nx = self->Nx;
    YAFL_CHECK(NX == mx,            YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->sp_info,  YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->sp_meth,  YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->sigmas_x, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->wm,       YAFL_ST_INV_ARG_1);
    YAFL_CHECK(result,              YAFL_ST_INV_ARG_2);
    YAFL_CHECK(sigmas_x,            YAFL_ST_INV_ARG_3);

    np = _UKF_SELF->sp_info->np; /*Get number of sigmas*/
    /*
    Here we must do:

    result = self.xmf(sigmas_x, self.wm)

    for linear case:

    result = self.wm.T.dot(sigmas_x)
    */

    return YAFL_ST_OK;
}
```

Measurement mean functions should look like this:
```C
yaflStatusEn zmf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *sigmas_z)
{
    yaflInt nz;
    yaflInt np;

    YAFL_CHECK(self,                YAFL_ST_INV_ARG_1);
    nz = self->Nz;
    YAFL_CHECK(NX == nz,            YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->sp_info,  YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->sp_meth,  YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->sigmas_z, YAFL_ST_INV_ARG_1);
    YAFL_CHECK(_UKF_SELF->wm,       YAFL_ST_INV_ARG_1);
    YAFL_CHECK(result,              YAFL_ST_INV_ARG_2);
    YAFL_CHECK(sigmas_z,            YAFL_ST_INV_ARG_3);

    np = _UKF_SELF->sp_info->np; /*Get number of sigmas*/
    /*
    Here we must do:

    result = self.zmf(sigmas_z, self.wm)

    for linear case:

    result = self.wm.T.dot(sigmas_z)
    */

    return YAFL_ST_OK;
}
```
#### Common UKF functions
Common UKF functions are:
```C
static inline yaflStatusEn yafl_ukf_post_init(yaflUKFBaseSt * self);
/* and */
static inline yaflStatusEn yafl_ukf_gen_sigmas(yaflUKFBaseSt * self);
```
These functions are applicable to all **UKF** variants.
The `yafl_ukf_post_init` must be called before any prediction or correction step.
The `yafl_ukf_gen_sigmas` may be called between prediction and correction steps to make the behaviour of filters more **"CLASSICAL"**.

#### Bierman UKF
The type of **Bierman UKF** control blocks is `yaflUKFBaseSt`.

The memory mixin used is `YAFL_UKF_BASE_MEMORY_MIXIN` which is similar to other memory mixins in YAFL.

The initializer used is `YAFL_UKF_BASE_INITIALIZER`.

Predict and update functions are simmilar to their **EKF** counterparts.
The predict function is `yafl_ukf_bierman_predict`.
The update function is `yafl_ukf_bierman_update`.

The filter bootstrap may look this:
```C

/*---------------------------------------------------------------------------*/
/*This is our filter memory structure*/
typedef struct {
    YAFL_UKF_BASE_MEMORY_MIXIN(NX, NZ);
    YAFL_UKF_MERWE_MEMORY_MIXIN(NX, NZ);
    /*Other fields*/
} myUKFMemorySt;

/*Now we will initiate our memory*/
/*Some constants*/
#define DP (0.1)
#define DX (1.0e-6)
#define DZ (400)

myUKFMemorySt ukf_memory =
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

/*Sigma points generator stuff*/
extern yaflStatusEn my_cool_addf(yaflUKFBaseSt * self, yaflFloat * delta, yaflFloat * pivot, yaflFloat mult);

/*                                             nx       addf     alpha beta kappa  memory*/
yaflUKFMerweSt sp = YAFL_UKF_MERWE_INITIALIZER(NX, my_cool_addf,  0.1,  2.,   0, ukf_memory);

/*Filter stuff*/
extern yaflStatusEn  fx(yaflKalmanBaseSt * self, yaflFloat * new_x, yaflFloat * old_x);
extern yaflStatusEn xmf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *sigmas_x);
extern yaflStatusEn xrf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *a, yaflFloat *b);

extern yaflStatusEn  hx(yaflKalmanBaseSt * self, yaflFloat * y, yaflFloat * x)
extern yaflStatusEn zmf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *sigmas_z);
extern yaflStatusEn zrf(yaflKalmanBaseSt * self, yaflFloat *result, yaflFloat *a, yaflFloat *b)

/*                                         points   points_metods    fx  xmf  xrf  hx  zmf  zrf  nx  nz  rff  memory */
yaflEKFBaseSt kf = YAFL_UKF_BASE_INITIALIZER(sp, yafl_ukf_merwe_spm, fx, xmf, xrf, hx, zmf, zrf, NX, NX, 0.0, ukf_memory);

```

The filter init sequence may look like this:
```C
yaflStatusEn status;

status = yafl_ukf_post_init(&kf);
if (YAFL_ST_OK != status)
{
    /*Handle errors here*/
    (void)status;
}
```

The filter run may look like this:
```C
extern volatile int some_condition;
extern yaflStatusEn get_some_measurement(yaflFloat * measurement_vector);

yaflFloat z[NZ]; /*Memory for measurement vectors*/

while (some_condition)
{
    /* This is Bierman filter predict step */
    status = yafl_ukf_bierman_predict(&kf);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }

    /*Now get one measurement vector*/
    status = get_some_measurement(&z[0]);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }

    /*OK! We have a correct measurement at this point. Let's process it...*/

    /*We have an option to generate sigma points manually at that point:*/
    //status = yafl_ukf_gen_sigmas(&kf);
    //if (YAFL_ST_OK != status)
    //{
    //    /*Handle errors here*/
    //    (void)status;
    //}

    /*The filter update*/
    status = yafl_ukf_bierman_update(&kf, &z[0]);
    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }
}
```

#### Adaptive Bierman UKF

The type of **Adaptive Bierman UKF** control blocks is `yaflUKFAdaptivedSt`. The basic type for `yaflUKFAdaptivedSt` is `yaflUKFBaseSt`.

The memory mixin used is `YAFL_UKF_BASE_MEMORY_MIXIN` which is similar to other memory mixins in YAFL.

The initializer used is `YAFL_UKF_ADAPTIVE_INITIALIZER` which is simmilar to `YAFL_UKF_BASE_INITIALIZER`.

Predict and update functions are simmilar to their **Adaptive EKF** counterparts.

The predict function is `yafl_ukf_adaptive_bierman_predict`.

The update function is `yafl_ukf_adaptive_bierman_update`.

The filter bootstap initialization and run are simmilar to **Bierman UKF**.

#### Robust Bierman UKF

The type of **Robust Bierman UKF** control blocks is `yaflUKFRobustSt`. The basic type for `yaflUKFRobustSt` is `yaflUKFBaseSt`.

The memory mixin used is `YAFL_UKF_BASE_MEMORY_MIXIN` which is similar to other memory mixins in YAFL.

The initializer used is:
```C
yaflEKFBaseSt kf = YAFL_UKF_ROBUST_INITIALIZER(points, points_methods, fx, xmf, xrf, hx, zmf, zrf, g, gdot, nx, nz, rff, memory);
```

where:
* `g` is influence function
* `gdot` is influence function first derivative. Yes, we are using derivatives here!

So **Robust Bierman UKF** is simmilar to **Robust Bierman EKF**.

Predict and update functions are simmilar to their **Adaptive EKF** counterparts.

The predict function is `yafl_ukf_robust_bierman_predict`.

The update function is `yafl_ukf_robust_bierman_update`.

The filter bootstap initialization and run are simmilar to **Bierman UKF** whith the following exception:

```C
extern yaflFloat psi(yaflKalmanBaseSt * self, yaflFloat normalized_error);
extern yaflFloat psi_dot(yaflKalmanBaseSt * self, yaflFloat normalized_error);

/*The filter is initiated differently from basic UKF:*/
/*                                           points   points_metods    fx  xmf  xrf  hx  zmf  zrf   g     gdot   nx  nz  rff  memory */
yaflEKFBaseSt kf = YAFL_UKF_ROBUST_INITIALIZER(sp, yafl_ukf_merwe_spm, fx, xmf, xrf, hx, zmf, zrf, psi, psi_dot, NX, NX, 0.0, ukf_memory);
```

#### Adaptive robust Bierman UKF

The type of **Adaptive robust Bierman UKF** control blocks is `yaflUKFAdaptiveRobustSt`. The basic type for `yaflUKFAdaptiveRobustSt` is `yaflUKFRobustSt`.

The memory mixin used is `YAFL_UKF_BASE_MEMORY_MIXIN` which is similar to other memory mixins in YAFL.

The initializer used is `YAFL_UKF_ADAPTIVE_ROBUST_INITIALIZER` which is simmilar to `YAFL_UKF_ROBUST_INITIALIZER`.

Predict and update functions are simmilar to their **Adaptive EKF** counterparts.

The predict function is `yafl_ukf_adaptive_robust_bierman_predict`.

The update function is `yafl_ukf_adaptive_robust_bierman_update`.

The filter bootstap initialization and run are simmilar to **Robust Bierman UKF**.

#### UD-factorized UKF

The filter does not use recusions simmilar to sequential filters.

The type of **UD-factorized UKF** control blocks is `yaflUKFSt`. The basic type for `yaflUKFSt` is `yaflUKFBaseSt`.

The memory mixin used is `YAFL_UKF_MEMORY_MIXIN` which is similar to other memory mixins in YAFL.

The initializer used is `YAFL_UKF_INITIALIZER` which is simmilar to `YAFL_UKF_BASE_INITIALIZER`.

Predict and update functions are simmilar to their **Bierman UKF** counterparts.

The predict function is `yafl_ukf_predict`.

The update function is `yafl_ukf_update`.

The filter bootstap initialization and run are simmilar to **Bierman UKF**.

#### UD-factorized adaptive UKF

The filter does not use recusions simmilar to sequential filters.

The type of **UD-factorized adaptive UKF** control blocks is `yaflUKFFullAdapiveSt`. The basic type for `yaflUKFFullAdapiveSt` is `yaflUKFSt`.

The memory mixin used is `YAFL_UKF_MEMORY_MIXIN` which is similar to other memory mixins in YAFL.

The initializer used is: `YAFL_UKF_INITIALIZER`
```C
/*                                                  points   points_metods    fx  xmf  xrf  hx  zmf  zrf  nx  nz  rff, chi2    memory */
yaflEKFBaseSt kf = YAFL_UKF_FULL_ADAPTIVE_INITIALIZER(sp, yafl_ukf_merwe_spm, fx, xmf, xrf, hx, zmf, zrf, NX, NX, rff, chi2, ukf_memory);
```
where `chi2` must be set to `scipy.stats.chi2.ppf(1.0 - some_small_number, NZ)` as we can not predict `chi2` for specific `nz` at compile time using C-preprocessor.

Predict and update functions are simmilar to their **Bierman UKF** counterparts.

The predict function is `yafl_ukf_adaptive_predict`.

The update function is `yafl_ukf_adaptive_update`.

The filter bootstap initialization and run are simmilar to **UD-factorized UKF** with exception of `YAFL_UKF_FULL_ADAPTIVE_INITIALIZER` which needs `chi2` value.

### Interracting multiple model stuff

IMM estimator used to track manuevering targets or other objects wich have multiple modes.
IMM uses separate KFs to estimate different object modes. Markov chain based algorithm is used to manage KFs in the estimator.

#### Data types

The IMM estimator control block is:
```C
typedef struct _yaflIMMCBSt {
    yaflFilterBankItemSt * bank;
    /*Markov chain parameters*/
    yaflFloat            * mu;
    yaflFloat            * M;
    /*Mixed state*/
    yaflFloat            * Up;
    yaflFloat            * Dp;
    yaflFloat            * x;
    /*Scratchpad memory*/
    yaflFloat            * cbar;
    yaflFloat            * omega;
    yaflFloat            * y;
    yaflFloat            * W;
    yaflFloat            * D;
    /*Number of filters in the bank*/
    yaflInt               Nb;
} yaflIMMCBSt;
```
where:
* `bank` is a pointer to the filter bank array of type `yaflFilterBankItemSt`;
* `mu` is a pointer to mode probability vector;
* `M` is a pointer to transition probability matrix;
* `Up` is a pointer to estimators U component of estimators covariance decomposition;
* `Dp` is a pointer to estimators D component of estimators covariance decomposition;
* `x` is a pointer to estimators state vector;
* `cbar` is a pointer scratchpad memory;
* `omega` is a pointer to mode mixing probabilities;
* `y` is a pointer scratchpad memory;
* `W` is a pointer scratchpad memory used for covariance uptdates;
* `D` is a pointer scratchpad memory used for covariance uptdates;
* `Nb` is a number of KFs in the filter bank;

THe filter bank consist of the folowing records:
```C
typedef struct _yaflFilterBankItemSt {
    yaflKalmanBaseSt   * filter;
    yaflKalmanUpdateCBP  predict;
    yaflKalmanUpdateCBP2 update;
    /*Scratchpad memory for mixed updates*/
    yaflFloat          * Us;
    yaflFloat          * Ds;
    yaflFloat          * Xs;
} yaflFilterBankItemSt;
```
where:
* `filter` is a pointer to a filter;
* `predict` is a pointer to filters predict function;
* `update` is a pointer to filters update function;
* `Us` is a pointer to scratchpad memory used to compute mixed covariance;
* `Ds` is a pointer to scratchpad memory used to compute mixed covariance;
* `Xs` is a pointer to scratchpad memory used to compute mixed state.

#### Memories
The filter bank items use filter memory pool for scratchpad memories so no extra memory needed.

The IMM estimator ise `YAFL_IMM_MEMORY_MIXIN(nb, nx)` to delacre estimators memory pool. The parameters are:
* `nb` is a number of filters in the bank.
* `nx` is a lenght of estimator and filter state vector.

#### Initializers
Filter bank items are initialized with:
* `YAFL_IMM_EKF_ITEM_INITIALIZER(kf, mem, pre, upd)` for EKF variants.
* `YAFL_IMM_UKF_ITEM_INITIALIZER(kf, mem, pre, upd)` for UKF variants.

The parameters are:
* `kf`is a pointer to the filter;
* `mem` is filters memory pool variable name;
* `pre` is a filters predict function name (pointer);
* `upd` is a filters update function name (pointer).

Different filter types can be used in one bank.

The IMM estimator is initialized with `YAFL_IMM_INITIALIZER(_bank, _nb, mem)`.

The parameters are:
* `_bank` is a filter bank array name (pointer);
* `_nb` is a filter bank size;
* `mem` is IMM estimator memory pool variable name.

#### Methods
IMM estimator have three methods:
* `yaflStatusEn yafl_imm_post_init(yaflIMMCBSt * self);` must be called to check filters integrity before use;
* `yaflStatusEn yafl_imm_predict(yaflIMMCBSt * self);` is ***predict*** function;
* `yaflStatusEn yafl_imm_update(yaflIMMCBSt * self, yaflFloat * z);` is ***update*** function.

#### Example:
Now let's put it all together:
```C
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <yafl.h>

/*=============================================================================
                            Kalman filter stuff
=============================================================================*/
#define NX 3
#define NZ 1
#define DT (0.1)

/*-----------------------------------------------------------------------------
                          Constant velocity model
-----------------------------------------------------------------------------*/
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

/*-----------------------------------------------------------------------------
                         Constant acceleration model
-----------------------------------------------------------------------------*/
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

/*-----------------------------------------------------------------------------
                              Measurement model
-----------------------------------------------------------------------------*/
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

/*-----------------------------------------------------------------------------
                           Constant velocity filter
-----------------------------------------------------------------------------*/
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

/*-----------------------------------------------------------------------------
                         Constant acceleration filter
-----------------------------------------------------------------------------*/
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

/*=============================================================================
                             IMM estimator stuff
=============================================================================*/
typedef struct {
    YAFL_IMM_MEMORY_MIXIN(2,3);
}immMemorySt;

//Yes we can use different filter types in one bank!
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

/*=============================================================================
                                 IMM usage
=============================================================================*/
extern volatile int some_condition;
extern yaflStatusEn get_some_measurement(yaflFloat * measurement_vector);

yaflFloat z[NZ]; /*Memory for measurement vectors*/

int main (void)
{
    yaflStatusEn status = yafl_imm_post_init(&imm);

    if (YAFL_ST_OK != status)
    {
        /*Handle errors here*/
        (void)status;
    }

    while (some_condition)
    {
        /* This is Bierman filter predict step */
        status = yafl_imm_predict(&imm);
        if (YAFL_ST_OK != status)
        {
            /*Handle errors here*/
            (void)status;
        }

        /*Now get one measurement vector*/
        status = get_some_measurement(&z[0]);
        if (YAFL_ST_OK != status)
        {
            /*Handle errors here*/
            (void)status;
        }

        /*OK! We have a correct measurement at this point. Let's process it...*/

        /*The filter update*/
        status = yafl_imm_update(&imm, &z[0]);
        if (YAFL_ST_OK != status)
        {
            /*Handle errors here*/
            (void)status;
        }
    }

    return 0;
}

```

## Good luck \%username\%
And we hope that you run YAFL on your MCUs and get some usefull results.
