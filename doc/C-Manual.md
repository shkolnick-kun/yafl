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
nx     - the dimension of state vector
nz     - the dimension of measurement vector
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
If this macro doesn't mean **never** and if some **printf-like** function oe macro is used in `YAFL_LOG` macro, then you can get something like **stacktrace** in your console og log file if error happens.

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
yaflEKFBaseSt kf = YAFL_EKF_BASE_INITIALIZER(fx, jfx, hx, jhx, zrf, nx, nz, memory);
```
where:
* `fx`     - state transition function pointer
* `jfx`    - state transition Jacobian function pointer
* `hx`     - measurement function pointer
* `jhx`    - measurement Jacobian function pointer
* `zrf`    - measurement Residual function pointer (needed to calculate the distance between forecast and measurement vectors)
* `nx`     - the dimension of state vector
* `nz`     - the dimension of measurement vector
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

### EKF stuff
The basic type for all **EKF** control blocks is `yaflEKFBaseSt`. It is initialized with `YAFL_EKF_BASE_INITIALIZER` macro.
In the example above we gave the initial explanation of its work. Now let's go deeper.

We will use numpy code to explain stuf so you should learn something about numpy.

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
yaflEKFBaseSt kf = YAFL_EKF_BASE_INITIALIZER(fx, jfx, hx, jhx, zrf, nx, nz, memory);
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
YAFL_EKF_ADAPTIVE_INITIALIZER(fx, jfx, hx, jhx, zrf, nx, nz, memory)
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
    if (3.0 >= normalized_error)
    {
        return normalized_error;
    }

    if (6.0 >= normalized_error)
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
    if (3.0 >= normalized_error)
    {
        return 1.0;
    }

    if (6.0 >= normalized_error)
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
YAFL_EKF_ROBUST_INITIALIZER(fx, jfx, hx, jhx, zrf, g, gdot,  nx, nz, memory)
```

Note that we have two additional arguments:
* `g` is influence function
* `gdot` is influence function first derivative

The example for `YAFL_EKF_ROBUST_INITIALIZER` call is:
```C
yaflEKFRobustSt kf = YAFL_EKF_ROBUST_INITIALIZER(fx, jfx, hx, jhx, zrf, psi, psi_dot,  nx, nz, memory);
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
YAFL_EKF_ADAPTIVE_ROBUST_INITIALIZER(fx, jfx, hx, jhx, zrf, g, gdot,  nx, nz, memory)
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

Work in progress...