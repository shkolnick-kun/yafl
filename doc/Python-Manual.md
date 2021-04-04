# YAFL Python manual

YAFL python extension is provided for prototyping and test purposes, to install it you need to:
* go to [Releases](https://github.com/shkolnick-kun/yafl/releases),
* download latest yaflpy-\<latest version\>.tar.gz,
* install it:
```bash
# Cython, numpy, scipy, setuptools, wheel
# are needed at this point
pip install path_to/yaflpy-\<latest version\>.tar.gz
```
* import the extension:
```Python
import yaflpy
```

If you don't want to install the extension then you can do this:
```Python
import numpy as np
import pyximport
import os
import scipy.stats
import sys

# Write the rigth path here!
SRC_PATH = 'path_to_yafl/src'
EXT_PATH = os.path.join(SRC_PATH, 'yaflpy')

sys.path.insert(0, EXT_PATH)

pyximport.install(
    pyimport=True,
    reload_support=True,
    language_level=3,
    setup_args={
        'include_dirs': [np.get_include(), SRC_PATH, EXT_PATH],
        }
    )

#Now we are ready to import yaflpy
import yaflpy
```

## YAFLPY API
Actually the API is quite simmilar to [FilterPy](https://filterpy.readthedocs.io/en/latest/index.html).
We will focus on YAFL features, goto FilterPy documentation for the rest.

### Constants
We have the following constants:
```Python
#Bit masks
yaflpy.ST_MSK_REGULARIZED  = YAFL_ST_MSK_REGULARIZED
yaflpy.ST_MSK_GLITCH_SMALL = YAFL_ST_MSK_GLITCH_SMALL
yaflpy.ST_MSK_GLITCH_LARGE = YAFL_ST_MSK_GLITCH_LARGE
yaflpy.ST_MSK_ANOMALY      = YAFL_ST_MSK_ANOMALY
yaflpy.ST_MSK_ERROR        = 0xFF0
yaflpy.ST_MSK_WARNING      = 0xF

#Everthing is OK
yaflpy.ST_OK = YAFL_ST_OK

#Error threshold
yaflpy.ST_ERR_THR = YAFL_ST_ERR_THR

#Errors:
# Invalid argument numer
yaflpy.ST_INV_ARG_1  = YAFL_ST_INV_ARG_1
yaflpy.ST_INV_ARG_2  = YAFL_ST_INV_ARG_2
yaflpy.ST_INV_ARG_3  = YAFL_ST_INV_ARG_3
yaflpy.ST_INV_ARG_4  = YAFL_ST_INV_ARG_4
yaflpy.ST_INV_ARG_5  = YAFL_ST_INV_ARG_5
yaflpy.ST_INV_ARG_6  = YAFL_ST_INV_ARG_6
yaflpy.ST_INV_ARG_7  = YAFL_ST_INV_ARG_7
yaflpy.ST_INV_ARG_8  = YAFL_ST_INV_ARG_8
yaflpy.ST_INV_ARG_9  = YAFL_ST_INV_ARG_9
yaflpy.ST_INV_ARG_10 = YAFL_ST_INV_ARG_10
yaflpy.ST_INV_ARG_11 = YAFL_ST_INV_ARG_11
```
These constants are used to determine statusses of `predict` and `update` methods calls.
See [C-Manual](./C-Manual.md) for details.

### Base class
The base class is
```Python
class yaflpy.yaflKalmanBase(dim_x, dim_z, dt, fx, hx, residual_z = None)
```
where:
* `dim_x :int` is state vector dimmension,
* `dim_z: int` is measurement vector dimmension,
* `dt :float` is timestep between measurements is seconds,
* `fx :function(x,dt)` is state transistion function,
* `hx :function(x)` is measurement function,
* `residual_z :function(x, y)`, optional, is measurement residual function.

All other classes extend `yaflKalmanBase`. This base class is not supposed to be used directly.
Its API is based on [**UKF**](https://filterpy.readthedocs.io/en/latest/kalman/UnscentedKalmanFilter.html) from FilterPy.

#### Attributes
The Kalman filter has these attrbutes:
* `P = Up.dot(Dp.dot(Up.T))` which is state covariance matrix
* `Q = Uq.dot(Dq.dot(Uq.T))` which is process noise matrix
* `R = Ur.dot(Dr.dot(Ur.T))` which is measurement noise matrix

Since YAFL implements UD-factorized filters we don't have these attributes directly, but have vectors to store UD-decomposition elements of these attributes:
* `Up :np.array((max(1, (dim_x * (dim_x - 1))//2), ))` is vector with upper triangular elements of P
* `Dp :np.array((dim_x, ))` is vector with diagonal elements of P
* `Uq :np.array((max(1, (dim_x * (dim_x - 1))//2), ))` is vector with upper triangular elements of Q
* `Dq :np.array((dim_x, ))` is vector with diagonal elements of Q
* `Ur :np.array((max(1, (dim_z * (dim_z - 1))//2), ))` is vector with upper triangular elements of R
* `Dr :np.array((dim_z, ))` is vector with diagonal elements of R

### EKF stuf
The base class for all **EKF** variants is
```Python
class yaflpy.yaflExtendedBase(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
where:
* `jfx :function(x,dt)` is function which computes the Jacobian of the `fx` function. Returns: `np.array((dim_x, dim_x))`.
* `jhx :function(x)` is function which computes the Jacobian of the `hx` function. Returns: `np.array((dim_x, dim_x))`.

This class extends `yaflpy.yaflKalmanBase`. This class is not supposet to be used directly.

#### Basic EKF variants
Basic UD-factorized **EKF** variants are:
```Python
class yaflpy.Bierman(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
and
```Python
class yaflpy.Joseph(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```

Both classes extend `yaflpy.yaflExtendedBase` by implemrnting `_update(self)` method with corresponding update procedures.

The example code is:
```Python
from yaflpy import Bierman as KF
#from yaflpy import Joseph as KF

# Measurement std
STD = 100.

# Bootstrap the filter
def _fx(x, dt, **fx_args):
    x = x.copy()
    x[0] += x[1] * dt
    x[2] += x[3] * dt
    return x

def _jfx(x, dt, **fx_args):
    F = np.array([
        [1., dt, 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., 1., dt],
        [0., 0., 0., 1.],
        ])
    return F

def _hx(x, **hx_args):
    if hx_args:
        print(hx_args)
    return np.array([x[0], x[2]])

def _jhx(x, **hx_args):
    H = np.array([
        [1., 0., 0., 0.],
        [0., 0., 1., 0.],
        ])
    return H

def _zrf(a,b):
    return a - b

kf = KF(4, 2, 1., _fx, _jfx, _hx, _jhx, residual_z=_zrf)

kf.x[0] = 0.
kf.x[1] = 0.3
kf.Dp *= .00001
kf.Dq *= 1.0e-8
kf.Dr *= STD*STD
kf.Dr[0] *= .75
kf.Ur += 0.5

# Generate some data
N = 6000

clean  = np.zeros((N, 2))
noisy  = np.zeros((N, 2))
kf_out = np.zeros((N, 2))
t      = np.zeros((N,), dtype=np.float)

for i in range(1, len(clean)):
    clean[i] = clean[i-1] + np.array([1.,1.])
    noisy[i] = clean[i]   + np.random.normal(scale=STD, size=2)
    t[i] = i

#Run the filter
for i, z in enumerate(noisy):
    kf.predict()
    kf.update(z)
    kf_out[i] = kf.x[::2]
```

#### Adaptive EKF variants

Basic class for adaptive UD-factorized **EKF** variants is
```Python
class yaflpy.yaflAdaptiveBase(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
This class is based on `yaflpy.yaflExtendedBase` and adds `chi2` attribute to it.
This class is not supposed to be used directly.

Adaptive UD-factorized **EKF** variants are:
```Python
class yaflpy.AdaptiveBierman(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
and
```Python
class yaflpy.AdaptiveJoseph(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```

Both classes extend `yaflpy.yaflAdaptiveBase` by implemrnting `_update(self)` method with corresponding update procedures.
Both classes APIs are simmilar to basic variants.

#### Robust EKF variants

Basic class for adaptive UD-factorized **EKF** variants is
```Python
class yaflpy.yaflRobustBase(dim_x, dim_z, dt, fx, jfx, hx, jhx, gz=None, gdotz=None, **kwargs)
```

This class is based on `yaflpy.yaflExtendedBase` and adds the following attrbutes to it:
* `gz :function(nu, **hx_args)` is influence function. Returns: `float`.
* `gdotz :function(nu, **hx_args)` is influence function first derivative. Returns: `float`.
Both functions take `hx_args` which is `hx` function **kwargs**

This class is not supposed to be used directly.

Adaptive UD-factorized **EKF** variants are:
```Python
class yaflpy.RobustBierman(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
and
```Python
class yaflpy.RobustJoseph(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```

Both classes extend `yaflpy.yaflRobustBase` by implemrnting `_update(self)` method with corresponding update procedures.
Both classes APIs are simmilar to basic variants.

#### Adaptive robust EKF variants

Basic class for adaptive UD-factorized **EKF** variants is
```Python
class yaflpy.yaflAdaptiveRobustBase(dim_x, dim_z, dt, fx, jfx, hx, jhx, gz=None, gdotz=None, **kwargs)
```

This class is based on `yaflpy.yaflRobustBase` and adds `chi2` attribute to it.
This class is not supposed to be used directly.

Adaptive UD-factorized **EKF** variants are:
```Python
class yaflpy.AdaptiveRobustBierman(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
and
```Python
class yaflpy.AdaptiveRobustJoseph(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```

Both classes extend `yaflpy.yaflAdaptiveRobustBase` by implemrnting `_update(self)` method with corresponding update procedures.
Both classes APIs are simmilar to basic variants.


Work in progress...
