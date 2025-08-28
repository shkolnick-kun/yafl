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
import sys

# Write the rigth path here!
SRC_PATH = 'path_to_yafl'
EXT_PATH = os.path.join(os.path.join(SRC_PATH, 'src'), 'yaflpy')

sys.path.insert(0, EXT_PATH)

pyximport.install(pyimport=True, reload_support=True)

#Now we are ready to import yaflpy
import yaflpy
```

## YAFL python extension API
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

**NOTE: Only `predict` and `update` methods are supported by filters in YAFL at the moment.**


#### Attributes
The Kalman filter has these attrbutes:
* `P = Up.dot(Dp.dot(Up.T))` which is state covariance matrix
* `Q = Uq.dot(Dq.dot(Uq.T))` which is process noise matrix
* `R = Ur.dot(Dr.dot(Ur.T))` which is measurement noise matrix
* `rff` which is `R` forgetting factor used to ajust measurement noice covariance at runtime

Since YAFL implements UD-factorized filters we don't have these attributes directly, but have vectors to store UD-decomposition elements of these attributes:
* `Up :np.array((max(1, (dim_x * (dim_x - 1))//2), ))` is vector with upper triangular elements of P
* `Dp :np.array((dim_x, ))` is vector with diagonal elements of P
* `Uq :np.array((max(1, (dim_x * (dim_x - 1))//2), ))` is vector with upper triangular elements of Q
* `Dq :np.array((dim_x, ))` is vector with diagonal elements of Q
* `Ur :np.array((max(1, (dim_z * (dim_z - 1))//2), ))` is vector with upper triangular elements of R
* `Dr :np.array((dim_z, ))` is vector with diagonal elements of R

### Notes on Q and R adjustments
We used [this paper](https://arxiv.org/pdf/1702.00884.pdf) to implement optional Q and R adaptive adjustments.
Here are som notes on our implementation:
* All filters in this lib have the optional measurement noice covariance adjustment which can be enabled by setting `rff` attribute to a small positive number e.g. `1e-4`.
* All EKF filters in this lib have the optional process noice covariance adjustment which can be enabled by setting `qff` attribute to a small positive number e.g. `1e-4`.
* None of UKF filters have the optional process noice covariance adjustment as it leads to filters instability.

### EKF stuf
The base class for all **EKF** variants is
```Python
class yaflpy.yaflExtendedBase(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
where:
* `jfx :function(x,dt)` is function which computes the Jacobian of the `fx` function. Returns: `np.array((dim_x, dim_x))`.
* `jhx :function(x)` is function which computes the Jacobian of the `hx` function. Returns: `np.array((dim_x, dim_x))`.

This class extends `yaflpy.yaflKalmanBase`. This class is not supposet to be used directly.

**NOTE: In addition to `rff` all EKF variants have `qff` attribute which is `Q` forgetting factor used to ajust process noice covariance at runtime**

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

Both classes extend `yaflpy.yaflAdaptiveBase` by implementing `_update(self)` method with corresponding update procedures.
Both classes APIs are simmilar to basic variants.

#### Robust EKF variants

Basic class for adaptive UD-factorized **EKF** variants is
```Python
class yaflpy.yaflRobustBase(dim_x, dim_z, dt, fx, jfx, hx, jhx, gz=None, gdotz=None, **kwargs)
```

This class is based on `yaflpy.yaflExtendedBase` and adds the following hidden attrbutes to it:
* `gz :function(nu, **hx_args)` is influence function. Returns: `float`.
* `gdotz :function(nu, **hx_args)` is influence function first derivative. Returns: `float`.
Both functions take `hx_args` which are `hx` function **kwargs**

This class is not supposed to be used directly.

Adaptive UD-factorized **EKF** variants are:
```Python
class yaflpy.RobustBierman(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```
and
```Python
class yaflpy.RobustJoseph(dim_x, dim_z, dt, fx, jfx, hx, jhx, residual_z = None)
```

Both classes extend `yaflpy.yaflRobustBase` by implementing `_update(self)` method with corresponding update procedures.
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

Both classes extend `yaflpy.yaflAdaptiveRobustBase` by implementing `_update(self)` method with corresponding update procedures.
Both classes APIs are simmilar to basic variants.

### UKF stuf
#### Sigma points generator base class
The base class for sigma points generators is:
```Python
class yaflpy.yaflSigmaBase(dim_x, addf=None)
```
where:
* `addf: function(mean, delta, mult)` which computes `sigma_point = mean + delta * mult`

This class is a wrapper around `yaflUKFBaseSt` type from [src/yafl.h](../src/yafl.h).
This class is not supposed to be used directly.

#### Van der Merwe sigma points generator
The generator class is:
```Python
class MerweSigmaPoints(dim_x, alpha, beta, kappa=0.0, **kwargs)
```
This class is a wrapper around `yaflUKFMerweSt`. It extends `yaflSigmaBase` and adds the following attributes:
* `alpha` is a scaling parameter that determines the width of the sigma point spread.
* `beta`  is used to account for prior knowledge about the distribution of the mean.
* `kappa` is a secondary scaling parameter usually set to 0 or 3-dim_x.

#### Julier sigma points generator
The generator class is:
```Python
class JulierSigmaPoints(dim_x, kappa=0.0, **kwargs)
```
This class is a wrapper around `yaflUKFJulierSt`. It extends `yaflSigmaBase` and adds the following attributes:
* `kappa` is a secondary scaling parameter usually set to `0` or `3-dim_x`.

#### Filter base class
The base class for all **UKF** variants is
```Python
class yaflpy.yaflUnscentedBase(dim_x, dim_z, dt, hx, fx, points, x_mean_fn=None, z_mean_fn=None, residual_x=None, residual_z=None)
```
where:
* `points :yaflSigmaBase` is a sigma points generator object.
* `x_mean_fn :function(state_sigmas)` is function which computes state sigma points mean. Takes `np.array((points_number, dim_x))`, returns `np.array((dim_x,))`.
* `z_mean_fn :function(measurement_sigmas)` is function which computes measurements sigma points mean. Takes `np.array((points_number, dim_z))`, returns `np.array((dim_z,))`.
* `residual_x :function(a, b)`, optional, is state residual function.

This class is very simmilar to `filterpy.kalman.UnscentedKalmanFilter`, see FilterPy [docs](https://filterpy.readthedocs.io/en/latest/kalman/UnscentedKalmanFilter.html).
This class extends `yaflpy.yaflKalmanBase`. This class is not supposet to be used directly.

#### Bierman UKF variants

##### Basic Bierman **UKF** class
```Python
class yaflpy.UnscentedBierman(dim_x, dim_z, dt, hx, fx, points, x_mean_fn=None, z_mean_fn=None, residual_x=None, residual_z=None)
```

This class extends `yaflpy.yaflUnscentedBase` by implementing `_update(self)` method.

The example code is:
```Python
from yaflpy import MerweSigmaPoints as SP
from yaflpy import UnscentedBierman as KF

# Bootstrap the filter
STD = 100.

def _fx(x, dt, **fx_args):
    x = x.copy()
    x[0] += x[1] * dt
    x[2] += x[3] * dt
    return x


def _hx(x, **hx_args):
    if hx_args:
        print(hx_args)
    return np.array([x[0], x[2]])


def _zrf(a,b):
    return a - b

sp = SP(4, 0.1, 2., 0)
kf = KF(4, 2, 1., _hx, _fx, sp, residual_z=_zrf)

kf.x[0] = 0.
kf.x[1] = 0.3
kf.Dp *= .00001
kf.Dq *= 1.0e-8
kf.Dr *= STD*STD
kf.Dr[0] *= .75
kf.Ur += .5

# Generate some data
N = 6000
clean = np.zeros((N, 2))
noisy = np.zeros((N, 2))
kf_out = np.zeros((N, 2))
t     = np.zeros((N,), dtype=np.float)

for i in range(1, len(clean)):
    clean[i] = clean[i-1] + np.array([1.,1.])
    noisy[i] = clean[i]   + np.random.normal(scale=STD, size=2)
    t[i] = i

# Run the filter
for i, z in enumerate(noisy):
    kf.predict()
    kf.update(z)
    kf_out[i] = kf.zp
```

##### Adaptive Bierman **UKF** class
```Python
class yaflpy.UnscentedAdaptiveBierman(dim_x, dim_z, yaflFloat dt, hx, fx, points, **kwargs)
```
This class extends `yaflpy.yaflUnscentedBase` by addition of `chi2` attribute and implementation of `_update(self)` method.

THe usage is similar to `yaflpy.UnscentedBierman`.

##### Robust Bierman **UKF** class
```Python
class yaflpy.UnscentedRobustBierman(dim_x, dim_z, dt, hx, fx, points, gz=None, gdotz=None, **kwargs)
```

This class extends `yaflpy.yaflRobustUKFBase` by implementation of `_update(self)` method.

THe usage is similar to `yaflpy.UnscentedBierman`.

The `yaflpy.yaflRobustUKFBase` class extends `yaflpy.yaflUnscentedBase` by addition of the following hidden attrbutes:
* `gz :function(nu, **hx_args)` is influence function. Returns: `float`.
* `gdotz :function(nu, **hx_args)` is influence function first derivative. Returns: `float`.

##### Adaptive robust Bierman **UKF** class
```Python
class yaflpy.UnscentedAdaptiveRobustBierman(dim_x, dim_z, dt, hx, fx, points, **kwargs)
```

This class extends `yaflpy.yaflAdaptiveRobustUKFBase` by implementation of `_update(self)` method.

THe usage is similar to `yaflpy.UnscentedBierman`.

The `yaflpy.yaflAdaptiveRobustUKFBase` class extends `yaflpy.yaflRobustUKFBase` by addition `chi2` attribute.

#### UD-factorized "Full" UKF variants
##### Unscented
```Python
class yaflpy.Unscented(dim_x, dim_z, dt, hx, fx, points, **kwargs)
```

This class extends `yaflpy.yaflUnscentedBase` by implementation of `_update` method and by addition of the following attributes:
* `Us :np.array((max(1, (dim_x * (dim_x - 1))//2), ))` is vector with upper triangular elements of S
* `Ds :np.array((dim_x, ))` is vector with diagonal elements of S
See FilterPy [docs](https://filterpy.readthedocs.io/en/latest/kalman/UnscentedKalmanFilter.html) for details.

THe usage is similar to `yaflpy.UnscentedBierman`.

##### UnscendedAdaptive
```Python
class yaflpy.UnscentedAdaptive(dim_x, dim_z, dt, hx, fx, points, aplha = 0.000001, **kwargs)
```
where:
* `alpha: float` is probability margin which is use to initiate `chi2` attribute by `scipy.stats.chi2.ppf(1.0 - aplha, dim_z)`

This class extends `yaflpy.Unscented` by implementation of `_update` method and by addition of `chi2` attribute.

THe usage is similar to `yaflpy.UnscentedBierman`.

### Interracting multiple models estimator

The class is:
```Python
class yaflpy.IMMEstimator(filters, mu, M, dt)
```
where:
* `filters :list` is a list of `yaflpy` filters. Supported filter classes are:
  * `yaflpy.Bierman`
  * `yaflpy.Joseph`
  * `yaflpy.AdaptiveBierman`
  * `yaflpy.AdaptiveJoseph`
  * `yaflpy.RobustBierman`
  * `yaflpy.RobustJoseph`
  * `yaflpy.AdaptiveRobustBierman`
  * `yaflpy.AdaptiveRobustJoseph`
  * `yaflpy.UnscentedBierman`
  * `yaflpy.UnscentedAdaptiveBierman`
  * `yaflpy.UnscentedRobustBierman`
  * `yaflpy.UnscentedAdaptiveRobustBierman`
  * `yaflpy.Unscented`
  * `yaflpy.UnscentedAdaptive`
* `mu :np.array((len(filters),))` mode probability: `mu[i]` is the probability that filter `i` is the correct one.
* `M :np.array((len(filters), len(filters)))` is Markov chain transition matrix. `M[i,j]` is the probability of switching from filter `j` to filter `i`.
* `dt :float` is timestep between measurements is seconds.

The example is:
```Python
###############################################################################
# Import yaflpy classes

from yaflpy import Bierman as KF
from yaflpy import IMMEstimator

###############################################################################
# Define constants
_dt = 0.1
STD = 1.

###############################################################################
# Define filters
def _ca(x, dt, **fx_args):
    x = x.copy()
    x[0] += (x[1] + 0.5 * x[2] * dt) * dt
    x[1] += x[2] * dt
    return x    

def _jca(x, dt, **fx_args):
    F = np.array([
        [1., dt, dt*dt*0.5],
        [0., 1.,    dt    ],
        [0., 0.,    1.    ],
        ])
    return F

def _cv(x, dt, **fx_args):
    x = x.copy()
    x[0] += x[1] * dt
    return x    

def _jcv(x, dt, **fx_args):
    F = np.array([
        [1., dt, 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        ])
    return F

def _hx(x, **hx_args):
    return np.array([x[0]])

def _jhx(x, **hx_args):
    H = np.array([[1., 0., 0.]])
    return H

cv = KF(3, 1, _dt, _cv, _jcv, _hx, _jhx)
cv.Dp *= .1
cv.Dq *= .000001
cv.Dr = STD*STD
cv.x[0] = 0.
cv.x[1] = 0.
cv.x[2] = 0.

ca = KF(3, 1, _dt, _ca, _jca, _hx, _jhx)
ca.Dp *= .1
ca.Dq *= .000001
ca.Dr = STD*STD
ca.x[0] = 0.
ca.x[1] = 0.
ca.x[2] = 0.

###############################################################################
# Define IMM estimator
mu = np.array([0.5, 0.5])
M  = np.array([[0.95, 0.05],
               [0.05, 0.95]])

imm = IMMEstimator([ca, cv], mu, M, _dt)

###############################################################################
#Generate test data
N = 100
time  = []
clean = []
noisy = []

#CV
t = 0.
tgt = np.array([-0., 0., 0.])
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CA
tgt[2] = 0.1
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)
#CV
tgt[2] = 0.
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)
#CA
tgt[2] = -0.2
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CV
tgt[2] = 0.
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

#CA
tgt[2] = 0.1
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

tgt[2] = 0.
for i in range(N):
    time.append(t)
    clean.append(tgt[0])
    noisy.append(tgt[0] + np.random.normal(scale=STD, size=1))
    t += _dt
    tgt = _ca(tgt, _dt)

###############################################################################
#Process test data

out = []
ps  = []
mus = []

for i,z in enumerate(clean):
     imm.predict()             #dt and **fx_args are supported
     imm.update(np.array([z])) #**hx_args are supported

     out.append(imm.x.copy())
     ps.append(imm.P.copy())
     mus.append(imm.mu.copy())
```

## Good luck \%usename\%

Happy prototyping with YAFL python extension.
