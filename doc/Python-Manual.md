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

### Base class
The base class is
```Python
yaflpy.yaflKalmanBase(dim_x, dim_z, dt, fx, hx, residual_z = None)
```
where:
* `dim_x :int` is state vector dimmension,
* `dim_z: int` is measurement vector dimmension,
* `dt :float` is timestep between measurements is seconds,
* `fx :function(x,dt)` is state transistion function,
* `hx :function(x)` is measurement function,
* `residual_z :function(x, y)`, optional, is measurement residual function.

This base class is not supposed to be used directly.
Its API is based on [**UKF**](https://filterpy.readthedocs.io/en/latest/kalman/UnscentedKalmanFilter.html) form filtery.

#### Attributes
Kalman filter hase these attrbutes:
* `P = Up.dot(Dp.dot(Up.T))` which is state covariance matrix
* `Q = Uq.dot(Dq.dot(Uq.T))` which is process noise matrix
* `R = Ur.dot(Dr.dot(Ur.T))` which is measurement noise matrix

Since YAFL implements UD-factorized filters we don't have these attributes directly, but have vectors to store UD-decomposition elements of these attributes:
* `Up :np.array((max(1, (dim_x * (dim_x - 1))//2), ))` is vector with upper triangular elements of P
* `Dp :np.array((max(1, (dim_x, ))` is vector with diagonal elements of P
* `Uq :np.array((max(1, (dim_x * (dim_x - 1))//2), ))` is vector with upper triangular elements of Q
* `Dq :np.array((max(1, (dim_x, ))` is vector with diagonal elements of Q
* `Ur :np.array((max(1, (dim_z * (dim_z - 1))//2), ))` is vector with upper triangular elements of R
* `Dr :np.array((max(1, (dim_z, ))` is vector with diagonal elements of R

### EKF stuf


Work in progress...
