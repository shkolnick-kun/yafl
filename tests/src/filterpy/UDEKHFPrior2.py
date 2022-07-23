# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,too-many-instance-attributes, too-many-arguments
"""
Copyright 2022 Paul A Beltyukov
Copyright 2019 Paul A Beltyukov
Copyright 2015 Roger R Labbe Jr.

FilterPy library.
http://github.com/rlabbe/filterpy

Documentation at:
https://filterpy.readthedocs.org

Supporting book at:
https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python

This is licensed under an MIT license. See the readme.MD file
for more information.
"""
import numpy as np
from numpy import dot, outer
from scipy.stats import chi2
from UDEKF import mwgs, UDExtendedKalmanFilter


class UDExtendedKalmanHinfFilterPrior2(UDExtendedKalmanFilter):

    """
    Implements an UD modification of extended Kalman/Hinfinity filter (EKHF)
    with prior resuduals used for Hinfinity correction. You are responsible
    for setting the various state variables to reasonable values;
    the defaults will  not give you a functional filter.

    You will have to set the following attributes after constructing this
    object for the filter to perform properly. Please note that there are
    various checks in place to ensure that you have made everything the
    'correct' size. However, it is possible to provide incorrectly sized
    arrays such that the linear algebra can not perform an operation.
    It can also fail silently - you can end up with matrices of a size that
    allows the linear algebra to work, but are the wrong shape for the problem
    you are trying to solve.

    Parameters
    ----------

    dim_x : int
        Number of state variables for the Kalman filter. For example, if
        you are tracking the position and velocity of an object in two
        dimensions, dim_x would be 4.

        This is used to set the default size of P, Q, and u

    dim_z : int
        Number of of measurement inputs. For example, if the sensor
        provides you with position in (x,y), dim_z would be 2.

    Attributes
    ----------
    x : numpy.array(dim_x, 1)
        State estimate vector

    P : numpy.array(dim_x, dim_x)
        Covariance matrix

    x_prior : numpy.array(dim_x, 1)
        Prior (predicted) state estimate. The *_prior and *_post attributes
        are for convienence; they store the  prior and posterior of the
        current epoch. Read Only.

    P_prior : numpy.array(dim_x, dim_x)
        Prior (predicted) state covariance matrix. Read Only.

    x_post : numpy.array(dim_x, 1)
        Posterior (updated) state estimate. Read Only.

    P_post : numpy.array(dim_x, dim_x)
        Posterior (updated) state covariance matrix. Read Only.

    R : numpy.array(dim_z, dim_z)
        Measurement noise matrix

    Q : numpy.array(dim_x, dim_x)
        Process noise matrix

    F : numpy.array()
        State Transition matrix

    H : numpy.array(dim_x, dim_x)
        Measurement function

    y : numpy.array
        Residual of the update step. Read only.

    K : numpy.array(dim_x, dim_z)
        Kalman gain of the update step. Read only.

    S :  numpy.array
        Systen uncertaintly projected to measurement space. Read only.

    z : ndarray
        Last measurement used in update(). Read only.

    log_likelihood : float
        log-likelihood of the last measurement. Read only.

    likelihood : float
        likelihood of last measurment. Read only.

        Computed from the log-likelihood. The log-likelihood can be very
        small,  meaning a large negative value such as -28000. Taking the
        exp() of that results in 0.0, which can break typical algorithms
        which multiply by this value, so by default we always return a
        number >= sys.float_info.min.

    mahalanobis : float
        mahalanobis distance of the innovation. E.g. 3 means measurement
        was 3 standard deviations away from the predicted value.

        Read only.
    """

    def __init__(self, dim_x, dim_z, dim_u=0, alpha = 0.001):
        UDExtendedKalmanFilter.__init__(self, dim_x, dim_z, dim_u)
        self.beta_1 = chi2.ppf(1.0 - alpha, 1)

    def _scalar_update(self, nu, h, r):
        """Joseph/Hinfinity scalar update using prior errors

        Parameters
        ----------

        axis_residual : function which returns current axis residual
                        returns scalar, float.

        axis_hjacobian : function which returns current axis HJacobian row
                         returns np.array, float.

        r : scalar, float, current axis state disp
        """
        u, d, n = self.U, self.D, self.dim_x

        f = h.dot(u)
        v = d * f
        c = f.dot(v)
        s = c + r

        a = (nu * nu) / self.beta_1 - s
        if a > 0.0:
            #Divergence detected, H-infinity correction needed
            k = a / c + 1.0
            d *= k
            v *= k
            s = k * c + r

        K = u.dot(v / s).reshape((n, 1))

        WW = np.concatenate((outer(K, f) - u, K), axis = 1)
        DD = np.concatenate((d, np.array([r])))

        self.U, self.D = mwgs(WW, DD)
        self.x += (K*nu).reshape(self.x.shape)