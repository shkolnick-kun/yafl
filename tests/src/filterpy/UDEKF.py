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
from copy import deepcopy
from math import log, exp, sqrt
import sys
import numpy as np
from numpy import dot, outer, eye, zeros, ones, diag
import scipy.linalg as linalg
from filterpy.stats import logpdf
from filterpy.common import pretty_str, reshape_z
'''
UDU decomposition:
P = U * diag(D) * U^T
'''
def udu(p):

    if 2 != len(p.shape):
        return None

    if p.shape[0] != p.shape[1]:
        return None

    n = p.shape[0]

    u = zeros((n, n))
    d = zeros((n))

    d[n-1] = p[n-1,n-1]
    u[:,n-1]   = p[:,n-1] / d[n-1]

    for j in range(n-2, -1, -1):
        dd = d[j+1:]
        c = dd * u[j,j+1:] #dd is meant to be diag(d[j+1:])
        d[j] = p[j,j] - np.dot(u[j,j+1:].T, c)

        if d[j] == 0:
            return None

        for i in range(j, -1, -1):
            c = dd * u[j,j+1:]
            u[i,j] = (p[i,j] - np.dot(u[i,j+1:].T, c))/d[j]

    return u, d
'''
MWGS update:
U * diag(D) * U^T = w * diag(d) * w^T

Params:
w - is n*k float full rank
d - is k*None float
where k>n
return:
u - is n*n float upper triangular
D - id n*None float
'''
def mwgs(w,d):

    if 1 != len(d.shape):
        return None

    if 2 != len(w.shape):
        return None

    if w.shape[1] != d.shape[0]:
        return None

    if w.shape[0] >= d.shape[0]:
        return None

    n = w.shape[0]

    u = np.eye(n)
    D = np.zeros((n))

    for i in range(n-1, -1, -1):
        c = w[i,:] * d
        D[i] = np.dot(w[i,:], c)

        if D[i] <= 1e-15:
            # How about partial reset heuristics here?
            D[i] = 1e-15
            for j in range(0, i):
                u[j,i] = 0
            continue

        dd = c/D[i]

        for j in range(0, i):
            u[j,i]  = np.dot(dd, w[j,:])
            w[j,:] -= u[j,i] * w[i,:]

    return u, D


class UDExtendedKalmanFilter(object):

    """ Implements an UD modification of extended Kalman filter (EKF).
    You are responsible for setting the various state variables to
    reasonable values; the defaults will  not give you a functional filter.

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

    def __init__(self, dim_x, dim_z, dim_u=0):

        self.dim_x = dim_x
        self.dim_z = dim_z
        self.dim_u = dim_u

        self.x = zeros((dim_x, 1)) # state
        # uncertainty covariance
        self.U = eye(dim_x)
        self.D = ones((dim_x))
        self.B = 0                 # control transition matrix
        self.F = np.eye(dim_x)     # state transition matrix
        # state uncertainty
        self.Dm = eye(dim_z)       #Decorrelation matrix
        self.Ur = eye(dim_z)       #Decorrelation matrix
        self.Dr = ones((dim_z))
        # process uncertainty
        self.Uq = eye(dim_x)
        self.Dq = ones((dim_x))

        z = np.array([None]*self.dim_z)
        self.z = reshape_z(z, self.dim_z, self.x.ndim)

        # residual is computed during the innovation step. We
        # save them so that in case you want to inspect it for various
        # purposes
        self.y = zeros((dim_z, 1)) # residual
        self.S = np.zeros((dim_z, dim_z))   # system uncertainty
        self.SI = np.zeros((dim_z, dim_z))  # inverse system uncertainty

        self._log_likelihood = log(sys.float_info.min)
        self._likelihood = sys.float_info.min
        self._mahalanobis = None

        # these will always be a copy of x,P after predict() is called
        self.x_prior = self.x.copy()
        self.U_prior = self.U.copy()
        self.D_prior = self.D.copy()

        # these will always be a copy of x,P after update() is called
        self.x_post = self.x.copy()
        self.U_post = self.U.copy()
        self.D_post = self.D.copy()

    @property
    def Q(self):
        """ Process uncertainty"""
        return dot(self.Uq, dot(diag(self.Dq), self.Uq.T))

    @Q.setter
    def Q(self, value):
        """ Process uncertainty"""
        self.Uq, self.Dq = udu(value)

    @property
    def P(self):
        """ covariance matrix"""
        return dot(self.U, dot(diag(self.D), self.U.T))

    @property
    def P_prior(self):
        """ covariance matrix of the prior"""
        return dot(self.U_prior, dot(diag(self.D_prior), self.U_prior.T))

    @property
    def P_post(self):
        """ covariance matrix of the posterior"""
        return dot(self.U_post, dot(diag(self.D_post), self.U_post.T))

    @P.setter
    def P(self, value):
        """ covariance matrix"""
        self.U,self.D = udu(value);

    @property
    def R(self):
        """ measurement uncertainty"""
        return dot(self.Ur, dot(diag(self.Dr), self.Ur.T))

    @R.setter
    def R(self, value):
        """ measurement uncertainty"""
        self.Ur, self.Dr = udu(value)
        self.Dm = linalg.inv(self.Ur)

    def predict_x(self, u=0):
        """
        Predicts the next state of X. If you need to
        compute the next state yourself, override this function. You would
        need to do this, for example, if the usual Taylor expansion to
        generate F is not providing accurate results for you.
        """
        self.x = dot(self.F, self.x) + dot(self.B, u)

    def predict(self, u=0):
        """
        Predict next state (prior) using the Kalman filter state propagation
        equations.

        Parameters
        ----------

        u : np.array
            Optional control vector. If non-zero, it is multiplied by B
            to create the control input into the system.
        """

        self.predict_x(u)

        W = np.concatenate((self.Uq, dot(self.F, self.U)), axis=1)
        D = np.concatenate((self.Dq, self.D))
        self.U, self.D = mwgs(W, D)

        # save prior
        self.x_prior = np.copy(self.x)
        self.U_prior = np.copy(self.U)
        self.D_prior = np.copy(self.D)

    def _scalar_update(self, nu, h, r):
        """Joseph scalar update

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

        a = r + f.dot(v)
        K = u.dot(v / a).reshape((n, 1))

        WW = np.concatenate((outer(K, f) - u, K), axis = 1)
        DD = np.concatenate((d, np.array([r])))

        self.U, self.D = mwgs(WW, DD)
        self.x += (K*nu).reshape(self.x.shape)

    def update(self, z, HJacobian, Hx, R=None, args=(), hx_args=(),
               residual=np.subtract):
        """ Performs the update innovation of the extended Kalman filter.

        Parameters
        ----------

        z : np.array
            measurement for this step.
            If `None`, posterior is not computed

        HJacobian : function
           function which computes the Jacobian of the H matrix (measurement
           function). Takes state variable (self.x) as input, returns H.

        Hx : function
            function which takes as input the state variable (self.x) along
            with the optional arguments in hx_args, and returns the measurement
            that would correspond to that state.

        R : np.array, scalar, or None
            Optionally provide R to override the measurement noise for this
            one call, otherwise  self.R will be used.

        args : tuple, optional, default (,)
            arguments to be passed into HJacobian after the required state
            variable. for robot localization you might need to pass in
            information about the map and time of day, so you might have
            `args=(map_data, time)`, where the signature of HCacobian will
            be `def HJacobian(x, map, t)`

        hx_args : tuple, optional, default (,)
            arguments to be passed into Hx function after the required state
            variable.

        residual : function (z, z2), optional
            Optional function that computes the residual (difference) between
            the two measurement vectors. If you do not provide this, then the
            built in minus operator will be used. You will normally want to use
            the built in unless your residual computation is nonlinear (for
            example, if they are angles)
        """
        if z is None:
            self.z = np.array([[None]*self.dim_z]).T
            self.x_post = self.x.copy()
            self.P_post = self.P.copy()
            return

        if not isinstance(args, tuple):
            args = (args,)

        if not isinstance(hx_args, tuple):
            hx_args = (hx_args,)

        if R is None:
            Dm = self.Dm
            Dr = self.Dr
        elif np.isscalar(R):
            Dm = eye(self.dim_z)       #Decorrelation matrix
            Dr = ones((self.dim_z)) * R
        else:
            u,d = udu(R);
            Dm = linalg.inv(u)
            Dr = d

        if np.isscalar(z) and self.dim_z == 1:
            z = np.asarray([z], float)

        #The ExtendedKalmanFilter class has self.y, self.S, and self.SI
        #so we have to update them for [partial] compatibility.
        #And yes, this is completely ineffective!!!
        H = HJacobian(self.x, *args)
        hx = Hx(self.x, *hx_args)
        self.y = residual(z, hx)

        self.S = dot(H, dot(self.P, H.T)) + self.R
        self.SI = linalg.inv(self.S)

        #Scalar updates
        dy = dot(Dm, self.y)
        dh = dot(Dm, H)
        for j in range(self.dim_z):
            self._scalar_update(dy[j], dh[j], Dr[j])

        # set to None to force recompute
        self._log_likelihood = None
        self._likelihood = None
        self._mahalanobis = None

        # save measurement and posterior state
        self.z = deepcopy(z)
        self.x_post = self.x.copy()
        self.U_post = self.U.copy()
        self.D_post = self.D.copy()

    def predict_update(self, z, HJacobian, Hx, args=(), hx_args=(), u=0):
        """ Performs the predict/update innovation of the extended Kalman
        filter.

        Parameters
        ----------

        z : np.array
            measurement for this step.
            If `None`, only predict step is perfomed.

        HJacobian : function
           function which computes the Jacobian of the H matrix (measurement
           function). Takes state variable (self.x) as input, along with the
           optional arguments in args, and returns H.

        Hx : function
            function which takes as input the state variable (self.x) along
            with the optional arguments in hx_args, and returns the measurement
            that would correspond to that state.

        args : tuple, optional, default (,)
            arguments to be passed into HJacobian after the required state
            variable.

        hx_args : tuple, optional, default (,)
            arguments to be passed into Hx after the required state
            variable.

        u : np.array or scalar
            optional control vector input to the filter.
        """
        self.predict(u)
        self.update(z, HJacobian, Hx, self.R, args, hx_args, residual=np.subtract)

    @property
    def log_likelihood(self):
        """
        log-likelihood of the last measurement.
        """

        if self._log_likelihood is None:
            self._log_likelihood = logpdf(x=self.y, cov=self.S)
        return self._log_likelihood

    @property
    def likelihood(self):
        """
        Computed from the log-likelihood. The log-likelihood can be very
        small,  meaning a large negative value such as -28000. Taking the
        exp() of that results in 0.0, which can break typical algorithms
        which multiply by this value, so by default we always return a
        number >= sys.float_info.min.
        """
        if self._likelihood is None:
            self._likelihood = exp(self.log_likelihood)
            if self._likelihood == 0:
                self._likelihood = sys.float_info.min
        return self._likelihood

    @property
    def mahalanobis(self):
        """
        Mahalanobis distance of innovation. E.g. 3 means measurement
        was 3 standard deviations away from the predicted value.

        Returns
        -------
        mahalanobis : float
        """
        if self._mahalanobis is None:
            self._mahalanobis = sqrt(float(dot(dot(self.y.T, self.SI), self.y)))
        return self._mahalanobis

    def __repr__(self):
        return '\n'.join([
            'KalmanFilter object',
            pretty_str('x', self.x),
            pretty_str('P', self.P),
            pretty_str('x_prior', self.x_prior),
            pretty_str('P_prior', self.P_prior),
            pretty_str('F', self.F),
            pretty_str('Q', self.Q),
            pretty_str('R', self.R),
            pretty_str('K', self.K),
            pretty_str('y', self.y),
            pretty_str('S', self.S),
            pretty_str('likelihood', self.likelihood),
            pretty_str('log-likelihood', self.log_likelihood),
            pretty_str('mahalanobis', self.mahalanobis)
            ])
