import numpy as np
import numba


def simulateBOLD(Z, dt, voxelCounts, X=None, F=None, Q=None, V=None,
                 rho=None, alpha=None, V0=None, k1_mul=None, k2=None, k3_mul=None, gamma=None, k=None, tau=None):
    """ Adopted function from neurolib, added parameters for the shape of bold activation to the argument,
    the only difference
    https://github.com/neurolib-dev/neurolib/blob/master/neurolib/models/bold/timeIntegration.py

    Simulate BOLD activity using the Balloon-Windkessel model.
    See Friston 2000, Friston 2003 and Deco 2013 for reference on how the BOLD signal is simulated.
    The returned BOLD signal should be downsampled to be comparable to a recorded fMRI signal.

    :param Z: Synaptic activity
    :type Z: numpy.ndarray
    :param dt: dt of input activity in s
    :type dt: float
    :param voxelCounts: Number of voxels in each region (not used yet!)
    :type voxelCounts: numpy.ndarray
    :param X: Initial values of Vasodilatory signal, defaults to None
    :type X: numpy.ndarray, optional
    :param F: Initial values of Blood flow, defaults to None
    :type F: numpy.ndarray, optional
    :param Q: Initial values of Deoxyhemoglobin, defaults to None
    :type Q: numpy.ndarray, optional
    :param V: Initial values of Blood volume, defaults to None
    :type V: numpy.ndarray, optional

    :return: BOLD, X, F, Q, V
    :rtype: (numpy.ndarray,)
    """

    N = np.shape(Z)[0]

    if "voxelCounts" not in globals():
        voxelCounts = np.ones((N,))

    # Balloon-Windkessel model parameters (from Friston 2003):
    # Friston paper: Nonlinear responses in fMRI: The balloon model, Volterra kernels, and other hemodynamics
    # Note: the distribution of each Balloon-Windkessel models parameters are given per voxel
    # Since we usually average the empirical fMRI of each voxel for a given area, the standard
    # deviation of the gaussian distribution should be divided by the number of voxels in each area
    # voxelCountsSqrtInv = 1 / np.sqrt(voxelCounts)
    #
    # See Friston 2003, Table 1 mean values and variances:
    # rho     = np.random.normal(0.34, np.sqrt(0.0024) / np.sqrt( sum(voxelCounts) ) )    # Capillary resting net oxygen extraction
    # alpha   = np.random.normal(0.32, np.sqrt(0.0015) / np.sqrt( sum(voxelCounts) ) )    # Grubb's vessel stiffness exponent
    # V0      = 0.02
    # k1      = 7 * rho
    # k2      = 2.0
    # k3      = 2 * rho - 0.2
    # Gamma   = np.random.normal(0.41 * np.ones(N), np.sqrt(0.002) * voxelCountsSqrtInv)   # Rate constant for autoregulatory feedback by blood flow
    # K       = np.random.normal(0.65 * np.ones(N), np.sqrt(0.015) * voxelCountsSqrtInv)   # Vasodilatory signal decay
    # Tau     = np.random.normal(0.98 * np.ones(N), np.sqrt(0.0568) * voxelCountsSqrtInv)   # Transit time
    #
    # If no voxel counts are given, we can use scalar values for each region's parameter:

    # Capillary resting net oxygen extraction (dimensionless), E_0 in Friston2000
    rho = 0.34 if rho is None else rho
    voxelCountsSqrtInv = 1 / np.sqrt(voxelCounts)
    Rho = np.random.normal(rho, np.sqrt(0.0024)*voxelCountsSqrtInv, size=(N,))
    # Grubb's vessel stiffness exponent (dimensionless), \alpha in Friston2000
    alpha = 0.32 if alpha is None else alpha
    Alpha = np.random.normal(alpha, np.sqrt(0.0015)*voxelCountsSqrtInv, size=(N,))

    # Resting blood volume fraction (dimensionless)
    V0 = 0.02 if V0 is None else V0

    K1 = 7 * Rho if k1_mul is None else k1_mul*Rho # (dimensionless)
    k2 = 2.0 if k2 is None else k2  # (dimensionless)
    K3 = 2 * Rho - 0.2 if k3_mul is None else k3_mul * Rho - 0.2  # (dimensionless)
    gamma = 0.41 if gamma is None else gamma
    #Gamma = gamma * np.ones((N,))  # Rate constant for autoregulatory feedback by blood flow (1/s)
    Gamma = np.random.normal(gamma, np.sqrt(0.002)* voxelCountsSqrtInv, size=(N,))
    k = 0.65 if k is None else k
    #K = k * np.ones((N,))  # Vasodilatory signal decay (1/s)
    K = np.random.normal(k * np.ones(N), np.sqrt(0.015) * voxelCountsSqrtInv, size=(N,))
    #tau = 0.98 if tau is None else tau
    #Tau = tau * np.ones((N,))  # Transit time  (s)
    Tau = np.random.normal(0.98 * np.ones(N), np.sqrt(0.0568) * voxelCountsSqrtInv, size=(N,))
    # initialize state variables
    # NOTE: We need to use np.copy() because these variables
    # will be overwritten later and numba doesn't like to do that
    # with anything that was defined outside the scope of the @njit'ed function
    X = np.zeros((N,)) if X is None else np.copy(X)  # Vasso dilatory signal
    F = np.zeros((N,)) if F is None else np.copy(F)  # Blood flow
    Q = np.zeros((N,)) if Q is None else np.copy(Q)  # Deoxyhemoglobin
    V = np.zeros((N,)) if V is None else np.copy(V)  # Blood volume

    BOLD = np.zeros(np.shape(Z))
    # return integrateBOLD_numba(BOLD, X, Q, F, V, Z, dt, N, rho, alpha, V0, k1, k2, k3, Gamma, K, Tau)
    BOLD, X, F, Q, V = integrateBOLD_numba(BOLD, X, Q, F, V, Z, dt, N, Rho, Alpha, V0, K1, k2, K3, Gamma, K, Tau)
    return BOLD, X, F, Q, V


@numba.njit
def integrateBOLD_numba(BOLD, X, Q, F, V, Z, dt, N, Rho, Alpha, V0, K1, k2, K3, Gamma, K, Tau):
    """Integrate the Balloon-Windkessel model.

    Reference:

    Friston et al. (2000), Nonlinear responses in fMRI: The balloon model, Volterra kernels, and other hemodynamics.
    Friston et al. (2003), Dynamic causal modeling

    Variable names in Friston2000:
    X = x1, Q = x4, V = x3, F = x2

    Friston2003: see Equation (3)

    NOTE: A very small constant EPS is added to F to avoid F become too small / negative
    and cause a floating point error in EQ. Q due to the exponent **(1 / F[j])

    """

    EPS = 1e-120  # epsilon for softening
    #X[0] = 0
    #F[0] = 1
    #V[0] = 1
    #Q[0] = 1

    for i in range(len(Z[0, :])):  # loop over all timesteps
        # component-wise loop for compatibilty with numba
        for j in range(N):  # loop over all areas
            X[j] = X[j] + dt * (Z[j, i] - K[j] * X[j] - Gamma[j] * (F[j] - 1))
            Q[j] = Q[j] + dt / Tau[j] * (
                    F[j] / Rho[j] * (1 - (1 - Rho[j]) ** (1 / F[j])) - Q[j] * V[j] ** (1 / Alpha[j] - 1))
            V[j] = V[j] + dt / Tau[j] * (F[j] - V[j] ** (1 / Alpha[j]))
            F[j] = F[j] + dt * X[j]

            F[j] = max(F[j], EPS)

            BOLD[j, i] = V0 * (K1[j] * (1 - Q[j]) + k2 * (1 - Q[j] / V[j]) + K3[j] * (1 - V[j]))
    return BOLD, X, F, Q, V
