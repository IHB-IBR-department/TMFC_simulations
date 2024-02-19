import numpy as np
import scipy.signal

def getPowerSpectrum(activity, dt, maxfr=70, spectrum_windowsize=1.0, normalize=False):
    """Returns a power spectrum using Welch's method.

    :param activity: One-dimensional timeseries
    :type activity: np.ndarray
    :param dt: Simulation time step
    :type dt: float
    :param maxfr: Maximum frequency in Hz to cutoff from return, defaults to 70
    :type maxfr: int, optional
    :param spectrum_windowsize: Length of the window used in Welch's method (in seconds), defaults to 1.0
    :type spectrum_windowsize: float, optional
    :param normalize: Maximum power is normalized to 1 if True, defaults to False
    :type normalize: bool, optional

    :return: Frquencies and the power of each frequency
    :rtype: [np.ndarray, np.ndarray]
    """
    # convert to one-dimensional array if it is an (1xn)-D array
    if activity.shape[0] == 1 and activity.shape[1] > 1:
        activity = activity[0]
    assert len(activity.shape) == 1, "activity is not one-dimensional!"

    f, Pxx_spec = scipy.signal.welch(
        activity,
        1000 / dt,
        window="hann",
        nperseg=int(spectrum_windowsize * 1000 / dt),
        scaling="spectrum",
    )
    f = f[f < maxfr]
    Pxx_spec = Pxx_spec[0 : len(f)]
    if normalize:
        Pxx_spec /= np.max(Pxx_spec)
    return f, Pxx_spec


def getMeanPowerSpectrum(activities, dt, maxfr=70, spectrum_windowsize=1.0, normalize=False):
    """Returns the mean power spectrum of multiple timeseries.

    :param activities: N-dimensional timeseries
    :type activities: np.ndarray
    :param dt: Simulation time step
    :type dt: float
    :param maxfr: Maximum frequency in Hz to cutoff from return, defaults to 70
    :type maxfr: int, optional
    :param spectrum_windowsize: Length of the window used in Welch's method (in seconds), defaults to 1.0
    :type spectrum_windowsize: float, optional
    :param normalize: Maximum power is normalized to 1 if True, defaults to False
    :type normalize: bool, optional

    :return: Frquencies and the power of each frequency
    :rtype: [np.ndarray, np.ndarray]
    """

    powers = np.zeros(getPowerSpectrum(activities[0], dt, maxfr, spectrum_windowsize)[0].shape)
    ps = []
    for rate in activities:
        f, Pxx_spec = getPowerSpectrum(rate, dt, maxfr, spectrum_windowsize)
        ps.append(Pxx_spec)
        powers += Pxx_spec
    powers /= len(ps)
    if normalize:
        powers /= np.max(powers)
    return f, powers