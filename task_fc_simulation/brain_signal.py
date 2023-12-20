import numpy as np
import xarray as xr
from scipy import signal


class BrainSignal:
    """
    special class for processing brain signals from different sources
    """
    name = ""
    type = "neuronal activity"

    def __init__(self, data, time=None, dt=5, time_in_ms=True):
        """

        Args:
            data (np.ndarray): data array with shape N regions to
            dt (float): sampling rate
            time_in_ms: if time in ms
        """
        data_vars = {"neural_activity": (['region', 'time'], data)}
        if time_in_ms:
            units = 'm/s'
            self.nyquist = 1000 / dt / 2
        else:
            units = 's'
            self.nyquist = 1 / dt / 2
        if time is None:
            time = np.arange(data.shape[1]) * dt

        coords = {'time': (['time'], time,
                           {'units': units,
                            'sampling_rate': dt})
                  }
        self.ds = xr.Dataset(data_vars=data_vars, coords=coords)

    @classmethod
    def read_from_netcdf(cls, path_to_data):
        ds = xr.load_dataset(path_to_data)
        assert 'neural_activity' in list(ds.keys()), "Variable with the name 'neural_activity' should be in list"
        data = ds['neural_activity'].to_numpy()
        assert 'time' in list(ds.coords), "time should be in ds coords"
        time = ds['time'].to_numpy()
        assert 'sampling_rate' in ds.time.attrs.keys(), 'sampling rate should be in time attrs'
        dt = ds['time'].sampling_rate
        if ds['time'].units == 'm/s':
            time_in_ms = True
        else:
            time_in_ms = False
        return cls(data, time=time, dt=dt, time_in_ms=time_in_ms)

    def compute_angles(self, low_f, high_f):
        high_band = high_f / self.nyquist
        low_band = low_f / self.nyquist
        b1, a1 = signal.butter(4, [high_band, low_band], btype='bandpass')
        data = self.ds.neural_activity.to_numpy()
        filtered_data = signal.filtfilt(b1, a1, data)
        analytic_data = signal.hilbert(filtered_data)
        angles = np.angle(analytic_data)
        phase_diffs = np.zeros((self.ds.dims['region'], self.ds.dims['region'], self.ds.dims['time']))
        for i in range(self.ds.dims['region']):
            for j in range(i, self.ds.dims['region']):
                phase_diffs[i, j, :] = angles[i, :] - angles[j, :]
