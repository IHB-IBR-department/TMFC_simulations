from scipy import signal, stats, io
import xarray as xr
import matplotlib.pyplot as plt
from neurolib.models.wc import WCModel
from neurolib.utils.collections import dotdict
from .boldIntegration import simulateBOLD
from .read_utils import generate_sw_matrices_from_mat
from .read_utils import read_onsets_from_mat
import numpy as np
from typing import Optional, Union
import numpy.typing as npt
from .task_utils import (create_task_design_activation,
                         create_reg_activations,
                        create_activations_per_module)


class WCTaskSim:
    """Class for simulation block design fMRI with WC model,
    where used the own matrix synaptic matrix for each state
    with predefined onset times
    """

    def __init__(self,
                 Wij_task_dict: dict[str, npt.NDArray],
                 Wij_rest: npt.NDArray, D: npt.NDArray,
                 onset_time_list: Optional[list[float]] = None,
                 task_name_list: Optional[list[str]] = None,
                 duration_list: Union[list, float] = 3,
                 rest_before: bool = True,
                 first_duration: Union[list, float] = 6,
                 last_duration: Union[list, float] = 8,
                 append_outputs: bool = False,
                 bold: bool = False,
                 chunkwise: bool = False,
                 num_modules: int = 3,
                 act_sampling_rate: float = 20,
                 exc_ext: float = 0.75,
                 K_gl: float = 2.85,
                 sigma_ou: float = 5 * 1e-3,
                 **kwargs):
        """ The most arguments chosen in accordance to neurolib inputs. The key difference
        is changing weight matrix during simulation

        Args:
            Wij_task_dict (dict):
                dictionary with the matrices for all tasks used in experiment
            Wij_rest (np.ndarray):
                matrix of synaptic weight for rest
            D (np.ndarray):
                delay matrix
            onset_time_list (list):
                list with the start time for the tasks
            task_name_list (list):
                name of the task to start
            duration_list (list):
                duration for each task or duration for all tasks if single float
            rest_before (bool):
                if need to generate rest before the start of experiment
            first_duration (float):
                duration of first rest signal, in sec
            last_duration (float):
                duration of last rest signal, in sec
            append_outputs (bool):
                if need to append raw outputs after each integration, default False,
                because it is very memory consuming
            bold (bool):
                True if use internal bold from neurolib, resampling will be with 0.5Hz
            chunkwise (bool):
                True if use internal chunkwise from neurolib, default False
            num_modules (int):
                number of modules in matrix
            act_sampling_rate (float):
                sampling rate of neuronal activity to store, in ms
                original neural activity simulated with sampling rate 0.1 ms
            normalize_max (float):
                normalizing constant for bold input transformation
            output_activation (str):
                type of output activation to pass to bold signal, possible params: 'syn_act', 'sum'
            exc_ext (float):
                excitation parameter
            K_gl (float):
                coupling parameter
            sigma_ou (float):
                sigma parameter
            **kwargs:
                other parameters for wc model
        """

        self.Dmat = D
        self.hrf = None
        self.BOLD = None
        self.normalize_max = None
        self.t_BOLD = None
        self.output_activation = None
        self.rest_before = rest_before
        self.Wij_rest = Wij_rest
        self.wc = WCModel(Cmat=self.Wij_rest, Dmat=self.Dmat)
        self.Wij_task_dict = Wij_task_dict
        self.num_tasks = len(Wij_task_dict.keys())
        if task_name_list is None:
            task_name_list = list(Wij_task_dict.keys())
        self.task_name_list = task_name_list
        assert sum([task_name in Wij_task_dict.keys()
                    for task_name in task_name_list]) == len(
            task_name_list), "Name in task names not from the dict"
        if onset_time_list is None:
            self.onset_time_list = list(np.arange(len(self.task_name_list)))
        else:
            self.onset_time_list = onset_time_list
        assert len(self.task_name_list) == len(
            self.onset_time_list), (
            "Number of tasks should be equal to number of onset times")
        if len(self.onset_time_list) > 1:
            assert (np.diff(self.onset_time_list) > 0).any(), (
                "Next onset time should be more than previous ")

        if isinstance(duration_list, (list, tuple)):
            self.duration_list = duration_list
            assert len(self.duration_list) == len(
                self.onset_time_list), (
                "Lenght of duration list should be equal to onset_time list")
        else:
            assert isinstance(duration_list, (int, float))
            self.duration_list = [duration_list] * len(self.onset_time_list)

        self.append_outputs = append_outputs
        self.chunkwise = chunkwise
        self.exc_rest = None
        # activity sampling rate in ms
        self.activity = {"exc_series": None, "inh_series": None, "sa_series": None,
                         "sampling_rate": act_sampling_rate, "idx_last_t": 0, "t": None}
        self.bold = bold
        self.bold_input_ready = False
        self.TR = None
        self.last_duration = last_duration
        self.first_duration = first_duration
        time_idxs_dict = {"Rest": []}
        time_idxs_dict.update({key: [] for key in self.Wij_task_dict.keys()})
        self.time_idxs_dict = time_idxs_dict
        self.num_modules = num_modules
        self.num_regions_per_modules = None

        kw_defaults = {'inh_ext': 0,
                       'tau_ou': 5,
                       'a_exc': 1.5,
                       'a_inh': 1.5,
                       'c_excexc': 16,
                       'c_excinh': 15,
                       'c_inhexc': 12,
                       'c_inhinh': 3,
                       'mu_exc': 3,
                       'mu_inh': 3,
                       'tau_exc': 2.5,
                       'tau_inh': 3.75,
                       'signalV': 10}
        for kw in kw_defaults.keys():
            kwargs.setdefault(kw, kw_defaults[kw])

        self.kwargs = kwargs
        self.exc_ext = exc_ext
        self.K_gl = K_gl
        self.sigma_ou = sigma_ou
        self.init_wc_model()

    def init_wc_model(self):
        self.wc.params['duration'] = self.duration_list[0] * 1000
        # this default values optimised for 30 rois
        self.wc.params['exc_ext'] = self.exc_ext
        self.wc.params['K_gl'] = self.K_gl
        self.wc.params['sigma_ou'] = self.sigma_ou

        # this default values common for all
        self.wc.params['inh_ext'] = self.kwargs['inh_ext']
        self.wc.params['tau_ou'] = self.kwargs['tau_ou']
        self.wc.params['a_exc'] = self.kwargs['a_exc']
        self.wc.params['a_inh'] = self.kwargs['a_inh']
        self.wc.params['mu_exc'] = self.kwargs['mu_exc']
        self.wc.params['mu_inh'] = self.kwargs['mu_inh']
        self.wc.params['tau_exc'] = self.kwargs['tau_exc']
        self.wc.params['tau_inh'] = self.kwargs['tau_inh']
        self.wc.params['signalV'] = self.kwargs['signalV']
        self.wc.params['c_excexc'] = self.kwargs['c_excexc']
        self.wc.params['c_excinh'] = self.kwargs['c_excinh']
        self.wc.params['c_inhinh'] = self.kwargs['c_inhinh']
        self.wc.params['c_inhexc'] = self.kwargs['c_inhexc']

    @classmethod
    def from_matlab_structure(cls, mat_path: str,
                              sigma: float = 0.1,
                              norm_type: str = 'cols',
                              num_regions: int = 30,
                              num_modules: int = 3,
                              num_regions_per_modules: Optional[list[int]] =None,
                              gen_type: str = 'simple_prod',
                              bold: bool = False,
                              chunkwise: bool = False,
                              delay: float = 250,
                              append_outputs: bool = False,
                              normalize_max: float = 50,
                              output_activation: str = 'exc',
                              rest_before: bool = True,
                              first_duration: float = 12,
                              last_duration: float = 8,
                              exc_ext=0.75,
                              K_gl=2.85,
                              sigma_ou=5 * 1e-3,
                              **kwargs):
        """Generate simulation from matlab structure, see in description in
        read_onsets_from_mat from read_utils. sigma, norm_type, gen_type, num_regions and
        num modules  - parameters from generate_synaptic_weights_matrices

        :param mat_path: path to matfile
        :param sigma: noise parameter
        :param norm_type: type of normalization
        :param num_regions:  number of regions
        :param num_modules: number of modules
        :param gen_type: [simple_prod, equal_var]
        :param bold: if use original bold
        :param chunkwise: if use original chunkwise
        :param delay:
        :param append_outputs: if save raw outputs (memory consuming)
        :param normalize_max:
        :param output_activation:
        :param rest_before:
        :param first_duration:
        :param last_duration:
        :param exc_ext:
        :param K_gl:
        :param sigma_ou:
        :param kwargs:
        :return:

        Args:
            num_regions_per_modules:
        """
        cls.num_modules = num_modules
        cls.num_regions_per_modules = num_regions_per_modules
        Wij_rest, Wij_task_dict = generate_sw_matrices_from_mat(mat_path,
                                                                num_regions,
                                                                num_modules,
                                                                num_regions_per_modules,
                                                                sigma=sigma,
                                                                gen_type=gen_type,
                                                                norm_type=norm_type)
        D = np.ones((num_regions, num_regions)) * delay

        onset_time_list, task_names_list, duration_list = read_onsets_from_mat(mat_path)

        return cls(Wij_task_dict,
                   Wij_rest,
                   D,
                   onset_time_list=onset_time_list,
                   task_name_list=task_names_list,
                   duration_list=duration_list,
                   rest_before=rest_before,
                   first_duration=first_duration,
                   last_duration=last_duration,
                   normalize_max=normalize_max,
                   output_activation=output_activation,
                   append_outputs=append_outputs,
                   bold=bold,
                   chunkwise=chunkwise,
                   num_modules=num_modules,
                   exc_ext=exc_ext,
                   K_gl=K_gl,
                   sigma_ou=sigma_ou,
                   **kwargs)

    def _generate_single_block(self,
                               Cmat: npt.NDArray,
                               duration: float = 10,
                               activity: bool = True,
                               a_s_rate: float = 0.02,
                               syn_act: bool = True):
        """
        Function for generation neural oscillations based on Cmat data.
        On top of activations two options for handling bold signal is available: built in neurolib,
        in this case attribute self.chunkwise should be set True and self.bold which will be further used in wc.run.
        In this case bold will be generated along with excitation and inhibition activity
        Args:
            Cmat (np.ndarray): matrix of synaptic weight connection
            duration (int): generation length in seconds
            activity (bool): if activity should additionally saved with resampling with a_s_rate
            a_s_rate (float): sampling rate for activity saving, in seconds
            syn_act (bool): True if need to generate synaptic activity (weighted
            sum of inhibition and excitation)

        Returns:

        """
        # buit-in nativ neurolib chunkwise integration
        if self.chunkwise:
            assert duration % 2 == 0, (
                "For faster integration time duration is chunkwise"
                " duration should be divisible by two ")
        self.wc.params['Cmat'] = Cmat
        np.fill_diagonal(self.wc.params['Cmat'], 0)
        self.wc.params['duration'] = duration * 1000  # duration in ms
        self.wc.run(append_outputs=self.append_outputs,
                    bold=self.bold,
                    continue_run=True,
                    chunkwise=self.chunkwise)
        # in-place generation of synaptic activity
        if syn_act:
            syn_act_series = self.generate_neuronal_oscill()
        # for saving downsampled neural activity
        if activity:
            idx_last_t = self.activity["idx_last_t"]
            self.activity["sampling_rate"] = a_s_rate
            sampling_rate_dt = int(round(1000 * a_s_rate / self.wc.params["dt"]))
            new_exc_activity = self.wc.exc[:,
                               sampling_rate_dt -
                               np.mod(idx_last_t - 1,
                                      sampling_rate_dt)::sampling_rate_dt
                               ]
            new_inh_activity = self.wc.inh[:,
                               sampling_rate_dt -
                               np.mod(idx_last_t - 1,
                                      sampling_rate_dt)::sampling_rate_dt
                               ]
            if syn_act:
                new_sa_activity = syn_act_series[:,
                                  sampling_rate_dt -
                                  np.mod(idx_last_t - 1,
                                         sampling_rate_dt)::sampling_rate_dt
                                  ]
            if self.activity["exc_series"] is None:
                self.activity["exc_series"] = new_exc_activity
                self.activity["inh_series"] = new_inh_activity
                if syn_act:
                    self.activity["sa_series"] = new_sa_activity
            else:
                # noinspection PyTypedDict
                self.activity["exc_series"] = np.hstack([self.activity["exc_series"],
                                                         new_exc_activity])
                # noinspection PyTypedDict
                self.activity["inh_series"] = np.hstack([self.activity["inh_series"],
                                                         new_inh_activity])
                if syn_act:
                    # noinspection PyTypedDict
                    self.activity["sa_series"] = (
                        np.hstack([self.activity["sa_series"],
                                   new_sa_activity]))

            new_idx_t = idx_last_t + np.arange(self.wc.exc.shape[1])
            t_activity = self.wc.params["dt"] * (
                new_idx_t[sampling_rate_dt -
                          np.mod(idx_last_t - 1,
                                 sampling_rate_dt)::sampling_rate_dt])
            if self.activity["t"] is None:
                self.activity["t"] = t_activity
            else:
                self.activity["t"] = np.hstack([self.activity["t"], t_activity])
            self.activity["idx_last_t"] = idx_last_t + self.wc.exc.shape[1]

    def _generate_first_rest(self,
                             activity: bool = True,
                             a_s_rate: float = 0.02,
                             syn_act: bool = True):
        """ Wrapper for first resting state block """
        # first block in design always started with resting state
        if self.rest_before:
            start_time_rest = -self.first_duration
            duration = -start_time_rest + self.onset_time_list[0]
            end_time_rest = self.onset_time_list[0]
        else:
            start_time_rest = 0
            duration = self.onset_time_list[0]
            end_time_rest = self.onset_time_list[0]
        Cmat = self.Wij_rest
        self._generate_single_block(Cmat,
                                    duration=duration,
                                    activity=activity,
                                    a_s_rate=a_s_rate,
                                    syn_act=syn_act)
        self.time_idxs_dict["Rest"].append([round(start_time_rest, 3),
                                            round(end_time_rest, 3)])

    def _generate_last_rest(self,
                            activity: bool = True,
                            a_s_rate: float = 0.02,
                            syn_act: bool = True):
        """ Wrapper for last resting state block """
        start_time_rest = self.onset_time_list[-1] + self.duration_list[-1]
        Cmat = self.Wij_rest
        # set last rest duration equal to previous gap between onset times
        duration = self.last_duration
        end_time_rest = start_time_rest + duration
        self._generate_single_block(Cmat,
                                    duration=duration,
                                    activity=activity,
                                    a_s_rate=a_s_rate, syn_act=syn_act)
        self.time_idxs_dict["Rest"].append([round(start_time_rest, 3),
                                            round(end_time_rest, 3)])

    def generate_full_series(self,
                             bold_chunkwise: bool = True,
                             TR: float = 2,
                             activity: bool = True,
                             a_s_rate: float = 0.02,
                             normalize_max: float = 2,
                             output_activation: str = 'syn_act',
                             clear_raw: bool = True,
                             fix_bold=True,
                             **kwargs):
        """
        Function for generation full length activity based on experiment design
        and matrices set up

        Args:
            bold_chunkwise (bool): in this case at each block or event bold signal
            will be computed on top of raw activity
            TR (float): time repetition for bold signal generation
            activity (bool): if resampled activity need to be saved
            a_s_rate (float): sampling rate for activity resampling, im s,
            typical values 5*1e-3 equal to 5 ms
            clear_raw (bool): True if clear raw excitation activity from wc model
            **kwargs:

        Returns:
        object with BOLD and activity parameters filled
        """

        # set up bold parameters
        self.output_activation = output_activation
        self.normalize_max = normalize_max
        self.fix_bold = fix_bold
        self.TR = TR
        if bold_chunkwise:
            self.bold_input_ready = True
            self.append_outputs = False
            if self.bold is True:
                print("Options bold_chunkwise incompatible with bold True, "
                      "set to False")
            self.bold = False
            self.chunkwise = False
            chunksize = TR * 1000 / self.wc.params["dt"]
        if self.output_activation != 'exc':
            syn_act = True
        else:
            syn_act = False
        self._generate_first_rest(activity=activity, a_s_rate=a_s_rate, syn_act=syn_act)

        if bold_chunkwise:
            assert self.wc['exc'].shape[1] >= chunksize, \
                "First rest series should be longer than TR"
            self.generate_bold_chunkwise(TR=TR,
                                         input_type=self.output_activation,
                                         normalize_max=normalize_max,
                                         is_first=True,
                                         fix=self.fix_bold,
                                         **kwargs)

        for i in range(len(self.onset_time_list)):
            task_name = self.task_name_list[i]
            Cmat = self.Wij_task_dict[task_name]
            onset_time = self.onset_time_list[i]
            duration = self.duration_list[i]
            start_time_block = onset_time
            self._generate_single_block(Cmat,
                                        duration=duration,
                                        activity=activity,
                                        a_s_rate=a_s_rate,
                                        syn_act=syn_act)
            if bold_chunkwise:
                self.generate_bold_chunkwise(TR=TR,
                                             input_type=self.output_activation,
                                             normalize_max=normalize_max,
                                             is_first=False,
                                             fix=self.fix_bold,
                                             **kwargs)
            end_time_block = onset_time + duration
            self.time_idxs_dict[task_name].append([round(start_time_block, 3),
                                                   round(end_time_block, 3)])
            if i < len(self.onset_time_list) - 1:
                duration = self.onset_time_list[i + 1] - \
                           self.onset_time_list[i] - self.duration_list[i]
                if duration > 0:
                    Cmat = self.Wij_rest
                    start_time_rest = self.onset_time_list[i] + self.duration_list[i]
                    end_time_rest = self.onset_time_list[i + 1]
                    self._generate_single_block(Cmat,
                                                duration=duration,
                                                activity=activity,
                                                a_s_rate=a_s_rate,
                                                syn_act=syn_act)
                    if bold_chunkwise:
                        self.generate_bold_chunkwise(TR=TR,
                                                     input_type=self.output_activation,
                                                     normalize_max=normalize_max,
                                                     is_first=False,
                                                     fix=self.fix_bold,
                                                     **kwargs)
                    self.time_idxs_dict["Rest"].append([round(start_time_rest, 3),
                                                        round(end_time_rest, 3)])
            else:
                self._generate_last_rest(activity=activity,
                                         a_s_rate=a_s_rate,
                                         syn_act=syn_act)
                if bold_chunkwise:
                    self.generate_bold_chunkwise(TR=TR,
                                                 input_type=self.output_activation,
                                                 normalize_max=normalize_max,
                                                 is_first=False,
                                                 fix=self.fix_bold,
                                                 **kwargs)

            if clear_raw:
                self.wc.exc = []
                self.wc.inh = []
            self.wc.t = []
            self.wc.outputs = dotdict()
            self.wc.state = dotdict()

    def generate_neuronal_oscill(self):
        """ function for computing integrated synaptic
        activity (without integration over 50ms)
        as described in Horwitz and Tagamets (1999)
        IN = wEE*E + WEI*E+WIE*IE +C*E

        Args:
            source (string): basis for computing isa series,
             possible values raw or act

        Returns:
            isa_series: np.array with the same shape as input activity
        """
        ee = self.wc.params['c_excexc']
        ei = self.wc.params['c_excinh']
        ie = self.wc.params['c_inhexc']
        ii = self.wc.params['c_inhinh']

        exc = self.wc.exc
        inh = self.wc.inh
        sa = ee * exc + ei * exc + ie * inh + ii * inh + self.wc.Cmat @ exc
        return sa

    def generate_coactivation_by_mat(self,
                                     mat_path,
                                     act_scaling,
                                     fix_bold=True,
                                     **kwargs):
        """
        Generate outer activation for each node defined with a tasks and activation info, where
        written which model sensitive to which task. All information defined in mat file
        with the description of  tasks, onsets, durations and activation info for each modules

        Args:
            mat_path (str): path to matfile
            act_scaling (float): scaling factor for hrf function
            **kwargs : kwargs for HRF class

        Returns:

        """

        num_regions = self.wc.params["N"]
        num_regions_per_modules = self.num_regions_per_modules
        first_rest = self.first_duration
        last_rest = self.last_duration
        TR = self.TR
        dt = self.wc.params["dt"]
        input_data = io.loadmat(mat_path)
        num_tasks = input_data['onsets'].shape[1]
        onsets_list = []
        activations = []
        durations_list = []
        for i in range(num_tasks):
            onsets_list.append(list(input_data['onsets'][0, i].squeeze()))
            activations.append(input_data['activations'][0, i].squeeze())
            durations_list.append(input_data["durations"][0, i].squeeze())

        box_car_activations = create_task_design_activation(onsets_list,
                                                            durations_list,
                                                            dt,
                                                            first_rest=first_rest,
                                                            last_rest=last_rest)
        activations_by_module = create_activations_per_module(activations,
                                                              box_car_activations)
        hrf = HRF(self.num_modules, dt=dt, TR=TR, normalize_max=act_scaling, fix=fix_bold)
        hrf.bw_convolve(activations_by_module, append=False, **kwargs)
        t_res_activ, res_activ = hrf.resample_to_TR(activations_by_module)

        res_activ = create_reg_activations(res_activ,
                                           num_regions,
                                           num_regions_per_modules)
        convolved_activ = create_reg_activations(hrf.BOLD, num_regions)
        return t_res_activ, res_activ, convolved_activ



    def generate_coactivations(self,
                               mat_path,
                               act_scaling=0.5,
                               all_rois=True,
                               fix_bold=True,
                               **kwargs):
        N = self.wc.params["N"]
        first_rest = self.first_duration
        last_rest = self.last_duration
        TR = self.TR
        dt = self.wc.params["dt"]
        input_data = io.loadmat(mat_path)
        task_duration = float(input_data["durations"][0, 0].squeeze())

        onset_taskA = list(input_data['onsets'][0, 0].squeeze().round(2))
        onset_taskB = list(input_data['onsets'][0, 1].squeeze().round(2))
        onset_taskAB, _, _ = read_onsets_from_mat(mat_path)
        onset_taskAB = list(onset_taskAB)
        activations_A = input_data['activations'][0][0].squeeze()
        activations_B = input_data['activations'][0][1].squeeze()
        assert len(activations_A) == len(activations_B), "Length of activations arrays should be equal for all task"
        assert len(activations_A) == self.num_modules, "Length of activations should be equal to number of modules"
        onsets = []

        for i in range(self.num_modules):
            if (activations_A[i] == 1) and (activations_B[i] == 1):
                # onsets.extend([onset_taskAB]*(N // self.num_modules))
                onsets.extend([onset_taskAB])
            elif (activations_A[i] == 1) and (activations_B[i] == 0):
                # onsets.extend([onset_taskA] * (N // self.num_modules))
                onsets.extend([onset_taskA])
            else:
                # onsets.extend([onset_taskB] * (N // self.num_modules))
                onsets.extend([onset_taskB])
        # hrf = HRF(N, dt=dt, TR=TR, normalize_max=act_scaling)
        hrf = HRF(self.num_modules, dt=dt, TR=TR, normalize_max=act_scaling, fix=fix_bold)
        local_activation = hrf.create_task_design_activation(onsets,
                                                             duration=task_duration,
                                                             first_rest=first_rest,
                                                             last_rest=last_rest)
        hrf.bw_convolve(local_activation, append=False, **kwargs)
        t_res_activ, res_activ = hrf.resample_to_TR(local_activation)

        if all_rois:
            return t_res_activ, np.repeat(res_activ,
                                          repeats=N // self.num_modules, axis=0), np.repeat(hrf.BOLD,
                                                                                            repeats=N // self.num_modules,
                                                                                            axis=0)
        else:
            return t_res_activ, res_activ, hrf.BOLD

    def generate_bold_chunkwise(self,
                                TR: float = 2,
                                input_type: str = "syn_act",
                                normalize_max: float = 2,
                                is_first: bool = True,
                                fix:bool = True,
                                **kwargs):
        """Custom bold signal generation on top of neuronal oscillations with
        user-defined TR and type of bold signal

        Args:
            fix: (bool): is all parametetrs for BW model os fixed with mean, or variable for each node
            is_first (bool): if it the first chunk for generation
            (do not need to init with previous values)
            normalize_max: normalization constant
            input_type(): what signal to convolve, only excitation or both or
            synaptic activity by formula
            TR (float): time repetition for bold signal resampling, default 2c
         """

        self.TR = TR
        N = self.wc.params["N"]
        dt = self.wc.params["dt"]
        chunksize = TR * 1000 / dt

        assert input_type in ['exc', 'sum', 'syn_act'], "input type could be 'exc', 'sum' or 'syn_act'"
        if is_first:
            input_rest = np.array([], dtype=np.float64).reshape(N, 0)
        else:
            input_rest = self.input_rest
        if input_type == 'exc':
            new_exc = np.hstack((input_rest, self.wc['exc']))
        elif input_type == 'sum':
            new_exc = np.hstack((input_rest, self.wc['exc'] + self.wc['inh']))
        elif input_type == 'syn_act':
            new_exc = np.hstack((input_rest, self.generate_neuronal_oscill()))

        if is_first:
            self.hrf = HRF(N,
                           dt=dt,
                           TR=TR,
                           normalize_input=True,
                           normalize_max=normalize_max,
                           fix=fix)

        if new_exc.shape[1] > chunksize:
            used_last_idxs = int(new_exc.shape[1] - new_exc.shape[1] % chunksize)
            self.input_rest = new_exc[:, used_last_idxs:]
            bold_input = new_exc[:, :used_last_idxs]

            if is_first:
                self.hrf.bw_convolve(bold_input, append=False, **kwargs)
            else:
                self.hrf.bw_convolve(bold_input, append=True, **kwargs)
            self.BOLD = self.hrf.BOLD
            self.t_BOLD = self.hrf.t_BOLD
        else:
            self.input_rest = new_exc

    def draw_envelope_bold_compare(self, node_id=2, series_name='sa_series',
                                   low_f=10, high_f=50, low_pass=None,
                                   drop_first_sec=7, shift_sec=4, plot_first=1, to_plot=True):

        a_s_rate = self.activity["sampling_rate"]
        TR = self.TR
        nyquist = 1 / a_s_rate / 2
        plot_first_dt = int(plot_first / a_s_rate)
        raw_signal = self.activity[series_name][node_id, :]
        high_band = high_f / nyquist
        low_band = low_f / nyquist
        b1, a1 = signal.butter(4, [low_band, high_band], btype='bandpass')

        emg_filtered = signal.filtfilt(b1, a1, raw_signal)

        if low_pass is not None:
            low_pass = low_pass / nyquist
            b2, a2 = signal.butter(4, low_pass, btype='lowpass')
            hilbert_envelope = signal.filtfilt(b2, a2, abs(emg_filtered))
        else:
            hilbert_envelope = np.abs(signal.hilbert(emg_filtered))
            # plt.plot(emg_envelope)

        drop_first_sec_TR = int(drop_first_sec / TR)
        drop_first_sec_dt = int(drop_first_sec / a_s_rate)
        step = int(TR / a_s_rate)
        bold_scaled = normalize(self.BOLD[node_id, drop_first_sec_TR:]).flatten()
        # compute all shifts
        shift_list_sec = list(np.linspace(-3, 7, 41))
        rcoeff_list = []
        for shift in shift_list_sec:
            shift = int(shift / a_s_rate)
            env_scaled_shifted = normalize(hilbert_envelope[drop_first_sec_dt - shift::step]).flatten()
            sig_len = min(len(env_scaled_shifted), len(bold_scaled))
            rcoeff, p_val = stats.pearsonr(env_scaled_shifted[:sig_len], bold_scaled[:sig_len])
            rcoeff_list.append(rcoeff)

        shift = int(shift_sec / a_s_rate)
        env_scaled_shifted = normalize(hilbert_envelope[drop_first_sec_dt - shift::step]).flatten()
        sig_len = min(len(env_scaled_shifted), len(bold_scaled))
        rcoeff, p_val = stats.pearsonr(env_scaled_shifted[:sig_len], bold_scaled[:sig_len])
        time = np.arange(sig_len) * TR

        if to_plot:
            fig = plt.figure(figsize=(12, 8))
            gs = fig.add_gridspec(3, 3)
            ax1 = fig.add_subplot(gs[0, :])
            # ax1.set_title('gs[0, :]')
            ax2 = fig.add_subplot(gs[1, :])
            # ax2.set_title('gs[1, :]')
            ax3 = fig.add_subplot(gs[-1, :-1])
            # ax3.set_title('gs[2, :-1]')
            ax4 = fig.add_subplot(gs[-1, -1])
            # ax4.set_title('gs[-1, -1]')
            ax1.plot(raw_signal[:plot_first_dt]);
            ax1.set_title("Raw neuronal activity")
            ax2.plot(emg_filtered[:plot_first_dt]);
            ax2.set_title(f"Filtered neuronal with low:{low_f} high:{high_f}")

            ax2.plot(hilbert_envelope[:plot_first_dt]);

            ax3.plot(time, env_scaled_shifted[:sig_len], label="Envelope");
            ax3.plot(time, bold_scaled[:sig_len], 'orange', label="BOLD");
            ax3.legend();
            ax3.set_title(f"Shifted envelope with {shift_sec} s and BOLD, rcoeff {rcoeff:.2f}, p_val {(p_val):.3f} ")
            ax4.plot(shift_list_sec, rcoeff_list)
            ax4.set_title("Bold-Envelope correlation with different time lag")
            ax4.set_xlabel("Time lag (seconds) ")
            ax4.set_ylabel("Pearson r")

            fig.tight_layout()

        return {'rs': np.array([shift_list_sec, rcoeff_list]),
                'raw_first': raw_signal[:plot_first_dt],
                'bold_envelope': np.array([time, env_scaled_shifted[:sig_len], bold_scaled[:sig_len]]),
                'filtered': np.array([emg_filtered[:plot_first_dt], hilbert_envelope[:plot_first_dt]])}

    def compute_phase_diff(self, series_name='sa_series', low_f=30, high_f=40, return_xr=True):
        activity = self.activity[series_name]
        act_time = self.activity['t']
        N_ROIs = activity.shape[0]
        s_rate = self.activity["sampling_rate"]
        coeff = 1 / s_rate
        zero_shift = self.time_idxs_dict["Rest"][0][0]
        len_tasks = int(np.ceil(coeff * (self.time_idxs_dict["Rest"][-1][1] - zero_shift)))
        assert len_tasks <= activity.shape[1], 'Computed length and series len should be less or equal'
        if len_tasks < activity.shape[1]:
            activity = activity[:,:len_tasks]
            act_time = self.activity['t'][:len_tasks]
        task_type = np.array(['Rest'] * int(len_tasks))
        trial_time_point = -1 + np.zeros(int(len_tasks), dtype=int)
        tasks = ["Task_A", "Task_B"]
        trial_number = -1 + np.zeros(int(len_tasks), dtype=int)
        for task in tasks:
            for i in range(len(self.time_idxs_dict[task])):
                idx_start = int((self.time_idxs_dict[task][i][0] - zero_shift) * coeff)
                idx_end = int((self.time_idxs_dict[task][i][1] - zero_shift) * coeff)
                task_type[idx_start:idx_end] = task.split('_')[-1]
                trial_number[idx_start:idx_end] = int(i)
                trial_time_point[idx_start:idx_end] = np.arange(idx_end - idx_start)

        roi_idx1, roi_idx2 = np.triu_indices(N_ROIs, k=1)
        # phase_diffs = np.zeros((roi_idx1.shape[0], int(len_tasks)), dtype=complex)
        phase_diffs = np.zeros((roi_idx1.shape[0], int(len_tasks)))
        nyquist = 1 / s_rate / 2
        high_band = high_f / nyquist
        low_band = low_f / nyquist
        b1, a1 = signal.butter(4, [low_band, high_band], btype='bandpass')
        filtered_data = signal.filtfilt(b1, a1, activity)
        analytic_data = signal.hilbert(filtered_data)
        angles = np.angle(analytic_data)
        for r in range(roi_idx1.shape[0]):
            # phase_diffs[r, :] = np.exp(1j * (angles[roi_idx1[r], :] - angles[roi_idx2[r], :]))
            phase_diffs[r, :] = angles[roi_idx1[r], :] - angles[roi_idx2[r], :]
        act_dict = {'activity': activity, 'phase_diff': phase_diffs,
                    'time': act_time, 's_rate': s_rate, 'task_type': task_type,
                    'trial_time': trial_time_point, 'trial_number': trial_number}
        if return_xr:
            act_vars = {'neural_activity': (['region', 'time'], act_dict['activity'],
                                            {'long_name': 'wc simulated wc neural activity'}),
                        'phase_diff': (['reg_reg', 'time'], act_dict['phase_diff'])
                        }

            coords = {'time': (['time'], act_dict['time'], {'units': 'm/s',
                                                            'sampling_rate': act_dict['s_rate']}),
                      'task_type': ('time', act_dict['task_type']),
                      'trial_number': ('time', act_dict['trial_number']),
                      'trial_time_point': ('time', act_dict['trial_time'])
                      }
            attrs = {"model": "WC"}
            ds = xr.Dataset(data_vars=act_vars,
                            coords=coords, attrs=attrs)

            return ds
        else:
            return act_dict

    def generate_bold(self,
                      bold_input,
                      dt=5,
                      TR=2,
                      drop_first=12,
                      normalize_max=2,
                      conv_type='BW',
                      fix=True,
                      **kwargs):
        """
        Generate bold signal on top of bold_input signal
        Args:
            bold_input (np.ndarry): with shape N rois vs time length, default excitation activity
            TR (float): time repetition
            drop_first (float): number of seconds droped after bold convolving
            conv_type (string): how to convolve with neurolib or with gamma, possible values "BW" or "Gamma"

        Returns:
            Convolved Bold signal

        """

        N = bold_input.shape[0]
        hrf = HRF(N, dt=dt, TR=TR, normalize_input=True, normalize_max=normalize_max, fix=fix)
        if conv_type == "BW":
            hrf.bw_convolve(bold_input, append=False, **kwargs)
        elif conv_type == "Gamma":
            hrf.gamma_convolve(bold_input, append=False, **kwargs)
        else:
            NotImplementedError
        drop_idxs = int(drop_first / TR)
        return (hrf.t_BOLD[drop_idxs:], hrf.BOLD[:, drop_idxs:])


def normalize(signal):
    if signal.ndim == 1:
        signal = signal.reshape(1, -1)
    mean_ = np.mean(signal, axis=1).reshape(-1, 1)
    std_ = np.std(signal, axis=1).reshape(-1, 1)
    return (signal - mean_) / std_


class HRF:
    """
        Balloon-Windkessel BOLD simulator class.
        BOLD activity is downsampled according to TR.
        BOLD simulation results are saved in t_BOLD, BOLD instance attributes.
        parameter fix responsible for the variability parameters
    """

    def __init__(self,
                 N: int,
                 dt: float = 10,
                 TR: float = 2,
                 normalize_input: bool = True,
                 normalize_max: float = 50,
                 fix: bool = True):
        self.N = N
        self.dt = dt  # in ms
        self.TR = TR  # in seconds
        self.normalize_input = normalize_input
        self.normalize_max = normalize_max
        self.fix = fix
        self.samplingRate_NDt = int(round(TR * 1000 / dt))

        # return arrays
        self.t_BOLD = np.array([], dtype="f", ndmin=2)
        self.BOLD = np.array([], dtype="f", ndmin=2)
        self.BOLD_chunk = np.array([], dtype="f", ndmin=2)

        self.idxLastT = 0  # Index of the last computed t

        # initialize BOLD model variables
        # self.X_BOLD = np.ones((N,))
        self.X_BOLD = np.zeros((N,))
        # Vasso dilatory signal
        self.F_BOLD = np.ones((N,))
        # Blood flow
        self.Q_BOLD = np.ones((N,))
        # Deoxyhemoglobin
        self.V_BOLD = np.ones((N,))
        # Blood volume

    def create_task_design_activation(self,
                                      onsets,
                                      duration,
                                      first_rest=5,
                                      last_rest=5):
        """
        Create external activation separately for each region
        Args:
            onsets (list of list of int or list): onset list for each region, for example [10, 12, 15], N lists
            duration (float or list of lists): duration of each task
            last_rest (float): duration of the last rest part

        Returns:

        """
        if isinstance(onsets[0], (int, float)):
            onsets = self.N * [onsets]  # just duplicate for all regions
        max_onset = np.max([np.max(onset) for onset in onsets])
        length = int((max_onset + duration + last_rest) * 1000 / self.dt)
        length_first_rest = int(first_rest * 1000 / self.dt)
        activation = np.zeros((self.N, length))
        assert isinstance(onsets, list), "Onsets should be a list or list of lists"

        for i in range(self.N):
            for onset in onsets[i]:
                start = int(round((1000 / self.dt) * onset))
                end = int(round((1000 / self.dt) * (onset + duration)))
                activation[i, start:end] = 1
        return np.hstack((np.zeros((self.N, length_first_rest)), activation))

    def resample_to_TR(self,
                       signal,
                       idxLastT=0):
        """ Resampling made with accordance to neurolib
        Args:
            signal (np.ndaray): numpy nd array

        Returns:
            resampled to TR signal
        """
        signal_resampled = signal[:,
                           self.samplingRate_NDt - np.mod(idxLastT - 1,
                                                          self.samplingRate_NDt):: self.samplingRate_NDt]
        t_new_idx = idxLastT + np.arange(signal.shape[1])
        t_resampled = (
                t_new_idx[self.samplingRate_NDt - np.mod(idxLastT - 1,
                                                         self.samplingRate_NDt):: self.samplingRate_NDt]
                * self.dt
        )
        return t_resampled, signal_resampled

    def bw_convolve(self,
                    activity,
                    append=False,
                    **kwargs):
        assert activity.shape[0] == self.N, "Input shape must be equal to Number of activations to times"
        if self.normalize_input:
            activity = self.normalize_max * activity

            # Compute the BOLD signal for the chunk
        BOLD_chunk, self.X_BOLD, self.F_BOLD, self.Q_BOLD, self.V_BOLD = simulateBOLD(activity,
                                                                                      self.dt * 1e-3,
                                                                                      X=self.X_BOLD,
                                                                                      F=self.F_BOLD,
                                                                                      Q=self.Q_BOLD,
                                                                                      V=self.V_BOLD,
                                                                                      fix=self.fix,
                                                                                      **kwargs)

        t_BOLD_resampled, BOLD_resampled = self.resample_to_TR(BOLD_chunk, idxLastT=self.idxLastT)

        if self.BOLD.shape[1] == 0:
            # add new data
            self.t_BOLD = t_BOLD_resampled
            self.BOLD = BOLD_resampled
        elif append is True:
            # append new data to old data
            self.t_BOLD = np.hstack((self.t_BOLD, t_BOLD_resampled))
            self.BOLD = np.hstack((self.BOLD, BOLD_resampled))
        else:
            # overwrite old data
            self.t_BOLD = t_BOLD_resampled
            self.BOLD = BOLD_resampled

        self.BOLD_chunk = BOLD_resampled

        self.idxLastT = self.idxLastT + activity.shape[1]

    def gamma_convolve(self, activity, append=False, **kwargs):
        assert activity.shape[0] == self.N, "Input shape must be equal to Number of activations to times"
        if self.normalize_input:
            activity = self.normalize_max * activity
        hrf_at_dt = self._gamma_hrf(**kwargs)
        BOLD_chunk = np.zeros_like(activity)
        for i in range(self.N):
            BOLD_chunk[i, :] = np.convolve(activity[i], hrf_at_dt)[:-(len(hrf_at_dt) - 1)]

        BOLD_resampled = BOLD_chunk[
                         :, self.samplingRate_NDt - np.mod(self.idxLastT - 1,
                                                           self.samplingRate_NDt):: self.samplingRate_NDt
                         ]
        t_new_idx = self.idxLastT + np.arange(activity.shape[1])
        t_BOLD_resampled = (
                t_new_idx[
                self.samplingRate_NDt - np.mod(self.idxLastT - 1,
                                               self.samplingRate_NDt):: self.samplingRate_NDt]
                * self.dt
        )

        if self.BOLD.shape[1] == 0:
            # add new data
            self.t_BOLD = t_BOLD_resampled
            self.BOLD = BOLD_resampled
        elif append is True:
            # append new data to old data
            self.t_BOLD = np.hstack((self.t_BOLD, t_BOLD_resampled))
            self.BOLD = np.hstack((self.BOLD, BOLD_resampled))
        else:
            # overwrite old data
            self.t_BOLD = t_BOLD_resampled
            self.BOLD = BOLD_resampled

        self.BOLD_chunk = BOLD_resampled

        self.idxLastT = self.idxLastT + activity.shape[1]

    def _gamma_hrf(self, length=32, peak=6, undershoot=16, beta=0.1667, scaling=0.6):
        """ Return values for HRF at given times

        Args:
            peak (float):  time to peak (in seconds)
            undershoot:
            beta:
            scaling:

        Returns:

        """

        # Gamma pdf for the peak
        from scipy.stats import gamma
        times = np.arange(0, length, self.dt * 1e-3)
        peak_values = gamma.pdf(times, peak)
        # Gamma pdf for the undershoot
        undershoot_values = gamma.pdf(times, undershoot)
        # Combine them
        values = peak_values - beta * undershoot_values
        # Scale max to 0.6
        return values / np.max(values) * scaling
