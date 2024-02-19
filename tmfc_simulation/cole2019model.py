import numpy as np
from scipy.special import expit as sigmoid
from typing import Optional


class NeuralMassModel:
    """
    Refactored code for running the neural mass simulation from:
     Cole MW, Ito T, Schultz D, Mill R, Chen R, Cocuzza C (2019).
     "Task activations produce spurious but systematic inflation
     of task functional connectivity estimates".
     NeuroImage. doi:10.1016/j.neuroimage.2018.12.054
     https://github.com/ColeLab/TaskFCRemoveMeanActivity/blob/master/neuralmassmodel/NeuralMassModel.ipynb
    """

    def __init__(self, num_regions: int, num_modules: int = 3,
                 num_regions_per_modules: Optional[list[int]] = None,
                 struct_conn_probs: Optional[dict] = None,
                 syn_com_mult: Optional[dict] = None,
                 syn_weight_std: float = 0.001):
        """

        Args:
            num_regions (int):
                number of interacted regions
            num_modules (int):
                number of modules in structural connection
            struct_conn_probs (dict):
                dictionary with structural connections probabilities
                    between the modules
            num_regions_per_module (list of int):
                len should be equal to num_modules, and sum of
                    regions should be equal to num_regions
            syn_com_mult (dict):
                dictionary with synaptic weights to existing
                    structural connections
            syn_weight_std (float):
                standard deviation for synaptic weight
        """
        self.num_regions = num_regions
        self.num_modules = num_modules
        self.syn_weight_std = syn_weight_std
        if struct_conn_probs is None:
            self.struct_conn_probs = {'in': 0.8, 'out': 0.2}
        else:
            self.struct_conn_probs = struct_conn_probs

        if num_regions_per_modules is None:
            num_equal = int(round(num_regions / num_modules))
            self.num_regions_per_modules = (num_modules - 1) * [num_equal] + [
                num_regions - (num_modules - 1) * num_equal]
        else:
            self.num_regions_per_modules = num_regions_per_modules
        assert num_regions == np.sum(
            self.num_regions_per_modules), """Sum number in each regions (num_regions_per_modules)) 
                                            should be equal to num_regions"""
        self.module_borders = ([0] +
                               list(np.cumsum(self.num_regions_per_modules))
                               )

        self._init_struct_matrix()

        if syn_com_mult is None:
            self.syn_com_mult = {
                0: [1, 1, 0],
                1: [1, 1, 0],
                2: [0, 0, 1]
            }
        else:
            self.syn_com_mult = syn_com_mult

        assert (len(self.syn_com_mult.keys()) == self.num_modules,
                "Number of keys should corresponds to number of modules ")
        self._init_synaptic_weights()

    def _init_struct_matrix(self):
        self.struct_matrix = (
                np.random.uniform(0, 1, (self.num_regions, self.num_regions))
                > 1 - self.struct_conn_probs['out'])
        for i in range(self.num_modules):
            self.struct_matrix[self.module_borders[i]:self.module_borders[i + 1],
            self.module_borders[i]:self.module_borders[i + 1]] = (
                    np.random.uniform(0, 1, (self.num_regions_per_modules[i],
                                             self.num_regions_per_modules[i]))
                    > 1 - self.struct_conn_probs['in']
            )

        np.fill_diagonal(self.struct_matrix, 1)

    def _init_synaptic_weights(self):
        """ User defined synaptic weight matrix with predefined structural matrix

        Returns: init synaptic weight matrix

        """
        self.synaptic_weight = (self.struct_matrix * (
                1 + np.random.standard_normal(
            (self.num_regions, self.num_regions))
                * self.syn_weight_std))
        for row in range(self.num_modules):
            for col in range(self.num_modules):
                self.synaptic_weight[self.module_borders[row]:self.module_borders[row + 1],
                self.module_borders[col]:self.module_borders[col + 1]] *= (
                    self.syn_com_mult[row][col])
        np.fill_diagonal(self.synaptic_weight, 0)
        for node in range(self.num_regions):
            k = np.sum(self.synaptic_weight[node, :])
            if k > 0:
                self.synaptic_weight[node, :] = (
                    np.divide(self.synaptic_weight[node, :], k))

    def _init_synpatic_cole(self):
        """
        Synaptic weights suggested in Cole
        https://github.com/ColeLab/
            TaskFCRemoveMeanActivity/blob/master/
                neuralmassmodel/NeuralMassModel.ipynb

        Cole MW, Ito T, Schultz D, Mill R, Chen R, Cocuzza C (2019).
        "Task activations produce spurious
         but systematic inflation of task functional
         connectivity estimates". NeuroImage.
        Returns: init synaptic weigths, unweighted matrix
        (with 1s indicating edges and 0s otherwise)

        """
        struct_conn_vector = np.random.uniform(0, 1, (self.num_regions, self.num_regions)) > .90
        # Add self-connections (important if adding autocorrelation later)
        np.fill_diagonal(struct_conn_vector, 10)
        self.struct_matrix = struct_conn_vector
        # Create modular structural network (3 modules)
        num_modules = 3
        numr_per_module = int(round(self.num_regions / num_modules))
        lastModuleNode = -1
        for moduleNum in range(0, num_modules):
            for thisNodeNum in range(lastModuleNode + 1,
                                     lastModuleNode +
                                     numr_per_module + 1):
                # Set this node to connect to 10 random other nodes in module
                for i in range(1,
                               numr_per_module // 2):
                    randNodeInModule = int(
                        np.random.uniform(lastModuleNode + 1,
                                          lastModuleNode + numr_per_module
                                          + 1, (1, 1)))
                    struct_conn_vector[thisNodeNum, randNodeInModule] = 1
            lastModuleNode = lastModuleNode + numr_per_module

        # Adding synaptic weights to existing structural connections
        # (small random synapse strength variation)
        synaptic_weight_vector = struct_conn_vector * (
                1 + np.random.standard_normal((self.num_regions,
                                               self.num_regions)) * .001)

        # Adding synaptic mini-communities (within community 1)
        synaptic_weight_vector[0:numr_per_module // 2,
        numr_per_module //
        2:numr_per_module] = (
                synaptic_weight_vector[0:numr_per_module // 2,
                numr_per_module // 2:numr_per_module] * -0.2)
        synaptic_weight_vector[numr_per_module // 2:
                               numr_per_module,
        0:numr_per_module // 2] = (
                synaptic_weight_vector[
                numr_per_module // 2:numr_per_module,
                0:numr_per_module // 2] * -0.2)
        synaptic_weight_vector[0:numr_per_module // 2,
        0:numr_per_module // 2] = (synaptic_weight_vector[
                                   0:numr_per_module // 2,
                                   0:numr_per_module // 2] * 1.2)
        synaptic_weight_vector[numr_per_module // 2:numr_per_module,
        numr_per_module // 2:
        numr_per_module] = (synaptic_weight_vector[
                            numr_per_module // 2:numr_per_module,
                            numr_per_module // 2:numr_per_module] * 1.2)

        # MODIFICATION: 0 connectivity between structural community 1 and 3
        synaptic_weight_vector[0:numr_per_module,
        2 * numr_per_module:
        3 * numr_per_module] = (synaptic_weight_vector[0:numr_per_module,
                                2 * numr_per_module:3 * numr_per_module] * 0)
        synaptic_weight_vector[2 * numr_per_module:
                               3 * numr_per_module,
                               0:numr_per_module] = synaptic_weight_vector[
                                                    2 * numr_per_module:3 * numr_per_module,
                                                    0:numr_per_module] * 0
        synaptic_weight_vector[numr_per_module:2 * numr_per_module,
        2 * numr_per_module:3 * numr_per_module] = synaptic_weight_vector[numr_per_module:2 * numr_per_module,
                                                   2 * numr_per_module:3 * numr_per_module] * 0
        synaptic_weight_vector[2 * numr_per_module:3 * numr_per_module,
        numr_per_module:2 * numr_per_module] = synaptic_weight_vector[2 * numr_per_module:3 * numr_per_module,
                                               numr_per_module:2 * numr_per_module] * 0

        # Normalize each region's inputs to have a mean of 1/k, where k is the number of incoming connections (in degree)
        # Based on Barral, J. & Reyes, A. D. Synaptic scaling rule preserves excitatory-inhibitory balance and salient neuronal network dynamics. Nat. Neurosci. (2016)
        # This ensures that all inputs into each node sum to 1
        for nodeNum in range(self.num_regions):
            k = np.sum(synaptic_weight_vector[nodeNum, :])
            if k > 0:
                synaptic_weight_vector[nodeNum, :] = np.divide(synaptic_weight_vector[nodeNum, :], k)

        self.synaptic_weight = synaptic_weight_vector

    def compute_network_model_cole(self, num_time_points: int, bias_param: float = -20,
                                   spont_act_level: float = 3,
                                   stim_regions: np.ndarray = np.array([0]),
                                   stim_times: np.ndarray = np.array([0]),
                                   g: float = 5.0, indep: float = 1.0, k: float = 1,
                                   stim_act_mult: float = 0.3, ind: int = 1):
        """
        Adapted to python 3 version of WC model from Cole (2019)
        https://github.com/ColeLab/TaskFCRemoveMeanActivity/blob/master/neuralmassmodel/NeuralMassModel.ipynb

        Args:
            num_time_points (int): number of time points to generate
            bias_param (float): (population resting potential, or excitability).
            spont_act_level (float): standard deviation for noise described spontaneous level of activity
            stim_regions (np.array, int):
            stim_times (np.array, int):
            g (float): global coupling parameter
            indep: parameter for self-connection modulation
            k (float):
            ind (int): random number for seed initialization

        Returns:
            dict with input and output activity

        """

        np.random.seed(np.random.randint(100 + ind))
        # np.random.seed(10)

        outputdata = {'input_activity': np.zeros((num_time_points, self.num_regions)),
                      'output_activity': np.zeros((num_time_points, self.num_regions))}

        bias = np.zeros(shape=(self.num_regions,))
        bias[range(self.num_regions)] = bias_param
        autocorr = 0.0
        global_coupling_mat = self.synaptic_weight.copy() * g
        np.fill_diagonal(global_coupling_mat, 0)
        indep_var_mat = np.identity(self.num_regions) * indep
        synaptic_weight = global_coupling_mat + indep_var_mat
        outputvect = np.zeros(self.num_regions)

        for this_time_point in range(0, num_time_points):
            stim_act_vector = np.zeros(self.num_regions)
            if this_time_point == 0:
                # Generate spontaneous activity for initial state
                act_vector = sigmoid(k * (bias + spont_act_level * np.random.normal(0, 1, (self.num_regions,))))
                input_activity = np.zeros(self.num_regions)
            else:
                # Bring outputs from previous time point as inputs for this timepoint
                act_vector = outputvect

            input_activity = synaptic_weight @ act_vector

            if np.any(stim_regions != 0):
                if this_time_point in stim_times:
                    stim_act = np.ones(len(stim_regions)) * stim_act_mult
                    stim_act_vector[stim_regions] = stim_act
            input_activity = input_activity + stim_act_vector + spont_act_level * np.random.normal(0, 1,
                                                                                                   (self.num_regions,))
            outputvect = sigmoid(k * (bias + input_activity))
            outputdata['input_activity'][this_time_point, :] = input_activity
            outputdata['output_activity'][this_time_point, :] = outputvect
        return outputdata
