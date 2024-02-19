from scipy import io
from .synaptic_weights_matrices import normalize
from .synaptic_weights_matrices import generate_synaptic_weights_matrices
from typing import Optional


def read_onsets_from_mat(mat_path: str) -> tuple[list, list, list]:
    """Extract onset moment from matlab files with the matlab file with all
     the information about task. Matfile represented with a structure with the
     next field:
        activations  - N arrays, where N corresponds to number of tasks,
        with the length equal to number of time moments with 0 and 1,
        where 1 means that at this onset task is active. Used only for outer activations
        durations : array with duration of each onset
        names: array with the task names
        onsets: time moments when new task started, number of tasks corresponds
        to the shape of array, number of vectors corresponds to number of tasks,
        one array per task
        rest matrix: matrix with synaptic connections between modules
        task_matrices: matrices with synaptic connections between modules,
        number of matrix should correponds to number of tasks


    :param mat_path:
    :return: 3 lists with onsets, task_names and duration
    """

    input_data = io.loadmat(mat_path)
    assert "onsets" in input_data.keys(), 'onset key should be in structure'
    assert "names" in input_data.keys(), 'names key should be in structure'
    assert "durations" in input_data.keys(), 'durations key should be in structure'

    num_tasks = input_data['onsets'].shape[1]
    onset_tasks = []
    names = []
    durations = []
    for i in range(num_tasks):
        cur_onset = list(input_data['onsets'][0, i].squeeze().round(2))
        onset_tasks.extend(cur_onset)
        names.extend([input_data['names'][0, i][0]] * len(cur_onset))
        durations.extend(
            [float(input_data['durations'][0, i].squeeze())] * len(cur_onset))
    onset_time_list, task_names_list, duration_list = (
        list(zip(*sorted(list(zip(onset_tasks, names, durations))))))

    return onset_time_list, task_names_list, duration_list


def generate_sw_matrices_from_mat(mat_path: str,
                                  num_regions: int,
                                  num_modules: int = 3,
                                  num_regions_per_modules: Optional[list] = None,
                                  sigma: float = 0.01,
                                  norm_type: str = "cols",
                                  gen_type: str = 'simple_prod'):
    """
    Generate task and rest matrices from mat file
    Matfile represented with a structure with the
     next field:
        activations  - N arrays, where N corresponds to number of tasks,
        with the length equal to number of time moments with 0 and 1,
        where 1 means that at this onset task is active. Used only for outer activations
        durations : array with duration of each onset
        names: array with the task names
        onsets: time moments when new task started, number of tasks corresponds
        to the shape of array, number of vectors corresponds to number of tasks,
        one array per task
        rest matrix: matrix with synaptic connections between modules
        task_matrices: matrices with synaptic connections between modules,
        number of matrix should correponds to number of tasks

       num_regions (int):
            number of regions generated during the simulation
        num_modules (int):
            number of modules (or connected block)
        num_regions_per_modules (list of int):
            number of regions in each module (should sum to num_regions)
        gen_type (str):
            if simple_prod  -
            generation is equal to scaling normal distribution
            with factors, else - scaling with equal variance,
            possible values [simple_prod, equal_var]

    return: rest matrix and list of task matrices
    """
    # TODO add num regions per module
    input_data = io.loadmat(mat_path)
    coeff_rest_matrix = input_data["rest_matrix"]
    coeff_task_matrices = input_data["task_matrices"]
    num_tasks = coeff_task_matrices.shape[1]
    names = []
    for i in range(num_tasks):
        names.append(input_data["names"][0, i][0])
    rest_factors = coeff_rest_matrix[0, 0]
    Wij_task_list = []
    for i in range(num_tasks):
        Wij_task = generate_synaptic_weights_matrices(num_regions,
                                                      num_modules,
                                                      num_regions_per_modules,
                                                      factors=coeff_task_matrices[0, i],
                                                      sigma=sigma,
                                                      gen_type=gen_type)
        Wij_task = normalize(Wij_task,
                             norm_type=norm_type)
        Wij_task_list.append(Wij_task)

    Wij_rest = generate_synaptic_weights_matrices(num_regions,
                                                  num_modules,
                                                  num_regions_per_modules,
                                                  factors=rest_factors,
                                                  sigma=sigma,
                                                  gen_type=gen_type)
    Wij_rest = normalize(Wij_rest, norm_type=norm_type)
    Wij_task_dict = dict(zip(names, Wij_task_list))

    return Wij_rest, Wij_task_dict
