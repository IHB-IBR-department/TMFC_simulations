from scipy import io
from task_fc_simulation.weight_matrix_utils import normalize, generate_modulars
import numpy as np


def read_onsets_from_input(mat_path):
    input_data = io.loadmat(mat_path)
    assert "onsets" in input_data.keys()
    num_tasks = input_data['onsets'].shape[1]
    onset_tasks = []
    names = []
    durations = []
    for i in range(num_tasks):
        cur_onset = list(input_data['onsets'][0, i].squeeze().round(2))
        onset_tasks.extend(cur_onset)
        names.extend([input_data['names'][0, i][0]] * len(cur_onset))
        durations.extend([float(input_data['durations'][0, i].squeeze())] * len(cur_onset))
    onset_time_list, task_names_list, duration_list = list(zip(*sorted(list(zip(onset_tasks, names, durations)))))

    return onset_time_list, task_names_list, duration_list


def read_generate_task_matrices(mat_path, num_regions, num_modules=3,
                                sigma=0.01, norm_type="cols", gen_type='simple_prod'):
    """
    Generate task and rest matrices from mat file
    todo: description of mat file
    return: rest matrix and list of task matrices
    """
    input_data = io.loadmat(mat_path)
    coeff_rest_matrix = input_data["rest_matrix"]
    coeff_task_matrices = input_data["task_matrices"]
    num_tasks = coeff_task_matrices.shape[1]
    names = []
    for i in range(num_tasks):
        names.append(input_data["names"][0, i][0])
    rest_factors = coeff_rest_matrix[0, 0]
    C_task_list = []
    for i in range(num_tasks):
        C_task = generate_modulars(num_regions, num_modules, factors=coeff_task_matrices[0, i],
                                   sigma=sigma, gen_type=gen_type)
        C_task = normalize(C_task, norm_type=norm_type)
        C_task_list.append(C_task)

    C_rest = generate_modulars(num_regions, num_modules, factors=rest_factors,
                               sigma=sigma, gen_type=gen_type)
    C_rest = normalize(C_rest, norm_type=norm_type)
    C_task_dict = dict(zip(names, C_task_list))

    return C_rest, C_task_dict