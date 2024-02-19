import numpy as np
from typing import Optional, Union
import numpy.typing as npt


def create_task_design_activation(onsets_list: list[list],
                                  durations_list: Union[list[float], list[list]],
                                  dt: float = 0.1,
                                  first_rest: float = 5,
                                  last_rest: float = 5) -> npt.NDArray:
    """
    Create external activation array separately for each task, return box car with the same size as task

    Args:
        onsets_list (list of list of float or list): onset list for each task,
                    for example [10, 12, 15], N lists equal to number of tasks, onsets in seconds
        durations_list (list of float or list of lists): duration of each task,
        one number for each task or list of durations corresponds to each onset
        dt (float): sampling time in ms, i.e. 0.1 means 0.1ms
        last_rest (float): duration of the last rest part
        first_rest (float): duration of the first part


    Returns:
        box-car function corresponded to design
    """
    n_tasks = len(onsets_list)
    if (n_tasks == 1) & ~isinstance(onsets_list[0], list):
        onsets_list = [onsets_list]
    assert all([isinstance(onsets, list) for onsets in onsets_list]), \
        'For each task should be onset list'
    assert len(durations_list) == len(onsets_list), \
        "Duration should be specified for each task"
    # find max onset and task corresponded
    max_onsets = [np.max(onsets) for onsets in onsets_list]
    task_max_onset, max_onset = np.argmax(max_onsets), np.max(max_onsets)

    for i, durations in enumerate(durations_list):
        if isinstance(durations, (list, tuple)):
            assert len(durations) == len(onsets_list[i]), \
                "if list of duration provided, length should corresponds to onsets list"
        else:
            durations_list[i] = [durations] * len(onsets_list[i])
    max_onset = np.max([np.max(onsets) for onsets in onsets_list])
    max_duration = durations_list[task_max_onset][-1]
    length = int((max_onset + max_duration + last_rest) * 1000 / dt)
    length_first_rest = int(first_rest * 1000 / dt)
    activation = np.zeros((n_tasks, length))

    for i, onsets in enumerate(onsets_list):
        for onset, duration in zip(onsets, durations_list[i]):
            start = int(round((1000 / dt) * onset))
            end = int(round((1000 / dt) * (onset + duration)))
            activation[i, start:end] = 1
    return np.hstack((np.zeros((n_tasks, length_first_rest)), activation))


def module_activation(tasks_responded: Union[int, list[int]],
                      box_car_response: npt.NDArray):
    """

    Args:
        tasks_responded (int or list of int): task for which there is activation in module
        box_car_response (npt.NDArray):  numpy array with the box-car responses for each task

    Returns:

    """
    assert isinstance(tasks_responded, (int, list)), \
        "Variable should contain task indexes corresponded  with box car size"

    result = 0
    if isinstance(tasks_responded, int):
        result = box_car_response[tasks_responded]
    else:

        max_task_value = max(tasks_responded)
        assert max_task_value < len(box_car_response), \
            'Task number should corresponds to box_car_responce'
        for task_idx in tasks_responded:
            result += box_car_response[task_idx]
    return result


def create_activations_per_module(activations: list[Union[list[bool], list[int]]],
                                  box_car_response: npt.NDArray):
    """

    Args:
        activations: list of list with length equal of the number of tasks,
                    where in each list inside indicator function for each modules
        box_car_response: array with the shape equal to number of tasks

    Returns: numpy array with activations for each module

    """

    task_numbers = len(activations)
    assert task_numbers == box_car_response.shape[0], \
        "Number of tasks should be equal to box-car series numbers "

    num_modules = len(activations[0])
    activations_by_module = np.zeros((num_modules, box_car_response.shape[1]))
    for module in range(num_modules):
        tasks_responded = [i for i in range(task_numbers) if activations[i][module] > 0]
        activations_by_module[module] = module_activation(tasks_responded, box_car_response)

    return activations_by_module


def create_reg_activations(activations_by_module,
                           num_regions: int,
                           num_regions_per_modules: Optional[Union[int, list]] = None):
    num_modules = activations_by_module.shape[0]
    if num_regions_per_modules is None:
        num_equal = int(round(num_regions / num_modules))
        num_regions_per_modules = (
                (num_modules - 1) * [num_equal]
                + [num_regions - (num_modules - 1) * num_equal])
    assert np.sum(num_regions_per_modules) == num_regions
    reg_activation = np.zeros((num_regions, activations_by_module.shape[1]))
    last_filled_region = 0
    for m in range(num_modules):
        reg_activation[last_filled_region:last_filled_region+num_regions_per_modules[m], :] = activations_by_module[m]
        last_filled_region += num_regions_per_modules[m]
    return reg_activation
