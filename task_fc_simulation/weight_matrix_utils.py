import numpy as np


def generate_modulars(num_regions, num_modules, num_regions_per_modules=None,
                      factors=None, sigma=0.01, return_stats=False, gen_type='simple_prod'):
    """Function for generation of matrix with different module structure. Construction of the synaptic weight matrices
     involved three steps. First, synaptic weights (wji) were drawn from a Gaussian distribution (mean of 1,
     standard deviation of sigma) for each subject.
     Second, synaptic weights within and between functional modules were multiplied by weighting factors,
      that determined the network structure

    Args:
        return_stats (bool): if True also return mean and std by blocks
        sigma: standard deviation for gaussian distribution
        num_regions (int): number of regions generated during the simulation
        num_modules (int): number of modules
        num_regions_per_modules (list of int): number of regions in each module (should sum to num_regions)
        factors (list of list or np.ndarray): coefficient to multiply each factor
        gen_type (str): if simple_prod  - generation is equal to scaling normal dostribution with factors,
                        else - scaling with equal variance, possible values [simple_prod, equal_var]


    Returns:
        weight_matrix(np.ndarray of float): resulted weight matrix
    """
    assert gen_type in ["simple_prod", "equal_var"], "gen_type variable should be simple_prod or equal_var "

    if num_regions_per_modules == None:
        num_equal = int(round(num_regions / num_modules))
        num_regions_per_modules = (num_modules - 1) * [num_equal] + [num_regions - (num_modules - 1) * num_equal]
    module_borders = [0] + list(np.cumsum(num_regions_per_modules))
    if factors is None:
        factors = np.array([[0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]])
    if return_stats:
        stats = {'mean': np.zeros((num_modules, num_modules)),
                 'std': np.zeros((num_modules, num_modules))}
    assert np.array(factors).shape[0] == num_modules, "Number of modules should be compatible with the factors"
    if gen_type == "simple_prod":
        weight_matrix = 1 + np.random.normal(0, sigma, size=(num_regions, num_regions))
    else:
        weight_matrix = np.zeros((num_regions, num_regions))
    for row in range(num_modules):
        for col in range(num_modules):
            if gen_type == "simple_prod":
                weight_matrix[module_borders[row]:module_borders[row + 1],
                    module_borders[col]:module_borders[col + 1]] *= factors[row, col]

            else:

                block = weight_matrix[module_borders[row]:module_borders[row + 1],
                        module_borders[col]:module_borders[col + 1]]

                block = np.abs(np.random.normal(factors[row, col], sigma, size=block.shape))
                weight_matrix[module_borders[row]:module_borders[row + 1],
                module_borders[col]:module_borders[col + 1]] = block
            if return_stats:
                stats['mean'][row, col] = np.mean(weight_matrix[module_borders[row]:module_borders[row + 1],
                    module_borders[col]:module_borders[col + 1]])
                stats['std'][row, col] = np.std(weight_matrix[module_borders[row]:module_borders[row + 1],
                    module_borders[col]:module_borders[col + 1]])
    if return_stats:
        return weight_matrix, stats
    else:
        return weight_matrix


def normalize(weight_matrix, norm_type='sum'):
    norm_weight_matrix = weight_matrix.copy()

    if norm_type == 'cols':
        norm_weight_matrix = norm_weight_matrix / np.sum(norm_weight_matrix, axis=1)[:, None]
    elif norm_type == 'sym_cols':
        N = weight_matrix.shape[0]
        norm_weight_matrix = N * norm_weight_matrix / np.sum(weight_matrix)
    elif norm_type == 'max':
        norm_weight_matrix = norm_weight_matrix / np.max(weight_matrix)
    elif norm_type == 'sum':
        norm_weight_matrix = norm_weight_matrix / np.sum(weight_matrix)
    elif norm_type == 'raw':
        norm_weight_matrix = norm_weight_matrix

    return norm_weight_matrix
