import numpy as np


def mean_conf_interval(data, n_it, cover_percentage=95):
    """
    Compute the mean confidence interval using bootstrap resampling

    Parameters
    ----------
        data : np.ndarray
            1d or 2d array for which to compute the mean confidence intervals
            the LAST axis of `data``should correspond the axis along which the
            bootstrapping should take place. (See test_bootstrap.py)
            (i.e. one row of the array contains values
                for which the bootstrap interval is computed)
        n_it : the number of resamplings to
        cover_percentage : the total coverage in percentages (symmetric)

    Returns
    -------
    percentiles : np.ndarray
        for each

        The percentiles of the data [prop, percentile]
    """
    data = np.array(data)
    tail = (100 - cover_percentage) / 2.
    percentiles = [tail, 100 - tail]
    indices = np.random.randint(0, data.shape[-1], (n_it, data.shape[-1]))
    avgs = np.average(data[..., indices], axis=-1)
    return np.array(np.percentile(avgs, percentiles, axis=-1)).T


def mean_groupwise_conf_intervals_from_sim_matrix(matrix, n1, n_it,
                                                  cover_percentage):
    """
    Computes the estimates for the mean bootstrap confidence intervals
    'within-group' values of the matrix.

    With matrix of size (7,7) and n1 =3, the averages would be computed over
    the indices marked by '1' and '2'::

        0 1 1 0 0 0 0
        0 0 1 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 2 2 2
        0 0 0 0 0 2 2
        0 0 0 0 0 0 2
        0 0 0 0 0 0 0

    Parameters
    ----------
        matrix : np.ndarray
            1d or 2d array for which to compute the mean confidence intervals
            the LAST axis of `data``should correspond the axis along which the
            bootstrapping should take place. (See test_bootstrap.py)
            (i.e. one row of the array contains values
                for which the bootstrap interval is computed)
        n_it : the number of resamplings to
        cover_percentage : the total coverage in percentages (symmetric)

    Returns
    -------
    means: np.array of shape (2,)
        contains the 'within-group means' for the first and second group
    percentiles : np.array of shape (2,2)
        contains the corresponding bootstrap percentiles for the first and
        second group
    """
    n2 = matrix.shape[0] - n1
    indices1 = np.triu_indices(n1, k=1)
    indices2base = np.triu_indices(n2, k=1)
    indices2I = indices2base[0].copy() + n1
    indices2J = indices2base[1].copy() + n1
    indices2 = (indices2I, indices2J)
    values1 = matrix[indices1]
    values2 = matrix[indices2]
    mean1 = np.average(values1)
    mean2 = np.average(values2)
    mean_confs1 = mean_conf_interval(values1, n_it, cover_percentage)
    mean_confs2 = mean_conf_interval(values2, n_it, cover_percentage)
    return np.array([mean1, mean2]), np.vstack((mean_confs1, mean_confs2))
