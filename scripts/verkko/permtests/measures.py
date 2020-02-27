import numpy as np


def mean_difference(data_array, n1):
    """
    Computes the mean difference

    Parameters
    ----------
    data_array: numpy array
        1 or 2 dimensional numpy array, the n1 first elements in
        the zeroth axis of the array, should correspond to the
        values of the first group
    n1: int
        the number of elements in the first group

    Returns
    -------
    md: numpy array
        the mean difference
    """
    return (np.average(data_array[:n1], axis=0) -
            np.average(data_array[n1:], axis=0))


def paired_t_value(data_array, n1):
    """
    Computes the paired t-value for a one or
    two dimensional numpy array
    See wikipedia

    Parameters
    ----------
    data_array: numpy array
        1 or 2 dimensional numpy array, the n1 first elements in
        the zeroth axis of the array, should correspond to the
        values of the first group
    n1: int
        the number of elements in the first group

    Returns
    -------
    t-value(s) : float/numpy array
        the corresponding t-value(s)
    """
    assert len(data_array) / 2 == n1, "The data array is not "
    differences = data_array[:n1] - data_array[n1:]
    stds = np.std(differences, axis=0, ddof=1)
    avgs = np.average(differences, axis=0)
    return avgs / stds * np.sqrt(n1)


def unpaired_t_value(data_array, n1):
    """
    Computes the t-value (variance normalized mean difference) for a one or
    two dimensional numpy array

    Parameters
    ----------
    data_array: numpy array
        1 or 2 dimensional numpy array, the n1 first elements in
        the zeroth axis of the array, should correspond to the
        values of the first group
    n1: int
        the number of elements in the first group

    Returns
    -------
    t-value(s) : float/numpy array
        the corresponding t-value(s)
    """
    n2 = data_array.shape[0] - n1
    var1 = np.var(data_array[:n1], axis=0, ddof=1)
    var2 = np.var(data_array[n1:], axis=0, ddof=1)
    mu1 = np.average(data_array[:n1], axis=0)
    mu2 = np.average(data_array[n1:], axis=0)
    return (mu1 - mu2) / np.sqrt(var1 / n1 + var2 / n2)


def sim_matrix_within_group_means(matrix, n1):
    """
    Computes the mean of the upper triangle (k=1) for the blocks
    (0,n-1)*(0,n-1) and (n,2n-1)*(n,2n-1), and their difference
    (for convenience).

    Parameters
    ----------
    matrix : 2D symmetric numpy array
        1 or 2 dimensional numpy array, the n1 first indices in
        the zeroth axis of the array, should correspond to the
        values of the first group.
        The value of ``matrix[i][j]`` should correspond to
    n1 : int
        the number of elements in the first group

    Returns
    -------
    mean1 : float
        the average similarity between members in the first group
    mean2: float
        the average similarity between members in the second group
    mean1-mean2: float
        just mean1-mean2 (as a convenience for stat. testing)
    """
    n2 = matrix.shape[0] - n1
    indices1 = np.triu_indices(n1, k=1)
    indices2base = np.triu_indices(n2, k=1)
    indices2I = indices2base[0].copy() + n1
    indices2J = indices2base[1].copy() + n1
    indices2 = (indices2I, indices2J)
    mean1 = np.average(matrix[indices1])
    mean2 = np.average(matrix[indices2])
    return mean1, mean2, mean1 - mean2


def sim_matrix_mean_inter_group_similarity(mat, n1):
    """ Computes the average distance/similarity between groups

    Parameters
    ----------
    mat : 2D symmetric numpy array
        matrix : 2D symmetric numpy array
        1 or 2 dimensional numpy array, the n1 first indices in
        the zeroth axis of the array, should correspond to the
        values of the first group.
        The value of ``matrix[i][j]`` should correspond to
    n1 : int
        the number of elements in the first group

    Returns
    -------
    inter_group_mean : float
        the average distance between the two groups
    """
    n2 = np.shape(mat)[0] - n1
    # between group similarities:
    incidence_matrix12 = np.ones((n1, n2))
    between_I = incidence_matrix12.nonzero()[0].copy()  # read-only
    between_J = incidence_matrix12.nonzero()[1].copy()  # read-only
    between_J += n1
    between_indices = (between_I, between_J)
    inter_group_mean = np.average(mat[between_indices])
    return inter_group_mean


def paired_sim_matrix_inter_group_means(mat, n1=None):
    """
    Computes the inter-group average (and separately for the same subjects!)

    Parameters
    ----------
    mat : 2D symmetric numpy array
        matrix : 2D symmetric numpy array
        1 or 2 dimensional numpy array, the n1 first indices in
        the zeroth axis of the array, should correspond to the
        values of the first group.
        The value of ``matrix[i][j]`` should correspond to
    n1 : int, optional
        should equal to the number of members in the first group
        given as an input parameter only due to consistency issues within
        the permtests conventions

    Returns
    -------
    inter_group_mean : float
        the mean value of the inter-group area of the mat
    semidiag_mean : float
        the mean value of the `half-diagonal' corresponding
        to the same/paired subject in different conditions.
    semidiag_mean-inter_group_mean : float
        (convenience for statistical testing)
    """
    assert mat.shape[0] == mat.shape[1]
    assert len(mat.shape) == 2
    assert mat.shape
    if n1 is None:
        n1 = np.shape(mat)[0] / 2
    assert n1 is np.shape(mat)[0] / 2
    assert len(mat) == 2 * n1, "Matrix can not be a paired sim matrix, " \
        "should have even number of subjects"

    # between group similarities:
    incidence_matrix12 = np.ones((n1, n1))
    incidence_matrix12 -= np.eye(n1)
    between_I = incidence_matrix12.nonzero()[0].copy()
    between_I += n1
    between_J = incidence_matrix12.nonzero()[1].copy()
    between_indices = (between_I, between_J)
    inter_group_mean = np.average(mat[between_indices])

    semidiag_indices_I = np.eye(n1).nonzero()[0].copy()
    semidiag_indices_J = np.eye(n1).nonzero()[1].copy()
    semidiag_indices_J += n1
    semidiag_indices = (semidiag_indices_I, semidiag_indices_J)
    semidiag_mean = np.average(mat[semidiag_indices])
    return inter_group_mean, semidiag_mean, semidiag_mean - inter_group_mean


def sim_matrix_within_groups_mean_minus_inter_group_mean(mat, paired, n1=None):
    """
    Computes the difference between the average similarity within the (two)
    groups.

    Parameters
    ----------
    mat : 2D symmetric numpy array
        matrix : 2D symmetric numpy array
        1 or 2 dimensional numpy array, the n1 first indices in
        the zeroth axis of the array, should correspond to the
        values of the first group.
        The value of ``matrix[i][j]`` should correspond to
    paired : bool
        if paired=True, the data for the pairs are not taken into
        account. I.e. the diagonals and 'semi-diagonals' are not taken
        into account.

    Returns
    -------
    within_group_means - inter_group_mean : float
    """
    if paired:
        n1 = np.shape(mat)[0] / 2
    else:
        assert n1 is not None, ("give out n1 in sim_matrix_within_group_" +
                                "means_minus_inter_group_mean")
    within_group_mean = np.mean(
        sim_matrix_within_group_means(mat, n1)[0:2])
    if paired:
        inter_group_mean = paired_sim_matrix_inter_group_means(mat)[0]
    else:
        inter_group_mean = sim_matrix_mean_inter_group_similarity(mat, n1)
    return within_group_mean - inter_group_mean
