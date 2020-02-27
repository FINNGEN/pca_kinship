"""
This module contains functions and methods for permuting arrays.
Intended to use by the 'ptests' module, but may be suitable also
for other use.

Author: Rainer.Kujala@gmail.com
"""

import numpy as np
from numpy.random import RandomState


def get_random_state(seed=None):
    """
    Get a numpy.random.RandomState object with the given seed.
    If no seed is given or seed = None, a RandomState object with default
    numpy seeding method (from time) is returned.

    Parameters
    ----------
    seed : int
        Seed for the random number generator.

    Returns
    -------
    rng : numpy.random.RandomState
        A random number generator from numpy
    """
    if seed is not None:
        return RandomState(seed)
    else:
        return RandomState()


def get_permutation(paired, n1, n2, rng=None, i=None):
    """
    Helper function to compute a permutation with all the switches.

    Parameters
    ----------
    paired : bool
        True or False, whether to return a paired permutation
    n1 : int
        number of instances in the first group
    n2 : int
        number of instances in the second group
    rng : np.random.RandomState
        A random number generator (None for deterministic permutations)
    i : int
        The "index" of the permutation, when permSamples = 'all' is used

    Returns:
    --------
    The permutation as a 1D numpy array
    """
    if paired is True:
        assert n1 == n2, "With paired statistics setup, the number of" \
                         "elements in each group should be the same"
        if rng is None:  # go through full permutation dist.
            assert i is not None
            return _get_paired_permutation(n1, i)
        else:  # a paired setting but a random permutation
            return _get_random_paired_permutation(n1 + n2, rng)
    else:  # not a paired setting
        return rng.permutation(n1 + n2)


def _get_random_paired_permutation(n_tot, rng):
    """
    Get a random paired permutation.

    Parameters
    ----------
    n_tot : int
        the number of pairs * 2
    rng : numpy.random.RandomState
        the random number generator in use

    Returns
    -------
    A random permutation as a 1D numpy array.
    """
    n = n_tot / 2
    perm = np.arange(0, n_tot, dtype=np.int32)
    rands = rng.randint(0, 2, n)
    perm[:n] += rands * n
    perm[n:] -= rands * n
    return perm


def _get_paired_permutation(n, i):
    """Get a paired permutation corresponding to the index i

    Parameters
    ----------
    i : int
        index of the permutation (in range 0<= i < 2**n)
    n : int
        number of *pairs* /(half of the total number of elements
        to be permuted)
    """

    assert i >= 0
    assert i < 2 ** n
    permBase = np.zeros(n, dtype=np.int32)
    # compute the bin repr.
    for j in np.linspace(n - 1, 0, n).astype(np.int32):
        c = i // 2 ** j
        permBase[j] = c
        i -= 2 ** j * c
    perm = np.arange(0, 2 * n, dtype=np.int32)
    perm[:n] += permBase * n
    perm[n:] -= permBase * n
    return perm


def permute_array(data_array, perm):
    """
    Simply permute the array.

    Parameters
    ----------
    data_array : numpy array
        the numpy.ndarray to be permuted along the first axis
    perm : numpy array
        the permutation

    Returns
    -------
    data_array: numpy array
        A permuted view of the data_array
    """
    return data_array[perm]


def permute_matrix(matrix, perm):
    """
    Permute the rows and columns of the (similarity) matrix.
    (Destroys any 'grouping' effects)

    Parameters
    ----------
    matrix : 2D numpy array
    perm : 1D numpy array
        the permutation

    Returns
    -------
    A permuted view to the matrix.
    """
    assert matrix.shape[0] == matrix.shape[1], "should be a square matrix"
    return matrix[perm, :][:, perm]


def half_permute_paired_matrix(matrix, perm):
    """
    Permute only the latter half of the matrix, which
    destroys the paired structure of the matrix.
    The upper left corner of the matrix stays untouched.

    Parameters
    ----------
    matrix : numpy array
        the matrix to be permuted
    perm : numpy array
        The permutation (list of indices)

    Returns
    -------
    A view to the `half-permuted' matrix
    """
    assert len(matrix) % 2 == 0, "matrix rank should be paired"
    assert matrix.shape[0] == matrix.shape[1], "should be a square matrix"
    perm = np.array(perm)
    n1 = np.shape(matrix)[0] / 2
    perm = perm[perm >= n1]
    perm = np.hstack((np.arange(0, n1, dtype=np.int32), perm))
    return permute_matrix(matrix, perm)
