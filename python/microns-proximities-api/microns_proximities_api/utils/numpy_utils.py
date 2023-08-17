import numpy as np


def intersect2d(ar1, ar2, assume_unique=False, return_indices=False):
    """
    Use numpy intersect1d with 2-dimensional arrays
    
    https://stackoverflow.com/a/8317403
    """
    def get_dtype_for_view(arr):
        _, ncols = arr.shape
        return {
            'names':['f{}'.format(i) for i in range(ncols)],
            'formats': ncols * [arr.dtype]
        }
    
    arr1 = np.array(ar1)
    arr2 = np.array(ar2)
    
    return np.intersect1d(
        arr1.view(get_dtype_for_view(arr1)), 
        arr2.view(get_dtype_for_view(arr2)), 
        assume_unique=assume_unique, 
        return_indices=return_indices
    )


def unique_row_view(data, return_inverse=False):
    """
    Faster implementation of np.unique with axis=0 argument.

    Credit to nschloe:
    https://github.com/numpy/numpy/issues/11136
    """
    if data.shape[0] == 0:
        return data

    b = np.ascontiguousarray(data).view(
        np.dtype((np.void, data.dtype.itemsize * data.shape[1]))
        )
    
    u = np.unique(b, return_inverse=return_inverse)
    if return_inverse:
        u, inverse_idx = u
        return u.view(data.dtype).reshape(-1, data.shape[1]), inverse_idx
    else:
        return u.view(data.dtype).reshape(-1, data.shape[1])