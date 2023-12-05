import gc

import numpy as np

from scipy.sparse import csr_matrix, csc_matrix

# code in this file is a copy of code in scipy.sparse._construct
# somewhat memoy optimized for concatenating csr_matrices
def get_index_dtype(arrays=(), maxval=None, check_contents=False):
    """
    Based on input (integer) arrays `a`, determine a suitable index data
    type that can hold the data in the arrays.

    Parameters
    ----------
    arrays : tuple of array_like
        Input arrays whose types/contents to check
    maxval : float, optional
        Maximum value needed
    check_contents : bool, optional
        Whether to check the values in the arrays and not just their types.
        Default: False (check only the types)

    Returns
    -------
    dtype : dtype
        Suitable index data type (int32 or int64)

    """

    int32min = np.int32(np.iinfo(np.int32).min)
    int32max = np.int32(np.iinfo(np.int32).max)

    # not using intc directly due to misinteractions with pythran
    dtype = np.int32 if np.intc().itemsize == 4 else np.int64
    if maxval is not None:
        maxval = np.int64(maxval)
        if maxval > int32max:
            dtype = np.int64

    if isinstance(arrays, np.ndarray):
        arrays = (arrays,)

    for arr in arrays:
        arr = np.asarray(arr)
        if not np.can_cast(arr.dtype, np.int32):
            if check_contents:
                if arr.size == 0:
                    # a bigger type not needed
                    continue
                elif np.issubdtype(arr.dtype, np.integer):
                    maxval = arr.max()
                    minval = arr.min()
                    if minval >= int32min and maxval <= int32max:
                        # a bigger type not needed
                        continue

            dtype = np.int64
            break

    return dtype


def _compressed_sparse_stack(blocks, axis):
    """
    Stacking fast path for CSR/CSC matrices
    (i) vstack for CSR, (ii) hstack for CSC.
    """
    other_axis = 1 if axis == 0 else 0
    data = np.concatenate([b.data for b in blocks])
    constant_dim = blocks[0].shape[other_axis]
    idx_dtype = get_index_dtype(arrays=[b.indptr for b in blocks],
                                maxval=max(data.size, constant_dim))
    indices = np.empty(data.size, dtype=idx_dtype)
    indptr = np.empty(sum(b.shape[axis] for b in blocks) + 1, dtype=idx_dtype)
    last_indptr = idx_dtype(0)
    sum_dim = 0
    sum_indices = 0
    for b in blocks:
        if b.shape[other_axis] != constant_dim:
            raise ValueError(f'incompatible dimensions for axis {other_axis}')
        indices[sum_indices:sum_indices+b.indices.size] = b.indices
        sum_indices += b.indices.size
        idxs = slice(sum_dim, sum_dim + b.shape[axis])
        indptr[idxs] = b.indptr[:-1]
        indptr[idxs] += last_indptr
        sum_dim += b.shape[axis]
        last_indptr += b.indptr[-1]
    
    del blocks
    gc.collect()

    indptr[-1] = last_indptr
    if axis == 0:
        return csr_matrix((data, indices, indptr),
                          shape=(sum_dim, constant_dim))
    else:
        return csc_matrix((data, indices, indptr),
                          shape=(constant_dim, sum_dim))
    

def bmat(matrices, dtype=None):
    # stack along rows (axis 0):
    A = _compressed_sparse_stack(matrices, 0)
    if dtype is not None:
        A = A.astype(dtype)
    return A
