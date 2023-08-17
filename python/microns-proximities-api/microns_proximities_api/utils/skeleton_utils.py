import numpy as np
from . import numpy_utils as npu


def discretize_skeleton(full_edges, maximum_length, return_as_int=True, reshape=True):
    """
    from Christos
    """
    if full_edges.shape[0] == 0:
        return full_edges

    p0s = full_edges[:, 0]
    p1s = full_edges[:, 1]

    diffs = p1s - p0s
    distances = np.linalg.norm(diffs, axis=1)
    inc_nums = np.ceil(distances / maximum_length).astype(int)
    inc_nums[inc_nums<2] = 2
    diffs_inc = np.repeat(diffs / inc_nums[:, None], inc_nums, axis=0)

    p0s_stack = np.repeat(p0s, inc_nums, axis=0)
    max_arange = np.arange(inc_nums.max())
    multiplicative_incrementer = np.hstack([max_arange[0:i] for i in inc_nums.tolist()])
    evenly_spaced = p0s_stack + (multiplicative_incrementer[:, None] * diffs_inc)

    total = 0
    incremented_edges = list()
    for increment, p1 in zip(inc_nums, p1s):
        temp_total = total + increment
        inc_edge = evenly_spaced[total:temp_total]
        inc_edge = np.concatenate((inc_edge, p1[None]))
        incremented_edges.append(inc_edge)
        total = temp_total
    new_full_edges = np.vstack([np.array((inc_edge[:-1], inc_edge[1:])).transpose(1, 0, 2) for inc_edge in incremented_edges])
    
    if return_as_int:
        new_full_edges = np.round(new_full_edges).astype(int)
    
    if not reshape:
        return new_full_edges
    else:
        return npu.unique_row_view(new_full_edges.reshape(-1, 3))