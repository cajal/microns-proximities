import itertools
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from microns_utils.cluster_utils import cluster_point_cloud
from collections import Counter

from ..utils import numpy_utils as npu, skeleton_utils as sku

values = (-1, 0, 1)
neighbor_offsets = np.array(list(itertools.product(values, values, values)))

vx_to_grid = lambda pts, res: np.floor(pts / res).astype(int)


def create_mapping_from_unique_to_original(unique_inv):
    mapping = {}
    for original_idx, unique_idx in enumerate(unique_inv):
        if unique_idx not in mapping:
            mapping[unique_idx] = []
        mapping[unique_idx].append(original_idx)
    return mapping


def is_adjacent(pt1, pt2):
    """
    Determines if two points are adjacent in a 3D grid.

    Parameters:
    pt1 (tuple): A tuple of three integers representing the coordinates of the first point.
    pt2 (tuple): A tuple of three integers representing the coordinates of the second point.

    Returns:
    bool: True if the points are adjacent, False otherwise.
    """
    diff = np.array(pt1) - np.array(pt2)
    return np.all(np.abs(diff) <= 1) and not np.all(diff == 0)


def get_adjacent_voxels(voxel_list, test_voxel):
    """
    Finds all voxels in a list that are adjacent to a specified voxel.

    Parameters:
    voxel_list (list): A list of tuples, each representing the coordinates of a voxel.
    test_voxel (tuple): A tuple of three integers representing the coordinates of the voxel to test.

    Returns:
    list: A list of tuples representing the coordinates of all voxels adjacent to the test voxel.
    """
    matches = []
    indices_to_remove = []
    for i, vx in enumerate(voxel_list):
        if is_adjacent(vx, test_voxel):
            matches.append(vx)
            indices_to_remove.append(i)
    # Remove matched voxels from the list
    for i in sorted(indices_to_remove, reverse=True):
        voxel_list.pop(i)
    return matches


def group_adjacent_voxels(voxel_list):
    """
    Groups adjacent voxels into sets and returns a dictionary that maps voxel coordinates to group indices.

    Parameters:
    voxel_list (list): A list of tuples, each representing the coordinates of a voxel.

    Returns:
    dict: A dictionary that maps tuples representing voxel coordinates to integers representing group indices.
    """
    assert isinstance(voxel_list, list), 'voxel_list must be a list'
    voxel_list = voxel_list.copy()
    
    # make sets of adjacent voxels
    voxel_sets = []
    while len(voxel_list) > 0:
        vx_set = [voxel_list.pop()]
        i = 0
        while i < len(vx_set):
            vx = vx_set[i]
            adjacent_voxels = get_adjacent_voxels(voxel_list, vx)
            vx_set.extend(adjacent_voxels)
            i += 1
        voxel_sets.append(vx_set)

    # make voxel group lookup dict
    vx_lookup = {}
    for i, vx_set in enumerate(voxel_sets):
        for vx in vx_set:
            vx_lookup[tuple(vx)] = i
    return vx_lookup


def compute_proximities_voxel_set_intersection(verts1, verts2, resolution):
    """
    Calculates the intersection proximity between two sets of voxels.

    Parameters:
    - verts1 (np.array): An array of vertex coordinates for the first voxel set.
    - verts2 (np.array): An array of vertex coordinates for the second voxel set.
    - resolution (int/float): The grid resolution for voxelization.

    Returns:
    dict: A dictionary containing various data about the intersection, including
    - grids for both voxel sets
    - unique grids for both voxel sets
    - indices for intersecting unique grids for both voxel sets
    - proximal vertices for both voxel sets
    - proximal vertex indices for both voxel sets
    - distances of the proximities
    """
    verts1_grid = vx_to_grid(verts1, resolution)
    verts2_grid = vx_to_grid(verts2, resolution)
    
    verts1_grid_un, verts1_grid_un_inv = npu.unique_row_view(verts1_grid, return_inverse=True)
    verts2_grid_un, verts2_grid_un_inv = npu.unique_row_view(verts2_grid, return_inverse=True)
    
    # test 27 intersections between grids

    verts1_grid_un_inds_prx = []
    verts2_grid_un_inds_prx = []
    verts1_inds_prx = []
    verts2_inds_prx = []
    dists_prx = []
    for offset in neighbor_offsets:
        vx, verts1_grid_un_inds_sub, verts2_grid_un_inds_sub = npu.intersect2d(verts1_grid_un + offset, verts2_grid_un, assume_unique=True, return_indices=True)
        # TODO : consider concatenating grids
        # Calculate per grid-grid 
        for i, j in zip(verts1_grid_un_inds_sub, verts2_grid_un_inds_sub):
            verts1_inds_sub = np.where(verts1_grid_un_inv == i)[0]
            verts2_inds_sub = np.where(verts2_grid_un_inv == j)[0]
            
            verts1_sub = verts1[verts1_inds_sub]
            verts2_sub = verts2[verts2_inds_sub]
            
            dm = cdist(verts1_sub, verts2_sub)
            verts1_inds_sub_prx, verts2_inds_sub_prx = np.stack(np.where(dm <= resolution), -1).T
            
            if len(verts1_inds_sub_prx) > 0:
                for ind1, ind2, dist in zip(
                    verts1_inds_sub[verts1_inds_sub_prx], 
                    verts2_inds_sub[verts2_inds_sub_prx],
                    dm[verts1_inds_sub_prx, verts2_inds_sub_prx].round(decimals=0).astype(int)
                ):
                    verts1_grid_un_inds_prx.append(i)
                    verts2_grid_un_inds_prx.append(j)
                    verts1_inds_prx.append(ind1)
                    verts2_inds_prx.append(ind2)
                    dists_prx.append(dist)
    if len(dists_prx) == 0:
        return {}
    return {
        'verts1_prx': verts1[verts1_inds_prx],
        'verts2_prx': verts2[verts2_inds_prx],
        'verts1_inds_prx': verts1_inds_prx,
        'verts2_inds_prx': verts2_inds_prx,
        'dists_prx': np.array(dists_prx),
        'verts1_grid_un': verts1_grid_un,
        'verts2_grid_un': verts2_grid_un,
        'verts1_grid_un_inds_prx': verts1_grid_un_inds_prx,
        'verts2_grid_un_inds_prx': verts2_grid_un_inds_prx,
    }


def cluster_proximities_adjacent_voxels(verts1_prx, verts2_prx, verts1_inds_prx, verts2_inds_prx, dists_prx, verts1_grid_un, verts2_grid_un, verts1_grid_un_inds_prx, verts2_grid_un_inds_prx):
    """
    TODO fix compatibility with compute_proximity outputs

    Clusters proximities between adjacent voxels.

    Parameters:
    - verts1_prx (np.array): Proximal vertices from the first voxel set.
    - verts2_prx (np.array): Proximal vertices from the second voxel set.
    - verts1_inds_prx (list[int]): Indices of the proximal vertices from the first voxel set.
    - verts2_inds_prx (list[int]): Indices of the proximal vertices from the second voxel set.
    - dists_prx (list[int/float]): Distances between the proximal vertices.
    - verts1_grid_un (np.array): Unique grid values for the first voxel set.
    - verts2_grid_un (np.array): Unique grid values for the second voxel set.
    - verts1_grid_un_inds_prx (list[int]): Indices for intersecting unique grids from the first voxel set.
    - verts2_grid_un_inds_prx (list[int]): Indices for intersecting unique grids from the second voxel set.

    Returns:
    list[dict]: A list of dictionaries, where each dictionary represents a group of proximal voxels. Each dictionary contains:
    - Proximity ID
    - Number of voxels in the group
    - Proximal vertices for both voxel sets
    - Proximal vertex indices for both voxel sets
    - Distances of the proximities
    """
    
    verts1_voxels = verts1_grid_un[verts1_grid_un_inds_prx]
    verts2_voxels = verts2_grid_un[verts2_grid_un_inds_prx]
    all_voxels = npu.unique_row_view(np.vstack([verts1_voxels, verts2_voxels]))
    vx_group_lookup = group_adjacent_voxels(all_voxels.tolist())
    
    verts1_group = []
    for vx in verts1_voxels:
        verts1_group.append(vx_group_lookup[tuple(vx)])
    
    verts2_group = []
    for vx in verts2_voxels:
        verts2_group.append(vx_group_lookup[tuple(vx)])
    
    assert np.array_equal(verts1_group, verts2_group), 'sanity check failed. voxel groups 1 and 2 are not equal.'
    
    voxel_groups = verts1_group
    group_counts = Counter(voxel_groups)
    
    results = []
    for group_id in range(np.max(voxel_groups) + 1):
        group_result = {
            'prx_grp_id': group_id
        }
        group_verts1 = []
        group_verts2 = []
        group_inds1 = []
        group_inds2 = []
        group_dists = []
        for vert1, vert2, ind1, ind2, dist, vx_group in zip(verts1_prx, verts2_prx, verts1_inds_prx, verts2_inds_prx, dists_prx, voxel_groups):
            if group_id == vx_group:
                group_verts1.append(vert1)
                group_verts2.append(vert2)
                group_inds1.append(ind1)
                group_inds2.append(ind2)
                group_dists.append(dist)
        
        group_result['n_voxels'] = group_counts[group_id]
        group_result['verts1_prx'] = np.stack(group_verts1)
        group_result['verts2_prx'] = np.stack(group_verts2)
        group_result['verts1_inds_prx'] = group_inds1
        group_result['verts2_inds_prx'] = group_inds2
        group_result['dists_prx'] = group_dists
        results.append(group_result)

    return results


def compute_proximities_kdtree(verts1, verts2, radius):
    """
    """
    # initialize result variables
    verts1_prx = []
    verts2_prx = []
    verts1_inds_prx = []
    verts2_inds_prx = []
    dists_prx = []

    # compute proximities
    tree1 = KDTree(verts1)
    tree2 = KDTree(verts2)
    indices2 = tree1.query_ball_tree(tree2, r=radius)

    # extract proximities
    for ind1, inds2_group in enumerate(indices2):
        if len(inds2_group) > 0:
            for ind2 in inds2_group:
                verts1_inds_prx.append(ind1)
                verts2_inds_prx.append(ind2)

    if len(verts1_inds_prx) == 0:
        return {}

    verts1_inds_prx = np.array(verts1_inds_prx)
    verts2_inds_prx = np.array(verts2_inds_prx)
    verts1_prx = verts1[verts1_inds_prx]
    verts2_prx = verts2[verts2_inds_prx]
    dists_prx = np.linalg.norm(verts1_prx - verts2_prx, axis=1).round().astype(int)

    return {
            'verts1_prx': verts1_prx,
            'verts2_prx': verts2_prx,
            'verts1_inds_prx': verts1_inds_prx,
            'verts2_inds_prx': verts2_inds_prx,
            'dists_prx': dists_prx,
    }


def cluster_proximities_point_cloud(verts1_prx, verts2_prx, verts1_inds_prx, verts2_inds_prx, dists_prx, eps, algorithm='ball_tree', min_samples=2, return_indices=True):
    """
    TODO fix compatibility with compute_proximity outputs
    """
    verts1_prx_un, verts1_prx_un_inv = npu.unique_row_view(verts1_prx, return_inverse=True)
    verts2_prx_un, verts2_prx_un_inv = npu.unique_row_view(verts2_prx, return_inverse=True)

    verts1_prx_un_to_original_map = create_mapping_from_unique_to_original(verts1_prx_un_inv)
    
    # cluster the point cloud
    _, cluster_indices = cluster_point_cloud(
        np.vstack([verts1_prx_un, verts2_prx_un]), 
        eps=eps, 
        algorithm=algorithm, 
        min_samples=min_samples, 
        return_indices=return_indices
    )
    
    # format results
    results = []
    for i, cluster_group in enumerate(cluster_indices):
        verts1_prx_cluster_inds = []
        for idx in cluster_group:
            if idx < len(verts1_prx_un):
                verts1_prx_cluster_inds.extend(verts1_prx_un_to_original_map[idx])

        group_result = {
            'prx_grp_id': i,
            'verts1_prx': verts1_prx[verts1_prx_cluster_inds],
            'verts2_prx': verts2_prx[verts1_prx_cluster_inds],
            'verts1_inds_prx': np.array(verts1_inds_prx)[verts1_prx_cluster_inds],
            'verts2_inds_prx': np.array(verts2_inds_prx)[verts1_prx_cluster_inds],
            'dists_prx': np.array(dists_prx)[verts1_prx_cluster_inds],
            
        }
        results.append(group_result)

    return results