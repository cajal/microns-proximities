import numpy as np
from . import numpy_utils as npu


def convert_skeleton_to_nodes_edges(
    skeleton,
    verbose = False,):
    """
    from BCelli
    """
    
    all_skeleton_vertices = skeleton.reshape(-1,3)
    unique_rows,indices = np.unique(all_skeleton_vertices,return_inverse=True,axis=0)

    #need to merge unique indices so if within a certain range of each other then merge them together
    reshaped_indices = indices.reshape(-1,2)
    
    return unique_rows,reshaped_indices


def discretize_skeleton(skeleton, max_length, return_mapping=False):
    """
    Discretize the provided skeleton based on the specified maximum length.
    
    Parameters:
    - skeleton (ndarray): Array of shape (N, 2, 3) representing N edges with start and end points.
    - max_length (float): The maximum allowed length for any segment in the discretized skeleton.
    - return_mapping (bool): If True, a mapping array will be returned that would allow the original skeleton to be recovered. Default is False.
    
    Returns:
    - ndarray: Discretized skeleton based on the given parameters.
    """
    start_points, end_points = skeleton[:, 0], skeleton[:, 1]
    segment_diffs = end_points - start_points
    segment_lengths = np.linalg.norm(segment_diffs, axis=1)
    
    num_segments = np.ceil(segment_lengths / max_length).astype(int)
    num_segments[segment_lengths <= max_length] = 1

    if np.all(num_segments == 1):
        if return_mapping:
            return skeleton, np.arange(len(skeleton))
        return skeleton
    
    adjusted_diffs = segment_diffs / num_segments[:, None]
    segmented_diffs = np.repeat(adjusted_diffs, num_segments, axis=0)

    repeated_starts = np.repeat(start_points, num_segments, axis=0)
    increments = np.hstack([np.arange(n) for n in num_segments])
    interpolated_points = repeated_starts + increments[:, None] * segmented_diffs
    
    total_points = 0
    discretized_segments = []
    for idx, splits in enumerate(num_segments):
        segment_points = interpolated_points[total_points:total_points + splits]
        
        if splits == 1:
            segment_points = skeleton[idx]
        else:
            segment_points = np.vstack([segment_points, end_points[idx]])
        
        discretized_segments.append(segment_points)
        total_points += splits

    mapping = []
    for idx, splits in enumerate(num_segments):
        mapping.extend([idx] * splits)
    mapping = np.array(mapping)

    discretized_skeleton = np.vstack([np.array((segment[:-1], segment[1:])).transpose(1, 0, 2) for segment in discretized_segments])

    if return_mapping:
        return discretized_skeleton, mapping
    return discretized_skeleton


def discretized_to_original_skeleton(discretized_skeleton, mapping):
    """
    Convert a discretized skeleton back to its original form using the provided mapping.
    
    Parameters:
    - discretized_skeleton (ndarray): Array of shape (M, 2, 3) representing edges of a discretized skeleton.
    - mapping (ndarray): Mapping array of shape (M,) representing original edge index for each discretized edge.
    
    Returns:
    - ndarray: Original skeleton edges.
    """
    
    _, segment_counts = np.unique(mapping, return_counts=True)
    
    original_edges = []
    idx = 0
    for segment_count in segment_counts:
        start_point = discretized_skeleton[idx][0]
        end_point = discretized_skeleton[idx + segment_count - 1][1]
        idx += segment_count
        original_edges.append((start_point, end_point))
    
    return np.stack(original_edges)


def filter_edges(edges, vertices_inds_subset):
    """
    Find and return the subset of edges where at least one node is present in vertices_inds_subset.

    :param edges: List of lists or array of shape (N, 2) representing edges by their vertex indices.
    :param vertices_inds_subset: List or array of shape (M,) containing vertex indices to filter edges against.
    :return: Array of shape (P, 2) where P ≤ N, representing the filtered edges.
    """
    
    inds_set = set(vertices_inds_subset)
    filtered_edges = [edge for edge in edges if edge[0] in inds_set or edge[1] in inds_set]
    if not filtered_edges:
        return np.array([], dtype=int).reshape(0, 2)
    else:
        return np.stack(filtered_edges)



def compute_skeletal_length(vertices, edges, vertices_inds_subset=None):
    """
    Compute the total length of the skeletal structure given vertices and edges. 
    If a subset of vertex indices is provided, it will calculate the length only for 
    the edges related to that subset.

    :param vertices: N x 3 array-like object containing the coordinates of each vertex.
    :param edges: N x 2 array-like object of edge pairs, where each edge is represented by the indices of its endpoints.
    :param vertices_inds_subset: (Optional) A list of vertex indices for which to filter edges with. 
                                 If not provided, the function computes the skeletal length for all edges.
    :return: The total length of the skeletal structure (or the subset if vertices_inds_subset is provided).
    """
    vertices = np.array(vertices)
    edges = np.array(edges)
    
    if vertices_inds_subset is not None:
        edges = filter_edges(edges, vertices_inds_subset)

    assert vertices.shape[1] == 3, f"Invalid vertices shape: {vertices.shape}"
    assert edges.shape[1] == 2, f"Invalid edge shape: {edges.shape}"
    
    skeleton = vertices[edges]

    # Calculate the differences for each edge's endpoints
    diffs = skeleton[:, 1, :] - skeleton[:, 0, :]
    
    # Compute the norm of each difference (edge length)
    edge_lengths = np.linalg.norm(diffs, axis=1)
    
    return edge_lengths.sum()