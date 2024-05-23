import numpy as np

def rotate_point(point, origin, angle):
    """
    Rotates a point around a given origin by a specified angle.

    Args:
    point (tuple): The (x, y) coordinates of the point to rotate.
    origin (tuple): The (x, y) coordinates of the origin around which to rotate.
    angle (float): The angle in radians to rotate the point.

    Returns:
    tuple: The (x, y) coordinates of the rotated point.
    """
    px, py = point
    ox, oy = origin

    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)], 
                                [np.sin(angle),  np.cos(angle)]])
    
    result = np.dot(rotation_matrix, np.array([px - ox, py - oy])) + np.array([ox, oy])
    
    return result[0], result[1]

def calculate_deltas(original_junction, subpart, center_x, center_y, new_radius, xarray, yarray, get_first_nucleotides, get_last_nucleotides, get_junction_radius, parent=False):
    first_node, last_node = (get_last_nucleotides(subpart.strands) if parent else get_first_nucleotides(subpart.strands))
    midpoint_x = (xarray[first_node] + xarray[last_node]) / 2
    midpoint_y = (yarray[first_node] + yarray[last_node]) / 2
    dx = midpoint_x - center_x
    dy = midpoint_y - center_y

    old_radius = get_junction_radius(original_junction)
    radius_difference = new_radius - old_radius

    unit_vector_x = dx / old_radius
    unit_vector_y = dy / old_radius

    delta_x = radius_difference * unit_vector_x
    delta_y = radius_difference * unit_vector_y

    return delta_x, delta_y

def update_positions(xarray, yarray, junction, delta_x, delta_y, visited_junctions, visited_nodes):
    if junction in visited_junctions:
        return

    visited_junctions.add(junction)

    for node in junction.positions:
        if node not in visited_nodes:
            visited_nodes.add(node)
            xarray[node] += delta_x
            yarray[node] += delta_y

    if junction.is_junction:
        if hasattr(junction, "center_x"):
            junction.center_x += delta_x
            junction.center_y += delta_y

    for child in junction.children:
        update_positions(xarray, yarray, child, delta_x, delta_y, visited_junctions, visited_nodes)

    if junction.has_parent():
        update_positions(xarray, yarray, junction.parent, delta_x, delta_y, visited_junctions, visited_nodes)
