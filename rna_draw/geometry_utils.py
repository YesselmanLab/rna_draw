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

def calculate_translation_vectors(junction, element, center_x, center_y, new_radius, xarray, yarray, get_first_nucleotides, get_last_nucleotides, get_junction_radius, is_parent=False):
    """
    Calculates the translation vectors (delta_x, delta_y) needed to reposition an element
    relative to its parent junction when the junction radius changes.

    Args:
    junction (object): The original junction object.
    element (object): The child element or parent element of the junction.
    center_x (float): The x-coordinate of the junction center.
    center_y (float): The y-coordinate of the junction center.
    new_radius (float): The new radius of the junction.
    xarray (list): The x-coordinates of all nodes.
    yarray (list): The y-coordinates of all nodes.
    get_first_nucleotides (function): Function to get the first nucleotides of a strand.
    get_last_nucleotides (function): Function to get the last nucleotides of a strand.
    get_junction_radius (function): Function to get the radius of a junction.
    is_parent (bool): Flag indicating whether the element is a parent of the junction.

    Returns:
    tuple: The translation vectors (delta_x, delta_y).
    """
    if is_parent:
        first_node, last_node = get_last_nucleotides(element.strands)
    else:
        first_node, last_node = get_first_nucleotides(element.strands)

    midpoint_x = (xarray[first_node] + xarray[last_node]) / 2
    midpoint_y = (yarray[first_node] + yarray[last_node]) / 2
    dx = midpoint_x - center_x
    dy = midpoint_y - center_y

    old_radius = get_junction_radius(junction)
    radius_difference = new_radius - old_radius

    unit_vector_x = dx / old_radius
    unit_vector_y = dy / old_radius

    delta_x = radius_difference * unit_vector_x
    delta_y = radius_difference * unit_vector_y

    return delta_x, delta_y

def translate_element_positions(xarray, yarray, junction, delta_x, delta_y, visited_junctions, visited_nodes):
    """
    Translates the positions of a junction and its elements (nodes and child junctions)
    by the given translation vectors (delta_x, delta_y).

    Args:
    xarray (list): The x-coordinates of all nodes.
    yarray (list): The y-coordinates of all nodes.
    junction (object): The junction object to translate.
    delta_x (float): The translation vector along the x-axis.
    delta_y (float): The translation vector along the y-axis.
    visited_junctions (set): Set of visited junctions to avoid processing the same junction multiple times.
    visited_nodes (set): Set of visited nodes to avoid processing the same node multiple times.

    Returns:
    None
    """
    if junction in visited_junctions:
        return

    visited_junctions.add(junction)

    # Translate node positions
    for node in junction.positions:
        if node not in visited_nodes:
            visited_nodes.add(node)
            xarray[node] += delta_x
            yarray[node] += delta_y

    # Translate junction center if it has center coordinates
    if hasattr(junction, "center_x"):
        junction.center_x += delta_x
        junction.center_y += delta_y

    # Recursively translate child junctions
    for child in junction.children:
        translate_element_positions(xarray, yarray, child, delta_x, delta_y, visited_junctions, visited_nodes)

    # Recursively translate parent junction if applicable
    if junction.has_parent():
        translate_element_positions(xarray, yarray, junction.parent, delta_x, delta_y, visited_junctions, visited_nodes)
