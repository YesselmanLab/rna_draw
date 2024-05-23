import numpy as np
import math
import random
from itertools import combinations
from scipy.optimize import curve_fit
import pandas as pd

def estimate_inches_using_curve_fit(x_range, y_range):
    """
    Estimate the required inches for the x and y dimensions using curve fitting based on RNA area and figure size data.
    
    Args:
    x_range (float): The range of x coordinates.
    y_range (float): The range of y coordinates.
    
    Returns:
    tuple: Estimated inches for x and y dimensions.
    """
    # Data for RNA identity, area, and figure size variables
    data = {
        "rna_identity": ["hairpin", "t-RNA", "CO-VID19 5' UTR", "50S Ribosome"],
        "rna_area": [3781, 126207, 1150472, 4286761],
        "rna_figsize_variable": [25, 30, 35, 40],
    }
    df = pd.DataFrame(data)
    
    x = df["rna_area"]
    y = df["rna_figsize_variable"]
    
    # Define the function for curve fitting
    def test(x, a, b, c):
        return a * (x - b) ** c

    # Perform curve fitting
    param, param_cov = curve_fit(test, x, y)
    
    area = x_range * y_range
    estimated_inches_x = x_range / (param[0] * (area - param[1]) ** param[2])
    estimated_inches_y = y_range / (param[0] * (area - param[1]) ** param[2])
    
    return estimated_inches_x, estimated_inches_y

def get_between_strand(junction, between_strands, xarray):
    """
    Get the strand positions between the given junction.
    
    Args:
    junction (object): The current junction object.
    between_strands (list): List of tuples indicating positions of strands between junctions.
    xarray (list): List of x coordinates for the structure.
    
    Returns:
    list: Start and end positions of the strand.
    """
    if between_strands:
        for between_strand in between_strands:
            if between_strand[0] <= junction.positions[0] <= between_strand[1] and between_strand[0] <= junction.positions[1] <= between_strand[1]:
                return between_strand
    return [0, len(xarray)]

def explore_paths(renderer):
    """
    Explore possible paths to optimize the structure by minimizing overlap.
    
    Args:
    renderer (object): The renderer object containing the RNA structure and methods.
    
    Returns:
    int: The overlap count after exploring paths.
    """
    total_acc = 0
    total_tried = 0
    between_strands = renderer.between_strands

    for index, junction in enumerate(renderer.struct.get_junctions()):
        # Set possible angles for branches based on the number of children
        if len(junction.children) == 1:
            angles = [180, 270, 90]
        elif len(junction.children) <= 3:
            angles = [270, 180, 90]
        else:
            angles = [315, 270, 225, 180, 135, 90, 45]
        
        if renderer.global_best_combo.get(junction) is None:
            renderer.global_best_combo[junction] = None
        if renderer.global_best_overlap.get(junction) is None:
            renderer.global_best_overlap[junction] = None

        best_overlap = None
        best_combo = None

        # Try all combinations of angles for the children
        for comb in combinations(angles, len(junction.children)):
            for child_index, child in enumerate(junction.children):
                renderer.set_branch_angle(junction, child, comb[child_index])

            renderer.update_unpaired_strands_positions(junction)

            overlap_count, nodes_below_straight_strand = None, 0
            between = get_between_strand(junction, between_strands, renderer.xarray)

            overlap_count = renderer.update_overlap_count(node_array=list(range(between[0], between[1])), only_inside=True)

            # Count nodes below the straight strand
            for node_pos in renderer.yarray:
                if node_pos < renderer.yarray[0] + renderer.NODE_R:
                    nodes_below_straight_strand += 1
            
            # Update global best overlap and combination
            if renderer.global_best_overlap.get(junction) is None or overlap_count < renderer.global_best_overlap[junction]:
                renderer.global_best_overlap[junction] = overlap_count
                renderer.global_best_combo[junction] = comb

            # Update best local overlap and combination
            if best_overlap is None or overlap_count < best_overlap:
                best_overlap = overlap_count
                best_combo = comb
            elif best_overlap is not None and overlap_count - best_overlap > 0:
                # Simulated annealing approach
                T = np.sqrt(len(renderer.xarray)) / 8
                diff = overlap_count - best_overlap
                val = np.exp(-diff / T)
                r = random.uniform(0, 1)
                if val > r:
                    best_overlap = overlap_count
                    best_combo = comb
                    total_acc += 1
                total_tried += 1

        # Set the best angles found for the children
        for child, angle in zip(junction.children, renderer.global_best_combo[junction]):
            renderer.set_branch_angle(junction, child, angle)
            renderer.update_unpaired_strands_positions(junction)

    renderer.add_horizontal_distance_between(between_strands)

    overlap_count = renderer.update_overlap_count()
    return overlap_count

def straighten_branches(renderer, global_best_overlap):
    """
    Straighten branches to minimize overlap in the RNA structure.
    Note: This is due to cases in which 90 degree angles being prioritized in the RNA structure. If the algorithm is setup more effectively, this method won't be needed.
    
    Args:
    renderer (object): The renderer object containing the RNA structure and methods.
    global_best_overlap (int): The best overlap count found globally.
    """
    for junction in renderer.struct.get_junctions():
        if len(junction.children) == 1:
            parent_pos_x, parent_pos_y, parent_angle_deg = renderer.get_junction_parent_data(junction, junction.center_x, junction.center_y)
            # Note get_junction_branch_angle error. It is not always perfectly accurate. 
            original_angle = renderer.get_junction_branch_angle(junction, junction.children[0], parent_pos_x, parent_pos_y, parent_angle_deg)
            if original_angle > 185:
                original_angle = 270
            elif original_angle < 175:
                original_angle = 90
            else:
                original_angle = 180

            renderer.set_branch_angle(junction, junction.children[0], 180)
            renderer.update_unpaired_strands_positions(junction)

            between = get_between_strand(junction, renderer.between_strands, renderer.xarray)
            new_overlap_count = renderer.update_overlap_count(node_array=list(range(between[0], between[1])), only_inside=True)

            if new_overlap_count > global_best_overlap:
                renderer.set_branch_angle(junction, junction.children[0], original_angle)
                renderer.update_unpaired_strands_positions(junction)
