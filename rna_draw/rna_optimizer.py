import numpy as np
import math
import random
from itertools import combinations
from scipy.optimize import curve_fit
import pandas as pd
from rna_draw.geometry_utils import rotate_point, calculate_translation_vectors, translate_element_positions
from rna_draw.chunk_controller import ChunkContainer

class RNAOptimizer:
    def __init__(self, renderer):
        self.renderer = renderer

    def get_first_nucleotides(self, strands):
        """
        Returns the first nucleotide pair in the junction used to determine angle.

        Args:
        strands (list): List of strands in the junction.

        Returns:
        tuple: First and last nucleotide in the junction.
        """
        first_nucleotide = strands[0][0]
        last_nucleotide = strands[-1][-1]
        return first_nucleotide, last_nucleotide

    def get_last_nucleotides(self, strands):
        """
        Returns the last nucleotide pair in the junction.

        Args:
        strands (list): List of strands in the junction.

        Returns:
        tuple: Last nucleotide pair in the junction.
        """
        first_nucleotide = strands[0][-1]
        last_nucleotide = strands[-1][0]
        return first_nucleotide, last_nucleotide

    def get_junction_center(self, junction):
        """
        Returns the center x and y coordinates of the junction.

        Args:
        junction (object): The junction object.

        Returns:
        tuple: Center x and y coordinates of the junction.
        """
        if hasattr(junction, 'center_x') and hasattr(junction, 'center_y'):
            return junction.center_x, junction.center_y
        else:
            return self.update_junction_center(junction)

    def get_junction_radius(self, junction):
        """
        Returns the radius of the junction.

        Args:
        junction (object): The junction object.

        Returns:
        float: Radius of the junction.
        """
        if hasattr(junction, "radius"):
            return junction.radius
        else:
            junction.radius = self.calculate_junction_radius(junction)
            return junction.radius

    def calculate_junction_radius(self, junction):
        """
        Calculates the radius of the junction.

        Args:
        junction (object): The junction object.

        Returns:
        float: Calculated radius of the junction.
        """
        center_x, center_y = self.get_junction_center(junction)
        if junction.has_parent():
            parent_strand = junction.parent
        elif junction.has_children():
            parent_strand = junction.children[0]
        else:
            print("Junction has no parent or child.")
            return None

        first, last = self.get_last_nucleotides(parent_strand.strands)
        mid_x = (self.renderer.xarray[first] + self.renderer.xarray[last]) / 2
        mid_y = (self.renderer.yarray[first] + self.renderer.yarray[last]) / 2
        radius = np.sqrt((mid_x - center_x) ** 2 + (mid_y - center_y) ** 2)
        return radius

    def update_junction_center(self, junction, inversion_check=False):
        """
        Calculates the junction center and updates its coordinates.

        Args:
        junction (object): The junction object.
        inversion_check (bool): Flag to check for inversion.

        Returns:
        tuple: Updated center x and y coordinates of the junction.
        """
        if hasattr(junction, "radius") and junction.has_parent():
            parent_strand = junction.parent
            first, last = self.get_last_nucleotides(parent_strand.strands)
            mid_x = (self.renderer.xarray[first] + self.renderer.xarray[last]) / 2
            mid_y = (self.renderer.yarray[first] + self.renderer.yarray[last]) / 2

            direction_x = self.renderer.xarray[last] - self.renderer.xarray[first]
            direction_y = self.renderer.yarray[last] - self.renderer.yarray[first]

            magnitude = np.sqrt(direction_x**2 + direction_y**2)
            unit_x = direction_x / magnitude
            unit_y = direction_y / magnitude

            perp_x = -unit_y
            perp_y = unit_x

            centroid_x = sum(self.renderer.xarray[node] for node in junction.positions) / len(junction.positions)
            centroid_y = sum(self.renderer.yarray[node] for node in junction.positions) / len(junction.positions)

            direction_to_centroid_x = centroid_x - self.renderer.xarray[first]
            direction_to_centroid_y = centroid_y - self.renderer.yarray[first]

            dot_product = perp_x * direction_to_centroid_x + perp_y * direction_to_centroid_y

            if dot_product < 0 and inversion_check == True:
                perp_x, perp_y = -perp_x, -perp_y

            center_x = mid_x + junction.radius * perp_x
            center_y = mid_y + junction.radius * perp_y
        else:
            x_coords = [self.renderer.xarray[value] for value in junction.positions]
            y_coords = [self.renderer.yarray[value] for value in junction.positions]
            center_x = sum(x_coords) / len(x_coords)
            center_y = sum(y_coords) / len(y_coords)

        junction.center_x = center_x
        junction.center_y = center_y
        return center_x, center_y

    def get_junction_parent_data(self, junction, center_x, center_y):
        """
        Returns the parent data of the junction, including the coordinates and angle.

        Args:
        junction (object): The junction object.
        center_x (float): Center x-coordinate of the junction.
        center_y (float): Center y-coordinate of the junction.

        Returns:
        tuple: Parent x and y coordinates, and the angle in degrees.
        """
        if junction.has_parent():
            left, right = self.get_last_nucleotides(junction.parent.strands)
            parent_pos_x = (self.renderer.xarray[left] + self.renderer.xarray[right]) / 2
            parent_pos_y = (self.renderer.yarray[left] + self.renderer.yarray[right]) / 2

            parent_angle_rad = np.arctan2(parent_pos_y - center_y, parent_pos_x - center_x)
            parent_angle_deg = np.rad2deg(parent_angle_rad)

            return parent_pos_x, parent_pos_y, parent_angle_deg

    def get_junction_branch_angle(self, junction, child, center_x, center_y, parent_angle_deg):
        """
        Calculates the angle of a single branch with respect to the parent helix.

        Args:
        junction (object): The junction object.
        child (object): The child node.
        center_x (float): Center x-coordinate of the junction.
        center_y (float): Center y-coordinate of the junction.
        parent_angle_deg (float): Angle of the parent in degrees.

        Returns:
        float: Angle between the parent and the child in degrees.
        """
        data_new = self.renderer.struct[child.m_id]
        left, right = self.get_first_nucleotides(data_new.strands)
        child_pos_x = (self.renderer.xarray[left] + self.renderer.xarray[right]) / 2
        child_pos_y = (self.renderer.yarray[left] + self.renderer.yarray[right]) / 2

        child_angle_rad = math.atan2(child_pos_y - center_y, child_pos_x - center_x)
        child_angle_deg = math.degrees(child_angle_rad)

        angle_between = (child_angle_deg - parent_angle_deg + 360) % 360
        return angle_between

    def set_branch_angle(self, junction, child, target_angle):
        """
        Sets the angle of a child branch relative to its parent junction to a target angle. (Counter-Clockwise)

        Args:
        junction (object): The junction object which serves as the reference point.
        child (object): The child branch to be rotated.
        target_angle (float): The target angle in degrees to set for the child branch.

        Returns:
        None
        """
        center_x, center_y = self.get_junction_center(junction)
        x, y, parent_angle_deg = self.get_junction_parent_data(junction, center_x, center_y)

        current_angle = self.get_junction_branch_angle(junction, child, center_x, center_y, parent_angle_deg)
        angle_diff = math.radians(target_angle) - math.radians(current_angle)

        queue = [child]
        processed_nucleotides = set()

        # Process each branch in the queue
        while queue:
            current_branch = queue.pop(0)
            data_new = self.renderer.struct[current_branch.m_id]

            # Rotate each nucleotide in the current branch
            for strand in data_new.strands:
                for nucleotide in strand:
                    if nucleotide not in processed_nucleotides:
                        new_x, new_y = rotate_point((self.renderer.xarray[nucleotide], self.renderer.yarray[nucleotide]), (center_x, center_y), angle_diff)
                        self.renderer.xarray[nucleotide] = new_x
                        self.renderer.yarray[nucleotide] = new_y
                        processed_nucleotides.add(nucleotide)

            # Rotate the junction center if the current branch is a junction.
            if current_branch.is_junction():
                if hasattr(current_branch, 'center_x') and hasattr(current_branch, 'center_y'):
                    new_center_x, new_center_y = rotate_point((current_branch.center_x, current_branch.center_y), (center_x, center_y), angle_diff)
                    current_branch.center_x = new_center_x
                    current_branch.center_y = new_center_y

            # Add child branches of the current branch to the queue
            if current_branch.has_children():
                for child_branch in current_branch.children:
                    queue.append(child_branch)

    def update_unpaired_strands_positions(self, junction):
        """
        Places the unpaired nucleotides in a circle around the junction center,
        spaced evenly. Increases the radius if necessary to fit the nucleotides.

        Args:
        junction (object): The junction object containing unpaired strands.

        Returns:
        None
        """
        center_x, center_y = self.get_junction_center(junction)
        radius = self.get_junction_radius(junction)

        _, _, parent_angle_deg = self.get_junction_parent_data(junction, center_x, center_y)

        radius_updated = False

        # Process each strand in the junction
        for strand in junction.strands:
            num_nodes = len(strand[1:-1])
            space_needed = (num_nodes) * self.renderer.NODE_R * 2
            arc_length = -math.inf

            shifted_start_rad = None
            shifted_end_rad = None
            angle_shift_rad = None

            # Adjust the radius until the arc length can accommodate the nucleotides
            while arc_length < space_needed and space_needed != 0:
                shifted_start_rad = np.arctan2(self.renderer.yarray[strand[0]] - center_y, self.renderer.xarray[strand[0]] - center_x)
                shifted_end_rad = np.arctan2(self.renderer.yarray[strand[-1]] - center_y, self.renderer.xarray[strand[-1]] - center_x)

                shifted_start_rad %= (2 * np.pi)
                shifted_end_rad %= (2 * np.pi)

                angle_shift_rad = self.renderer.NODE_R * 2 / radius

                if num_nodes > 1:
                    shifted_start_rad -= angle_shift_rad
                    shifted_end_rad += angle_shift_rad

                if shifted_end_rad > shifted_start_rad:
                    shifted_end_rad -= 2 * np.pi

                angular_span_rad = ((shifted_start_rad + angle_shift_rad / 2) - (shifted_end_rad - angle_shift_rad / 2)) % (2 * np.pi)
                arc_length = radius * angular_span_rad

                if num_nodes == 1:
                    arc_length -= self.renderer.NODE_R * 4

                if space_needed > arc_length:
                    radius += 1
                    self.set_radius(junction, radius, auto_call=True)
                    radius_updated = True

            shifted_start_rad = np.arctan2(self.renderer.yarray[strand[0]] - center_y, self.renderer.xarray[strand[0]] - center_x)
            shifted_end_rad = np.arctan2(self.renderer.yarray[strand[-1]] - center_y, self.renderer.xarray[strand[-1]] - center_x)

            shifted_start_rad %= (2 * np.pi)
            shifted_end_rad %= (2 * np.pi)

            angle_shift_rad = self.renderer.NODE_R * 2 / radius
            shifted_start_rad -= angle_shift_rad
            shifted_end_rad += angle_shift_rad

            if shifted_end_rad > shifted_start_rad:
                shifted_end_rad -= 2 * np.pi

            if space_needed + self.renderer.NODE_R < arc_length or num_nodes == 1:
                start_rad = np.arctan2(self.renderer.yarray[strand[0]] - center_y, self.renderer.xarray[strand[0]] - center_x)
                end_rad = np.arctan2(self.renderer.yarray[strand[-1]] - center_y, self.renderer.xarray[strand[-1]] - center_x)

                start_rad %= (2 * np.pi)
                end_rad %= (2 * np.pi)

                if end_rad > start_rad:
                    end_rad -= 2 * np.pi

                total_angular_span_rad = start_rad - end_rad
                angle_between_nodes_rad = total_angular_span_rad / (num_nodes + 1)
                angles = [start_rad - angle_between_nodes_rad - i * angle_between_nodes_rad for i in range(num_nodes)]
            else:
                angles = np.linspace(shifted_start_rad, shifted_end_rad, num_nodes)

            for i, node in enumerate(strand[1:-1]):
                self.renderer.xarray[node] = center_x + radius * np.cos(angles[i])
                self.renderer.yarray[node] = center_y + radius * np.sin(angles[i])

        if radius_updated:
            self.update_unpaired_strands_positions(junction)

    def set_radius(self, junction, new_radius, auto_call=False):
        """
        Sets the radius of the junction and updates the positions of related nodes accordingly.

        Args:
        junction (object): The junction object whose radius is being set.
        new_radius (float): The new radius to set.
        auto_call (bool): Flag to indicate if the function is called automatically as part of an adjustment process.
                        If True, it prevents further recursive calls to `update_unpaired_strands_positions` to avoid infinite recursion.

        Returns:
        None
        """
        center_x, center_y = self.get_junction_center(junction)

        visited_junctions = set()
        visited_nodes = set()

        visited_junctions.add(junction)

        for child in junction.children:
            delta_x, delta_y = calculate_translation_vectors(
                junction,
                child,
                center_x,
                center_y,
                new_radius,
                self.renderer.xarray,
                self.renderer.yarray,
                self.get_first_nucleotides,
                self.get_last_nucleotides,
                self.get_junction_radius
            )
            translate_element_positions(
                self.renderer.xarray,
                self.renderer.yarray,
                child,
                delta_x,
                delta_y,
                visited_junctions,
                visited_nodes
            )

        if junction.has_parent():
            delta_x, delta_y = calculate_translation_vectors(
                junction,
                junction.parent,
                center_x,
                center_y,
                new_radius,
                self.renderer.xarray,
                self.renderer.yarray,
                self.get_first_nucleotides,
                self.get_last_nucleotides,
                self.get_junction_radius,
                True
            )
            translate_element_positions(
                self.renderer.xarray,
                self.renderer.yarray,
                junction.parent,
                delta_x,
                delta_y,
                visited_junctions,
                visited_nodes
            )

        if auto_call is not True:
            self.update_unpaired_strands_positions(junction)

        junction.radius = self.calculate_junction_radius(junction)

    def straighten_unpaired_strands(self):
        """
        Straightens unpaired strands of nucleotides in the RNA structure.
        Forms a straight horizontal line of nucleotides.
        Between strands are paired nucleotides that form junctions or helices.

        Args:
        renderer (object): The RNARenderer instance.

        Returns:
        list: A list of tuples representing the start and end positions of between strands.
        """
        if len(self.renderer.struct.get_single_strands()) > 0:
            single_strand_chunks = self.renderer.struct.get_single_strands()
            first_id = single_strand_chunks[0].positions[0]
            direction_vector = [1, 0]
            self.renderer.direction_vector = direction_vector

            number = 0
            between_strands = []
            between_strand_positions = {}

            total_offset = [0, 0]

            # Add first strand
            if first_id > 0:
                between_strands.append((0, first_id - 1))

                original_dx = self.renderer.xarray[first_id - 1] - self.renderer.xarray[0]
                original_dy = self.renderer.yarray[first_id - 1] - self.renderer.yarray[0]
                magnitude = math.sqrt(original_dx**2 + original_dy**2)
                extra_x = magnitude * direction_vector[0]
                extra_y = magnitude * direction_vector[1]

                between_strand_positions[first_id - 1] = (
                    self.renderer.xarray[first_id] - self.renderer.NODE_R * 2 * direction_vector[0],
                    self.renderer.yarray[first_id] - self.renderer.NODE_R * 2 * direction_vector[1]
                )

                between_strand_positions[0] = (
                    between_strand_positions[first_id - 1][0] - extra_x,
                    between_strand_positions[first_id - 1][1] - extra_y
                )

            for i, strand in enumerate(single_strand_chunks):
                positions = strand.positions
                if i != 0:
                    start, end = last_end+1, positions[0]-1
                    between_strands.append((start, end))

                    original_dx = self.renderer.xarray[end] - self.renderer.xarray[start]
                    original_dy = self.renderer.yarray[end] - self.renderer.yarray[start]
                    magnitude = math.sqrt(original_dx**2 + original_dy**2)
                    extra_x = magnitude * direction_vector[0]
                    extra_y = magnitude * direction_vector[1]
                    
                    between_strand_positions[start] = (self.renderer.xarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[0] * (number-2) + total_offset[0],
                                                    self.renderer.yarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[1] * (number-2) + total_offset[1])
                    
                    total_offset[0] += extra_x
                    total_offset[1] += extra_y
                    number -= 1

                    between_strand_positions[end] = (self.renderer.xarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[0] * (number-1) + total_offset[0],
                                                    self.renderer.yarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[1] * (number-1) + total_offset[1])

                for pos in positions:
                    self.renderer.xarray[pos] = self.renderer.xarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[0] * number + total_offset[0]
                    self.renderer.yarray[pos] = self.renderer.yarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[1] * number + total_offset[1]
                    number += 1

                last_end = positions[-1]
                number += 2

            if last_end < len(self.renderer.xarray) - 1:
                start, end = last_end+1, len(self.renderer.xarray) - 1
                between_strands.append((start, end))

                original_dx = self.renderer.xarray[end] - self.renderer.xarray[start]
                original_dy = self.renderer.yarray[end] - self.renderer.yarray[start]
                magnitude = math.sqrt(original_dx**2 + original_dy**2)
                extra_x = magnitude * direction_vector[0]
                extra_y = magnitude * direction_vector[1]
                
                between_strand_positions[start] = (self.renderer.xarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[0] * (number-2) + total_offset[0],
                                                self.renderer.yarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[1] * (number-2) + total_offset[1])
                
                total_offset[0] += extra_x
                total_offset[1] += extra_y
                number -= 1
            
                between_strand_positions[end] = (self.renderer.xarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[0] * (number-1) + total_offset[0],
                                                self.renderer.yarray[first_id] + self.renderer.NODE_R * 2 * direction_vector[1] * (number-1) + total_offset[1])
                
            # Validate between strands
            for i, between in enumerate(between_strands):
                valid = False
                valid_positions = None
                for group in self.renderer.struct.get_junctions() + self.renderer.struct.get_helices():
                    if between[0] in group.positions and between[1] in group.positions:
                        valid = True
                    elif between[0] in group.positions:
                        valid_positions = group.positions
                if not valid:
                    for group in self.renderer.struct.get_junctions() + self.renderer.struct.get_helices():
                        if valid_positions[-1] + 1 in group.positions:
                            between_strands[i] = (between[0], valid_positions[-1])
                            between_strands.append((group.positions[0], group.positions[-1]))

            between_strands.sort(key=lambda x: x[0])

            # Adjust between strands
            for i in range(len(between_strands) - 1):
                current_strand = between_strands[i]
                next_strand = between_strands[i+1]
                
                if current_strand[1] + 1 == next_strand[0]:
                    position = between_strand_positions[current_strand[0]]

                    original_dx = self.renderer.xarray[current_strand[0]] - self.renderer.xarray[current_strand[1]]
                    original_dy = self.renderer.yarray[current_strand[0]] - self.renderer.yarray[current_strand[1]]
                    magnitude = math.sqrt(original_dx**2 + original_dy**2)
                    extra_x = magnitude * direction_vector[0]
                    extra_y = magnitude * direction_vector[1]

                    between_strand_positions[current_strand[1]] = (position[0] + extra_x,
                                                                position[1])
                    between_strand_positions[next_strand[0]] = (position[0] + 2 * 2 * self.renderer.NODE_R,
                                                                position[1])

            # Update positions
            for set in between_strands:
                saved_original = {}
                for node in set:
                    saved_original[node] = (self.renderer.xarray[node], self.renderer.yarray[node])
                    self.renderer.xarray[node] = between_strand_positions[node][0]
                    self.renderer.yarray[node] = between_strand_positions[node][1]

                start, end = set
                old_dx = saved_original[end][0] - saved_original[start][0]
                old_dy = saved_original[end][1] - saved_original[start][1]

                new_dx = self.renderer.xarray[end] - self.renderer.xarray[start]
                new_dy = self.renderer.yarray[end] - self.renderer.yarray[start]

                for i in range(start + 1, end):
                    old_pos = (self.renderer.xarray[i], self.renderer.yarray[i])
                    dx = old_pos[0] - saved_original[start][0]
                    dy = old_pos[1] - saved_original[start][1]

                    angle_old = math.atan2(old_dy, old_dx)
                    angle_new = math.atan2(new_dy, new_dx)
                    rotation_angle = angle_new - angle_old

                    rotated_dx = dx * math.cos(rotation_angle) - dy * math.sin(rotation_angle)
                    rotated_dy = dx * math.sin(rotation_angle) + dy * math.cos(rotation_angle)

                    self.renderer.xarray[i] = self.renderer.xarray[start] + rotated_dx
                    self.renderer.yarray[i] = self.renderer.yarray[start] + rotated_dy

                if hasattr(self.renderer, 'center_x') and hasattr(self.renderer, 'center_y'):
                    dx = self.renderer.center_x - saved_original[start][0]
                    dy = self.renderer.center_y - saved_original[start][1]

                    rotated_dx = dx * math.cos(rotation_angle) - dy * math.sin(rotation_angle)
                    rotated_dy = dx * math.sin(rotation_angle) + dy * math.cos(rotation_angle)

                    self.renderer.center_x = self.renderer.xarray[start] + rotated_dx
                    self.renderer.center_y = self.renderer.yarray[start] + rotated_dy

            return between_strands

    def add_horizontal_distance_between(self, between_strands):
        """
        Prevents overlap by adding spacing along the x-axis between "between_nodes".

        Args:
        renderer (object): The RNARenderer instance.
        between_strands (list): List of tuples representing the start and end positions of between strands.

        Returns:
        None
        """
        distance = 0
        n = 0

        if not between_strands:
            # No between strands to add space between.
            return
                        
        def shift_nodes_x_axis(start, value):
            """
            Shifts nodes along the x-axis by a given value starting from a specific node.

            Args:
            start (int): Starting node position.
            value (float): Distance to shift nodes along the x-axis.

            Returns:
            None
            """
            for node in range(start, self.renderer.length):
                self.renderer.xarray[node] += value

            for junction in self.renderer.struct.get_junctions():
                if junction.positions[0] >= start:
                    junction.center_x += value

        # Process each between strand to add horizontal distance
        for between in between_strands[1:]:
            best_overlap = float('inf')
            distance = 0
            n = 0
            
            overlap = self.update_overlap_count(node_array=list(range(between[0], between[1])), specific_check=True)

            while overlap > 0 and n < 15:
                magnitude = 1 + n ** 2
                shift_nodes_x_axis(start=between[0], value=magnitude)
                distance += magnitude

                overlap = self.update_overlap_count(node_array=list(range(between[0], between[1])), specific_check=True)

                if overlap < best_overlap:
                    best_overlap = overlap

                n += 1

        # Update junction centers after adjusting the distances
        for junction in self.renderer.struct.get_junctions():
            self.update_junction_center(junction)

    def update_overlap_count(self, node_array=None, draw=False, specific_check=False, only_inside=False):
        """
        Updates the overlap count for nodes in the RNA structure. Determines if there is any overlap between nodes.

        Args:
        renderer (object): The RNARenderer instance.
        node_array (list, optional): List of node positions to check for overlap. If None, all nodes are checked.
        draw (bool, optional): Flag to indicate if the drawing should be updated. Defaults to False.
        specific_check (bool, optional): Flag to check for overlap specifically within the node_array. Defaults to False.
        only_inside (bool, optional): Flag to check for overlap only within the node_array. Defaults to False.

        Returns:
        int: The number of overlaps found.
        """
        container = ChunkContainer(self.renderer.struct, self.renderer.xarray, self.renderer.yarray, node_array=node_array, draw=draw, draw_num=self.renderer.structure_draw_num)
        self.renderer.structure_draw_num += 1

        if specific_check:
            # Check all chunks (in node_array) against all other chunks (after node_array) for overlap.
            return container.check_specific_overlap()
        
        if only_inside:
            # Check chunks (in node_array) against themselves and all other chunks (in node_array) for overlap.
            return container.check_any_overlap_in_node_array()

        # Check all chunks against all other chunks for overlap.
        return container.check_any_overlap()

    def estimate_inches_using_curve_fit(self, x_range, y_range):
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

    def get_between_strand(self, junction, between_strands, xarray):
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

    def explore_paths(self):
        """
        Explore possible paths to optimize the structure by minimizing overlap.
        
        Args:
        renderer (object): The renderer object containing the RNA structure and methods.
        
        Returns:
        int: The overlap count after exploring paths.
        """
        total_acc = 0
        total_tried = 0
        between_strands = self.renderer.between_strands

        for index, junction in enumerate(self.renderer.struct.get_junctions()):
            # Set possible angles for branches based on the number of children
            if len(junction.children) == 1:
                angles = [180, 270, 90]
            elif len(junction.children) <= 3:
                angles = [270, 180, 90]
            else:
                angles = [315, 270, 225, 180, 135, 90, 45]
            
            if self.renderer.global_best_combo.get(junction) is None:
                self.renderer.global_best_combo[junction] = None
            if self.renderer.global_best_overlap.get(junction) is None:
                self.renderer.global_best_overlap[junction] = None

            best_overlap = None
            best_combo = None

            # Try all combinations of angles for the children
            for comb in combinations(angles, len(junction.children)):
                for child_index, child in enumerate(junction.children):
                    self.set_branch_angle(junction, child, comb[child_index])

                self.update_unpaired_strands_positions(junction)

                overlap_count, nodes_below_straight_strand = None, 0
                between = self.get_between_strand(junction, between_strands, self.renderer.xarray)

                overlap_count = self.update_overlap_count(node_array=list(range(between[0], between[1])), only_inside=True)

                # Count nodes below the straight strand
                for node_pos in self.renderer.yarray:
                    if node_pos < self.renderer.yarray[0] + self.renderer.NODE_R:
                        nodes_below_straight_strand += 1
                
                # Update global best overlap and combination
                if self.renderer.global_best_overlap.get(junction) is None or overlap_count < self.renderer.global_best_overlap[junction]:
                    self.renderer.global_best_overlap[junction] = overlap_count
                    self.renderer.global_best_combo[junction] = comb

                # Update best local overlap and combination
                if best_overlap is None or overlap_count < best_overlap:
                    best_overlap = overlap_count
                    best_combo = comb
                elif best_overlap is not None and overlap_count - best_overlap > 0:
                    # Simulated annealing approach
                    T = np.sqrt(len(self.renderer.xarray)) / 8
                    diff = overlap_count - best_overlap
                    val = np.exp(-diff / T)
                    r = random.uniform(0, 1)
                    if val > r:
                        best_overlap = overlap_count
                        best_combo = comb
                        total_acc += 1
                    total_tried += 1

            # Set the best angles found for the children
            for child, angle in zip(junction.children, self.renderer.global_best_combo[junction]):
                self.set_branch_angle(junction, child, angle)
                self.update_unpaired_strands_positions(junction)

        self.add_horizontal_distance_between(between_strands)

        overlap_count = self.update_overlap_count()
        return overlap_count

    def straighten_branches(self, global_best_overlap):
        """
        Straighten branches to minimize overlap in the RNA structure.
        Note: This is due to cases in which 90 degree angles being prioritized in the RNA structure. If the algorithm is setup more effectively, this method won't be needed.
        
        Args:
        renderer (object): The renderer object containing the RNA structure and methods.
        global_best_overlap (int): The best overlap count found globally.
        """
        for junction in self.renderer.struct.get_junctions():
            if len(junction.children) == 1:
                parent_pos_x, parent_pos_y, parent_angle_deg = self.get_junction_parent_data(junction, junction.center_x, junction.center_y)
                # Note get_junction_branch_angle error. It is not always perfectly accurate. 
                original_angle = self.get_junction_branch_angle(junction, junction.children[0], parent_pos_x, parent_pos_y, parent_angle_deg)
                if original_angle > 185:
                    original_angle = 270
                elif original_angle < 175:
                    original_angle = 90
                else:
                    original_angle = 180

                self.set_branch_angle(junction, junction.children[0], 180)
                self.update_unpaired_strands_positions(junction)

                between = self.get_between_strand(junction, self.renderer.between_strands, self.renderer.xarray)
                new_overlap_count = self.update_overlap_count(node_array=list(range(between[0], between[1])), only_inside=True)

                if new_overlap_count > global_best_overlap:
                    self.set_branch_angle(junction, junction.children[0], original_angle)
                    self.update_unpaired_strands_positions(junction)