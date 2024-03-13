import os
import shutil
import sys
import re
import random
import math
import time
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch
from matplotlib.axes import Axes
from rna_secstruct import SecStruct, MotifSearchParams
import numpy as np
from scipy.interpolate import CubicSpline
import scipy.interpolate
from matplotlib.patches import Arc
from itertools import combinations
from scipy.spatial import cKDTree
import pandas as pd
from scipy.optimize import curve_fit
from rna_draw.chunk_controller import ChunkContainer
#import pygame
#from pygame.locals import *


class Test:
    def __init__(self, something):
        self.something = something


test_1 = Test("something")


class RNATreeNode:
    def __init__(self):
        self.children_ = []
        self.is_pair_ = False
        self.index_a_ = -1
        self.index_b_ = -1
        self.x_ = 0
        self.y_ = 0
        self.go_x_ = 0
        self.go_y_ = 0


def get_pairmap_from_secstruct(secstruct):
    """
    generates dictionary containing pair mappings
    args:
    secstruct contains secondary structure string
    returns:
    dictionary with pair mappings
    """
    pair_stack = []
    end_stack = []
    pairs_array = []
    i_range = range(0, len(secstruct))

    # initialize all values to -1, meaning no pair
    for ii in i_range:
        pairs_array.append(-1)

    # assign pairs based on secstruct
    for ii in i_range:
        if secstruct[ii] == "(":
            pair_stack.append(ii)
        elif secstruct[ii] == ")":
            if not pair_stack:
                end_stack.append(ii)
            else:
                index = pair_stack.pop()
                pairs_array[index] = ii
                pairs_array[ii] = index
    if len(pair_stack) == len(end_stack):
        n = len(pair_stack)
        for ii in range(n):
            pairs_array[pair_stack[ii]] = end_stack[-ii]
            pairs_array[end_stack[-ii]] = pair_stack[ii]
    else:
        print("ERROR: pairing incorrect %s" % secstruct)

    return pairs_array


def add_nodes_recursive(bi_pairs, rootnode, start_index, end_index):
    #if start_index > end_index:
        #print("Error occured while drawing RNA %d %d" % (start_index, end_index))
        #sys.exit(0)

    if bi_pairs[start_index] == end_index:
        newnode = RNATreeNode()
        newnode.is_pair_ = True
        newnode.index_a_ = start_index
        newnode.index_b_ = end_index

        add_nodes_recursive(bi_pairs, newnode, start_index + 1, end_index - 1)

    else:
        newnode = RNATreeNode()
        jj = start_index
        while jj <= end_index:
            if bi_pairs[jj] > jj:
                add_nodes_recursive(bi_pairs, newnode, jj, bi_pairs[jj])
                jj = bi_pairs[jj] + 1
            else:
                newsubnode = RNATreeNode()
                newsubnode.is_pair_ = False
                newsubnode.index_a_ = jj
                newnode.children_.append(newsubnode)
                jj += 1

    rootnode.children_.append(newnode)


def setup_coords_recursive(
    rootnode,
    parentnode,
    start_x,
    start_y,
    go_x,
    go_y,
    NODE_R,
    PRIMARY_SPACE,
    PAIR_SPACE,
    radii
):
    cross_x = -go_y
    cross_y = go_x

    children_width = len(rootnode.children_) * NODE_R * 2

    rootnode.go_x_ = go_x
    rootnode.go_y_ = go_y

    if len(rootnode.children_) == 1:
        rootnode.x_ = start_x
        rootnode.y_ = start_y

        if rootnode.children_[0].is_pair_:
            setup_coords_recursive(
                rootnode.children_[0],
                rootnode,
                start_x + go_x * PRIMARY_SPACE,
                start_y + go_y * PRIMARY_SPACE,
                go_x,
                go_y,
                NODE_R,
                PRIMARY_SPACE,
                PAIR_SPACE,
                radii
            )
        elif (
            rootnode.children_[0].is_pair_ == False
            and rootnode.children_[0].index_a_ < 0
        ):
            setup_coords_recursive(
                rootnode.children_[0],
                rootnode,
                start_x,
                start_y,
                go_x,
                go_y,
                NODE_R,
                PRIMARY_SPACE,
                PAIR_SPACE,
                radii
            )
        else:
            setup_coords_recursive(
                rootnode.children_[0],
                rootnode,
                start_x + go_x * PRIMARY_SPACE,
                start_y + go_y * PRIMARY_SPACE,
                go_x,
                go_y,
                NODE_R,
                PRIMARY_SPACE,
                PAIR_SPACE,
                radii
            )

    elif len(rootnode.children_) > 1:
        npairs = 0
        for ii in range(0, len(rootnode.children_)):
            if rootnode.children_[ii].is_pair_:
                npairs += 1

        circle_length = (len(rootnode.children_) + 1) * PRIMARY_SPACE + (
            npairs + 1
        ) * PAIR_SPACE
        circle_radius = circle_length / (2 * math.pi)
        length_walker = PAIR_SPACE / 2.0

        if parentnode == None:
            rootnode.x_ = go_x * circle_radius
            rootnode.y_ = go_y * circle_radius
        else:
            rootnode.x_ = parentnode.x_ + go_x * circle_radius
            rootnode.y_ = parentnode.y_ + go_y * circle_radius

        for ii in range(0, len(rootnode.children_)):
            length_walker += PRIMARY_SPACE

            if rootnode.children_[ii].is_pair_:
                length_walker += PAIR_SPACE / 2.0

            radii.insert(rootnode.children_[ii].index_a_, circle_radius)

            rad_angle = length_walker / circle_length * 2 * math.pi - math.pi / 2.0
            child_x = (
                rootnode.x_
                + math.cos(rad_angle) * cross_x * circle_radius
                + math.sin(rad_angle) * go_x * circle_radius
            )
            child_y = (
                rootnode.y_
                + math.cos(rad_angle) * cross_y * circle_radius
                + math.sin(rad_angle) * go_y * circle_radius
            )

            child_go_x = child_x - rootnode.x_
            child_go_y = child_y - rootnode.y_
            child_go_len = math.sqrt(child_go_x * child_go_x + child_go_y * child_go_y)

            setup_coords_recursive(
                rootnode.children_[ii],
                rootnode,
                child_x,
                child_y,
                child_go_x / child_go_len,
                child_go_y / child_go_len,
                NODE_R,
                PRIMARY_SPACE,
                PAIR_SPACE,
                radii
            )

            if rootnode.children_[ii].is_pair_:
                length_walker += PAIR_SPACE / 2.0

    else:
        rootnode.x_ = start_x
        rootnode.y_ = start_y

    return radii


def get_coords_recursive(rootnode, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE):
    if rootnode.is_pair_:
        cross_x = -rootnode.go_y_
        cross_y = rootnode.go_x_

        xarray[rootnode.index_a_] = rootnode.x_ + cross_x * PAIR_SPACE / 2.0
        xarray[rootnode.index_b_] = rootnode.x_ - cross_x * PAIR_SPACE / 2.0

        yarray[rootnode.index_a_] = rootnode.y_ + cross_y * PAIR_SPACE / 2.0
        yarray[rootnode.index_b_] = rootnode.y_ - cross_y * PAIR_SPACE / 2.0
    elif rootnode.index_a_ >= 0:
        xarray[rootnode.index_a_] = rootnode.x_
        yarray[rootnode.index_a_] = rootnode.y_

    for ii in range(0, len(rootnode.children_)):
        get_coords_recursive(
            rootnode.children_[ii], xarray, yarray, PRIMARY_SPACE, PAIR_SPACE
        )

class RNARenderer:
    def __init__(self):
        self.root_ = None
        self.xarray_ = None
        self.yarray_ = None
        self.size_ = None
        self.fig = plt.Figure()
        self.ax = self.fig.add_subplot(111, aspect="equal")
        self.chunks = []
        self.struct = None
        self.structure_draw_num = 0

    # Returns Nucleotide Pair in Junction used to Determine Angle
    def get_first_nucleotides(self, strands):
        first_nucleotide = strands[0][0]
        last_nucleotide = strands[-1][-1]
        return first_nucleotide, last_nucleotide
    
    def get_last_nucleotides(self, strands):
        first_nucleotide = strands[0][-1]
        last_nucleotide = strands[-1][0]
        return first_nucleotide, last_nucleotide
    
    def get_junction_center(self, junction):
        if hasattr(junction, 'center_x') and hasattr(junction, 'center_y'):
            return junction.center_x, junction.center_y
        else:
            return self.update_junction_center(junction)
        
    def junction_substituents(self, junction):
        return len(junction.children) + (1 if junction.has_parent() else 0)
    
    def junction_substituents_from_strands(self, strands):
        count = 0

        for strand in strands:
            if len(strand[1:-1]) > 0:
                count += 1

        return count
    
    def get_junction_radius(self, junction):
        if hasattr(junction, "radius"):
            return junction.radius
        else:
            junction.radius = self.calculate_junction_radius(junction)
            return junction.radius
    
    def calculate_junction_radius(self, junction):
        center_x, center_y = self.get_junction_center(junction)

        if junction.has_parent():
            parent_strand = junction.parent
        elif junction.has_children():
            parent_strand = junction.children[0]
        else:
            print("Junction has no parent or child.")
            return None

        first, last = self.get_last_nucleotides(parent_strand.strands)

        mid_x = (self.xarray[first] + self.xarray[last]) / 2
        mid_y = (self.yarray[first] + self.yarray[last]) / 2

        radius = np.sqrt((mid_x - center_x) ** 2 + (mid_y - center_y) ** 2)

        return radius
            
    def update_junction_center(self, junction, inversion_check=False):
        if (hasattr(junction, "radius") and junction.has_parent()):
            parent_strand = junction.parent
            first, last = self.get_last_nucleotides(parent_strand.strands)
            mid_x = (self.xarray[first] + self.xarray[last]) / 2
            mid_y = (self.yarray[first] + self.yarray[last]) / 2

            direction_x = self.xarray[last] - self.xarray[first]
            direction_y = self.yarray[last] - self.yarray[first]

            magnitude = np.sqrt(direction_x**2 + direction_y**2)
            unit_x = direction_x / magnitude
            unit_y = direction_y / magnitude

            perp_x = -unit_y
            perp_y = unit_x

            # Fixes issue when flipping across x axis.
            centroid_x = sum(self.xarray[node] for node in junction.positions) / len(junction.positions)
            centroid_y = sum(self.yarray[node] for node in junction.positions) / len(junction.positions)
            
            direction_to_centroid_x = centroid_x - self.xarray[first]
            direction_to_centroid_y = centroid_y - self.yarray[first]
            
            dot_product = perp_x * direction_to_centroid_x + perp_y * direction_to_centroid_y
            
            if dot_product < 0 and inversion_check == True:
                perp_x, perp_y = -perp_x, -perp_y

            center_x = mid_x + junction.radius * perp_x
            center_y = mid_y + junction.radius * perp_y
        else:
            x_coords = [self.xarray[value] for value in junction.positions]
            y_coords = [self.yarray[value] for value in junction.positions]
            center_x = sum(x_coords) / len(x_coords)
            center_y = sum(y_coords) / len(y_coords)

        junction.center_x = center_x
        junction.center_y = center_y

        return center_x, center_y

    # Gives Information About a Junction, such as the X,Y cords and the angle relative to the origin.
    def get_junction_parent_data(self, junction, center_x, center_y):
        if junction.has_parent():
            left, right = self.get_last_nucleotides(junction.parent.strands)
            parent_pos_x = (self.xarray[left] + self.xarray[right]) / 2
            parent_pos_y = (self.yarray[left] + self.yarray[right]) / 2

            parent_angle_rad = np.arctan2(parent_pos_y - center_y, parent_pos_x - center_x)
            parent_angle_deg = np.rad2deg(parent_angle_rad)

            return parent_pos_x, parent_pos_y, parent_angle_deg
    
    # Calculates the Angle of a Single Branch with Respect to the Parent Helix. (Angle is Counter-Clockwise)
    def get_junction_branch_angle(self, junction, child, center_x, center_y, parent_angle_deg):
        data_new = self.struct[child.m_id]
        left, right = self.get_first_nucleotides(data_new.strands)
        child_pos_x = (self.xarray[left] + self.xarray[right]) / 2
        child_pos_y = (self.yarray[left] + self.yarray[right]) / 2

        child_angle_rad = math.atan2(child_pos_y - center_y, child_pos_x - center_x)
        child_angle_deg = math.degrees(child_angle_rad)

        angle_between = (child_angle_deg - parent_angle_deg + 360) % 360

        return angle_between

    def get_children_angles(self, junction):
        angles = []

        center_x, center_y = self.get_junction_center(junction)
        x, y, parent_angle_deg = self.get_junction_parent_data(junction, center_x, center_y)

        if junction.has_children():
            for branch in junction.children:
                angles.append(self.get_junction_branch_angle(junction, branch, center_x, center_y, parent_angle_deg))
        
        return angles
    
    def rotate_point(self, point, origin, angle):
        px, py = point
        ox, oy = origin

        rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)], 
                                    [np.sin(angle),  np.cos(angle)]])
        
        result = np.dot(rotation_matrix, np.array([px - ox, py - oy])) + np.array([ox, oy])
        
        return result[0], result[1]
    
    def set_branch_angle(self, junction, child, target_angle):
        center_x, center_y = self.get_junction_center(junction)
        x, y, parent_angle_deg = self.get_junction_parent_data(junction, center_x, center_y)

        current_angle = self.get_junction_branch_angle(junction, child, center_x, center_y, parent_angle_deg)
        angle_diff = math.radians(target_angle) - math.radians(current_angle)

        queue = [child]
        processed_nucleotides = set()

        while queue:
            current_branch = queue.pop(0)
            data_new = self.struct[current_branch.m_id]

            for strand in data_new.strands:
                for nucleotide in strand:
                    if nucleotide not in processed_nucleotides:
                        new_x, new_y = self.rotate_point((self.xarray[nucleotide], self.yarray[nucleotide]), (center_x, center_y), angle_diff)
                        self.xarray[nucleotide] = new_x
                        self.yarray[nucleotide] = new_y
                        processed_nucleotides.add(nucleotide)

            if current_branch.is_junction():
                if (hasattr(current_branch, 'center_x') and hasattr(current_branch, 'center_y')):
                    new_center_x, new_center_y = self.rotate_point((current_branch.center_x, current_branch.center_y), (center_x, center_y), angle_diff)
                    current_branch.center_x = new_center_x
                    current_branch.center_y = new_center_y

            if current_branch.has_children():
                for child_branch in current_branch.children:
                    queue.append(child_branch)
        
        #left, right = self.get_adjacent_unpaired_strands(junction, junction.children.index(child))

        #self.update_unpaired_strands_positions(junction)

    def get_adjacent_unpaired_strands(self, junction, index):
        if index < 0 or index >= len(junction.strands):
            raise ValueError("Index out of range in junction strands")

        left_strand = junction.strands[index]

        right_strand = junction.strands[index + 1] if index + 1 < len(junction.strands) else junction.parent

        return left_strand, right_strand
    
    def update_unpaired_strands_positions(self, junction):
        center_x, center_y = self.get_junction_center(junction)
        radius = self.get_junction_radius(junction)

        _, _, parent_angle_deg = self.get_junction_parent_data(junction, center_x, center_y)

        radius_updated = False

        for strand in junction.strands:
            num_nodes = len(strand[1:-1])
            space_needed = (num_nodes) * self.NODE_R * 2
            arc_length = -math.inf

            shifted_start_rad = None
            shifted_end_rad = None
            angle_shift_rad = None

            tolerance = self.NODE_R

            while arc_length < space_needed and space_needed != 0:
            #while space_needed != 0 and abs(arc_length - space_needed) < tolerance:
                shifted_start_rad = np.arctan2(self.yarray[strand[0]] - center_y, self.xarray[strand[0]] - center_x)
                shifted_end_rad = np.arctan2(self.yarray[strand[-1]] - center_y, self.xarray[strand[-1]] - center_x)

                shifted_start_rad %= (2 * np.pi)
                shifted_end_rad %= (2 * np.pi)

                angle_shift_rad = self.NODE_R * 2 / radius

                if num_nodes > 1:
                    shifted_start_rad -= angle_shift_rad
                    shifted_end_rad += angle_shift_rad

                if shifted_end_rad > shifted_start_rad:
                    shifted_end_rad -= 2 * np.pi

                angular_span_rad = ((shifted_start_rad + angle_shift_rad/2) - (shifted_end_rad - angle_shift_rad/2)) % (2 * np.pi)
                arc_length = radius * angular_span_rad

                if num_nodes == 1:
                    arc_length -= self.NODE_R * 4

                tolerance = self.NODE_R

                #if abs(arc_length - space_needed) > tolerance:
                if space_needed > arc_length:
                    radius += 1
                    self.set_radius(junction, radius, auto_call=True)
                    radius_updated = True
                    #elif space_needed < arc_length:
                        #radius *= 0.9
                        #self.set_radius(junction, radius, auto_call=True)
                        #radius_updated = True

            shifted_start_rad = np.arctan2(self.yarray[strand[0]] - center_y, self.xarray[strand[0]] - center_x)
            shifted_end_rad = np.arctan2(self.yarray[strand[-1]] - center_y, self.xarray[strand[-1]] - center_x)

            shifted_start_rad %= (2 * np.pi)
            shifted_end_rad %= (2 * np.pi)

            angle_shift_rad = self.NODE_R * 2 / radius
            shifted_start_rad -= angle_shift_rad
            shifted_end_rad += angle_shift_rad

            if shifted_end_rad > shifted_start_rad:
                shifted_end_rad -= 2 * np.pi

            if space_needed + self.NODE_R < arc_length or num_nodes == 1:
                start_rad = np.arctan2(self.yarray[strand[0]] - center_y, self.xarray[strand[0]] - center_x)
                end_rad = np.arctan2(self.yarray[strand[-1]] - center_y, self.xarray[strand[-1]] - center_x)

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
                self.xarray[node] = center_x + radius * np.cos(angles[i])
                self.yarray[node] = center_y + radius * np.sin(angles[i])

        if radius_updated:
            self.update_unpaired_strands_positions(junction)

    def set_radius(self, junction, new_radius, auto_call=False):
        center_x, center_y = self.get_junction_center(junction)

        visited_junctions = set()
        visited_nodes = set()

        visited_junctions.add(junction)

        for child in junction.children:
            delta_x, delta_y = self._calculate_deltas(original_junction=junction, subpart=child, center_x=center_x, center_y=center_y, new_radius=new_radius)
            self._update_positions(child, delta_x, delta_y, visited_junctions, visited_nodes, upward_recursion=False)

        if junction.has_parent():
            delta_x, delta_y = self._calculate_deltas(original_junction=junction, subpart=junction.parent, center_x=center_x, center_y=center_y, new_radius=new_radius, parent=True)
            self._update_positions(junction.parent, delta_x, delta_y, visited_junctions, visited_nodes, upward_recursion=True)

        if auto_call is not True:
            self.update_unpaired_strands_positions(junction)

        junction.radius = self.calculate_junction_radius(junction)

    def _calculate_deltas(self, original_junction=None, subpart=None, center_x=None, center_y=None, new_radius=None, parent=False):
        first_node = None
        last_node = None

        if parent == False:
            first_node,last_node = self.get_first_nucleotides(subpart.strands)
        else:
            first_node,last_node = self.get_last_nucleotides(subpart.strands)

        midpoint_x = (self.xarray[first_node] + self.xarray[last_node]) / 2
        midpoint_y = (self.yarray[first_node] + self.yarray[last_node]) / 2
        dx = midpoint_x - center_x
        dy = midpoint_y - center_y

        old_radius = self.get_junction_radius(original_junction)

        radius_difference = new_radius - old_radius

        unit_vector_x = dx / old_radius
        unit_vector_y = dy / old_radius

        delta_x = radius_difference * unit_vector_x
        delta_y = radius_difference * unit_vector_y

        return delta_x, delta_y

    def _update_positions(self, junction, delta_x, delta_y, visited_junctions, visited_nodes, upward_recursion=None):
        if junction in visited_junctions:
            return

        visited_junctions.add(junction)

        for node in junction.positions:
            if node in visited_nodes:
                continue
            visited_nodes.add(node)
            self.xarray[node] += delta_x
            self.yarray[node] += delta_y

        if junction.is_junction:
            if hasattr(junction, "center_x"):
                junction.center_x += delta_x
                junction.center_y += delta_y

        if junction.has_children():
            for child in junction.children:
                self._update_positions(child, delta_x, delta_y, visited_junctions, visited_nodes)

        if upward_recursion == True:
            if junction.has_parent():
                self._update_positions(junction.parent, delta_x, delta_y, visited_junctions, visited_nodes, upward_recursion)

    def straighten_unpaired_strands(self):
        if len(self.struct.get_single_strands()) > 0:
            single_strand_chunks = self.struct.get_single_strands()
            first_id = single_strand_chunks[0].positions[0]
            direction_vector = [1, 0]
            self.direction_vector = direction_vector

            number = 0
            between_strands = []
            between_strand_positions = {}

            total_offset = [0, 0]

            if first_id > 0:
                between_strands.append((0, first_id - 1))

                original_dx = self.xarray[first_id - 1] - self.xarray[0]
                original_dy = self.yarray[first_id - 1] - self.yarray[0]
                magnitude = math.sqrt(original_dx**2 + original_dy**2)
                extra_x = magnitude * direction_vector[0]
                extra_y = magnitude * direction_vector[1]

                between_strand_positions[first_id - 1] = (
                    self.xarray[first_id] - self.NODE_R * 2 * direction_vector[0],
                    self.yarray[first_id] - self.NODE_R * 2 * direction_vector[1]
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

                    original_dx = self.xarray[end] - self.xarray[start]
                    original_dy = self.yarray[end] - self.yarray[start]
                    magnitude = math.sqrt(original_dx**2 + original_dy**2)
                    extra_x = magnitude * direction_vector[0]
                    extra_y = magnitude * direction_vector[1]
                    
                    between_strand_positions[start] = (self.xarray[first_id] + self.NODE_R * 2 * direction_vector[0] * (number-2) + total_offset[0],
                                                        self.yarray[first_id] + self.NODE_R * 2 * direction_vector[1] * (number-2) + total_offset[1])
                    
                    total_offset[0] += extra_x
                    total_offset[1] += extra_y
                    number -= 1

                    between_strand_positions[end] = (self.xarray[first_id] + self.NODE_R * 2 * direction_vector[0] * (number-1) + total_offset[0],
                                                    self.yarray[first_id] + self.NODE_R * 2 * direction_vector[1] * (number-1) + total_offset[1])

                for pos in positions:
                    self.xarray[pos] = self.xarray[first_id] + self.NODE_R * 2 * direction_vector[0] * number + total_offset[0]
                    self.yarray[pos] = self.yarray[first_id] + self.NODE_R * 2 * direction_vector[1] * number + total_offset[1]
                    number += 1

                last_end = positions[-1]
                number += 2

            if last_end < len(self.xarray) - 1:
                start, end = last_end+1, len(self.xarray) - 1
                between_strands.append((start, end))

                original_dx = self.xarray[end] - self.xarray[start]
                original_dy = self.yarray[end] - self.yarray[start]
                magnitude = math.sqrt(original_dx**2 + original_dy**2)
                extra_x = magnitude * direction_vector[0]
                extra_y = magnitude * direction_vector[1]
                
                between_strand_positions[start] = (self.xarray[first_id] + self.NODE_R * 2 * direction_vector[0] * (number-2) + total_offset[0],
                                                    self.yarray[first_id] + self.NODE_R * 2 * direction_vector[1] * (number-2) + total_offset[1])
                
                total_offset[0] += extra_x
                total_offset[1] += extra_y
                number -= 1
            
                between_strand_positions[end] = (self.xarray[first_id] + self.NODE_R * 2 * direction_vector[0] * (number-1) + total_offset[0],
                                                self.yarray[first_id] + self.NODE_R * 2 * direction_vector[1] * (number-1) + total_offset[1])
                
            for i, between in enumerate(between_strands):
                valid = False
                valid_positions = None
                for group in self.struct.get_junctions() + self.struct.get_helices():
                    if between[0] in group.positions and between[1] in group.positions:
                        valid = True
                    elif between[0] in group.positions:
                        valid_positions = group.positions
                if not valid:
                    for group in self.struct.get_junctions() + self.struct.get_helices():
                        if valid_positions[-1] + 1 in group.positions:
                            between_strands[i] = (between[0], valid_positions[-1])
                            between_strands.append((group.positions[0], group.positions[-1]))

            for i in range(len(between_strands) - 1):
                current_strand = between_strands[i]
                next_strand = between_strands[i+1]
                
                if current_strand[1] + 1 == next_strand[0]:
                    position = between_strand_positions[current_strand[0]]

                    original_dx = self.xarray[current_strand[0]] - self.xarray[current_strand[1]]
                    original_dy = self.yarray[current_strand[0]] - self.yarray[current_strand[1]]
                    magnitude = math.sqrt(original_dx**2 + original_dy**2)
                    extra_x = magnitude * direction_vector[0]
                    extra_y = magnitude * direction_vector[1]

                    between_strand_positions[current_strand[1]] = (position[0] + extra_x,
                                                                   position[1])
                    between_strand_positions[next_strand[0]] = (position[0] + 2 * 2 * self.NODE_R,
                                                                   position[1])

            for set in between_strands:
                saved_original = {}
                for node in set:
                    saved_original[node] = (self.xarray[node], self.yarray[node])
                    self.xarray[node] = between_strand_positions[node][0]
                    self.yarray[node] = between_strand_positions[node][1]

                start, end = set
                old_dx = saved_original[end][0] - saved_original[start][0]
                old_dy = saved_original[end][1] - saved_original[start][1]

                new_dx = self.xarray[end] - self.xarray[start]
                new_dy = self.yarray[end] - self.yarray[start]

                for i in range(start + 1, end):
                    old_pos = (self.xarray[i], self.yarray[i])
                    dx = old_pos[0] - saved_original[start][0]
                    dy = old_pos[1] - saved_original[start][1]

                    angle_old = math.atan2(old_dy, old_dx)
                    angle_new = math.atan2(new_dy, new_dx)
                    rotation_angle = angle_new - angle_old

                    rotated_dx = dx * math.cos(rotation_angle) - dy * math.sin(rotation_angle)
                    rotated_dy = dx * math.sin(rotation_angle) + dy * math.cos(rotation_angle)

                    self.xarray[i] = self.xarray[start] + rotated_dx
                    self.yarray[i] = self.yarray[start] + rotated_dy

                if hasattr(self, 'center_x') and hasattr(self, 'center_y'):
                    dx = self.center_x - saved_original[start][0]
                    dy = self.center_y - saved_original[start][1]

                    rotated_dx = dx * math.cos(rotation_angle) - dy * math.sin(rotation_angle)
                    rotated_dy = dx * math.sin(rotation_angle) + dy * math.cos(rotation_angle)

                    self.center_x = self.xarray[start] + rotated_dx
                    self.center_y = self.yarray[start] + rotated_dy

            return between_strands
        
    # Checks for Node Overlap between all nodes.
    # Only checks overlap between node_array and other nodes if node_array is specified.

    def in_junction(self, node_number: int) -> bool:
        for junction in self.struct.get_junctions():
            if (node_number in junction.positions):
                return True
        return False    

    # Goes back through and minimizes the distance between.
    def minimize_horizontal_distance_between(self, between_strands):
        return 
    
        for between in between_strands[1:]:
            best_distance = float('inf')
            best_overlap = 0
            distance = 0
            n = 0
            
            def shift_nodes_x_axis(start, value):
                for node in range(start, self.length):
                    self.xarray[node] += value

            while True:
                magnitude = -(1 + n ** 2)  # Negative magnitude to move left
                shift_nodes_x_axis(start=between[0], value=magnitude)
                distance += magnitude
                
                overlap = 0
                #overlap, nodes_below_straight_strand = self.check_node_overlap(node_array=list(range(0, between[0])))
                
                if overlap == 0 and abs(distance) < best_distance:
                    best_distance = abs(distance)
                    
                if overlap > 0:
                    shift_nodes_x_axis(start=between[0], value=-magnitude)  # Move back to the previous position if overlap is detected
                    break
                    
                n += 1

        for junction in self.struct.get_junctions():
            self.update_junction_center(junction)

    # Prevents overlap by adding spacing along the x axis between nodes.
    def add_horizontal_distance_between(self, between_strands):
        distance = 0
        n = 0
        
        def shift_nodes_x_axis(start, value):
            for node in range(start, self.length):
                self.xarray[node] += value

        for between in between_strands[1:]:
            best_distance = 0
            best_overlap = float('inf')
            distance = 0
            n = 0
                        
            while True:
                magnitude = 1 + n ** 2
                shift_nodes_x_axis(start=between[0], value=magnitude)
                distance += magnitude

                overlap = self.update_overlap_count(node_array=list(range(0, between[0])))

                #print("Overlaps:", overlap, best_overlap)
                #print(" ")
                
                #overlap,nodes_below_straight_strand = self.check_node_overlap(node_array=list(range(0, between[0])))
                
                if overlap < best_overlap:
                    print("Updated Best Overlap from", best_overlap, " to ", overlap)
                    best_overlap = overlap
                    best_distance = distance
                    
                if overlap == 0 or n == 10:
                    break
                
                n += 1

        for junction in self.struct.get_junctions():
            self.update_junction_center(junction)

    def get_bounds(self):
        min_x = min(self.xarray)
        max_x = max(self.xarray)
        min_y = min(self.yarray)
        max_y = max(self.yarray)

        return min_x, max_x, min_y, max_y
    
    def update_overlap_count(self, node_array=None, draw=False):
        container = None

        container = ChunkContainer(self.struct, self.xarray, self.yarray, node_array=node_array, draw=draw, draw_num=self.structure_draw_num)
        self.structure_draw_num += 1

        return container.check_any_overlap()

    def setup_tree(self, secstruct, NODE_R, PRIMARY_SPACE, PAIR_SPACE, seq):
        dangling_start = 0
        dangling_end = 0
        bi_pairs = get_pairmap_from_secstruct(secstruct)
        self.length = len(bi_pairs)

        directory = "DrawVisualizer"

        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)

        self.NODE_R = NODE_R
        self.root_ = None

        for ii in range(0, len(bi_pairs)):
            if bi_pairs[ii] < 0:
                dangling_start += 1
            else:
                break

        for ii in (len(bi_pairs) - 1, -1, -1):
            if bi_pairs[ii] < 0:
                dangling_end += 1
            else:
                break

        self.root_ = RNATreeNode()

        # for jj in range(0,len(bi_pairs)):
        jj = 0
        while jj < len(bi_pairs):
            if bi_pairs[jj] > jj:
                add_nodes_recursive(bi_pairs, self.root_, jj, bi_pairs[jj])
                jj = bi_pairs[jj] + 1
            else:
                newsubnode = RNATreeNode()
                newsubnode.is_pair_ = False
                newsubnode.index_a_ = jj
                self.root_.children_.append(newsubnode)
                jj += 1
        xarray = []
        yarray = []

        for ii in range(0, len(secstruct)):
            xarray.append(0.0)
            yarray.append(0.0)

        radii = self.setup_coords(NODE_R, PRIMARY_SPACE, PAIR_SPACE)
        self.get_coords(xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)

        ## Everything Below is New to Modify Angles

        self.xarray = xarray
        self.yarray = yarray

        self.struct = SecStruct(seq, secstruct) 

        for junction in self.struct.get_junctions():
            for strand in junction.strands:
                for node in strand[1:-1]:
                    junction.radius = radii[node]
                    self.get_junction_center(junction)

        between_strands = self.straighten_unpaired_strands()
        #between_strands = None

        if between_strands:
            self.add_horizontal_distance_between(between_strands)

        def estimate_inches_using_curve_fit(x_range, y_range):
            data = {
                "rna_identity": ["hairpin", "t-RNA", "CO-VID19 5' UTR", "50S Ribosome"],
                "rna_area": [3781, 126207, 1150472, 4286761],
                "rna_figsize_variable": [25, 30, 35, 40],
            }
            df = pd.DataFrame(data)
            
            x = df["rna_area"]
            y = df["rna_figsize_variable"]
            
            def test(x, a, b, c):
                return a * (x - b) ** c

            param, param_cov = curve_fit(test, x, y)
            
            area = x_range * y_range
            estimated_inches_x = x_range / (param[0] * (area - param[1]) ** param[2])
            estimated_inches_y = y_range / (param[0] * (area - param[1]) ** param[2])
            
            return estimated_inches_x, estimated_inches_y
        
        x_range = max(self.xarray) - min(self.xarray)
        y_range = max(self.yarray) - min(self.yarray)

        estimated_inches_x, estimated_inches_y = estimate_inches_using_curve_fit(x_range, y_range)

        if estimated_inches_x > 650 or estimated_inches_y > 650:
            print("Structure is too large to output using matplotlib.")
            print('Estimated Inches Required (X):', estimated_inches_x)
            print('Estimated Inches Required (Y):', estimated_inches_y)
            print('Structure BP Length', len(self.xarray))
            raise Exception('Too Large')
        
        #self.setup_pygame_objects()

        start_time = time.time()

        self.global_best_overlap = {}
        self.global_best_combo = {}
        
        def explore_paths():
            total_acc = 0
            total_tried = 0

            # first node is too close sometimes, increase the radius of that first node to prevent closeness.

            # a* or dijstras

            # panda have a weighted distribution for picking different junctions, prefering those with 5 possible solutions.
            # random.choices and supply weights

            # pick random junctions
            for index, junction in enumerate(self.struct.get_junctions()):
                if len(junction.children) == 1:
                    angles = [180]
                    angles = [180, 270, 90]
                elif len(junction.children) <= 3:
                    angles = [270, 180, 90]
                else:
                    angles = [315, 270, 225, 180, 135, 90, 45]

                #global_best_overlap,global_best_combo=None,None
                    
                if self.global_best_combo.get(junction) is None:
                    self.global_best_combo[junction] = None
                if self.global_best_overlap.get(junction) is None:
                    self.global_best_overlap[junction] = None

                best_overlap = None
                best_combo = None

                for comb in combinations(angles, len(junction.children)):
                    for child_index, child in enumerate(junction.children):
                        self.set_branch_angle(junction, child, comb[child_index])

                    self.update_unpaired_strands_positions(junction)

                    overlap_count, nodes_below_straight_strand = None, 0

                    overlap_count = self.update_overlap_count(draw=False)

                    for node_pos in self.yarray:
                        if node_pos < self.yarray[0] + self.NODE_R:
                            nodes_below_straight_strand += 1

                    #if between_strands and nodes_below_straight_strand > 0:
                        #self.add_horizontal_distance_between(between_strands)
                            
                    self.add_horizontal_distance_between(between_strands)

                    if self.global_best_overlap.get(junction) is None or overlap_count < self.global_best_overlap[junction]:
                        self.global_best_overlap[junction] = overlap_count
                        self.global_best_combo[junction] = comb

                    if best_overlap is None or overlap_count < best_overlap:
                        best_overlap = overlap_count
                        best_combo = comb
                    elif best_overlap is not None and overlap_count - best_overlap > 0:
                        T = np.sqrt(len(self.xarray)) / 8 # Temperature
                        # update tempertature based on acceptance percentage
                        # simulated annealing, decrease temperature by a fixed amount
                        diff = overlap_count - best_overlap
                        val = np.exp(-diff / T)
                        r = random.uniform(0,1)
                        if val > r:
                            best_overlap = overlap_count
                            best_combo = comb
                            #print("OVERWRITTEN", diff, val, r)
                            total_acc += 1
                        else:
                            pass
                            #print("NOT OVERWRITTEN", diff, val, r)
                        total_tried += 1
                        #print("Percentage Acceptance:", (total_acc/total_tried) * 100)

                for child, angle in zip(junction.children, self.global_best_combo[junction]):
                    self.set_branch_angle(junction, child, angle)
                    self.update_unpaired_strands_positions(junction)

            overlap_count = self.update_overlap_count()
            #print('Elapsed Seconds:', time.time() - start_time)
            #print('Overlap Count:', overlap_count)
        
            return overlap_count
        
        def straighten_branches(global_best_overlap):
            for junction in self.struct.get_junctions():
                if len(junction.children) == 1:
                    # Retrieve parent data for the junction
                    parent_pos_x, parent_pos_y, parent_angle_deg = self.get_junction_parent_data(junction, junction.center_x, junction.center_y)
                    
                    # Retrieve the original angle of the branch
                    original_angle = self.get_junction_branch_angle(junction, junction.children[0], parent_pos_x, parent_pos_y, parent_angle_deg)

                    if original_angle > 185:
                        original_angle = 270
                    elif original_angle < 175:
                        original_angle = 90
                    else:
                        original_angle = 180

                    # Try setting the angle to 180 degrees
                    self.set_branch_angle(junction, junction.children[0], 180)
                    self.update_unpaired_strands_positions(junction)

                    # Check if this change increases the overlap count
                    new_overlap_count = self.update_overlap_count()
                    if new_overlap_count > global_best_overlap:
                        # Revert to original angle if overlap increases
                        self.set_branch_angle(junction, junction.children[0], original_angle)
                        self.update_unpaired_strands_positions(junction)
                        
        last_best = None
        ovp = explore_paths()

        while ovp > 0 and last_best is not ovp:
            last_best = ovp
            ovp = explore_paths()

        for i in range(0,5):
            straighten_branches(ovp)

        self.add_horizontal_distance_between(between_strands)

        overlap_count = self.update_overlap_count(draw=True)

        print("Final Overlap:", overlap_count)
        print(last_best)

        # afterwards, go through paths with overlap still occuring
        # increase radius for those, attempting to see if that fixes the issue with the angle checking?

        #self.set_branch_angle(self.struct.get_junctions()[1], self.struct.get_junctions()[1].children[3], 90)
        #self.set_branch_angle(self.struct.get_junctions()[1], self.struct.get_junctions()[1].children[2], 135)
        #self.set_branch_angle(self.struct.get_junctions()[1], self.struct.get_junctions()[1].children[1], 225)
        #self.set_branch_angle(self.struct.get_junctions()[1], self.struct.get_junctions()[1].children[0], 270)

        #self.update_unpaired_strands_positions(self.struct.get_junctions()[1])

        if False:
            self.minimize_horizontal_distance_between(between_strands)
        #print("overlap", overlap, nodes_below_straight_strand)

        # Sets up the Information to Connect Junction Nodes
        self.junction_data = []
        self.between_strands = between_strands

        for junction in self.struct.get_junctions():
            center_x, center_y = self.get_junction_center(junction)
            radius = self.get_junction_radius(junction)
            
            for strand in junction.strands:
                nodes_data = []
                start_node, end_node = strand[0], strand[-1]

                start_angle = np.arctan2(self.yarray[start_node] - center_y, self.xarray[start_node] - center_x)
                end_angle = np.arctan2(self.yarray[end_node] - center_y, self.xarray[end_node] - center_x)

                def normalize_angle(start_angle, end_angle):
                    if start_angle < 0:
                        start_angle += 2 * np.pi
                    if end_angle < 0:
                        end_angle += 2 * np.pi
                    if start_angle > end_angle:
                        end_angle += 2 * np.pi

                    start_angle = start_angle % (2 * np.pi)
                    end_angle = end_angle % (2 * np.pi)

                    # Choose the shorter angle
                    if end_angle - start_angle > np.pi:
                        end_angle -= 2 * np.pi
                    elif end_angle - start_angle < -np.pi:
                        end_angle += 2 * np.pi

                    # Reverse the angle for Matplotlib's clockwise drawing
                    start_angle, end_angle = end_angle, start_angle

                    return start_angle, end_angle

                start_angle, end_angle = normalize_angle(start_angle, end_angle)

                a,b = self.get_last_nucleotides(junction.parent.strands)
                
                for node in strand:
                    x = self.xarray[node]
                    y = self.yarray[node]
                    nodes_data.append({'Node': node, 'X': x, 'Y': y})
                    
                self.junction_data.append({
                    'Junction': junction,
                    'Strand': strand,
                    'Center': (center_x, center_y),
                    'First': a,
                    'Last': b,
                    'Radius': radius,
                    'Nodes': nodes_data,
                    'Start Angle': start_angle,
                    'End Angle': end_angle
                })

        xarray = self.xarray
        yarray = self.yarray

        ## Below Is From Prior

        min_x = xarray[0] - NODE_R
        min_y = yarray[0] - NODE_R
        max_x = xarray[0] + NODE_R
        max_y = xarray[0] + NODE_R

        for x in xarray:
            if x - NODE_R < min_x:
                min_x = x - NODE_R
            if x + NODE_R > max_x:
                max_x = x + NODE_R

        for y in yarray:
            if y - NODE_R < min_y:
                min_y = y - NODE_R
            if y + NODE_R > max_y:
                max_y = y + NODE_R

        for ii in range(0, len(xarray)):
            xarray[ii] -= min_x
            yarray[ii] -= min_y

        self.size_ = [max_x - min_x, max_y - min_y]
        self.xarray_ = xarray
        self.yarray_ = yarray

        return overlap_count

    def get_size(self):
        return self.size_

    def draw(
        self, offset_x, offset_y, colors, pairs, sequence, render_in_letter, line=False
    ):
        for strand_info in self.junction_data:
            radius = strand_info['Radius']
            start_angle = np.degrees(strand_info['Start Angle'])
            end_angle = np.degrees(strand_info['End Angle'])

            if len(strand_info['Strand']) > 2 or len(strand_info['Junction'].children) > 1:
                first = strand_info['First']
                last = strand_info['Last']

                mid_x = (self.xarray[first] + self.xarray[last]) / 2
                mid_y = (self.yarray[first] + self.yarray[last]) / 2

                direction_x = self.xarray[last] - self.xarray[first]
                direction_y = self.yarray[last] - self.yarray[first]

                magnitude = np.sqrt(direction_x**2 + direction_y**2)
                unit_x = direction_x / magnitude
                unit_y = direction_y / magnitude

                perp_x = -unit_y
                perp_y = unit_x

                # Fixes issue when flipping across x axis.
                junction = strand_info['Junction']
                centroid_x = sum(self.xarray[node] for node in junction.positions) / len(junction.positions)
                centroid_y = sum(self.yarray[node] for node in junction.positions) / len(junction.positions)
                
                direction_to_centroid_x = centroid_x - self.xarray[first]
                direction_to_centroid_y = centroid_y - self.yarray[first]
                
                dot_product = perp_x * direction_to_centroid_x + perp_y * direction_to_centroid_y
                
                if dot_product < 0:
                    perp_x, perp_y = -perp_x, -perp_y

                center_x = mid_x + radius * perp_x + offset_x
                center_y = mid_y + radius * perp_y + offset_y

                arc = Arc(
                    (center_x, center_y),
                    2 * radius,
                    2 * radius,
                    angle=0,
                    theta1=start_angle,
                    theta2=end_angle,
                    linewidth=7.5,
                    edgecolor="#969696",
                )
                self.ax.add_patch(arc)
            else:
                node_1 = strand_info['Strand'][0]
                node_2 = strand_info['Strand'][1]

                x1, y1 = self.xarray_[node_1] + offset_x, self.yarray_[node_1] + offset_y
                x2, y2 = self.xarray_[node_2] + offset_x, self.yarray_[node_2] + offset_y

                rec = ConnectionPatch(
                    (x1, y1),
                    (x2, y2),
                    coordsA="data",
                    linewidth=7.5,
                    edgecolor="#969696",
                )
                self.ax.add_patch(rec)

        if self.between_strands:
            for between_strand in self.between_strands[0:]:
                rec = ConnectionPatch(
                    (self.xarray[between_strand[0]-1] + offset_x, self.yarray[between_strand[0]-1] + offset_y),
                    (self.xarray[between_strand[0]] + offset_x, self.yarray[between_strand[0]] + offset_y),
                    coordsA="data",
                    linewidth=7.5,
                    edgecolor="#969696",
                )
                #Draws Horizontal Line between. This is where extra check should occur.
                self.ax.add_patch(rec)
        '''
        for between_strand in self.between_strands[0:]:
            # We are drawing from start_node to end_node.
            start_node = between_strand[0] - 1
            end_node = start_node + 1

            # Check if Straight Line Should be Drawn
            #not self.between_strands[0] == between_strand and
            if min(self.yarray[between_strand[0] + 1:between_strand[1] - 1]) < self.yarray[0] + self.NODE_R:
                print('There is something below, needs to draw different')
                line_distance_offset_y = self.NODE_R * 4

                #Left Line Vertical
                rec = ConnectionPatch(
                    (self.xarray[start_node] + offset_x, self.yarray[start_node] + offset_y),
                    (self.xarray[start_node] + offset_x, min(self.yarray[between_strand[0] + 1:between_strand[1] - 1]) + offset_y - line_distance_offset_y),
                    coordsA="data",
                    linewidth=7.5,
                    edgecolor="#969696",
                )
                self.ax.add_patch(rec)

                #Right Line Vertical
                rec = ConnectionPatch(
                    (self.xarray[end_node] + offset_x, self.yarray[end_node] + offset_y),
                    (self.xarray[end_node] + offset_x, min(self.yarray[between_strand[0] + 1:between_strand[1] - 1]) + offset_y - line_distance_offset_y),
                    coordsA="data",
                    linewidth=7.5,
                    edgecolor="#969696",
                )
                self.ax.add_patch(rec)

                #Horizontal Connection Between Vertical Lines
                rec = ConnectionPatch(
                    (self.xarray[start_node] + offset_x, min(self.yarray[between_strand[0] + 1:between_strand[1] - 1]) + offset_y - line_distance_offset_y),
                    (self.xarray[end_node] + offset_x, min(self.yarray[between_strand[0] + 1:between_strand[1] - 1]) + offset_y - line_distance_offset_y),
                    coordsA="data",
                    linewidth=7.5,
                    edgecolor="#969696",
                )
                self.ax.add_patch(rec)
            else:
                rec = ConnectionPatch(
                    (self.xarray[between_strand[0]-1] + offset_x, self.yarray[between_strand[0]-1] + offset_y),
                    (self.xarray[between_strand[0]] + offset_x, self.yarray[between_strand[0]] + offset_y),
                    coordsA="data",
                    linewidth=7.5,
                    edgecolor="#969696",
                )
                #Draws Horizontal Line between. This is where extra check should occur.
                self.ax.add_patch(rec)

            rec = ConnectionPatch(
                (self.xarray[between_strand[0]] + offset_x, self.yarray[between_strand[0]] + offset_y),
                (self.xarray[between_strand[0] + 1] + offset_x, self.yarray[between_strand[0] + 1] + offset_y),
                coordsA="data",
                linewidth=7.5,
                edgecolor="#969696",
            )
            # I forget what this was for.
            #self.ax.add_patch(rec)

            rec = ConnectionPatch(
                (self.xarray[between_strand[1]] + offset_x, self.yarray[between_strand[1]] + offset_y),
                (self.xarray[between_strand[1] - 1] + offset_x, self.yarray[between_strand[1] - 1] + offset_y),
                coordsA="data",
                linewidth=7.5,
                edgecolor="#969696",
            )
            # I forget what this was for.
            #self.ax.add_patch(rec)
        '''
        if self.xarray_ != None:
            if line:
                for ii in range(len(self.xarray_) - 1):
                    if colors == None:
                        pass
                    else:
                        pass
            else:
                if pairs:
                    for pair in pairs:
                        x1, y1 = (
                            [
                                offset_x + self.xarray_[pair["from"]],
                                offset_y + self.yarray_[pair["from"]],
                            ],
                            [
                                offset_x + self.xarray_[pair["to"]],
                                offset_y + self.yarray_[pair["to"]],
                            ],
                        )
                        rec = ConnectionPatch(
                            (x1[0], x1[1]),
                            (y1[0], y1[1]),
                            coordsA="data",
                            linewidth=15,
                            edgecolor="#969696",
                        )
                        self.ax.add_patch(rec)

                if not render_in_letter:
                    for ii in range(0, len(self.xarray_)):
                        if colors == None:
                            x = self.xarray_[ii] + offset_x
                            y = self.yarray_[ii] + offset_y
                            cir = Circle(
                                (x, y), radius=self.NODE_R, facecolor="k", edgecolor="k"
                            )
                            self.ax.add_patch(cir)
                        else:
                            x = self.xarray_[ii] + offset_x
                            y = self.yarray_[ii] + offset_y
                            cir = Circle(
                                (x, y),
                                radius=self.NODE_R,
                                facecolor=(colors[ii][0], colors[ii][1], colors[ii][2]),
                                edgecolor="k",
                            )
                            self.ax.add_patch(cir)
                if sequence:
                    for ii in range(0, len(self.xarray_)):
                        if not render_in_letter:
                            text_size = 20
                            if colors[ii] == [0, 0, 0]:
                                color = "w"
                            else:
                                color = "k"
                            text_offset_x = 0
                            text_offset_y = 0
                        else:
                            if colors == None:
                                color = "k"
                            else:
                                color = colors[ii]
                            text_size = 20
                            text_offset_x = 0
                            text_offset_y = 0
                        
                        self.ax.text(
                            self.xarray_[ii] + offset_x + text_offset_x,
                            self.yarray_[ii] + offset_y + text_offset_y,
                            sequence[ii],
                            family="monospace",
                            fontsize=text_size,
                            ha="center",
                            va="center",
                        )

    def get_coords(self, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE):
        if self.root_ != None:
            get_coords_recursive(self.root_, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)
        else:
            for ii in range(0, len(xarray)):
                xarray[ii] = 0
                yarray[ii] = ii * PRIMARY_SPACE

    def setup_coords(self, NODE_R, PRIMARY_SPACE, PAIR_SPACE):
        if self.root_ != None:
            radii = [0] * self.length

            return setup_coords_recursive(
                self.root_, None, 0, 0, 0, 1, NODE_R, PRIMARY_SPACE, PAIR_SPACE, radii
            )