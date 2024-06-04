import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch
import numpy as np
from matplotlib.patches import Arc
from rna_draw.rna_optimizer import RNAOptimizer
from rna_secstruct import SecStruct, MotifSearchParams

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
    Generates dictionary containing pair mappings.
    
    Args:
    secstruct (str): Contains secondary structure string.

    Returns:
    list: List with pair mappings.
    """
    pair_stack = []
    end_stack = []
    pairs_array = []
    i_range = range(0, len(secstruct))

    # Initialize all values to -1, meaning no pair.
    for ii in i_range:
        pairs_array.append(-1)

    # Assign pairs based on secstruct.
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
    """
    Recursively adds nodes to the RNA tree structure.

    Args:
    bi_pairs (list): List of paired indices.
    rootnode (RNATreeNode): The root node of the tree.
    start_index (int): Starting index of the current segment.
    end_index (int): Ending index of the current segment.

    Returns:
    None
    """
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
    """
    Recursively sets up the coordinates for the nodes in the RNA tree.

    Args:
    rootnode (RNATreeNode): The current root node.
    parentnode (RNATreeNode): The parent node.
    start_x (float): Starting x-coordinate.
    start_y (float): Starting y-coordinate.
    go_x (float): Direction vector x-component.
    go_y (float): Direction vector y-component.
    NODE_R (float): Radius of each node.
    PRIMARY_SPACE (float): Primary spacing between nodes.
    PAIR_SPACE (float): Spacing between paired nodes.
    radii (list): List to store radii of nodes.

    Returns:
    list: Updated list of radii.
    """
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
    """
    Recursively retrieves coordinates for the nodes in the RNA tree.

    Args:
    rootnode (RNATreeNode): The current root node.
    xarray (list): List to store x-coordinates.
    yarray (list): List to store y-coordinates.
    PRIMARY_SPACE (float): Primary spacing between nodes.
    PAIR_SPACE (float): Spacing between paired nodes.

    Returns:
    None
    """
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
        self.optimizer = RNAOptimizer(self)

    def setup_tree(self, secstruct, NODE_R, PRIMARY_SPACE, PAIR_SPACE, seq):
        """
        Sets up the RNA tree for rendering, initializing all necessary structures and optimizing their layout.

        Args:
        secstruct (str): The secondary structure string of the RNA.
        NODE_R (float): Radius of each node in the tree.
        PRIMARY_SPACE (float): Primary spacing between nodes.
        PAIR_SPACE (float): Spacing between paired nodes.
        seq (str): The RNA sequence.

        Returns:
        int: The count of overlapping nodes.
        """
        print("Initializing tree...")
        self.initialize_tree(secstruct, NODE_R)
        print("Setting up coordinates...")
        self.setup_coordinates(NODE_R, PRIMARY_SPACE, PAIR_SPACE)
        print("Initializing junctions...")
        self.initialize_junctions(seq, secstruct, PRIMARY_SPACE, PAIR_SPACE)
        print("Optimizing structure...")
        self.optimize_structure()
        final_overlap_count = self.optimizer.update_overlap_count()
        print("Preparing drawing...")
        self.prepare_drawing(NODE_R)
        return final_overlap_count

    def initialize_tree(self, secstruct, NODE_R):
        """
        Initializes the RNA tree based on the given secondary structure.

        Args:
        secstruct (str): The secondary structure string of the RNA.
        NODE_R (float): Radius of each node in the tree.

        Returns:
        None
        """
        self.NODE_R = NODE_R
        bi_pairs = get_pairmap_from_secstruct(secstruct)
        self.length = len(bi_pairs)
        self.root_ = RNATreeNode()

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

    def setup_coordinates(self, NODE_R, PRIMARY_SPACE, PAIR_SPACE):
        """
        Sets up the coordinates for the nodes in the RNA tree.

        Args:
        NODE_R (float): Radius of each node in the tree.
        PRIMARY_SPACE (float): Primary spacing between nodes.
        PAIR_SPACE (float): Spacing between paired nodes.

        Returns:
        list: List of radii for the nodes.
        """
        xarray = [0.0] * self.length
        yarray = [0.0] * self.length

        radii = self.setup_coords(NODE_R, PRIMARY_SPACE, PAIR_SPACE)
        self.get_coords(xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)

        self.xarray = xarray
        self.yarray = yarray
        return radii

    def initialize_junctions(self, seq, secstruct, PRIMARY_SPACE, PAIR_SPACE):
        """
        Initializes the junctions in the RNA structure.

        Args:
        seq (str): The RNA sequence.
        secstruct (str): The secondary structure string of the RNA.
        PRIMARY_SPACE (float): Primary spacing between nodes.
        PAIR_SPACE (float): Spacing between paired nodes.

        Returns:
        None
        """
        self.struct = SecStruct(seq, secstruct)
        radii = self.setup_coordinates(self.NODE_R, PRIMARY_SPACE, PAIR_SPACE)
        for junction in self.struct.get_junctions():
            for strand in junction.strands:
                for node in strand[1:-1]:
                    junction.radius = radii[node]
                    self.optimizer.get_junction_center(junction)

    def optimize_structure(self):
        """
        Optimizes the RNA structure to minimize overlap and prepare it for drawing.

        Args:
        None

        Returns:
        None
        """
        self.between_strands = self.optimizer.straighten_unpaired_strands()
        if self.between_strands:
            self.optimizer.add_horizontal_distance_between(self.between_strands)

        x_range = max(self.xarray) - min(self.xarray)
        y_range = max(self.yarray) - min(self.yarray)

        estimated_inches_x, estimated_inches_y = self.optimizer.estimate_inches_using_curve_fit(x_range, y_range)

        if estimated_inches_x > 650 or estimated_inches_y > 650:
            print("Structure is too large to output using matplotlib.")
            print(f'Estimated Inches Required (X): {estimated_inches_x}')
            print(f'Estimated Inches Required (Y): {estimated_inches_y}')
            print(f'Structure BP Length: {len(self.xarray)}')
            raise Exception('Too Large')

        self.global_best_overlap = {}
        self.global_best_combo = {}

        last_best = None
        ovp = self.optimizer.explore_paths()
        print(f'overlap check #1: {ovp}')

        while ovp > 0 and last_best != ovp:
            last_best = ovp
            ovp = self.optimizer.explore_paths()
            print(f'last overlap best: {ovp}')

        for _ in range(5):
            self.optimizer.straighten_branches(ovp)

        self.optimizer.add_horizontal_distance_between(self.between_strands)

    def prepare_drawing(self, NODE_R):
        """
        Prepares the RNA structure for drawing by calculating necessary coordinates and angles.

        Args:
        NODE_R (float): Radius of each node in the tree.

        Returns:
        None
        """
        self.junction_data = []
        for junction in self.struct.get_junctions():
            center_x, center_y = self.optimizer.get_junction_center(junction)
            radius = self.optimizer.get_junction_radius(junction)

            for strand in junction.strands:
                nodes_data = []
                start_node, end_node = strand[0], strand[-1]
                start_angle = np.arctan2(self.yarray[start_node] - center_y, self.xarray[start_node] - center_x)
                end_angle = np.arctan2(self.yarray[end_node] - center_y, self.xarray[end_node] - center_x)

                start_angle, end_angle = self.normalize_angle(start_angle, end_angle)
                first, last = self.optimizer.get_last_nucleotides(junction.parent.strands)

                for node in strand:
                    x = self.xarray[node]
                    y = self.yarray[node]
                    nodes_data.append({'Node': node, 'X': x, 'Y': y})

                self.junction_data.append({
                    'Junction': junction,
                    'Strand': strand,
                    'Center': (center_x, center_y),
                    'First': first,
                    'Last': last,
                    'Radius': radius,
                    'Nodes': nodes_data,
                    'Start Angle': start_angle,
                    'End Angle': end_angle
                })

        min_x, max_x = min(self.xarray) - NODE_R, max(self.xarray) + NODE_R
        min_y, max_y = min(self.yarray) - NODE_R, max(self.yarray) + NODE_R

        self.xarray = [x - min_x for x in self.xarray]
        self.yarray = [y - min_y for y in self.yarray]
        self.size_ = [max_x - min_x, max_y - min_y]

        self.xarray_ = self.xarray
        self.yarray_ = self.yarray

    def normalize_angle(self, start_angle, end_angle):
        """
        Normalizes the start and end angles to ensure they are within the correct range for drawing.

        Args:
        start_angle (float): The starting angle in radians.
        end_angle (float): The ending angle in radians.

        Returns:
        tuple: The normalized end and start angles in radians.
        """
        if start_angle < 0:
            start_angle += 2 * np.pi
        if end_angle < 0:
            end_angle += 2 * np.pi
        if start_angle > end_angle:
            end_angle += 2 * np.pi

        start_angle %= (2 * np.pi)
        end_angle %= (2 * np.pi)

        if end_angle - start_angle > np.pi:
            end_angle -= 2 * np.pi
        elif end_angle - start_angle < -np.pi:
            end_angle += 2 * np.pi

        return end_angle, start_angle  # Reversed for Matplotlib's clockwise drawing

    def get_size(self):
        return self.size_

    def draw(
        self, offset_x, offset_y, colors, pairs, sequence, render_in_letter, line=False
    ):
        """
        Draws the RNA structure using Matplotlib.

        Args:
        offset_x (float): X-axis offset for the drawing.
        offset_y (float): Y-axis offset for the drawing.
        colors (list): List of colors for each nucleotide.
        pairs (list): List of base pairs.
        sequence (str): RNA sequence.
        render_in_letter (bool): Flag to render in letter format.
        line (bool): Flag to draw lines instead of nodes.

        Returns:
        None
        """
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
                # Draws Horizontal Line between.
                self.ax.add_patch(rec)
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
        """
        Retrieves coordinates for the nodes in the RNA tree.

        Args:
        xarray (list): List to store x-coordinates.
        yarray (list): List to store y-coordinates.
        PRIMARY_SPACE (float): Primary spacing between nodes.
        PAIR_SPACE (float): Spacing between paired nodes.

        Returns:
        None
        """
        if self.root_ != None:
            get_coords_recursive(self.root_, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)
        else:
            for ii in range(0, len(xarray)):
                xarray[ii] = 0
                yarray[ii] = ii * PRIMARY_SPACE

    def setup_coords(self, NODE_R, PRIMARY_SPACE, PAIR_SPACE):
        """
        Sets up the coordinates for the nodes in the RNA tree.

        Args:
        NODE_R (float): Radius of each node in the tree.
        PRIMARY_SPACE (float): Primary spacing between nodes.
        PAIR_SPACE (float): Spacing between paired nodes.

        Returns:
        list: List of radii for the nodes.
        """
        if self.root_ != None:
            radii = [0] * self.length

            return setup_coords_recursive(
                self.root_, None, 0, 0, 0, 1, NODE_R, PRIMARY_SPACE, PAIR_SPACE, radii
            )
