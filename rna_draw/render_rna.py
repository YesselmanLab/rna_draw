import sys
import re
import random
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch
from matplotlib.axes import Axes


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
    if start_index > end_index:
        print("Error occured while drawing RNA %d %d" % (start_index, end_index))
        sys.exit(0)

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
            )

            if rootnode.children_[ii].is_pair_:
                length_walker += PAIR_SPACE / 2.0

    else:
        rootnode.x_ = start_x
        rootnode.y_ = start_y


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

class Node:
    def __init__(self, index = None, pair = None, x = None, y = None, chunk = None) -> None:
        self.index = index
        self.pair = pair
        self.x = x
        self.y = y
        self.chunk = chunk

    def get_chunk(self):
        return self.chunk

class RNAChunk:
    def __init__(self, chunk_Type):
        self.chunkType = chunk_Type
        self.nodes = []
        self.Connected_Chunks = []

    def add_node(self, node):
        self.nodes.append(node)

    def connect_chunks(self, other_chunk):
        self.Connected_Chunks.append(other_chunk)

    def __str__(self):
        str_repr = "RNA Chunk\n"
        str_repr += f"Type: {self.chunkType}\n"
        str_repr += f"Number of Nodes: {len(self.nodes)}\n"
        str_repr += "Connected to Chunks: \n"

        for chunk in self.Connected_Chunks:
            str_repr += f"- {chunk.chunkType}\n"
        
        return str_repr
    
    def __len__(self) -> int:
        return len(self.nodes)


class RNARenderer:
    def __init__(self):
        self.root_ = None
        self.xarray_ = None
        self.yarray_ = None
        self.size_ = None
        self.fig = plt.Figure()
        self.ax = self.fig.add_subplot(111, aspect="equal")
        self.chunks = []

    def setup_node_structure(self, bi_pairs):
        chunk_type = "Linear"
        current_chunk = None
        chunks = []
        node_registry = [None] * len(bi_pairs)

        for i, pair_node in enumerate(bi_pairs):
            if pair_node == -1:
                # Unpaired
                if chunk_type == "Linear" or chunk_type == None:
                    if chunk_type is not None and current_chunk != None and len(current_chunk) > 0:
                        chunks.append(current_chunk)
                    chunk_type = "Circular"
                    current_chunk = RNAChunk(chunk_type)
            else:
                # Paired
                if chunk_type == "Circular" or chunk_type == None:
                    if chunk_type is not None and current_chunk != None and len(current_chunk) > 0:
                        chunks.append(current_chunk)
                    chunk_type = "Linear"
                    current_chunk = RNAChunk(chunk_type)

                if node_registry[bi_pairs[i]] == None:
                    new_node = Node(bi_pairs[i], chunk = current_chunk)
                    current_chunk.add_node(new_node)
                    node_registry[bi_pairs[i]] = new_node
            
            if node_registry[i] == None:
                new_node = Node(i, chunk = current_chunk)
                current_chunk.add_node(new_node)
                node_registry[i] = new_node
        
        if current_chunk != None and len(current_chunk) > 0:
            chunks.append(current_chunk)
            
        return chunks

    def setup_tree(self, secstruct, NODE_R, PRIMARY_SPACE, PAIR_SPACE):
        dangling_start = 0
        dangling_end = 0
        bi_pairs = get_pairmap_from_secstruct(secstruct)

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

        # New Data Structure Below

        chunks = self.setup_node_structure(bi_pairs)

        for chunk in chunks:
            print(chunk)
            print(" ")
            print(" ")

        # OLD DATA STRUCTURE BELOW

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

        self.setup_coords(NODE_R, PRIMARY_SPACE, PAIR_SPACE)
        self.get_coords(xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)

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

    def get_size(self):
        return self.size_

    def draw(
        self, offset_x, offset_y, colors, pairs, sequence, render_in_letter, line=False
    ):
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
            setup_coords_recursive(
                self.root_, None, 0, 0, 0, 1, NODE_R, PRIMARY_SPACE, PAIR_SPACE
            )
