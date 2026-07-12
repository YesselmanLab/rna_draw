import sys
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch


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
        elif rootnode.children_[0].is_pair_ == False and rootnode.children_[0].index_a_ < 0:
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

        circle_length = (len(rootnode.children_) + 1) * PRIMARY_SPACE + (npairs + 1) * PAIR_SPACE
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
        get_coords_recursive(rootnode.children_[ii], xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)


class RNARenderer:
    def __init__(self):
        self.root_ = None
        self.xarray_ = None
        self.yarray_ = None
        self.size_ = None
        self.fig = plt.Figure()
        self.ax = self.fig.add_subplot(111, aspect="equal")

    def setup_tree(self, secstruct, NODE_R, PRIMARY_SPACE, PAIR_SPACE):
        bi_pairs = get_pairmap_from_secstruct(secstruct)

        self.NODE_R = NODE_R
        self.root_ = None

        # NOTE: a dangling_start/dangling_end counting block previously lived
        # here. Those variables were written but never read anywhere (the
        # second loop also used a malformed `for ii in (n, -1, -1)` tuple
        # instead of `range(n, -1, -1)`), so it was deleted as dead code; this
        # is a provable no-op on the coordinate output (see
        # tests/test_layout_baseline.py).
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

        self.setup_coords(NODE_R, PRIMARY_SPACE, PAIR_SPACE)
        self.get_coords(xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)

        self.set_coords(xarray, yarray, NODE_R)

    def set_coords(self, xarray, yarray, NODE_R):
        """Inject externally computed coords; compute bounds/shift/size like setup_tree does.

        Lets an external `rna_draw.layout` engine feed its own
        coordinates through this renderer without going through
        `setup_tree`'s tree-recursion layout.

        Args:
            xarray: Nucleotide x-coordinates, one per nucleotide.
            yarray: Nucleotide y-coordinates, one per nucleotide.
            NODE_R: Nucleotide disk radius (also the padding used when
                computing the shifted bounding box).
        """
        self.NODE_R = NODE_R

        min_x = xarray[0] - NODE_R
        min_y = yarray[0] - NODE_R
        max_x = xarray[0] + NODE_R
        max_y = yarray[0] + NODE_R

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

    def draw(self, offset_x, offset_y, colors, pairs, sequence, render_in_letter, line=False):
        if self.xarray_ != None:
            if line:
                # TODO(M1): 'line' backbone-drawing mode is unimplemented and
                # currently unreachable (draw() is only ever called with the
                # default line=False in the shipped pipeline). Fail loudly
                # instead of silently drawing nothing.
                raise NotImplementedError("line rendering mode not implemented")
            else:
                if pairs:
                    # NOTE: endpoints verified from->to; geometry intentionally
                    # unchanged in M1 (no pixel baseline exists yet).
                    for pair in pairs:
                        from_xy = [
                            offset_x + self.xarray_[pair["from"]],
                            offset_y + self.yarray_[pair["from"]],
                        ]
                        to_xy = [
                            offset_x + self.xarray_[pair["to"]],
                            offset_y + self.yarray_[pair["to"]],
                        ]
                        rec = ConnectionPatch(
                            (from_xy[0], from_xy[1]),
                            (to_xy[0], to_xy[1]),
                            coordsA="data",
                            linewidth=15,
                            edgecolor=pair.get("color", "#969696"),
                        )
                        self.ax.add_patch(rec)

                if not render_in_letter:
                    for ii in range(0, len(self.xarray_)):
                        if colors == None:
                            x = self.xarray_[ii] + offset_x
                            y = self.yarray_[ii] + offset_y
                            cir = Circle((x, y), radius=self.NODE_R, facecolor="k", edgecolor="k")
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

    def draw_routed_lines(self, lines, offset_x, offset_y, color, linewidth=15):
        """Stroke each PK-A `RoutedLine` polyline (`rna_draw.layout.pseudoknot`).

        Reuses the same `ConnectionPatch` path `draw()` uses for a
        straight pair, one segment per consecutive pair of `line.points`,
        so a bowed/ring-routed crossing connector renders with the same
        stroke style as an ordinary base pair, just in `color` and
        possibly multi-segment.

        Args:
            lines: `RoutedLine`s to draw (points already in the same
                pre-offset coordinate space as `self.xarray_`/`yarray_`).
            offset_x: Same x offset `draw()`'s pairs/disks use.
            offset_y: Same y offset `draw()`'s pairs/disks use.
            color: Matplotlib edgecolor for every segment.
            linewidth: Stroke width; matches `draw()`'s pair connectors.
        """
        for line in lines:
            for (x0, y0), (x1, y1) in zip(line.points, line.points[1:]):
                patch = ConnectionPatch(
                    (offset_x + x0, offset_y + y0),
                    (offset_x + x1, offset_y + y1),
                    coordsA="data",
                    linewidth=linewidth,
                    edgecolor=color,
                )
                self.ax.add_patch(patch)

    def get_coords(self, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE):
        if self.root_ != None:
            get_coords_recursive(self.root_, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)
        else:
            for ii in range(0, len(xarray)):
                xarray[ii] = 0
                yarray[ii] = ii * PRIMARY_SPACE

    def setup_coords(self, NODE_R, PRIMARY_SPACE, PAIR_SPACE):
        if self.root_ != None:
            setup_coords_recursive(self.root_, None, 0, 0, 0, 1, NODE_R, PRIMARY_SPACE, PAIR_SPACE)
