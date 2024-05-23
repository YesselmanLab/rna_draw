import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from rna_draw import parameters

class Rectangle:
    """
    Class to represent a rectangle with methods to update its bounding box and node positions.
    """

    def __init__(self, nodes):
        """
        Initialize a Rectangle with a list of nodes.

        :param nodes: List of (x, y) coordinates representing the corners of the rectangle.
        """
        self.nodes = nodes
        self.update_bounding_box()

    def update_bounding_box(self):
        """
        Update the bounding box of the rectangle based on its nodes.
        """
        min_x = min(node[0] for node in self.nodes)
        min_y = min(node[1] for node in self.nodes)
        max_x = max(node[0] for node in self.nodes)
        max_y = max(node[1] for node in self.nodes)

        self.bounding_box = {
            'top_left': (min_x, min_y),
            'bottom_right': (max_x, max_y)
        }

    def update_node_position(self, node_index, new_x, new_y):
        """
        Update the position of a node and refresh the bounding box.

        :param node_index: Index of the node to be updated.
        :param new_x: New x-coordinate of the node.
        :param new_y: New y-coordinate of the node.
        """
        self.nodes[node_index][0] = new_x
        self.nodes[node_index][1] = new_y
        self.update_bounding_box()

class Circle:
    """
    Class to represent a circle with methods to update its position.
    """

    def __init__(self, x, y, radius):
        """
        Initialize a Circle with its center coordinates and radius.

        :param x: x-coordinate of the circle's center.
        :param y: y-coordinate of the circle's center.
        :param radius: Radius of the circle.
        """
        self.x = x
        self.y = y
        self.radius = radius

    def update_position(self, new_x, new_y):
        """
        Update the position of the circle's center.

        :param new_x: New x-coordinate of the circle's center.
        :param new_y: New y-coordinate of the circle's center.
        """
        self.x = new_x
        self.y = new_y

class ChunkContainer:
    """
    Container for handling chunks of RNA secondary structure elements, with methods for visualization and overlap checking.
    """

    def __init__(self, struct, xarray, yarray, node_array=None, draw=False, draw_num=0):
        """
        Initialize the ChunkContainer with RNA structure elements and their positions.

        :param struct: RNA secondary structure.
        :param xarray: Array of x-coordinates of the RNA nodes.
        :param yarray: Array of y-coordinates of the RNA nodes.
        :param node_array: Optional array of node indices to consider.
        :param draw: Flag to enable drawing of the chunks.
        :param draw_num: Drawing number for visualization filename.
        """
        self.chunks = {}
        self.struct = struct
        self.draw = draw
        self.node_array = node_array or [t for t in range(0, len(xarray))]
        self.draw_parms = parameters.DrawParameters()
        self.draw_num = draw_num

        self.chunks_before = {}
        self.chunks_in_node_array = {}

        for m in struct:
            if m in struct.get_junctions() and hasattr(m, 'center_x') and hasattr(m, 'center_y') and hasattr(m, 'radius'):
                nodes = m.positions
                circle = Circle(m.center_x, m.center_y, m.radius)

                if any(node in self.node_array for node in nodes):
                    self.add_chunk(circle, m)
                    
                if m.positions[0] < self.node_array[0]:
                    self.chunks_before[circle] = m
                elif m.positions[0] < self.node_array[-1]:
                    self.chunks_in_node_array[circle] = m
            else:
                nodes = []
                for node in m.positions:
                    nodes.append((xarray[node], yarray[node]))
                rect = Rectangle(nodes)
                self.add_chunk(rect, m)

                if any(node in self.node_array for node in nodes):
                    self.add_chunk(rect, m)

                if m.positions[0] < self.node_array[0]:
                    self.chunks_before[rect] = m
                elif m.positions[0] < self.node_array[-1]:
                    self.chunks_in_node_array[rect] = m

    def fitted_size_function(self, area):
        """
        Calculate the figure size based on the area of the chunks.

        :param area: Area of the bounding box containing all chunks.
        :return: Figure size for plotting.
        """
        base_size = 10
        return base_size + np.log1p(area) * 0.05

    def visualize_chunks(self, overlapping_chunks=[]):
        """
        Visualize the chunks, highlighting any overlapping chunks.

        :param overlapping_chunks: List of chunks that are overlapping.
        """
        filename = f'DrawVisualizer/chunks_visualization_{self.draw_num}.png'

        all_x = []
        all_y = []
        for chunk in self.chunks.keys():
            if isinstance(chunk, Circle):
                all_x.extend([chunk.x - chunk.radius, chunk.x + chunk.radius])
                all_y.extend([chunk.y - chunk.radius, chunk.y + chunk.radius])
            elif isinstance(chunk, Rectangle):
                all_x.extend([chunk.bounding_box['top_left'][0], chunk.bounding_box['bottom_right'][0]])
                all_y.extend([chunk.bounding_box['top_left'][1], chunk.bounding_box['bottom_right'][1]])
        
        if not all_x or not all_y:
            # No shapes to visualize.
            return

        min_x, max_x = min(all_x), max(all_x)
        min_y, max_y = min(all_y), max(all_y)
        area = (max_x - min_x) * (max_y - min_y)

        fig_size = self.fitted_size_function(area)

        linewidth = max(1, fig_size / 100)

        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        for chunk in self.chunks.keys():
            color = 'g'
            if chunk in overlapping_chunks:
                color = 'r'
            if isinstance(chunk, Circle):
                circle = patches.Circle((chunk.x, chunk.y), chunk.radius, edgecolor=color, facecolor='none', linewidth=linewidth)
                ax.add_patch(circle)
            elif isinstance(chunk, Rectangle):
                # Use the actual node positions to define the polygon vertices
                vertices = chunk.nodes  # Assuming `chunk.nodes` holds the vertices in order
                polygon = patches.Polygon(vertices, closed=True, edgecolor=color, facecolor='none', linewidth=linewidth)
                ax.add_patch(polygon)

        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
        ax.set_aspect('equal', 'box')
        plt.axis('off')
        ax.autoscale_view()

        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()

    def add_chunk(self, chunk, motif):
        """
        Add a chunk to the container.

        :param chunk: The chunk (Circle or Rectangle) to add.
        :param motif: The associated RNA motif.
        """
        self.chunks[chunk] = motif

    def circle_rectangle_overlap(self, circle, rectangle):
        """
        Check if a circle overlaps with a rectangle.

        :param circle: The Circle object.
        :param rectangle: The Rectangle object.
        :return: True if there is an overlap, False otherwise.
        """
        closest_x = max(rectangle.bounding_box['top_left'][0], min(circle.x, rectangle.bounding_box['bottom_right'][0]))
        closest_y = max(rectangle.bounding_box['top_left'][1], min(circle.y, rectangle.bounding_box['bottom_right'][1]))

        distance_x = circle.x - closest_x
        distance_y = circle.y - closest_y

        distance_squared = distance_x**2 + distance_y**2
        return distance_squared < circle.radius**2

    def circles_overlap(self, circle1, circle2):
        """
        Check if two circles overlap.

        :param circle1: The first Circle object.
        :param circle2: The second Circle object.
        :return: True if there is an overlap, False otherwise.
        """
        dx = circle1.x - circle2.x
        dy = circle1.y - circle2.y
        distance_squared = dx**2 + dy**2
        radius_sum = circle1.radius + circle2.radius + self.draw_parms.NODE_R * 2

        return round(distance_squared, 2) < round(radius_sum**2, 2)

    def rectangles_overlap(self, rect1, rect2):
        """
        Check if two rectangles overlap.

        :param rect1: The first Rectangle object.
        :param rect2: The second Rectangle object.
        :return: True if there is an overlap, False otherwise.
        """
        bb1 = rect1.bounding_box
        bb2 = rect2.bounding_box

        return not (bb1['top_left'][0] > bb2['bottom_right'][0] or 
                    bb1['bottom_right'][0] < bb2['top_left'][0] or 
                    bb1['top_left'][1] > bb2['bottom_right'][1] or 
                    bb1['bottom_right'][1] < bb2['top_left'][1])

    def check_specific_overlap(self):
        """
        Check chunks in given node_array for overlaps with chunks that come after node_array.

        :return: Number of overlaps found.
        """
        overlaps = 0
        overlapping_chunks = []

        for chunk1, chunk1_key in self.chunks_before.items():
            for chunk2, chunk2_key in self.chunks_in_node_array.items():

                if not ((chunk1_key in self.chunks and chunk2_key in self.chunks[chunk1_key].children) or 
                        (chunk2_key in self.chunks and chunk1_key in self.chunks[chunk2_key].children)):
                    if isinstance(chunk1, Rectangle) and isinstance(chunk2, Rectangle):
                        if self.rectangles_overlap(chunk1, chunk2):
                            overlapping_chunks.append(chunk1_key)
                            overlapping_chunks.append(chunk2_key)
                            overlaps += 1
                    elif isinstance(chunk1, Circle) and isinstance(chunk2, Circle):
                        if self.circles_overlap(chunk1, chunk2):
                            overlaps += 1
                            overlapping_chunks.append(chunk1_key)
                            overlapping_chunks.append(chunk2_key)
                    elif isinstance(chunk1, Rectangle) and isinstance(chunk2, Circle):
                        if self.circle_rectangle_overlap(chunk2, chunk1):
                            overlaps += 1
                            overlapping_chunks.append(chunk1_key)
                            overlapping_chunks.append(chunk2_key)
                    elif isinstance(chunk1, Circle) and isinstance(chunk2, Rectangle):
                        if self.circle_rectangle_overlap(chunk1, chunk2):
                            overlaps += 1
                            overlapping_chunks.append(chunk1_key)
                            overlapping_chunks.append(chunk2_key)

        if self.draw:
            self.visualize_chunks(overlapping_chunks=overlapping_chunks)

        return overlaps

    def check_any_overlap_in_node_array(self):
        """
        Check chunks in node_array for overlaps with all other chunks in node_array.

        :return: Number of overlaps found.
        """
        overlaps = 0

        chunk_keys = list(self.chunks_in_node_array.keys())

        overlapping_chunks = []

        for i in range(len(chunk_keys)):
            for j in range(i + 1, len(chunk_keys)):
                chunk1 = chunk_keys[i]
                chunk2 = chunk_keys[j]

                if not (self.chunks_in_node_array[chunk2] in self.chunks_in_node_array[chunk1].children or self.chunks_in_node_array[chunk1] in self.chunks_in_node_array[chunk2].children):
                    if len(self.chunks_in_node_array[chunk2].parent.positions) > 2 and len(self.chunks_in_node_array[chunk1].positions) + len(self.chunks_in_node_array[chunk2].positions) > 2:
                        if isinstance(chunk1, Rectangle) and isinstance(chunk2, Rectangle):
                            if self.rectangles_overlap(chunk1, chunk2):
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)
                                overlaps += 1
                        elif isinstance(chunk1, Circle) and isinstance(chunk2, Circle):
                            if self.circles_overlap(chunk1, chunk2):
                                overlaps += 1
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)
                        elif isinstance(chunk1, Rectangle) and isinstance(chunk2, Circle):
                            if self.circle_rectangle_overlap(chunk2, chunk1):
                                overlaps += 1
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)
                        elif isinstance(chunk1, Circle) and isinstance(chunk2, Rectangle):
                            if self.circle_rectangle_overlap(chunk1, chunk2):
                                overlaps += 1
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)

        if self.draw:
            self.visualize_chunks(overlapping_chunks=overlapping_chunks)

        return overlaps

    def check_any_overlap(self):
        """
        Check for any overlaps among all chunks.

        :return: Number of overlaps found.
        """
        overlaps = 0

        chunk_keys = list(self.chunks.keys())

        overlapping_chunks = []

        for i in range(len(chunk_keys)):
            for j in range(i + 2, len(chunk_keys)):
                chunk1 = chunk_keys[i]
                chunk2 = chunk_keys[j]

                if not (self.chunks[chunk2] in self.chunks[chunk1].children or self.chunks[chunk1] in self.chunks[chunk2].children):
                    if len(self.chunks[chunk2].parent.positions) > 2 and len(self.chunks[chunk1].positions) + len(self.chunks[chunk2].positions) > 2:
                        if isinstance(chunk1, Rectangle) and isinstance(chunk2, Rectangle):
                            if self.rectangles_overlap(chunk1, chunk2):
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)
                                overlaps += 1
                        elif isinstance(chunk1, Circle) and isinstance(chunk2, Circle):
                            if self.circles_overlap(chunk1, chunk2):
                                overlaps += 1
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)
                        elif isinstance(chunk1, Rectangle) and isinstance(chunk2, Circle):
                            if self.circle_rectangle_overlap(chunk2, chunk1):
                                overlaps += 1
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)
                        elif isinstance(chunk1, Circle) and isinstance(chunk2, Rectangle):
                            if self.circle_rectangle_overlap(chunk1, chunk2):
                                overlaps += 1
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)

        if self.draw:
            self.visualize_chunks(overlapping_chunks=overlapping_chunks)

        return overlaps