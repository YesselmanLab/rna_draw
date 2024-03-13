import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from rna_draw import parameters

class Rectangle:
    def __init__(self, nodes):
        self.nodes = nodes
        self.update_bounding_box()

    def update_bounding_box(self):
        min_x = min(node[0] for node in self.nodes)
        min_y = min(node[1] for node in self.nodes)
        max_x = max(node[0] for node in self.nodes)
        max_y = max(node[1] for node in self.nodes)

        self.bounding_box = {
            'top_left': (min_x, min_y),
            'bottom_right': (max_x, max_y)
        }

    def update_node_position(self, node_index, new_x, new_y):
        self.nodes[node_index][0] = new_x
        self.nodes[node_index][1] = new_y
        self.update_bounding_box()

class Circle:
    def __init__(self, x, y, radius):
        self.x = x
        self.y = y
        self.radius = radius

    def update_position(self, new_x, new_y):
        self.x = new_x
        self.y = new_y

class ChunkContainer:
    def __init__(self, struct, xarray, yarray, node_array=None, draw=False, draw_num=0):
        self.chunks = {}
        self.struct = struct
        self.draw = draw
        self.node_array = node_array or [t for t in range(0, len(xarray))]
        self.draw_parms = parameters.DrawParameters()
        self.draw_num = draw_num

        for m in struct:
            if m in struct.get_junctions() and hasattr(m, 'center_x') and hasattr(m, 'center_y') and hasattr(m, 'radius'):
                circle = Circle(m.center_x, m.center_y, m.radius)
                self.add_chunk(circle, m)
            else:
                nodes = []
                for node in m.positions:
                    if self.node_array is not None and node in self.node_array:
                        nodes.append((xarray[node], yarray[node]))
                if len(nodes) > 0:
                    rect = Rectangle(nodes)
                    self.add_chunk(rect, m)

    def fitted_size_function(self, area):
        base_size = 10
        return base_size + np.log1p(area) * 0.05

    def visualize_chunks(self, overlapping_chunks=[]):
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

        # get true error % as of right now.

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
            elif isinstance(chunk, Rectangle):
                rect = patches.Rectangle(chunk.bounding_box['top_left'], 
                                        chunk.bounding_box['bottom_right'][0] - chunk.bounding_box['top_left'][0], 
                                        chunk.bounding_box['bottom_right'][1] - chunk.bounding_box['top_left'][1], 
                                        edgecolor=color, facecolor='none', linewidth=linewidth)
                ax.add_patch(rect)

        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
        ax.set_aspect('equal', 'box')
        plt.axis('off')
        ax.autoscale_view()

        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()

    def add_chunk(self, chunk, motif):
        self.chunks[chunk] = motif

    def circle_rectangle_overlap(self, circle, rectangle):
        # Find the closest point to the circle within the rectangle
        closest_x = max(rectangle.bounding_box['top_left'][0], min(circle.x, rectangle.bounding_box['bottom_right'][0]))
        closest_y = max(rectangle.bounding_box['top_left'][1], min(circle.y, rectangle.bounding_box['bottom_right'][1]))

        # Calculate the distance between the circle's center and this closest point
        distance_x = circle.x - closest_x
        distance_y = circle.y - closest_y

        # If the distance is less than the circle's radius, an overlap occurs
        distance_squared = distance_x**2 + distance_y**2
        return distance_squared < circle.radius**2
    
    def circles_overlap(self, circle1, circle2):
        dx = circle1.x - circle2.x
        dy = circle1.y - circle2.y
        distance_squared = dx**2 + dy**2
        radius_sum = circle1.radius + circle2.radius + self.draw_parms.NODE_R * 2

        return round(distance_squared, 2) < round(radius_sum**2, 2)
    
    def rectangles_overlap(self, rect1, rect2):
        bb1 = rect1.bounding_box
        bb2 = rect2.bounding_box

        return not (bb1['top_left'][0] > bb2['bottom_right'][0] or 
                    bb1['bottom_right'][0] < bb2['top_left'][0] or 
                    bb1['top_left'][1] > bb2['bottom_right'][1] or 
                    bb1['bottom_right'][1] < bb2['top_left'][1])

    def check_any_overlap(self):
        overlaps = 0

        chunk_keys = list(self.chunks.keys())
        overlapping_chunks = []

        for i in range(len(chunk_keys)):
            for j in range(i + 2, len(chunk_keys)):
                chunk1 = chunk_keys[i]
                chunk2 = chunk_keys[j]

                if not (self.chunks[chunk2] in self.chunks[chunk1].children or self.chunks[chunk1] in self.chunks[chunk2].children):
                    if len(self.chunks[chunk2].parent.positions) > 2 and len(self.chunks[chunk1].positions) + len(self.chunks[chunk2].positions) > 2:
                        overlaps_before = overlaps
                        if isinstance(chunk1, Rectangle) and isinstance(chunk2, Rectangle):
                            if self.rectangles_overlap(chunk1, chunk2):
                                overlapping_chunks.append(chunk1)
                                overlapping_chunks.append(chunk2)
                                overlaps += 1
                        elif isinstance(chunk1, Circle) and isinstance(chunk2, Circle):
                            if self.circles_overlap(chunk1, chunk2, ):
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

