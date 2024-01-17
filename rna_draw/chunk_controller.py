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
    def __init__(self, struct, xarray, yarray):
        self.chunks = {}
        self.struct = struct

        for m in struct:
            if m in struct.get_junctions() and hasattr(m, 'center_x') and hasattr(m, 'center_y') and hasattr(m, 'radius'):
                circle = Circle(m.center_x, m.center_y, m.radius)
                self.add_chunk(circle, m)
            else:
                nodes = []
                for node in m.positions:
                    nodes.append((xarray[node], yarray[node]))
                rect = Rectangle(nodes)
                self.add_chunk(rect, m)

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
        radius_sum = circle1.radius + circle2.radius

        return distance_squared < radius_sum**2
    
    def rectangles_overlap(self, rect1, rect2):
        bb1 = rect1.bounding_box
        bb2 = rect2.bounding_box

        return not (bb1['top_left'][0] > bb2['bottom_right'][0] or 
                    bb1['bottom_right'][0] < bb2['top_left'][0] or 
                    bb1['top_left'][1] > bb2['bottom_right'][1] or 
                    bb1['bottom_right'][1] < bb2['top_left'][1])

    def check_any_overlap(self):
        overlaps = 0
        overlap_table = []

        chunk_keys = list(self.chunks.keys())

        for i in range(len(chunk_keys)):
            for j in range(i + 2, len(chunk_keys)):
                chunk1 = chunk_keys[i]
                chunk2 = chunk_keys[j]

                if not (self.chunks[chunk2] in self.chunks[chunk1].children or self.chunks[chunk1] in self.chunks[chunk2].children):
                    if len(self.chunks[chunk2].parent.positions) > 2 and len(self.chunks[chunk1].positions) + len(self.chunks[chunk2].positions) > 2:
                        if isinstance(chunk1, Rectangle) and isinstance(chunk2, Rectangle):
                            if self.rectangles_overlap(chunk1, chunk2):
                                overlaps += 1
                        elif isinstance(chunk1, Circle) and isinstance(chunk2, Circle):
                            if self.circles_overlap(chunk1, chunk2):
                                overlaps += 1
                        elif isinstance(chunk1, Rectangle) and isinstance(chunk2, Circle):
                            if self.circle_rectangle_overlap(chunk2, chunk1):
                                overlaps += 1
                        elif isinstance(chunk1, Circle) and isinstance(chunk2, Rectangle):
                            if self.circle_rectangle_overlap(chunk1, chunk2):
                                overlaps += 1

        return overlaps

