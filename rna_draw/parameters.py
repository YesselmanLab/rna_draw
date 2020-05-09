

class Parameters:
    def __init__(self, seq, ss, data=None, color_str=None, render_type=None):
        self.sequence = seq
        self.structure = ss
        self.data = data
        self.color_str = color_str
        self.render_type = render_type
        self.NODE_R = 10
        self.PRIMARY_SPACE = 20
        self.PAIR_SPACE = 23
        self.CELL_PADDING = 40
        self.TEXT_SIZE = 50
        self.RENDER_IN_LETTERS = False
