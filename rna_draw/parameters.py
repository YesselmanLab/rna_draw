class DrawParameters:
    def __init__(self):
        self.NODE_R = 10
        self.PRIMARY_SPACE = 20
        self.PAIR_SPACE = 23
        self.CELL_PADDING = 40
        self.TEXT_SIZE = 50
        self.RENDER_IN_LETTERS = False
        self.output_width = 1280
        self.output_height = 960


class RenderType:
    RES_TYPE = 0
    PAIRED = 1
    MOTIF = 2
    STRAND = 3
