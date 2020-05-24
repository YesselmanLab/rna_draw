import unittest

from rna_draw import draw, colorer

# only do this is something changed in the algorithm and that visually inspect
# images
def generate_correct_imgs():
    pass


class Unittest(unittest.TestCase):

    def test_basic(self):
        rd = draw.RNADrawer()
        rd.draw("(....)", "CUUCGG", "test_renders/test_1")

        rd.draw("(....)", "CUUCGG", "test_renders/test_2",
                color_str="1-6:r")

        rd.draw("(....)", "CUUCGG", "test_renders/test_3",
                render_type=colorer.RenderType.RES_TYPE)

        rd.draw("(....)", "CUUCGG", "test_renders/test_4",
                default_color=colorer.COLORS["g"])


def main():
    unittest.main()

if __name__ == '__main__':
    main()