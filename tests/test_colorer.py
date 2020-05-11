import unittest
import seaborn as sns
import matplotlib.colors


from rna_draw import colorer, parameters

class ColorerUnittest(unittest.TestCase):

    def test_init(self):
        col = colorer.Colorer()


    def test_basic(self):
        col = colorer.Colorer()
        seq, ss = 'AAA', '...'

        with self.subTest("default colors"):
            expected = [colorer.COLORS["e"] for i in range(3)]
            self.assertTrue(col.get_rgb_colors(seq, ss) == expected)

        with self.subTest("new default color"):
            expected = [colorer.COLORS["r"] for i in range(3)]
            c = colorer.COLORS["r"]
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, default_color=c) == expected)


    def test_basic_color_strs(self):
        col = colorer.Colorer()
        seq, ss = 'AAA', '...'

        with self.subTest("standard 1 character color code"):
            color_str = "rbg"
            expected = [colorer.COLORS["r"], colorer.COLORS["b"], colorer.COLORS["g"]]
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, color_str=color_str) == expected)

        with self.subTest("equivlent to previous test, using optional ; to divide color names"):
            color_str = "r;b;g"
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, color_str=color_str) == expected)

        with self.subTest("using xkcd color names"):
            color_str = "windows blue;windows blue;faded green"
            expected = [colorer.xkcd_color_name_to_rgb(x) for x in color_str.split(";")]
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, color_str=color_str) == expected)

        with self.subTest("using range arguments, specify color for group of residues"):
            color_str = "1-3:r"
            expected = [colorer.COLORS["r"], colorer.COLORS["r"], colorer.COLORS["r"]]
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, color_str=color_str) == expected)

        with self.subTest("using both range and xkcd colors"):
            color_str = "1-3:windows blue;"
            expected = [colorer.xkcd_color_name_to_rgb("windows blue") for x in range(3)]
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, color_str=color_str) == expected)

        with self.subTest("using multiple non overlapping range arguments"):
            color_str = "1:r;2-3:g"
            expected = [colorer.COLORS["r"], colorer.COLORS["g"], colorer.COLORS["g"]]
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, color_str=color_str) == expected)

        with self.subTest("using range check if default is still applied"):
            color_str = "1-2:r"
            expected = [colorer.COLORS["r"], colorer.COLORS["r"], colorer.COLORS["e"]]
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, color_str=color_str) == expected)


    def test_basic_color_str_errors(self):
        col = colorer.Colorer()
        seq, ss = 'AAA', '...'

        with self.subTest("too many colors specifed longer than seq"):
            with self.assertRaises(ValueError):
                color_str = "rgbb"
                col.get_rgb_colors(seq, ss, color_str=color_str)

        with self.subTest("color string too short"):
            with self.assertRaises(ValueError):
                color_str = "r;g"
                col.get_rgb_colors(seq, ss, color_str=color_str)

        with self.subTest("overlapping colors"):
            with self.assertRaises(ValueError):
                color_str = "1-3:r;2-3:g"
                col.get_rgb_colors(seq, ss, color_str=color_str)


    def test_render_type(self):
        col = colorer.Colorer()
        seq, ss = 'ACGU', '(..)'

        with self.subTest("test coloring by restype"):
            expected = colorer.parse_color_single_letter_codes("ygrb")
            self.assertTrue(
                    col.get_rgb_colors(seq, ss,
                                       render_type=colorer.RenderType.RES_TYPE) == expected)

        with self.subTest("test coloring by restype with override"):
            expected = colorer.parse_color_single_letter_codes("egrb")
            self.assertTrue(
                    col.get_rgb_colors(seq, ss,
                                       render_type=colorer.RenderType.RES_TYPE,
                                       color_str="1:e") == expected)

        with self.subTest("test coloring by pairing"):
            expected = colorer.parse_color_single_letter_codes("byyb")
            self.assertTrue(
                    col.get_rgb_colors(seq, ss,
                                       render_type=colorer.RenderType.PAIRED) == expected)


    def test_data_color(self):
        col = colorer.Colorer()
        seq, ss = 'ACGU', '(..)'
        data = [0, 1, 2, 3]
        colors = colorer.color_by_data(data)

        with self.subTest("test coloring by data"):
            self.assertTrue(
                    col.get_rgb_colors(seq, ss, data=data) == colors)

        with self.subTest("cannot set both render_type and data for the same res"):
            with self.assertRaises(ValueError):
                col.get_rgb_colors(seq, ss, render_type=colorer.RenderType.RES_TYPE,
                                   data=data)

def main():
    unittest.main()

if __name__ == '__main__':
    main()


























