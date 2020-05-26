import unittest

import rna_draw as rd
from rna_draw import colorer


class Unittest(unittest.TestCase):


    def test_basic(self):
        rd.rna_draw(ss="(....)", seq="CUUCGG", out="test_renders/test_1")

        rd.rna_draw(ss="(....)", seq="CUUCGG", out="test_renders/test_2",
                    color_str="1-6:r")

        rd.rna_draw(ss="(....)", seq="CUUCGG", out="test_renders/test_3",
                    render_type="res_type")

        rd.rna_draw(ss="(....)", seq="CUUCGG", out="test_renders/test_4",
                    default_color="g")

    def _test_w_data(self):
        ss  = ".(((.....)))."
        seq = "AUGAAAAAAUCAA"

        with self.subTest("simplest data use"):
            rd.rna_draw(ss, seq, "test_renders/test_5",
                        data="0;1;2;3;4;5;6;7;8;9;10;11;12")

        with self.subTest("ignore restypes"):
            rd.rna_draw(ss, seq, "test_renders/test_6",
                        data="0;1;2;3;4;5;6;7;8;9;10;11;12",
                        data_ignore_restype="A")

        with self.subTest("ignore restypes multiple"):
            rd.rna_draw(ss, seq, "test_renders/test_7",
                        data="0;1;2;3;4;5;6;7;8;9;10;11;12",
                        data_ignore_restype="GCU")


def main():
    unittest.main()

if __name__ == '__main__':
    main()