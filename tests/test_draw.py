import unittest

import rna_draw as rd
from rna_draw import settings

path = settings.Paths.UNITTEST_PATH


class Unittest(unittest.TestCase):
    def test_basic(self):
        rd.rna_draw(ss="(....)", seq="CUUCGG", out=path / "test_renders/test_1")

        rd.rna_draw(
            ss="(....)",
            seq="CUUCGG",
            out=path / "test_renders/test_2",
            color_str="1-6:r",
        )

        rd.rna_draw(
            ss="(....)",
            seq="CUUCGG",
            out=path / "test_renders/test_3",
            render_type="res_type",
        )

        rd.rna_draw(
            ss="(....)",
            seq="CUUCGG",
            out=path / "test_renders/test_4",
            default_color="g",
        )

    def test_w_data(self):
        ss = ".(((.....)))."
        seq = "AUGAAAAAAUCAA"
        data_str = "0;1;2;3;4;5;6;7;8;9;10;11;12"

        with self.subTest("simplest data use"):
            rd.rna_draw(
                ss=ss, seq=seq, out=path / "test_renders/test_5", data_str=data_str
            )

        with self.subTest("ignore restypes"):
            rd.rna_draw(
                ss=ss,
                seq=seq,
                out=path / "test_renders/test_6",
                data_str=data_str,
                data_ignore_restype="A",
            )

        with self.subTest("ignore restypes multiple"):
            rd.rna_draw(
                ss=ss,
                seq=seq,
                out=path / "test_renders/test_7",
                data_str=data_str,
                data_ignore_restype="GCU",
            )


def main():
    unittest.main()


if __name__ == "__main__":
    main()
