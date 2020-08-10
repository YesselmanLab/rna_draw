import unittest

from rna_draw import data


class DataUnittest(unittest.TestCase):
    def test_load_values(self):
        vals = [float(x) for x in range(1, 6)]
        with self.subTest("load data from str"):
            d = data.Data("1;2;3;4;5")
            self.assertTrue(d.vals == vals)

        with self.subTest("load data from file"):
            d = data.Data(data_file="resources/test_data_file.dat")
            self.assertTrue(d.vals == vals)

    def test_min_max(self):
        d = data.Data("1;2;3;4;5", vmin=0, vmax=10)


def main():
    unittest.main()


if __name__ == "__main__":
    main()
