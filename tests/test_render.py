from rna_draw import render_rna
import unittest

class RenderUnittest(unittest.TestCase):
    def test_structure1(self):
        ss1 = "((((....))))"
        pairmap1 = render_rna.get_pairmap_from_secstruct(ss1)
        expected_pairmap1=[11, 10, 9, 8, -1, -1, -1, -1, 3, 2, 1, 0]
        self.assertTrue(pairmap1==expected_pairmap1)
        
    def test_structure2(self):
        ss2 = "(((.(((..((((........))))..((((.......))))......(((((.......)))))))).)))...."
        pairmap2 = render_rna.get_pairmap_from_secstruct(ss2)
        expected_pairmap2=[71, 70, 69, -1, 67, 66, 65, -1, -1, 24, 23, 22, 21, -1, -1, -1, -1,
        -1, -1, -1, -1, 12, 11, 10, 9, -1, -1, 41, 40, 39, 38, -1, -1, -1, -1, -1, -1, -1, 30,
        29, 28, 27, -1, -1, -1, -1, -1, -1, 64, 63, 62, 61, 60, -1, -1, -1, -1, -1, -1, -1, 52,
        51, 50, 49, 48, 6, 5, 4, -1, 2, 1, 0, -1, -1, -1, -1]
        self.assertTrue(pairmap2==expected_pairmap2)
        
    def test_structure3(self):
        ss3 = "(((..(((...)))..)))"
        pairmap3 = render_rna.get_pairmap_from_secstruct(ss3)
        expected_pairmap3 = [18, 17, 16, -1, -1, 13, 12, 11, -1, -1, -1, 7, 6, 5, -1, -1, 2, 1, 0]
        self.assertTrue(pairmap3==expected_pairmap3)
        
        
if __name__ == '__main__':
    unittest.main()
