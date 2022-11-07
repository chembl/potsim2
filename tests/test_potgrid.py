# from numpy.testing import assert_equal
from potsim2 import PotGrid
import unittest
import os

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestPotGrid(unittest.TestCase):
    def test_skin_similarity(self):
        pdb1 = os.path.join(TESTS_DIR, "data/A__02_01.pdb")
        pdb2 = os.path.join(TESTS_DIR, "data/A__11_01.pdb")
        grid1 = os.path.join(TESTS_DIR, "data/A__02_01.pkl")
        grid2 = os.path.join(TESTS_DIR, "data/A__11_01.pkl")

        grid1 = PotGrid(pdb1, grid1)
        grid2 = PotGrid(pdb2, grid2)

        skin_mask1 = grid1.get_skin_mask()
        grid1.apply_mask(skin_mask1)

        skin_mask2 = grid2.get_skin_mask()
        grid2.apply_mask(skin_mask2)
        hsi, dis = grid1.score(grid2)

        self.assertAlmostEqual(hsi, 0.9081104462828233, 13)
        self.assertAlmostEqual(dis, 0.4286946552435117, 13)

    def test_residue_sphere_sim(self):

        resnum = 149

        pdb1 = os.path.join(TESTS_DIR, "data/A__02_01.pdb")
        pdb2 = os.path.join(TESTS_DIR, "data/A__11_01.pdb")
        grid1 = os.path.join(TESTS_DIR, "data/A__02_01.pkl")
        grid2 = os.path.join(TESTS_DIR, "data/A__11_01.pkl")

        grid1 = PotGrid(pdb1, grid1)
        grid2 = PotGrid(pdb2, grid2)

        res_mask = grid1.get_residue_sphere_mask(resnum)
        grid1.apply_mask(res_mask)

        res_mask2 = grid2.get_residue_sphere_mask(resnum)
        grid2.apply_mask(res_mask2)

        hsi, dis = grid1.score(grid2)

        self.assertAlmostEqual(hsi, 0.9895875459929891, 13)
        self.assertAlmostEqual(dis, 0.14430837818374131, 13)

    def test_resample(self):

        pdb1 = os.path.join(TESTS_DIR, "data/A__02_01.pdb")
        grid1 = os.path.join(TESTS_DIR, "data/A__02_01.pkl")

        grid1 = PotGrid(pdb1, grid1)
        grid1 = grid1.resample_factor(0.5)

        self.assertEqual(grid1.grid.shape, (18, 18, 18))


if __name__ == "__main__":
    unittest.main()
