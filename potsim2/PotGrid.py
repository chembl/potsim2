import numpy as np
from gridData import Grid
from .skinlib import fill_exclusion_mask
from prody import parsePDB, confProDy
from .vdw import get_vdw_radius

confProDy(verbosity='none')

class PotGrid(Grid):

    def __init__(self, grid=None, edges=None, origin=None, delta=None,
                 metadata=None, interpolation_spline_order=3, file_format=None):
        super().__init__(grid, edges, origin, delta,
                         metadata, interpolation_spline_order, file_format)

    def _calc_pdb_center(self, atom_coords):
        center = atom_coords.mean(axis=0)
        return center

    def get_skin_mask(self, pdb_fname, probe=3, skin=4):
        atoms = parsePDB(pdb_fname)

        vdw_radii = np.asarray([get_vdw_radius(a.getElement())
                                for a in atoms], dtype='float64')
        atom_coords = np.asarray([a.getCoords()
                                  for a in atoms], dtype='float64')
        pdb_center = self._calc_pdb_center(atom_coords)
        oe = self.origin - pdb_center

        mask1 = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(mask1, self.delta, oe, pdb_center,
                            atom_coords, vdw_radii, vdw_radii.max(), probe)
        mask1 = np.invert(mask1)

        mask2 = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(mask2, self.delta, oe, pdb_center,
                            atom_coords, vdw_radii, vdw_radii.max(), probe + skin)

        return mask1 * mask2

    def apply_mask(self, mask):
        self.grid *= mask

    def score(self, other, skin=True):
        """
        Calcs hodgkin SI on grid intersection.
        """
        g1arr = self.grid.ravel()
        g2arr = other.grid.ravel()

        if skin:
            mask = (g1arr == 0) | (g2arr == 0)
            g1arr[mask] = 0
            g2arr[mask] = 0

        p1p1 = g1arr.dot(g1arr)
        p2p2 = g2arr.dot(g2arr)
        p1p2 = g1arr.dot(g2arr)

        # hodgkin SI
        hodgkin_si = 2 * p1p2 / (p1p1 + p2p2)

        # PIPSA returns a distance based score when using mkdismx
        pipsa_diff_index = np.sqrt(2 - 2 * hodgkin_si)
        return hodgkin_si, pipsa_diff_index
