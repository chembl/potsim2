import numpy as np
from gridData import Grid
from .skinlib import fill_exclusion_mask
import prody as pr
from .vdw import get_vdw_radius

pr.confProDy(verbosity="none")


class PotGrid(Grid):
    def __init__(
        self,
        pdb_filename=None,
        grid=None,
        edges=None,
        origin=None,
        delta=None,
        metadata=None,
        interpolation_spline_order=3,
        file_format=None,
    ):

        super().__init__(
            grid,
            edges,
            origin,
            delta,
            metadata,
            interpolation_spline_order,
            file_format,
        )
        self._load_pdb(pdb_filename)

    def _load_pdb(self, pdb_filename):
        # read pdb file and calc protein center
        self.pdb = pr.parsePDB(pdb_filename)
        self.protein = self.pdb.select("protein")
        self.protein_center = pr.calcCenter(self.protein)
        self.oe = self.origin - self.protein_center

    def get_skin_mask(self, probe=3, skin_thickness=4):

        vdw_radii = np.asarray(
            [get_vdw_radius(element) for element in self.protein.getElements()],
            dtype="float64",
        )
        atom_coords = pr.getCoords(self.protein)

        mask1 = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask1,
            self.delta,
            self.oe,
            self.protein_center,
            atom_coords,
            vdw_radii,
            vdw_radii.max(),
            probe,
        )
        mask1 = np.invert(mask1)

        mask2 = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask2,
            self.delta,
            self.oe,
            self.protein_center,
            atom_coords,
            vdw_radii,
            vdw_radii.max(),
            probe + skin_thickness,
        )
        return mask1 * mask2

    def get_sphere_mask(self, sphere_center, radius):
        mask = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask,
            self.delta,
            self.oe,
            self.protein_center,
            np.asarray(sphere_center).reshape(1, -1),
            np.asarray([0]),
            0,
            radius,
        )
        return mask

    def get_ligand_mask(self, ligand, radius=4.5):
        atom_coords = pr.getCoords(ligand)
        mask = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask,
            self.delta,
            self.oe,
            self.protein_center,
            atom_coords,
            np.zeros(atom_coords.shape[0]),
            0,
            radius,
        )
        return mask

    def get_res_sphere_mask(self, resnum, radius=4.5):
        res_ca = self.pdb.select(f"ca resnum {resnum}")
        coords = pr.getCoords(res_ca)
        mask = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask,
            self.delta,
            self.oe,
            self.protein_center,
            coords,
            np.zeros(coords.shape[0]),
            0,
            radius,
        )
        return mask

    def apply_mask(self, mask):
        self.grid *= mask

    def score(self, other, intersection=True):
        """
        Calcs hodgkin SI on grid intersection.
        """
        g1arr = self.grid.ravel()
        g2arr = other.grid.ravel()

        if intersection:
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
