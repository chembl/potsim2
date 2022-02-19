from .skinlib import fill_exclusion_mask
from .vdw import get_vdw_radius
from gridData import Grid
import numpy as np
import prody as pr

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
        self.pdb_filename = pdb_filename
        self._load_pdb(pdb_filename)

    def _load_pdb(self, pdb_filename):
        # read pdb file and calc protein center
        self.pdb = pr.parsePDB(pdb_filename)
        self.protein = self.pdb.select("protein")
        self.protein_center = pr.calcCenter(self.protein)
        self.oe = self.origin - self.protein_center

    def __copy_with_new_grid(self, new_grid):
        """
        Creates a copy of the PotGrid object but with a new grid
        """
        new_grid = PotGrid(
            pdb_filename=self.pdb_filename,
            grid=new_grid,
            edges=self.edges,
            origin=self.origin,
            delta=self.delta,
            metadata=None,
            interpolation_spline_order=self.interpolation_spline_order,
            file_format=None,
        )
        return new_grid

    def __add__(self, other):
        return self.__copy_with_new_grid(self.grid + other.grid)

    def __sub__(self, other):
        return self.__copy_with_new_grid(self.grid - other.grid)

    def __mul__(self, other):
        return self.__copy_with_new_grid(self.grid * other.grid)

    def export(self, filename, file_format=None, grid_dtype=None, typequote='"'):
        if self.grid.dtype == np.dtype("bool"):
            self.grid = self.grid.astype(int)
        super().export(filename, file_format, grid_dtype, typequote)

    def apply_mask(self, mask):
        self.grid *= mask.grid

    def get_skin_mask(self, probe=3, skin_thickness=4):
        """
        Get a PIPSA like skin mask
        """
        vdw_radii = np.asarray(
            [get_vdw_radius(element) for element in self.protein.getElements()],
            dtype="float64",
        )
        atom_coords = pr.getCoords(self.protein)
        # center pdb structure in euclidean (0, 0, 0)
        atom_coords -= self.protein_center

        mask1 = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask1, self.delta, self.oe, atom_coords, vdw_radii, vdw_radii.max(), probe,
        )
        mask1 = np.invert(mask1)

        mask2 = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask2,
            self.delta,
            self.oe,
            atom_coords,
            vdw_radii,
            vdw_radii.max(),
            probe + skin_thickness,
        )
        return self.__copy_with_new_grid(mask1 * mask2)

    def get_sphere_mask(self, sphere_center, radius):
        """
        Get a new grid containing an sphere on arbitrary coords
        """
        mask = np.zeros(self.grid.shape, dtype=bool)
        sphere_center -= self.protein_center
        fill_exclusion_mask(
            mask,
            self.delta,
            self.oe,
            np.asarray(sphere_center).reshape(1, -1),
            np.asarray([0]),
            0,
            radius,
        )
        return self.__copy_with_new_grid(mask)

    def get_atom_list_mask(self, atom_list, radius=4.5):
        """
        Get a new grid containing an sphere on each atom in a ProDy atom selection
        """
        atom_coords = pr.getCoords(atom_list)
        atom_coords -= self.protein_center

        mask = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask,
            self.delta,
            self.oe,
            atom_coords,
            np.zeros(atom_coords.shape[0]),
            0,
            radius,
        )
        return self.__copy_with_new_grid(mask)

    def get_residue_sphere_mask(self, resnum, radius=10):
        """
        Get a new grid containing an sphere on the selected residue alpha carbon
        """
        res_ca = self.pdb.select(f"ca resnum {resnum}")
        res_coords = pr.getCoords(res_ca)
        res_coords -= self.protein_center

        mask = np.zeros(self.grid.shape, dtype=bool)
        fill_exclusion_mask(
            mask,
            self.delta,
            self.oe,
            res_coords,
            np.zeros(res_coords.shape[0]),
            0,
            radius,
        )
        return self.__copy_with_new_grid(mask)

    def score(self, other, intersection=True):
        """
        Calcs hodgkin similarity index on grids intersection.
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

        # hodgkin similarity index
        hodgkin_si = 2 * p1p2 / (p1p1 + p2p2)

        # mkdismx (PIPSA) like distance score
        pipsa_diff_index = np.sqrt(2 - 2 * hodgkin_si)
        return hodgkin_si, pipsa_diff_index
