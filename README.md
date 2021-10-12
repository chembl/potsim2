# PotSim2: Simple Python/C++ NumPy centric package to segment protein potential grids and to compare electrostatic potentials of proteins.

Using [gridDataFormats](https://github.com/MDAnalysis/GridDataFormats) to read the grids.

## Installation

```bash
git clone https://github.com/chembl/potsim2.git
pip install potsim2/
```

## Usage: Extract a volumetric ligand mask from a [ProdDy](https://github.com/prody/ProDy) set of atoms

Volumetric grids can be used to train DL models.

```python
from potsim2 import PotGrid
import prody as pr

# read grid
grid = PotGrid('B__40_01.pdb', 'B__40_01.dx.gz')

# read a pdb file (ideally one containing only the ligand atoms)
pdb = pr.parsePDB('B__40_01.pdb')

# chain C is an small peptide in this pdb file, select it
atoms = pdb.select('chain C')

# ligand mask is a boolean NumPy array, can be converted to int: ligand_mask.astype(int)
ligand_mask = grid.get_ligand_mask(atoms)

# export the (uncompressed) mask to visually verify it with PyMol/ChimeraX (hacky)
grid.grid = ligand_mask.astype(int)
grid.export('grid_c_chain.dx')
```

## Usage: Open two grids and calculate [PIPSA](https://pipsa.h-its.org/pipsa/) like scores

The protein electrostatic potential grids can generated with [APBS](https://github.com/Electrostatics/apbs).

```python
from potsim2 import PotGrid

# read grids
grid1 = PotGrid('B__40_01.pdb', 'B__40_01.dx.gz')
grid2 = PotGrid('B__40_02.pdb', 'B__40_02.dx.gz')

# calculate skin for grid1 
skin_mask1 = grid1.get_skin_mask()
grid1.apply_mask(skin_mask1)

# calculate skin for grid2 
skin_mask2 = grid2.get_skin_mask()
grid2.apply_mask(skin_mask2)

# calc the Hodgkin similarity index and PIPSA like distance 
hsi, dis = grid1.score(grid2)

# export (uncompressed) skins to be visualized in PyMol/ChimeraX
grid1.export('B__40_01_skin.dx')
grid2.export('B__40_02_skin.dx')
```
