# PotSim2: Simple Python/C++ package to compare electrostatic potentials of proteins.

Using [gridDataFormats](https://github.com/MDAnalysis/GridDataFormats) to read a variety of grid formats.

## Installation

```bash
git clone https://github.com/chembl/potsim2.git
pip install potsim2/
# GridDataFormats master branch needed since we submitted patches to make it work with gzipped dx files and fixed a bug
pip install https://github.com/MDAnalysis/GridDataFormats/archive/master.zip
```

## Usage: Open two grids and calculate Hodgkin similarity index and PIPSA like distance.

```python
from potsim2 import PotGrid

# read grids
grid1 = PotGrid('B__40_01.pdb', 'B__40_01.pkl')
grid2 = PotGrid('B__40_02.pdb', 'B__40_02.pkl')

# calculate skin for grid1 
skin_mask1 = grid1.get_skin_mask()
grid1.apply_mask(skin_mask1)

# calculate skin for grid2 
skin_mask2 = grid2.get_skin_mask()
grid2.apply_mask(skin_mask2)

# calc Hodgkin similarity index and PIPSA like distance 
hsi, dis = grid1.score(grid2)

# export skins so can be visualized in PyMol
grid1.export('B__40_01_skin.dx')
grid2.export('B__40_02_skin.dx')
```

## Usage: Extract ligand mask from a ProdDy set of atoms.


```python
from potsim2 import PotGrid
import prody as pr

# read grid
grid = PotGrid('B__40_01.pdb', 'B__40_01.pkl')

# read a pdb file (ideally one containing only the ligand atoms)
pdb = pr.parsePDB('B__40_01.pdb')

# chain C is a peptide in this pdb file, select it
peptide = pdb.select('chain C')

# ligand mask is a boolean NumPy array, can be converted to int: ligand_mask.astype(int)
ligand_mask = grid.get_ligand_mask(peptide)

# save it to dx format to open it with PyMol to verify we did what we wanted (a bit hacky way)
grid.grid = ligand_mask.astype(int)
grid.export('grid_c_chain.dx')
```