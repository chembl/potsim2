# PotSim2: Simple Python/C++ package to compare electrostatic potentials of proteins.

Using [gridDataFormats](https://github.com/MDAnalysis/GridDataFormats) to read a variety of grid formats.

## Installation (it requires current gridDataFormats github master branch to read/write .dx.gz files)

```bash
git clone https://github.com/chembl/potsim2.git
pip install potsim2/
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