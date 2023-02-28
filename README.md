# PotSim2: Simple package to segment and compare protein potential grids

Using [gridDataFormats](https://github.com/MDAnalysis/GridDataFormats) to read the grids and [MDAnalysis](https://www.mdanalysis.org/) to make pdb atom selections

## Installation

```bash
pip install potsim2
```

## Usage: Open two grids and calculate [PIPSA](https://pipsa.h-its.org/pipsa/) like scores

The protein electrostatic potential grids can generated with [APBS](https://github.com/Electrostatics/apbs)

```python
from potsim2 import PotGrid

# read grids
grid1 = PotGrid('A__02_01.pdb', 'A__02_01.pkl')
grid2 = PotGrid('A__11_01.pdb', 'A__11_01.pkl')

# calculate skin for grid1 
skin_mask1 = grid1.get_skin_mask()
grid1.apply_mask(skin_mask1)

# calculate skin for grid2 
skin_mask2 = grid2.get_skin_mask()
grid2.apply_mask(skin_mask2)

# calc the Hodgkin similarity index and PIPSA like distance 
hsi, dis = grid1.score(grid2)

# export the skins in uncompressed OpenDX format to be visualized in PyMol/ChimeraX
grid1.export('A__02_01.dx')
grid2.export('A__11_01.dx')
```
