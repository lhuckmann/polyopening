# polyopening
A python tool to determine the opening angle of coordination polyhedra of solids and molecular systems.

## Requirements

The following packages are required
* numpy >= 2.0,
* ase >= 3.25.0 ([Atomic Simulation Environment](https://docs.ase-lib.org/))
* qc-grid >= 0.0.9 

## Installation

## Usage
### Set-up
```
from polyopening import PolyhedralOpeningAngle

get_angles = PolyhedralOpeningAngle(
    cutoff=cutoff,
    lebedev_degree=131,
    atom_type=element_1,
    neighbor_type=element_2
)
```
* `cutoff`: cutoff radius in &#8491 for the first coordination shell of the reference atom. Can be a singe float number or a pair-wise element specific dictionary like ```cutoff = {('H', 'H'): 1.1, ('C', 'H'): 1.3, ('C', 'C'): 1.85}``` See also ASE [`neighbor_list`](https://docs.ase-lib.org/ase/neighborlist.html#ase.neighborlist.neighbor_list).
* `lebedev_degree`: Determines the number of grid points for the Lebedev Quadrature. Maximum degree is 131 ( = 5810 grid points).
* `atom_type`: The element type, for which the opening angles are to be calculated. Can be a single string (`"C"`) or a list of strings (`["C", "H"]`). If all atoms are to be determined, select `"all"`.
