# polyopening
A python tool to determine the opening angle of coordination polyhedra of solids and molecular systems.

## Requirements

The following packages are required
* numpy >= 2.0,
* ase >= 3.25.0 ([Atomic Simulation Environment](https://docs.ase-lib.org/))
* qc-grid >= 0.0.9
* scipy <= 1.16.2

## Installation

```bash
pip install polyopening
```

## Usage
### Set-up
```
from polyopening import PolyhedralOpeningAngle

get_angles = PolyhedralOpeningAngle(
    cutoff=cutoff,
    lebedev_degree=131,
    atom_type="C",
    neighbor_type=["C", "H"]
)
```
* `cutoff`: cutoff radius in &#8491; for the first coordination shell of the reference atom. Can be a single float number for a global cutoff or a pair-wise element specific dictionary like ```cutoff = {('H', 'H'): 1.1, ('C', 'H'): 1.3, ('C', 'C'): 1.85}``` See also ASE [`neighbor_list`](https://docs.ase-lib.org/ase/neighborlist.html#ase.neighborlist.neighbor_list).
* `lebedev_degree`: Determines the number of grid points for the [Lebedev Quadrature](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.lebedev_rule.html). Maximum degree is 131 ( = 5810 grid points).
* `atom_type`: The element type, for which the opening angles are to be calculated. Can be a single string (`"C"`) or a list of strings (`["C", "H"]`). If all atoms are to be determined, select `"all"`.
* `neighbor_type`: The element type(s) to be taken into account when determining the coordination polyhedron. Same syntax as `atom_type`.

### Run
```
angles = get_angles.compute(Atoms_object)
```
The calculation of the opening angles is executed by assigning an [ASE Atoms object](https://docs.ase-lib.org/ase/atoms.html#) to the `compute( )` function. Returns a 1D array with all the opening angles in [deg].

## References
When using the code, please cite

L. Hückmann, J. Cottom, J. Meyer *Adv. Phys. Res.* **2024**, *3*, 2300109. [https://doi.org/10.1002/apxr.202300109](https://doi.org/10.1002/apxr.202300109)
