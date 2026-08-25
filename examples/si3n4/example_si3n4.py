from ase.io import read
from polyopening import PolyhedralOpeningAngle
import matplotlib.pyplot as plt
import numpy as np

# set-up and load structure
a=14.6488

si3n4 = read("strc.xyz") 
si3n4.set_cell([a,a,a])
si3n4.set_pbc([True,True,True]) #works also for clusters
cutoff=2.25


# initiate
get_angles = PolyhedralOpeningAngle(
    cutoff=cutoff,
    lebedev_degree=131, #max 131
    atom_type="Si", # or ["Si"] / ["Si", "N"] / "all"
    neighbor_type="N" #see atom_type
)

# calculate phi for ASE Atoms object
angles = get_angles.compute(si3n4) #takes ASE Atoms object

# plot the result
plt.figure()
plt.hist(angles, bins=np.arange(60.5,180.5,1), label="data")
plt.axvline(70.53,0,1,ls=":",color='r',label="tetrahedron")
plt.axvline(90,0,1,ls="--",color='r',label="planar")
plt.axvline(180,0,1,ls="--",color='r',label="linear")
plt.xlabel(r"$\phi [\,$deg$\,]$")
plt.ylabel(r"$P(\phi)$")
plt.title(r"Opening angle $\phi$ of Si$_3$N$_4$")
plt.legend()
plt.show()
