import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from ase.neighborlist import neighbor_list
from grid.angular import AngularGrid
import sys

class PolyhedralOpeningAngle:
    """
    Analyzer for polyhedral opening angles.

    Notes
    -----
    The analyzer is stateless. Each call to ``compute()`` analyzes
    one ASE ``Atoms`` object independently.
    
    Parameters
    ----------
    cutoff : float
        Cutoff radius (Å) for identifying neighboring atoms.
    lebedev_degree : int
        Degree of the Lebedev angular grid.
    atom_type : str, optional
        Chemical symbol of the atoms for which the opening angle
        should be evaluated. Default is "all".
    """

    def __init__(
        self,
        cutoff: 2.0,
        lebedev_degree: int = 110,
        atom_type: str | list[str] = "all",
        neighbor_type: str | list[str] = "all"
    ) -> None:

        self.cutoff = cutoff
        self.atom_type = atom_type
        self.neighbor_type = neighbor_type
        self.lebedev_degree = lebedev_degree

        self.grid = AngularGrid(degree=lebedev_degree).points
    
    @staticmethod
    def angle(rij, rik):
        rad = np.arccos( ( np.dot( rij, rik ) ) / (np.sqrt(np.dot(rij, rij)) * np.sqrt(np.dot(rik, rik)) ) )
        return np.rad2deg(rad)
        
    def _analyze_cp(
            self,
            atoms: Atoms,
            i: np.ndarray,
            j: np.ndarray,
            refidx: int,
            ) -> float:
        neighbors = j[i == refidx]
        refatom=atoms.positions[refidx]
        if len(neighbors)==0:
            print("\n*** WARNING ***\nempty neighborlist at atom "+str(refidx))
            return np.nan
        else:
            neigh_coords=atoms.positions[neighbors]
            local_grid=self.grid+refatom
            dr = local_grid[:, None, :] - neigh_coords[None, :, :]
            dr = dr.reshape(-1, 3)

            dr_mic, lengths = find_mic(dr, atoms.cell, atoms.pbc)

            lengths = lengths.reshape(len(self.grid), len(neigh_coords))
        
            # Distance to closest neighbor
            nearest = np.min(lengths, axis=1)

            # Best grid point
            best_idx = np.argmax(nearest)
        
            closest_atom = np.argmin(lengths[best_idx])
            xyz1,_=find_mic(local_grid[best_idx]-refatom, atoms.cell, atoms.pbc)
            xyz2,_=find_mic(neigh_coords[closest_atom]-refatom, atoms.cell, atoms.pbc)
            return self.angle(xyz1,xyz2)
    
    def compute(
            self,
            atoms: Atoms,
            ) -> list[float]:
        elem=np.array(atoms.get_chemical_symbols())
        coords=atoms.get_positions()
        i,j=neighbor_list("ij",atoms,self.cutoff,self_interaction=False)
        
        if "all" in self.atom_type:
            i_idx=np.arange(0,len(i),1)
            idx_ref=np.arange(0,len(coords),1)
        else:
            idx_ref=np.where(np.isin(elem, self.atom_type))[0]
            i_idx=np.where(np.isin(i,idx_ref))[0]
            #print(i_idx) 
        if "all" in self.neighbor_type:
            j_idx=np.arange(0,len(j),1)
        else:
            idx_neigh=np.where(np.isin(elem, self.neighbor_type))[0]
            j_idx=np.where(np.isin(j,idx_neigh))[0]
        idx=np.intersect1d(i_idx, j_idx, return_indices=True)[0]

        PHI=[]
        for ii in idx_ref:
            phi=self._analyze_cp(atoms,i[idx],j[idx],ii)
            PHI.append(phi)
        return np.array(PHI)
    

            
            
            
