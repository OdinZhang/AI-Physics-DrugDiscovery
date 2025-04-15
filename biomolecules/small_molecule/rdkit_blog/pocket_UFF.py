import os
import os.path as osp
from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit import Chem
import scipy.spatial as spatial
import numpy as np

def get_rd_atom_res_id(rd_atom):
    '''
    Return an object that uniquely
    identifies the residue that the
    atom belongs to in a given PDB.
    '''
    res_info = rd_atom.GetPDBResidueInfo()
    return (
        res_info.GetChainId(),
        res_info.GetResidueNumber()
    )
def get_pocket(lig_mol, rec_mol, max_dist=8):
    lig_coords = lig_mol.GetConformer().GetPositions()
    rec_coords = rec_mol.GetConformer().GetPositions()
    dist = spatial.distance.cdist(lig_coords, rec_coords)

    # indexes of atoms in rec_mol that are
    #   within max_dist of an atom in lig_mol
    pocket_atom_idxs = set(np.nonzero((dist < max_dist))[1])

    # determine pocket residues
    pocket_res_ids = set()
    for i in pocket_atom_idxs:
        atom = rec_mol.GetAtomWithIdx(int(i))
        res_id = get_rd_atom_res_id(atom)
        pocket_res_ids.add(res_id)

    # copy mol and delete atoms
    pkt_mol = rec_mol
    pkt_mol = Chem.RWMol(pkt_mol)
    for atom in list(pkt_mol.GetAtoms()):
        res_id = get_rd_atom_res_id(atom)
        if res_id not in pocket_res_ids:
            pkt_mol.RemoveAtom(atom.GetIdx())

    Chem.SanitizeMol(pkt_mol)
    return pkt_mol

def read_sdf(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file)
    mols_list = [i for i in supp]
    return mols_list

def write_sdf(mol_list,file):
    writer = Chem.SDWriter(file)
    for i in mol_list:
        writer.write(i)
    writer.close()

def uff_geomopt(rd_mol,pkt_mol,lig_constraint=None,n_iters=200,n_tries=2, lig_h=True, pkt_h=False):
    if lig_h:
        rd_mol = Chem.AddHs(rd_mol,addCoords=True)
    if pkt_h:
        pkt_mol = Chem.AddHs(pkt_mol,addCoords=True)

    rd_mol = Chem.RWMol(rd_mol)
    uff_mol = Chem.CombineMols(pkt_mol, rd_mol)

    try:
        Chem.SanitizeMol(uff_mol)
    except Chem.AtomValenceException:
        print('Invalid valence')
    except (Chem.AtomKekulizeException, Chem.KekulizeException):
        print('Failed to kekulize')
    try:
        uff = AllChem.UFFGetMoleculeForceField(
                        uff_mol, confId=0, ignoreInterfragInteractions=False
                    )
        uff.Initialize()
        for i in range(pkt_mol.GetNumAtoms()): # Fix the rec atoms
            uff.AddFixedPoint(i)
        if lig_constraint is not None:
            for i in lig_constraint:
                uff.AddFixedPoint(pkt_mol.GetNumAtoms()+i) # Fix the specified lig atoms

        converged = False
        n_iters=n_iters
        n_tries=n_tries
        while n_tries > 0 and not converged:
            print('.', end='', flush=True)
            converged = not uff.Minimize(maxIts=n_iters)
            n_tries -= 1
        print(flush=True)
        print("Performed UFF with binding site...")
    except:
        print('Skip UFF...')

    return rd_mol
    
from copy import deepcopy
def set_rdmol_positions_(mol, pos):
    """
    Args:
        rdkit_mol:  An `rdkit.Chem.rdchem.Mol` object.
        pos: (N_atoms, 3)
    """
    for i in range(pos.shape[0]):
        mol.GetConformer(0).SetAtomPosition(i, pos[i].tolist())
    return mol
    
def set_rdmol_positions(rdkit_mol, pos):
    """
    Args:
        rdkit_mol:  An `rdkit.Chem.rdchem.Mol` object.
        pos: (N_atoms, 3)
    """
    mol = deepcopy(rdkit_mol)
    set_rdmol_positions_(mol, pos)
    return mol

if __name__ == '__main__':
    pkt_mol = Chem.MolFromPDBFile('./ligan_sampled/0/2z3h_A_rec.pdb')
    rd_mol = read_sdf('./ligan_sampled/0/ligan.sdf')[0]
    opt_mol = uff_geomopt(rd_mol,pkt_mol)