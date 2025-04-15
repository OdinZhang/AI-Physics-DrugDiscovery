from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import os.path as osp
import mdtraj as md
import numpy as np
import os
from collections import defaultdict


def hbd(mol):
    return rdMolDescriptors.CalcNumHBA(mol)
def hbd(mol):
    return rdMolDescriptors.CalcNumHBD(mol)
def tpsa(mol):
    return rdMolDescriptors.CalcTPSA(mol)

import mdtraj as md
import numpy as np

def sasa(mol):
    pdb_file = './tmp.pdb'
    Chem.MolToPDBFile(mol, pdb_file)
    sasa = compute_sasa(pdb_file)
    os.remove(pdb_file)
    return sasa

def compute_sasa(pdb_file):
    # Load the PDB file
    traj = md.load(pdb_file)

    # Compute the SASA
    sasa = md.shrake_rupley(traj, mode='residue')

    # Get the total SASA
    total_sasa = np.sum(sasa)

    return total_sasa

class FusedRingAnalyzer:
    def __init__(self, mol):
        if type(mol) == str:
            self.mol = Chem.MolFromSmiles(mol)
        else:
            self.mol = mol 

        self.ring_info = self.mol.GetRingInfo()
        self.bonds_in_rings = self.ring_info.BondRings()

    def _build_graph(self):
        graph = defaultdict(set)
        bond_to_ring = {}
        for idx, ring in enumerate(self.bonds_in_rings):
            for bond in ring:
                if bond in bond_to_ring:
                    graph[bond_to_ring[bond]].add(idx)
                    graph[idx].add(bond_to_ring[bond])
                bond_to_ring[bond] = idx
        return graph

    def find_max_fused_rings(self):
        graph = self._build_graph()
        visited = set()
        max_fused = 0

        def dfs(node):
            stack = [node]
            count = 0
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    count += 1
                    for neighbor in graph[current]:
                        if neighbor not in visited:
                            stack.append(neighbor)
            return count

        for node in graph:
            if node not in visited:
                fused_count = dfs(node)
                if fused_count > max_fused:
                    max_fused = fused_count

        return max_fused
    
def is_molecule_fragmented(mol):
    """
    Determine if a molecule is fragmented (i.e., consists of disconnected parts).

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to check.

    Returns:
    bool: True if the molecule is fragmented, False otherwise.
    """

    fragments = rdmolops.GetMolFrags(mol, asMols=False)

    return len(fragments) > 1

if __name__ == '__main__':
    smiles = "C1=CC2=C(C=C1)C=CC=C2"  # Naphthalene
    analyzer = FusedRingAnalyzer(smiles)
    print("Max fused rings:", analyzer.find_max_fused_rings())