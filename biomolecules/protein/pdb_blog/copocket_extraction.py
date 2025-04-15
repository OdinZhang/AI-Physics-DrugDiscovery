from Bio.PDB import PDBParser, Superimposer
import numpy as np
import pymol
from pymol import cmd
from Bio import pairwise2

from Bio import pairwise2
from Bio.PDB import PDBParser
import pymol
from pymol import cmd
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def get_sequence(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if residue.get_resname() in three_to_one:
                    seq += three_to_one[residue.get_resname()]
            return seq

def map_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]  # 假设第一个对齐是最佳对齐
    aligned_seq1, aligned_seq2, score, begin, end = best_alignment
    return aligned_seq1, aligned_seq2, begin, end


def get_pocket_residues_from_pymol(selection_name):
    pocket_residues = []
    pymol.cmd.iterate(selection_name, 'pocket_residues.append((resn, resi))', space=locals())
    return pocket_residues

def map_pocket_residues(crystal_pocket_residues, alignment):
    mapped_residues = []
    crystal_seq, alphafold_seq, start, end = alignment
    for res_name, res_id in crystal_pocket_residues:
        res_id_int = int(res_id)
        if crystal_seq[res_id_int - start] == alphafold_seq[res_id_int - start]:
            mapped_residues.append((res_name, res_id))
    return mapped_residues

def create_pocket_selection(structure, mapped_residues, selection_name):
    selection_query = " or ".join([f"resi {res_id} and resn {res_name}" for res_name, res_id in mapped_residues])
    cmd.select(selection_name, selection_query)


def copocket_extraction(crystal_file, alpha_file, ligand_file):
    # Load structures
    cmd.load(crystal_file, "crystal")
    cmd.load(alpha_file, "alphafold")
    cmd.load(ligand_file, "ligand")

    ligand_name = 'ligand'
    cmd.align("alphafold", "crystal")

    # Define the radius for the pocket selection
    radius = 10
    cmd.select("crystal_pocket", f"br. crystal within {radius} of {ligand_name}")

    # Sequence alignment
    crystal_seq = get_sequence(crystal_file)
    alphafold_seq = get_sequence(alpha_file)
    alignment = map_sequences(crystal_seq, alphafold_seq)

    # Pocket residue mapping
    crystal_pocket_residues = get_pocket_residues_from_pymol("crystal_pocket")
    mapped_residues = map_pocket_residues(crystal_pocket_residues, alignment)

    # Save pocket structures
    parser = PDBParser(QUIET=True)
    alphafold_structure = parser.get_structure('pdb', alpha_file)
    create_pocket_selection(alphafold_structure, mapped_residues, "aligned_pocket")
    cmd.save("alphafold_protein_pocket.pdb", "aligned_pocket")
    cmd.save("crystal_protein_pocket.pdb", "crystal_pocket")


if __name__ == '__main__':
    ligand_file = './4gv1_ori.sdf'
    crystal_file = './4gv1_protein.pdb'
    alpha_file = './4gv1_alpha.pdb'
    copocket_extraction(crystal_file, alpha_file, ligand_file)