"""
fix_cyclic_peptide.py
---------------------
A script to remove overlapping terminal atoms (N-term and OXT) from a cyclic peptide
and optionally add a LINK record to indicate a new peptide bond between the first
and last residue.
"""
import os
from Bio.PDB import PDBParser, PDBIO, Select
from protein_utils import parse_structure

class CyclicPeptideFixer:
    def __init__(self, pdb_path, output_path=None, chain_id='A', add_link_record=True):
        """
        :param pdb_path: Input PDB file (with a presumably cyclic peptide)
        :param output_path: Output PDB file after fixing
        :param chain_id: Which chain to fix (default 'A')
        :param add_link_record: Whether to append a LINK line to indicate the new bond
        """
        self.pdb_path = pdb_path
        self.output_path = output_path
        self.chain_id = chain_id
        self.add_link_record = add_link_record
        if not output_path:
            self.output_path = pdb_path 

    def fix_cyclic(self):
        # Parse structure
        structure = parse_structure(self.pdb_path, name="cyclic_pep")

        # get the chain (assume only one model)
        model = structure[0]
        if self.chain_id not in model:
            raise ValueError(f"Chain {self.chain_id} not found in {self.pdb_path}")
        chain = model[self.chain_id]

        # get the first and last residues
        residues = [res for res in chain]
        residues.sort(key=lambda r: r.id[1])  

        first_res = residues[0]
        last_res = residues[-1]

        # remove the OXT atom from the last residue, (if it exists) and detach it from the structure 
        # record the last residue's C atom (for LINK line)
        last_res_atoms = list(last_res.get_atoms())
        c_atom_last = None
        for atom in last_res_atoms:
            if atom.id == "OXT":
                last_res.detach_child("OXT")
            if atom.id == "C":
                c_atom_last = atom

        # Trim the first residue's N-terminus
        # and record the first residue's N atom coordinates
        first_res_atoms = list(first_res.get_atoms())
        n_atom_first = None
        n_side_h = []
        for atom in first_res_atoms:
            if atom.id == "N":
                n_atom_first = atom
            if atom.element == "H":
                connected_atoms = self._find_connected_heavy_atom(atom, first_res_atoms)
                if connected_atoms and any(a.id == "N" for a in connected_atoms):
                    n_side_h.append(atom)

        # If there are multiple N-terminal hydrogens, keep only one
        # (assuming at most 2~3 hydrogens at the N-terminus, if more precise logic is needed, modify it)
        if len(n_side_h) > 1:
            n_side_h.sort(key=lambda x: x.id)
            keep = n_side_h[0]
            for h in n_side_h[1:]:
                first_res.detach_child(h.id)

        io = PDBIO()
        io.set_structure(structure)
        temp_output = self.output_path
        io.save(temp_output, select=SelectAllButTER())  # 不写TER

        # if need link record, add it to the end of the PDB file
        #    "LINK         N   VAL A   1                C   VAL A  14     1555   1555  1.32"
        #    Common way：LINK  atom1 + resid1 + chain1 + atom2 + resid2 + chain2
        if self.add_link_record and c_atom_last and n_atom_first:
            self._append_link_record(temp_output, c_atom_last, n_atom_first)

        print(f"[Done] Saved fixed cyclic PDB to: {temp_output}")

    def _append_link_record(self, pdb_file, c_atom_last, n_atom_first):

        # residue info
        last_res = c_atom_last.get_parent()
        first_res = n_atom_first.get_parent()

        # construct LINK line
        # columns 13-16 = "LINK", 17=space, 18-20=atom name, 22=altLoc, 23-26=resName, 27=space, 28=chainID, 29-31=resSeq, 32=iCode ...)
        link_line = (
            f"LINK         {c_atom_last.id:<4}{last_res.resname:>4} {self.chain_id}"
            f"{last_res.id[1]:>4}                "
            f"{n_atom_first.id:<4}{first_res.resname:>4} {self.chain_id}"
            f"{first_res.id[1]:>4}      1.33\n"
        )

        with open(pdb_file, "a") as f:
            f.write(link_line)

    def _find_connected_heavy_atom(self, h_atom, residue_atoms):
        """
        English:
        Simple method to find the heavy atom connected to a hydrogen atom
        (Biopython does not store CONNECT information by default, so we use coordinates + threshold)
        """
        import numpy as np
        threshold = 1.2  # typical H-X bond length
        h_coord = h_atom.get_vector()
        connected = []
        for atom in residue_atoms:
            if atom.element != "H":
                dist = (h_coord - atom.get_vector()).norm()
                if dist < threshold:
                    connected.append(atom)
        return connected


class SelectAllButTER(Select):
    """English: Custom selector to write all atoms but not the internal TER lines."""
    def accept_residue(self, residue):
        return True


