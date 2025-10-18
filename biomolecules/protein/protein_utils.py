import numpy as np
import os
import warnings
from Bio.PDB import MMCIFParser, is_aa, Superimposer, PDBParser
from Bio.PDB import MMCIFParser, PDBIO, MMCIFIO, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Polypeptide import PPBuilder, three_to_one
from Bio.Data.IUPACData import protein_letters_3to1_extended
class SafeMMCIFParser(MMCIFParser):
    def get_structure(self, structure_id, filename):
        mmcif_dict = MMCIF2Dict(filename)

        def ensure(key, default):
            if key not in mmcif_dict:
                atom_count = len(mmcif_dict.get("_atom_site.label_atom_id", []))
                mmcif_dict[key] = [default] * atom_count

        ensure("_atom_site.B_iso_or_equiv", "0.0")
        ensure("_atom_site.occupancy", "1.0")
        ensure("_atom_site.pdbx_PDB_ins_code", "?")
        ensure("_atom_site.label_alt_id", ".")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self._mmcif_dict = mmcif_dict
            self._build_structure(structure_id)
            self._structure_builder.set_header(self._get_header())

        return self._structure_builder.get_structure()

def parse_structure(path, name="structure"):
    ext = os.path.splitext(path)[-1].lower()
    if ext == '.cif':
        parser = SafeMMCIFParser(QUIET=True)
    elif ext == '.pdb':
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Don't support: {ext}")
    return parser.get_structure(name, path)

def get_ca_atoms(structure):
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True) and 'CA' in residue:
                    ca_atoms.append(residue['CA'])
        break  
    return ca_atoms

def compute_ca_rmsd(file1, file2):
    structure1 = parse_structure(file1, "structure1")
    structure2 = parse_structure(file2, "structure2")

    ca_atoms1 = get_ca_atoms(structure1)
    ca_atoms2 = get_ca_atoms(structure2)

    if len(ca_atoms1) != len(ca_atoms2):
        raise ValueError(f"Cα numbers don't match: {len(ca_atoms1)} vs {len(ca_atoms2)}")

    sup = Superimposer()
    sup.set_atoms(ca_atoms1, ca_atoms2)
    return sup.rms

def get_sequence_from_chain(structure, chain_ids):
    """
    Extract 1-letter sequence(s) from specified chain(s) in a structure.
    Returns a dictionary: {chain_id: sequence}
    """
    seqs = {}
    model = structure[0]  # default: first model

    for chain_id in chain_ids:
        if chain_id not in model:
            raise ValueError(f"Chain {chain_id} not found in structure.")
        chain = model[chain_id]
        seq = ""
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            try:
                aa = three_to_one(residue.get_resname())
            except KeyError:
                # fallback for non-standard AAs
                aa = protein_letters_3to1_extended.get(residue.get_resname(), "X")
            seq += aa
        seqs[chain_id] = seq
    return seqs

def get_sequence_from_file(filepath, chain_ids):
    """
    Parse a structure file (.pdb or .cif) and return sequence(s) from specified chain(s).
    Example:
        filepath = "path/to/structure.pdb"
        chain_ids = ["A", "B"]
        sequences = get_sequence_from_file(filepath, chain_ids)
        print(sequences)
    Returns a dictionary: {chain_id: sequence}
    """
    structure = parse_structure(filepath)
    return get_sequence_from_chain(structure, chain_ids)

def extract_conformers(file, output_dir, file_type='pdb'):
    """
    Extract individual models (conformers) from a PDB or mmCIF file
    and save them as separate files in the specified format.

    Args:
        file (str): Input PDB or mmCIF file
        output_dir (str): Directory to save individual models
        file_type (str): 'pdb' or 'cif', default is 'pdb'
    """

    structure = parse_structure(file)
    os.makedirs(output_dir, exist_ok=True)

    for i, model in enumerate(structure):
        base_name = os.path.splitext(os.path.basename(file))[0]
        out_name = f"{base_name}_{i}.{file_type}"
        out_path = os.path.join(output_dir, out_name)

        if file_type == 'pdb':
            io = PDBIO()
            io.set_structure(model)
            io.save(out_path)
        elif file_type == 'cif':
            io = MMCIFIO()
            io.set_structure(model)
            io.save(out_path)
        else:
            raise ValueError("Unsupported file_type. Use 'pdb' or 'cif'.")
