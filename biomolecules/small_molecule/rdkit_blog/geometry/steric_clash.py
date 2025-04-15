from rdkit import Chem
from pdb_parser import PDBProtein
import numpy as np
import numpy as np
from scipy.spatial.distance import cdist

default_vdw_radii = {
    1: 1.2,  # Hydrogen
    6: 1.7,  # Carbon
    7: 1.55, # Nitrogen
    8: 1.52, # Oxygen
    9: 1.47, # Fluorine
    15: 1.8, # Phosphorus
    16: 1.8, # Sulfur
    17: 2.27, # Chlorine
    35: 1.85 # Bromine
}


def read_sdf(file_path):
    supplier = Chem.SDMolSupplier(file_path)
    molecules = []
    for mol in supplier:
        if mol is not None:
            molecules.append(mol)
    return molecules

def parse_sdf_file(input_mol):
    if type(input_mol) == str:
        mol = read_sdf(input_mol)[0]
    else:
        mol = input_mol
    mol_info = {}
    atomic_type = []
    atomic_number = []
    atomic_coords = []
    # Iterate through each atom in the molecule
    for atom in mol.GetAtoms():
        atomic_type.append(atom.GetSymbol())
        atomic_number.append(atom.GetAtomicNum())
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        atomic_coords.append((pos.x, pos.y, pos.z))

    mol_info['atom_name'] = atomic_type
    mol_info['element'] = np.array(atomic_number)
    mol_info['pos'] = np.array(atomic_coords)
    return mol_info

def steric_clash(mol, pdb_file, vdw_radii=default_vdw_radii, tolerance=0.1):
    '''
    This function detects steric clashes between a ligand and a protein based on VdW radii.
    Input:
        mol: RDKit molecule object or path to SDF file containing ligand
        pdb_file: path to PDB file containing protein
        vdw_radii: dictionary of atomic numbers and corresponding VdW radii
        tolerance: tolerance for steric clash detection
    Output:
        clash_detected: boolean indicating whether a steric clash is detected
        additional_info: dictionary containing additional information
    Example:
        clash_detected, additional_info = steric_clash('ligand.sdf', 'protein.pdb')
    '''
    ligand_info = parse_sdf_file(mol)
    protein_info = PDBProtein(pdb_file).to_dict_atom()
    protein_pos = np.array(protein_info['pos'])
    ligand_pos = np.array(ligand_info['pos'])
    protein_elements = protein_info['element']
    ligand_elements = ligand_info['element']

    # Compute pairwise distances between ligand and protein atoms
    distances = cdist(ligand_pos, protein_pos)

    # Retrieve VdW radii based on atomic numbers
    protein_vdw_radii = np.array([vdw_radii[a] for a in protein_elements])
    ligand_vdw_radii = np.array([vdw_radii[a] for a in ligand_elements])

    # Compute the sum of VdW radii for each atom pair
    vdw_sums_with_tolerance = ligand_vdw_radii[:, np.newaxis] + protein_vdw_radii - tolerance

    # Check for steric clashes (distance < sum of VdW radii)
    clashes = distances < vdw_sums_with_tolerance
    # Determine if any clashes are detected
    clash_detected = np.any(clashes)

    # Get indices where clashes occur
    clash_indices = np.where(clashes)
    clashed_distances = distances[clash_indices]
    clashed_vdw_sums = vdw_sums_with_tolerance[clash_indices]

    # Additional information to return
    additional_info = {
        'clashed_indices': clash_indices,
        'clashed_distances': clashed_distances,
        'clashed_vdw_sums': clashed_vdw_sums
    }

    return clash_detected, additional_info

if __name__ == '__main__':
    pdb_file = '1a2g_A_rec.pdb'
    sdf_file = '1a2g_A_rec_4jmv_1ly_lig_tt_min_0.sdf'
    clash_detected, additional_info = steric_clash(sdf_file, pdb_file)
    print(clash_detected)
    print(additional_info)