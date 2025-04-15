from rdkit import Chem
from rdkit.Chem import rdMolTransforms

def get_bond_type(bond):
    if bond is None:
        return None
    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
        return '-'
    elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
        return '='
    elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
        return '#'
    else:
        return '?'

def generate_angle_list(mol, reverse_count=True):
    '''
    Generates a list of all possible three-atom combination angles in the molecule, including the angle type, indices, and degree.
    
    Parameters:
        mol (rdkit.Chem.Mol): An RDKit molecule object.
        reverse_count (bool, optional): If True, for each found angle, its reverse angle type and indices will also be added. Defaults to False.
    
    Returns:
        list: A list containing tuples of the format ('angle type', (atom index1, atom index2, atom index3), degree), where the degree is rounded to two decimal places.
    
    Example:
        Returns a list of angles, e.g., [('C-C=C',(0,1,4),32.33), ('C=C-C',(4,1,0),32.33), ...].
    '''
    angle_list = []
    conf = mol.GetConformer()
    for idx1 in range(mol.GetNumAtoms()):
        for idx2 in range(mol.GetNumAtoms()):
            for idx3 in range(mol.GetNumAtoms()):
                if idx1 < idx2 and idx2 < idx3:
                    angle = round(rdMolTransforms.GetAngleDeg(conf, idx1, idx2, idx3), 2)
                    atom1, atom2, atom3 = mol.GetAtomWithIdx(idx1), mol.GetAtomWithIdx(idx2), mol.GetAtomWithIdx(idx3)
                    bond1_type = get_bond_type(mol.GetBondBetweenAtoms(idx1, idx2))
                    bond2_type = get_bond_type(mol.GetBondBetweenAtoms(idx2, idx3))
                    if bond1_type and bond2_type:
                        angle_type = f"{atom1.GetSymbol()}{bond1_type}{atom2.GetSymbol()}{bond2_type}{atom3.GetSymbol()}"
                        angle_list.append((angle_type, (idx1, idx2, idx3), angle))
                        if reverse_angle_type:
                            reverse_angle_type = f"{atom3.GetSymbol()}{bond2_type}{atom2.GetSymbol()}{bond1_type}{atom1.GetSymbol()}"
                            angle_list.append((reverse_angle_type, (idx3, idx2, idx1), angle))
    return angle_list

if __name__ == '__main__':

    mol = Chem.MolFromMolFile('./druglike_property/test.sdf')
    Chem.Kekulize(mol)

    angle_list = generate_angle_list(mol,reverse_count=True)
    print(angle_list)