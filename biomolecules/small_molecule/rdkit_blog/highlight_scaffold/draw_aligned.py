from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def draw_align_mols(mol_list, sharing_frag, align_mol=None,highlight_choice=None):
    '''
    Parameters:
    mol_list: list of mols to be aligned
    sharing_frag: the fragment to be aligned
    align_mol: the mol to be aligned to, could be None
    highlight_choice: None, 'match', 'remain'
    '''
    match_mols = [x for x in mol_list if x.HasSubstructMatch(sharing_frag)]
    if len(match_mols) == 0:
        print('No match')
        return
    if len(match_mols) < len(mol_list):
        print('some mols do not have the sharing fragment')
        
    remain_idxes = []
    aligned_mols = []
    match_idxes = []
    for mol in match_mols:
        match = mol.GetSubstructMatch(sharing_frag)
        AllChem.Compute2DCoords(mol)
        mol_copy = Chem.Mol(mol)
        if align_mol is not None:
            AllChem.AlignMol(mol_copy, align_mol, atomMap=list(zip(match, aligned_mols[0].GetSubstructMatch(sharing_frag))))
        else:
            if aligned_mols:
                AllChem.AlignMol(mol_copy, aligned_mols[0], atomMap=list(zip(match, aligned_mols[0].GetSubstructMatch(sharing_frag))))
        
        aligned_mols.append(mol_copy)
        match_idxes.append(match)
        match_set = set(match)
        remain_set = set(range(mol.GetNumAtoms())) - match_set
        remain_idxes.append(tuple(remain_set))

    if highlight_choice is None:
        img = Draw.MolsToGridImage(aligned_mols, molsPerRow=4, subImgSize=(200,200), legends=[mol.GetProp('_Name') for mol in aligned_mols])
    if highlight_choice == 'match':
        img = Draw.MolsToGridImage(aligned_mols, molsPerRow=4, subImgSize=(200,200), legends=[mol.GetProp('_Name') for mol in aligned_mols], highlightAtomLists=match_idxes)
    elif highlight_choice == 'remain':
        img = Draw.MolsToGridImage(aligned_mols, molsPerRow=4, subImgSize=(200,200), legends=[mol.GetProp('_Name') for mol in aligned_mols], highlightAtomLists=remain_idxes)
    return img 