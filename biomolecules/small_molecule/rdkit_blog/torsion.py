from rdkit import Chem
from rdkit.Chem import AllChem
import copy
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import AllChem, TorsionFingerprints


def GetEnergy(Mol, ConfID = None, FF_type='MMFF'):
    '''
    Calculate the molecular energy by empirical force field
    '''
    Status = True
    Energy = None

    if ConfID is None:
        ConfID = -1
    
    if FF_type == "UFF":
        UFFMoleculeForcefield = AllChem.UFFGetMoleculeForceField(Mol, confId = ConfID)
        if UFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = UFFMoleculeForcefield.CalcEnergy()
    elif FF_type == "MMFF":
        #In the RDKit, MMFF has two options, “MMFF94” / “MMFF94s”.
        MMFFMoleculeProperties = AllChem.MMFFGetMoleculeProperties(Mol)
        MMFFMoleculeForcefield = AllChem.MMFFGetMoleculeForceField(Mol, MMFFMoleculeProperties, confId = ConfID)
        if MMFFMoleculeForcefield is None:
            Status = False
        else:
            Energy = MMFFMoleculeForcefield.CalcEnergy()
    else:
        print('Only support the UFF or MMFF!')
    
    return (Status, Energy)


def rigid_torsion_scan(mol, start_angle, end_angle, step_angle, torsion_ids, addHs=True, \
    gen_init_conf=False, FF_option='MMFF'):
    '''
    Perform the rigid torsion scan
    torsion_id: a list which contains 4 int to specify the dihedral angle
    '''
    idx_1, idx_2, idx_3, idx_4 = torsion_ids
    dihedral_energies = []
    angles = [Angle for Angle in range(start_angle, end_angle, step_angle)]
    if addHs:
        mol = Chem.AddHs(mol)
    if gen_init_conf:
        AllChem.EmbedMolecule(mol)
        if FF_option == 'MMFF':
            AllChem.MMFFOptimizeMolecule(mol)
        elif FF_option == 'UFF':
            AllChem.UFFOptimizeMolecule(mol)
        else:
            print('No optimization since the FF option is not supported')
    mol_ = copy.deepcopy(mol)
    for dihedral_degree in angles:
        Torsionmol = Chem.Mol(mol_)
        conformer = Torsionmol.GetConformer(0)
        rdMolTransforms.SetDihedralDeg(conformer, idx_1, idx_2, idx_3, idx_4, dihedral_degree)
        CalcStatus, Energy = GetEnergy(Torsionmol,FF_type=FF_option)
        dihedral_energies.append(Energy)
    return angles, dihedral_energies


def print_torsions(mol):
    nonring, ring = TorsionFingerprints.CalculateTorsionLists(mol)
    conf = mol.GetConformer(id=0)
    tups = [atoms[0] for atoms, ang in nonring]

    degs = [Chem.rdMolTransforms.GetDihedralDeg(conf, *tup) for tup in tups]
    print(degs)

def get_torsion(mol):
    conf = mol.GetConformer(id=0)
    nonring, ring = TorsionFingerprints.CalculateTorsionLists(mol)
    tups = [atoms[0] for atoms, ang in nonring]
    degs = [Chem.rdMolTransforms.GetDihedralDeg(conf, *tup) for tup in tups]
    qudra_types = [''.join([mol.GetAtomWithIdx(i).GetSymbol() for i in tup]) for tup in tups]
    qudra_types_rever = [''.join([mol.GetAtomWithIdx(i).GetSymbol() for i in tup][::-1]) for tup in tups]
    return degs+degs, qudra_types+qudra_types_rever

def statis_torsion(mols):
    
    statis = {}
    for mol in mols:
        degs, qudra_types = get_torsion(mol)
        for i, qudra_type in enumerate(qudra_types):
            if qudra_type not in statis:
                statis[qudra_type] = [degs[i]]
            else:
                statis[qudra_type].append(degs[i])
    return statis

from itertools import combinations

def get_angles(mol):
    angles = []
    types = [] # double counting
    conf =  mol.GetConformer(id=0)
    for atom in mol.GetAtoms():
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(neighbors) < 2:  # At least two neighbors are required to form an angle
            continue
        for a1, a2 in combinations(neighbors, 2):
            angle = rdMolTransforms.GetAngleDeg(conf, a1, atom.GetIdx(), a2)
            tup = (a1, atom.GetIdx(), a2)
            type_ = [mol.GetAtomWithIdx(i).GetSymbol() for i in tup]
            type = ''.join(type_)
            type_rever = ''.join(type_[::-1])

            angles.append(angle)
            angles.append(angle)
            types.append(type)
            types.append(type_rever)
    return types, angles

def statis_angles(mols):

    statis = {}
    for mol in mols:
        types, angles = get_angles(mol)
        for i, tetra_type in enumerate(types):
            if tetra_type not in statis:
                statis[tetra_type] = [angles[i]]
            else:
                statis[tetra_type].append(angles[i])
    return statis