import numpy as np
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

atom_mapping = {"Hydrophobe": 0,
                "LumpedHydrophobe": 0,
                "Aromatic": 1,
                "Acceptor": 2,
                "Donor": 3,
                "PosIonizable": 4,
                "NegIonizable": 5}

radiidict = {'Ac': 2.0, 'Ag': 1.72, 'Al': 2.0, 'Am': 2.0, 'Ar': 1.88, 'As': 1.85, 'At': 2.0, 'Au': 1.66, 'B': 2.0,
             'Ba': 2.0, 'Be': 2.0, 'Bh': 2.0, 'Bi': 2.0, 'Bk': 2.0, 'Br': 1.85, 'C': 1.7, 'Ca': 1.37, 'Cd': 1.58,
             'Ce': 2.0, 'Cf': 2.0, 'Cl': 2.27, 'Cm': 2.0, 'Co': 2.0, 'Cr': 2.0, 'Cs': 2.1, 'Cu': 1.4, 'Db': 2.0,
             'Ds': 2.0, 'Dy': 2.0, 'Er': 2.0, 'Es': 2.0, 'Eu': 2.0, 'F': 1.47, 'Fe': 2.0, 'Fm': 2.0, 'Fr': 2.0,
             'Ga': 1.07, 'Gd': 2.0, 'Ge': 2.0, 'H': 1.2, 'He': 1.4, 'Hf': 2.0, 'Hg': 1.55, 'Ho': 2.0, 'Hs': 2.0,
             'I': 1.98, 'In': 1.93, 'Ir': 2.0, 'K': 1.76, 'Kr': 2.02, 'La': 2.0, 'Li': 1.82, 'Lr': 2.0, 'Lu': 2.0,
             'Md': 2.0, 'Mg': 1.18, 'Mn': 2.0, 'Mo': 2.0, 'Mt': 2.0, 'N': 1.55, 'Na': 1.36, 'Nb': 2.0, 'Nd': 2.0,
             'Ne': 1.54, 'Ni': 1.63, 'No': 2.0, 'Np': 2.0, 'O': 1.52, 'Os': 2.0, 'P': 1.8, 'Pa': 2.0, 'Pb': 2.02,
             'Pd': 1.63, 'Pm': 2.0, 'Po': 2.0, 'Pr': 2.0, 'Pt': 1.72, 'Pu': 2.0, 'Ra': 2.0, 'Rb': 2.0, 'Re': 2.0,
             'Rf': 2.0, 'Rg': 2.0, 'Rh': 2.0, 'Rn': 2.0, 'Ru': 2.0, 'S': 1.8, 'Sb': 2.0, 'Sc': 2.0, 'Se': 1.9,
             'Sg': 2.0, 'Si': 2.1, 'Sm': 2.0, 'Sn': 2.17, 'Sr': 2.0, 'Ta': 2.0, 'Tb': 2.0, 'Tc': 2.0, 'Te': 2.06,
             'Th': 2.0, 'Ti': 2.0, 'Tl': 1.96, 'Tm': 2.0, 'U': 1.86, 'V': 2.0, 'W': 2.0, 'X': 1.5, 'Xe': 2.16, 'Y': 2.0,
             'Yb': 2.0, 'Zn': 1.39, 'Zr': 2.0}

def _getAtomTypes(smallmol):
    _mol = smallmol._mol
    n_atoms = smallmol.numAtoms
    feats = factory.GetFeaturesForMol(_mol)
    properties = np.zeros((n_atoms, 8), dtype=bool)

    for feat in feats:
        fam = feat.GetFamily()
        if fam not in atom_mapping:  # Non relevant property
            continue
        properties[feat.GetAtomIds(), atom_mapping[fam]] = 1

    # Occupancy, ignoring hydrogens.
    properties[:, 7] = smallmol.get('element') != 'H'
    return properties


def _getChannelRadii(smallmol):
    radii = np.vectorize(radiidict.__getitem__)(smallmol.get('element')) * _getAtomTypes(smallmol).T
    return radii.T.copy()