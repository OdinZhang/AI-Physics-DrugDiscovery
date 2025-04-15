# Odin Zhang 20240422
from rdkit import Chem
from rdkit.Chem import AllChem 
import os.path as osp
from frag import remove_dummys_mol, transfer_coord, linkerize_mol, fragmentize_mol,read_sdf, frag2mols
from frag import check_frags, Murcko_decompose_anchor, check_linkers
from frag import check_frags, qsmis, find_anchor_indices_3d, get_atom_map_3d, get_mol_coord, rm_radical, remove_substructure

import torch

def linker_decomp(mol):
    '''
    Input: mol
    Output: decomp_infos
    '''
    Chem.SanitizeMol(mol)
    fragmentations = linkerize_mol(mol)
    fragmentations = check_linkers(fragmentations)
    decomp_infos = []
    for frag_idx, fragmentation in enumerate(fragmentations):
        try:
            linker, frags = fragmentation
            
            frag_mols = qsmis([frags])[0]
            frag_mols, _ = remove_dummys_mol(frag_mols)
            frag_mols, _ = remove_dummys_mol(frag_mols)
            frag_3d = transfer_coord(frag_mols, mol)
            anchor_in_frag, anchor_in_mol = find_anchor_indices_3d(mol, frag_3d)
            # frag1_id, frag2_id = anchor_in_frag 
            # frag_3d.SetProp('anchor_ids',f'{frag1_id}_{frag2_id}')
            frag_atom_mapping = get_atom_map_3d(mol, frag_3d)
            frag_in_mol_index = [i[1] for i in frag_atom_mapping]
            frag_atoms = torch.tensor([frag_3d.GetAtomWithIdx(i).GetAtomicNum() for i in range(frag_3d.GetNumAtoms())])
            frag_pos = get_mol_coord(frag_3d)

            linker_mol = qsmis([linker])[0]
            linker_mol, _ = remove_dummys_mol(linker_mol)
            linker_mol, _ = remove_dummys_mol(linker_mol)
            linker_mol = rm_radical(linker_mol)
            linker_3d = transfer_coord(linker_mol, mol)
            linker_atom_mapping = get_atom_map_3d(mol, linker_3d)
            linker_in_mol_index = [i[1] for i in linker_atom_mapping]
            linker_atoms = torch.tensor([linker_3d.GetAtomWithIdx(i).GetAtomicNum() for i in range(linker_3d.GetNumAtoms())])
            linker_pos = get_mol_coord(linker_3d)
            
            decomp_info = {
                'linker_mol': linker_3d,
                'linker_atoms': linker_atoms,
                'linker_pos': linker_pos,
                'linker_in_mol_index': linker_in_mol_index,
                'frag_mol': frag_3d,
                'frag_atoms': frag_atoms,
                'frag_pos': frag_pos,
                'frag_in_mol_index': frag_in_mol_index,
                'anchor_id': anchor_in_mol # fragment is the condition
            }
            decomp_infos.append(decomp_info)
        except Exception as e:
            print(e)
    return decomp_infos

def fragment_decomp(mol):
    Chem.SanitizeMol(mol)
    fragmentations = fragmentize_mol(mol)
    fragmentations = check_frags(fragmentations)
    decomp_infos = []
    for frag_idx, fragmentation in enumerate(fragmentations):
        try:
            frag_small, frag_large = frag2mols(fragmentation)
            frag_small_3d = transfer_coord(frag_small, mol)
            anchor_in_frag, anchor_in_mol = find_anchor_indices_3d(mol, frag_small_3d)
            frag_small_atom_mapping = get_atom_map_3d(mol, frag_small_3d)
            frag_small_in_mol_index = [i[1] for i in frag_small_atom_mapping]
            frag_small_atoms = torch.tensor([frag_small_3d.GetAtomWithIdx(i).GetAtomicNum() for i in range(frag_small_3d.GetNumAtoms())])
            frag_small_pos = get_mol_coord(frag_small_3d)

            frag_large_3d = transfer_coord(frag_large, mol)
            frag_large_atom_mapping = get_atom_map_3d(mol, frag_large_3d)
            frag_large_in_mol_index = [i[1] for i in frag_large_atom_mapping]
            frag_large_atoms = torch.tensor([frag_large_3d.GetAtomWithIdx(i).GetAtomicNum() for i in range(frag_large_3d.GetNumAtoms())])
            frag_large_pos = get_mol_coord(frag_large_3d)
            
            decomp_info = {
                'frag_small_mol': frag_small_3d,
                'frag_small_atoms': frag_small_atoms,
                'frag_small_pos': frag_small_pos,
                'frag_small_in_mol_index': frag_small_in_mol_index,

                'frag_large_mol': frag_large_3d,
                'frag_large_atoms': frag_large_atoms,
                'frag_large_pos': frag_large_pos,
                'frag_large_in_mol_index': frag_large_in_mol_index,
                'anchor_id': anchor_in_mol # large fragment is the condition
            }
            decomp_infos.append(decomp_info)
        except Exception as e:
            print(e)
    return decomp_infos


def scaffold_decompo(mol):
    Chem.SanitizeMol(mol)
    scaffold, _ = Murcko_decompose_anchor(mol)
    scaffold = rm_radical(scaffold)
    scaffold_mol_mapping_list = get_atom_map_3d(mol, scaffold)
    scaffold_in_mol_index = [i[1] for i in scaffold_mol_mapping_list]
    scaffold_atoms = torch.tensor([scaffold.GetAtomWithIdx(i).GetAtomicNum() for i in range(scaffold.GetNumAtoms())])
    scaffold_pos = get_mol_coord(scaffold)

    side_chains, anchor_idx = remove_substructure(mol, sub_mol=scaffold)
    anchor_in_frag, anchor_in_mol = find_anchor_indices_3d(mol, side_chains)
    side_chains_mol_mapping_list = get_atom_map_3d(mol, side_chains)
    side_chains_in_mol_index = [i[1] for i in side_chains_mol_mapping_list]
    side_chains_atoms = torch.tensor([side_chains.GetAtomWithIdx(i).GetAtomicNum() for i in range(side_chains.GetNumAtoms())])
    side_chains_pos = get_mol_coord(side_chains)
    decomp_info = {
        'scaffold_mol': scaffold,
        'scaffold_atoms': scaffold_atoms,
        'scaffold_pos': scaffold_pos,
        'scaffold_in_mol_index': scaffold_in_mol_index,

        'side_chains_mol': side_chains,
        'side_chains_atoms': side_chains_atoms,
        'side_chains_pos': side_chains_pos,
        'side_chains_in_mol_index': side_chains_in_mol_index,
        'anchor_id': anchor_in_mol # side chain is the condition
    }
    return [decomp_info] # only one scaffold decomposition


def side_chains_decomp(mol):
    Chem.SanitizeMol(mol)
    scaffold, _ = Murcko_decompose_anchor(mol)
    scaffold = rm_radical(scaffold)
    scaffold_mol_mapping_list = get_atom_map_3d(mol, scaffold)
    scaffold_in_mol_index = [i[1] for i in scaffold_mol_mapping_list]
    scaffold_atoms = torch.tensor([scaffold.GetAtomWithIdx(i).GetAtomicNum() for i in range(scaffold.GetNumAtoms())])
    scaffold_pos = get_mol_coord(scaffold)

    side_chains, anchor_idx = remove_substructure(mol, sub_mol=scaffold)
    anchor_in_frag, anchor_in_mol = find_anchor_indices_3d(mol, scaffold)
    side_chains_mol_mapping_list = get_atom_map_3d(mol, side_chains)
    side_chains_in_mol_index = [i[1] for i in side_chains_mol_mapping_list]
    side_chains_atoms = torch.tensor([side_chains.GetAtomWithIdx(i).GetAtomicNum() for i in range(side_chains.GetNumAtoms())])
    side_chains_pos = get_mol_coord(side_chains)
    decomp_info = {
        'scaffold_mol': scaffold,
        'scaffold_atoms': scaffold_atoms,
        'scaffold_pos': scaffold_pos,
        'scaffold_in_mol_index': scaffold_in_mol_index,

        'side_chains_mol': side_chains,
        'side_chains_atoms': side_chains_atoms,
        'side_chains_pos': side_chains_pos,
        'side_chains_in_mol_index': side_chains_in_mol_index,
        'anchor_id': anchor_in_mol # side chain is the condition
    }
    return [decomp_info] # only one scaffold decomposition

if __name__ == '__main__':
    # 3D fragmentation

    ligand_nm = './1djy_A_rec_1djz_ip2_lig_tt_min_0.sdf'
    mol = read_sdf(ligand_nm)[0]

    linker_decomp_infos = linker_decomp(mol)
    fragment_decomp_infos = fragment_decomp(mol)
    scaffold_decomp_infos = scaffold_decompo(mol)
    side_chains_decomp_infos = side_chains_decomp(mol)

    print('linker_decomp_infos:', linker_decomp_infos)
    
