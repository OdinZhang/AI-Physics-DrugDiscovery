from rdkit import Chem
from rdkit.Chem import AllChem
import torch
from collections import defaultdict
from collections import deque

def get_atom_map_3d(mol, frag, verbose=True):
    """
    Find the frag-mol atomic maps.
    """
    coords_mol = get_mol_coord(mol)
    coords_frag = get_mol_coord(frag)

    distance_matrix = torch.norm(coords_frag[:, None, :] - coords_mol[None, :, :], dim=2)
    epsilon = 0.01
    atom_mapping = torch.nonzero(distance_matrix < epsilon, as_tuple=True)
    atom_mapping_list = list(zip(atom_mapping[0].tolist(), atom_mapping[1].tolist()))
    if len(atom_mapping_list) != frag.GetNumAtoms() and verbose:
        print('frag do not exactly match the total molecule')
    return atom_mapping_list

def get_mol_coord(mol):
    conf = mol.GetConformer()
    return torch.tensor([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    tmp_mol = Chem.Mol(mol)
    for idx in range(atoms):
        tmp_mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(tmp_mol.GetAtomWithIdx(idx).GetIdx()))
    return tmp_mol

def read_sdf(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file)
    mols_list = [i for i in supp]
    return mols_list

def write_sdf(mol_list,file):
    writer = Chem.SDWriter(file)
    for i in mol_list:
        writer.write(i)
    writer.close()

def remove_atoms_by_idx(mol, atom_indices_to_remove):
    """
    Remove specified atoms from a molecule based on their indices.

    :param mol: Original RDKit molecule object.
    :param atom_indices_to_remove: List of atom indices to remove.
    :return: A new RDKit molecule object with specified atoms removed.
    """

    atom_indices_to_remove = set(atom_indices_to_remove)
    emol = Chem.EditableMol(mol)
    
    for idx in sorted(atom_indices_to_remove, reverse=True):
        emol.RemoveAtom(idx)
    
    new_mol = emol.GetMol()

    Chem.SanitizeMol(new_mol)
    
    return new_mol


def merge_fragments(results):
    
    fragment_dict = defaultdict(set)  

    for fragment, attachment in results:
        fragment_tuple = tuple(fragment)  
        fragment_dict[fragment_tuple].add(attachment)

    merged_results = [(list(frag), list(attachments)) for frag, attachments in fragment_dict.items()]
    return merged_results

def attribute_growing_frags(mol, frag):
    '''
    Find the growing fragment and the corresponding attachment points in the molecule.
    Return:
    - list of tuples: (fragment atoms (mol index), attachment atoms (frag index))
    '''
    atom_mapping = get_atom_map_3d(mol, frag)
    adjacency_matrix = Chem.GetAdjacencyMatrix(mol)
    atom_mapping_dict = {f: m for f, m in atom_mapping}
    reversed_atom_mapping_dict = {m: f for f, m in atom_mapping}
    visited = set()  # Global set to keep track of visited atoms
    result = []
    ori_pool_set = set(atom_mapping_dict.values())  # Original fragment atoms index

    def bfs(start_atom):
        if start_atom in ori_pool_set:
            return set()  # Skip BFS for original fragment atoms

        frag_pool = [start_atom]
        ori_pool = []
        queue = deque([start_atom])  # Use deque for efficient pop from left
        local_visited = set([start_atom])  # Track atoms visited in this BFS

        while queue:
            current_atom = queue.popleft()
            neighbors = list(adjacency_matrix[current_atom].nonzero()[0])

            for neighbor in neighbors:
                if neighbor not in local_visited:
                    local_visited.add(neighbor)
                    if neighbor in ori_pool_set:
                        # Only add to ori_pool if this is a boundary atom
                        if any(adjacency_matrix[neighbor, n] and n not in ori_pool_set for n in range(mol.GetNumAtoms())):
                            ori_pool.append(neighbor)
                    else:
                        frag_pool.append(neighbor)
                        queue.append(neighbor)

        return frag_pool, ori_pool, local_visited

    for start_atom in range(mol.GetNumAtoms()):
        if start_atom not in visited and start_atom not in ori_pool_set:
            frag_pool, ori_pool, new_visited = bfs(start_atom)
            visited.update(new_visited)
            if frag_pool and ori_pool:
                for attachment in ori_pool:
                    result.append((frag_pool, reversed_atom_mapping_dict[attachment]))

    result = merge_fragments(result)

    return result


def remove_atoms_by_idx(mol, atom_indices_to_remove):
    """
    Remove specified atoms from a molecule based on their indices.

    :param mol: Original RDKit molecule object.
    :param atom_indices_to_remove: List of atom indices to remove.
    :return: A new RDKit molecule object with specified atoms removed.
    """

    atom_indices_to_remove = set(atom_indices_to_remove)
    emol = Chem.EditableMol(mol)
    
    for idx in sorted(atom_indices_to_remove, reverse=True):
        emol.RemoveAtom(idx)
    
    new_mol = emol.GetMol()

    Chem.SanitizeMol(new_mol)
    
    return new_mol


def filter_genmol_ac_attach(gen_mol, frag, keep_attachs):
    '''
    Filter the gen_mol by removing the atoms that are not connected to the GIVEN attachment points.
    Input:
        gen_mol and frag should be the 3D mol object, otherwise we should use another way to determine the atom mapping
        keep_attachs is from the original fragment
    Output:
        filtered gen_mol
    '''
    frag_on_attachs = attribute_growing_frags(gen_mol, frag)

    remove_atoms = []
    for frag_on_attach in frag_on_attachs:
        for keep_attach in keep_attachs:
            if keep_attach in frag_on_attach[1]:
                continue
            else:
                remove_atoms.extend(frag_on_attach[0])
    remove_atoms = [int(i) for i in remove_atoms]
    filtered_genmol = remove_atoms_by_idx(gen_mol, remove_atoms)
    return filtered_genmol
if __name__ == '__main__':
    gen_mol = read_sdf('./5n2d_gen.sdf')[0]
    frag = read_sdf('./5n2d_frag.sdf')[0]
    filtered_mol = filter_genmol_ac_attach(gen_mol, frag, [6])