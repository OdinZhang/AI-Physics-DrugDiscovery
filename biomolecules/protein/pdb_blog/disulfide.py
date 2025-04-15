from Bio.PDB import PDBParser, NeighborSearch


pdb_file = './1g9i.pdb'
def find_diSContact(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb',pdb_file)
    disulfide_bonds = []

    for model in structure:
        atoms = [atom for atom in model.get_atoms()]
        ns = NeighborSearch(atoms)
        for chain in model:
            cysteines = [residue for residue in chain if residue.get_resname() == 'CYS']
            for cysteine in cysteines:
                sg_atom = cysteine['SG']
                neighbors = ns.search(sg_atom.get_coord(), 3.0)
                for neighbor in neighbors:
                    if neighbor.get_parent().get_resname() == 'CYS' and neighbor.get_name() == 'SG':
                        bonded_cys = neighbor.get_parent()
                        disulfide_bonds.append((cysteine, bonded_cys))

    disulfide_bonds = list(set(tuple(sorted(pair)) for pair in disulfide_bonds))
    disulfide_list = []
    for bond in disulfide_bonds:
        cys1, cys2 = bond
        cys1_id, cys1_sg_coord = cys1.get_id(), cys1['SG'].get_coord()
        cys2_id, cys2_sg_coord = cys2.get_id(), cys2['SG'].get_coord()
        if cys1_id[1]!=cys2_id[1]:
            disulfide_list.append((cys1_id[1], cys2_id[1]))
    return disulfide_list
