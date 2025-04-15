import subprocess
import os
import shutil
from .chem import read_sdf, write_sdf, set_mol_position, sdf2centroid, sdf2mol2, rmfiles
from rdkit import Chem
import numpy as np
import os.path as osp

lepro_bin = '/home/haotian/Molecule_Generation/Docking/classical/software/lepro_linux_x86'
ledock_bin = '/home/haotian/Molecule_Generation/Docking/classical/software/ledock_linux_x86'

def getbox(center, extending= 10.0):
    '''
    center: (x, y, z)
    extending: float
    '''
    cx, cy, cz = center
    minX = cx - float(extending)
    minY = cy - float(extending)
    minZ = cz - float(extending)
    maxX = cx + float(extending)
    maxY = cy + float(extending)
    maxZ = cz + float(extending)

    return {'minX':minX, 'maxX': maxX},{'minY':minY, 'maxY':maxY}, {'minZ':minZ,'maxZ':maxZ}

def lepro(pdb_file, outfile=None):
    '''
    lepro will automatically generate the pro.pdb at the cwd, so we need to relocate it to the outfile manually
    '''
    if outfile is None:
        outfile = pdb_file.replace('.pdb', '_pro.pdb')
    command = f'{lepro_bin} {pdb_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr)
        return None
    cwd = os.getcwd()
    out_pdb = os.path.join(cwd, 'pro.pdb')

    shutil.move(out_pdb, outfile)
    return outfile

def ledock_config_generate(protein_file, ligand_file, out_dir=None):
    if out_dir is None:
        out_dir = osp.dirname(ligand_file)
        
    protein_name = osp.basename(protein_file).split('.')[0]
    ligand_name = osp.basename(ligand_file).split('.')[0]
    ligand_list_name = f'{protein_name}_{ligand_name}.list'
    docking_config_name = f'{protein_name}_{ligand_name}.in'
    ligand_list_file = osp.join(out_dir, ligand_list_name)
    docking_config = osp.join(out_dir, docking_config_name)
    return ligand_list_file, docking_config

def generate_ledock_file(receptor,center, l_list=[], l_list_outfile='', out='dock.in', n_poses=10, rmsd=1.0):
    rmsd=str(rmsd)
    box_info = getbox(center)
    x = [str(box_info[0]['minX']), str(box_info[0]['maxX'])]
    y = [str(box_info[1]['minY']), str(box_info[1]['maxY'])]
    z = [str(box_info[2]['minZ']), str(box_info[2]['maxZ'])]

    n_poses=str(n_poses)

    with open(l_list_outfile,'w') as l_out:
        for element in l_list:
            l_out.write(element)
            l_out.write('\n')

    file=[
        'Receptor\n',
        receptor + '\n\n',
        'RMSD\n',
        rmsd +'\n\n',
        'Binding pocket\n',
        x[0],' ',x[1],'\n',
        y[0],' ',y[1],'\n',
        z[0],' ',z[1],'\n\n',
        'Number of binding poses\n',
        n_poses + '\n\n',
        'Ligands list\n',
        l_list_outfile + '\n\n',
        'END']
    
    with open(out,'w') as output:
        for line in file:
            output.write(line)
    output.close()

    return out

def docking_with_ledock(dock_config, verbose=True):
    '''
    ledock_bin: str
    dock_config: str
    '''
    command = f'{ledock_bin} {dock_config}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print('LeDock Docking finished')
        else:
            print('LeDock Docking failed: {}'.format(result.stderr))
    with open(dock_config, 'r') as file:
        for line in file:
            if line.startswith('Ligands list'):
                ligands_list_path = next(file).strip()
                break  
    docked_list = []
    with open(ligands_list_path, 'r') as file:
        for line in file:
            ligand_mol2 = line.strip()
            ligand_dok = ligand_mol2.replace('.mol2', '.dok')
            if os.path.exists(ligand_dok):
                docked_list.append(ligand_dok)
                continue
    if len(docked_list) == 0:
        raise ValueError('Docking seems failed, please check your input files')
            
    return docked_list

def dok2sdf(dok_file, ori_file, out_file=None):
    if ori_file.endswith('.sdf'):
        ori_mol = read_sdf(ori_file)[0]
    elif ori_file.endswith('.mol2'):
        ori_mol = Chem.MolFromMol2File(ori_file)
    else:
        raise ValueError('ori_file should be .sdf or .mol2 file')
    
    if out_file is None:
        out_file = dok_file.replace('.dok', '_ledock.sdf')
    
    conformations = []
    current_conformation = []
    affins = []
    with open(dok_file, 'r') as file:
        for line in file:
            if line.startswith("END"):
                conformations.append(current_conformation)
                current_conformation = []
            elif line.startswith("ATOM"):
                parts = line.split()
                x, y, z = map(float, parts[5:8])
                current_conformation.append((x, y, z))
            elif line.startswith("REMARK Cluster"):
                affin = line.split()[-2]
                affins.append(affin)

    mols = []
    for i, new_coord in enumerate(conformations):
        new_mol = set_mol_position(ori_mol, np.array(new_coord))
        affin = affins[i]
        ledock_remark = f' LeDock RESULT:      {affin}      0.000      0.000'
        new_mol.SetProp("REMARK", ledock_remark)
        mols.append(new_mol)
        
    write_sdf(mols, out_file)

    return out_file

import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, default='./test/4fny_protein.pdb')
    parser.add_argument('--sdf_file', type=str, default='./test/4fny_ori.sdf')
    parser.add_argument('--ori_file', type=str, default='./test/4fny_ori.sdf')
    args = parser.parse_args()
    
    # docking 
    protein_pro = lepro(args.pdb_file)
    ligand_mol2 = sdf2mol2(args.sdf_file)

    center = sdf2centroid(args.ori_file)
    # ledock is quite wired, it need user to prepare the docking config, and in the docking config, 
    # you need another ligand list file to store the ligands path. Thus I automatically generate their
    # name, and use them in the generate_ledock_file file. 
    ligand_list_file, docking_config = ledock_config_generate(protein_pro, ligand_mol2)
    docking_config = generate_ledock_file(protein_pro, center, l_list=[ligand_mol2], l_list_outfile=ligand_list_file, out=docking_config)
    
    dok_files = docking_with_ledock(docking_config)
    # since I only dock one files, the list only contains one file
    dok_file = dok_files[0]
    result_sdf = dok2sdf(dok_file, args.sdf_file)

    # remove the redundant files
    rmfiles(*dok_files, docking_config, ligand_list_file) 