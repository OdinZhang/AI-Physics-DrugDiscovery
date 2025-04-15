import subprocess
import os
import shutil
from .chem import read_sdf, write_sdf, set_mol_position
from rdkit import Chem
import numpy as np

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

def lepro(pdb_file, out_dir, outname=None):
    '''
    pdb_file: str
    out_path: str
    '''
    command = f'{lepro_bin} {pdb_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr)
        return None
    cwd = os.getcwd()
    out_pdb = os.path.join(cwd, 'pro.pdb')
    pdb_name = os.path.basename(pdb_file).split('.')[0]
    if outname is None:
        save_pdb = os.path.join(out_dir, f'{pdb_name}_pro.pdb')
    else:
        if not outname.endswith('.pdb'):
            outname = outname + '.pdb'
        save_pdb = os.path.join(out_dir, outname)
    shutil.move(out_pdb, save_pdb)
    return save_pdb

def sdf2mol2(sdf_file, mol2_file=None):
    '''
    sdf_file: str
    mol2_file: str
    '''
    mol = read_sdf(sdf_file)[0]
    if mol2_file is None:
        mol2_file = sdf_file.replace('.sdf', '.mol2')
    Chem.MolToMolFile(mol, mol2_file)
    return mol2_file


def generate_ledock_file(receptor='pro.pdb',rmsd=1.0,x=[0,0],y=[0,0],z=[0,0], n_poses=10, l_list=[],l_list_outfile='',out='dock.in'):
    rmsd=str(rmsd)
    x=[str(x) for x in x]
    y=[str(y) for y in y]
    z=[str(z) for z in z]
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

def docking_with_ledock(dock_config):
    '''
    ledock_bin: str
    dock_config: str
    '''
    command = f'{ledock_bin} {dock_config}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

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

def dok2sdf(dok_file, ori_file, sdf_file):
    if ori_file.endswith('.sdf'):
        ori_mol = read_sdf(ori_file)[0]
    elif ori_file.endswith('.mol2'):
        ori_mol = Chem.MolFromMol2File(ori_file)
    else:
        raise ValueError('ori_file should be .sdf or .mol2 file')
    
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
        
    write_sdf(mols, sdf_file)