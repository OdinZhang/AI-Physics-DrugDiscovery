import os.path as osp
import os
from rdkit import Chem
import numpy as np
from easydict import EasyDict
from rdkit import Chem
import subprocess
from rdkit.Chem.rdMolAlign import CalcRMS
import shutil
import re
import copy

def read_sdf(sdf_file, sanitize=False):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=sanitize)
    mols_list = [i for i in supp]
    return mols_list

def write_sdf(mol_list,file, voice=False):
    writer = Chem.SDWriter(file)
    mol_cnt = 0
    for i in mol_list:
        try:
            writer.write(i)
            mol_cnt+=1
        except:
            pass
    writer.close()
    if voice: 
        print('Write {} molecules to {}'.format(mol_cnt,file))

def set_mol_position(mol, pos):
    mol = copy.deepcopy(mol)
    for i in range(pos.shape[0]):
        mol.GetConformer(0).SetAtomPosition(i, pos[i].tolist())
    return mol 

def sdf2centroid(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=False)
    lig_xyz = supp[0].GetConformer().GetPositions()
    centroid_x = lig_xyz[:,0].mean()
    centroid_y = lig_xyz[:,1].mean()
    centroid_z = lig_xyz[:,2].mean()
    return centroid_x, centroid_y, centroid_z

def mol2centroid(mol2_file):
    mol = Chem.MolFromMol2File(mol2_file, sanitize=False)
    lig_xyz = mol.GetConformer().GetPositions()
    centroid_x, centroid_y, centroid_z = lig_xyz.mean(axis=0)
    return centroid_x, centroid_y, centroid_z

def get_result(docked_sdf, ref_mol=None):
    suppl = Chem.SDMolSupplier(docked_sdf,sanitize=False)
    results = []
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        line = mol.GetProp('REMARK').splitlines()[0].split()[2:]
        try:
            rmsd = CalcRMS(ref_mol, mol)
        except:
            rmsd = np.nan
        results.append(EasyDict({
            'rdmol': mol,
            'mode_id': i,
            'affinity': float(line[0]),
            'rmsd_lb': float(line[1]),
            'rmsd_ub': float(line[2]),
            'rmsd_ref': rmsd
        }))
    return results

def sdf2mol2(sdf_file, mol2_file=None, verbose=True):
    '''
    sdf_file: str
    mol2_file: str
    '''
    if mol2_file is None:
        mol2_file = sdf_file.replace('.sdf', '.mol2')
    if os.path.exists(mol2_file):
        return mol2_file
    
    command = f'obabel {sdf_file} -O {mol2_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print(f'Have been converted to mol2 file! {sdf_file}')
        else:
            print(result.stderr)
    return mol2_file

def pdb2mol2(pdb_file, out_file=None):
    '''
    SurfLex needs protein file in mol2 format
    '''
    if out_file is None:
        out_file = pdb_file.replace('.pdb', '.mol2')
    if os.path.exists(out_file):
        # print(f'{out_file}: Exists!')
        return out_file
    
    command = f'obabel {pdb_file} -O {out_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(f'{pdb_file}  Have been converted to mol2 file!')
    else:
        print(result.stderr)
    return out_file

def rmfiles(*files, verbose=True):
    for file in files:
        try:
            os.remove(file)
            if verbose:
                print(f"File {file} has been successfully removed.")
        except FileNotFoundError:
            print(f"File {file} not found.")
        except Exception as e:
            print(f"An error occurred while deleting {file}: {e}")

def compute_rmsd_obrms(true_sdf, docked_sdf):
    command = f'obrms {docked_sdf} {true_sdf}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    outputs = result.stdout.split('\n')
    rmsd_list = []
    for output in outputs[:-1]:
        rmsd_list.append(float(output.split()[-1]))
    return rmsd_list