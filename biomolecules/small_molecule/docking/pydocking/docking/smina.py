import subprocess
import os.path as osp
import os
import shutil
from .chem import read_sdf, write_sdf, sdf2centroid, sdf2mol2
from rdkit import Chem
import numpy as np
import re

smina_bin = '/home/haotian/molecules_confs/Protein_test/classsical/software/smina.static'

def process_sdf(sdf_file):
    mols = read_sdf(sdf_file)
    for mol in mols:
        affin = mol.GetProp('REMARK').splitlines()[0].split()[1]
        smina_remark = f'SMINA RESULT:      {affin}      0.000      0.000'
        mol.SetProp('REMARK', smina_remark)
    write_sdf(mols, sdf_file)

def prepare_target(protein_file, out_file=None, verbose=1):
    if out_file is None:
        out_file = protein_file + 'qt'
    if osp.exists(out_file):
        return out_file
    
    command = f'prepare_receptor -r {protein_file} -o {out_file}'
    if osp.exists(protein_file+'qt'):
        return protein_file+'qt'
        
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if verbose:
        if result.returncode == 0:
            print('prepare the target successfully:', out_file)
        else:
            print('prepare the target failed:', result.stderr)

    return out_file

def prepare_ligand(ligand_file, out_file=None, verbose=1):
    ligand_name = osp.basename(ligand_file)
    if ligand_name.endswith('.mol2'):
        ligand_mol2_file = ligand_file
    elif ligand_name.endswith('.sdf'):
        ligand_mol2_file = sdf2mol2(ligand_file, verbose=verbose)
    else:
        raise ValueError('Unsupported ligand file format: {}'.format(ligand_file))
    
    if out_file is None:
        out_file = ligand_mol2_file.replace('.mol2', '.pdbqt')
    if osp.exists(out_file):
        return out_file
    
    # prepare_target4.py cannot identify the absolute or relative path, so we need to cd to the ligand_mol2_dir and perform the prepare_liagnd
    ligand_dir = osp.dirname(ligand_mol2_file)
    ligand_mol2_name = osp.basename(ligand_mol2_file)
    out_file_name = osp.basename(out_file)

    command = f'cd {ligand_dir} && prepare_ligand -l {ligand_mol2_name} -A hydrogens -o {out_file_name}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print('prepare the ligand successfully:', out_file)
        else:
            print('prepare the ligand failed:', result.stderr)
    pdbqt_file = osp.join(ligand_dir, ligand_mol2_name)
    command = f'mv {pdbqt_file} {out_file}'
    return out_file


def docking_with_smina(protein_pdbqt, lig_pdbqt, centroid, verbose=1, out_lig_sdf=None, save_pdbqt=False):
    cx, cy, cz = centroid
    out_lig_pdbqt = lig_pdbqt.replace('.pdbqt', '_smina.pdbqt')
    if out_lig_sdf is None:
        out_lig_sdf = lig_pdbqt.replace('.pdbqt', '_smina.sdf')
    if osp.exists(out_lig_sdf):
        return out_lig_sdf
    if not osp.exists(protein_pdbqt) or not osp.exists(lig_pdbqt):
        raise NotImplementedError('no pdbqt file found')
    
    command = '''{smina_bin} \
        --receptor {receptor_pre} \
        --ligand {ligand_pre} \
        --center_x {centroid_x:.4f} \
        --center_y {centroid_y:.4f} \
        --center_z {centroid_z:.4f} \
        --size_x 20 --size_y 20 --size_z 20 \
        --out {out_lig_pdbqt} \
        --exhaustiveness {exhaust}
        obabel {out_lig_pdbqt} -O {out_lig_sdf} -h'''.format(smina_bin=smina_bin,
                                            receptor_pre = protein_pdbqt,
                                            ligand_pre = lig_pdbqt,
                                            centroid_x = cx,
                                            centroid_y = cy,
                                            centroid_z = cz,
                                            out_lig_pdbqt = out_lig_pdbqt,
                                            exhaust = 8,
                                            out_lig_sdf = out_lig_sdf)
    
    dock_result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if dock_result.returncode == 0:
            print('docking successfully:', out_lig_sdf)
        else:
            print('docking failed:', dock_result.stderr)
    
    if not save_pdbqt:
        os.remove(out_lig_pdbqt)
    process_sdf(out_lig_sdf)
    return out_lig_sdf

def docking_with_sminaflex(protein_pdbqt, lig_pdbqt, centroid, verbose=1, out_lig_sdf=None, save_pdbqt=False):
    cx, cy, cz = centroid
    out_lig_pdbqt = lig_pdbqt.replace('.pdbqt', '_smina.pdbqt')
    if out_lig_sdf is None:
        out_lig_sdf = lig_pdbqt.replace('.pdbqt', '_smina.sdf')
    if osp.exists(out_lig_sdf):
        return out_lig_sdf
    
    command = '''{smina_bin} \
        --receptor {receptor_pre} \
        --flexdist_ligand {ligand_pre} \
        --flexdist {flexdist} \
        --ligand {ligand_pre} \
        --center_x {centroid_x:.4f} \
        --center_y {centroid_y:.4f} \
        --center_z {centroid_z:.4f} \
        --autobox_add 8.0 \
        --size_x 20 --size_y 20 --size_z 20 \
        --out {out_lig_pdbqt} \
        --exhaustiveness {exhaust}
        obabel {out_lig_pdbqt} -O {out_lig_sdf} -h'''.format(smina_bin=smina_bin,
                                            receptor_pre = protein_pdbqt,
                                            flexdist = 6.0,
                                            ligand_pre = lig_pdbqt,
                                            centroid_x = cx,
                                            centroid_y = cy,
                                            centroid_z = cz,
                                            out_lig_pdbqt = out_lig_pdbqt,
                                            exhaust = 8,
                                            out_lig_sdf = out_lig_sdf)
    dock_result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if dock_result.returncode == 0:
            print('docking successfully:', out_lig_sdf)
        else:
            print('docking failed:', dock_result.stderr)
    
    if not save_pdbqt:
        os.remove(out_lig_pdbqt)
    process_sdf(out_lig_sdf)
    return out_lig_sdf

import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, default='./test/4fny_protein.pdb')
    parser.add_argument('--sdf_file', type=str, default='./test/4fny_ori.sdf')
    parser.add_argument('--ori_file', type=str, default='./test/4fny_ori.sdf')
    args = parser.parse_args()
    
    # docking 
    protein_pdbqt = prepare_target(args.pdb_file)
    ligand_pdbqt = prepare_ligand(args.sdf_file)
    centroid = sdf2centroid(args.ori_file)
    result_sdf = docking_with_smina(protein_pdbqt, ligand_pdbqt, centroid)

    # flexible docking 
    protein_pdbqt = prepare_target(args.pdb_file)
    ligand_pdbqt = prepare_ligand(args.sdf_file)
    centroid = sdf2centroid(args.ori_file)
    result_sdf = docking_with_sminaflex(protein_pdbqt, ligand_pdbqt, centroid)