import subprocess
import os.path as osp
import os
import shutil
from .chem import read_sdf, write_sdf
from rdkit import Chem
import numpy as np
import re

smina_bin = '/home/haotian/Molecule_Generation/Docking/classical/software/smina.static'

def prepare_target(work_dir, protein_file_name, verbose=1):
    '''
    work_dir is the dir which .pdb locates
    protein_file_name: .pdb file which contains the protein data
    '''
    protein_file = osp.join(work_dir, protein_file_name)
    command = 'prepare_receptor -r {protein} -o {protein_pdbqt}'.format(protein=protein_file,
                                                            protein_pdbqt = protein_file+'qt')
    if osp.exists(protein_file+'qt'):
        return protein_file+'qt'
        
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if verbose:
        if result.returncode == 0:
            print('prepare the target successfully')
        else:
            print(result.stderr)

    return protein_file+'qt'

def prepare_ligand(work_dir, lig_sdf, verbose=1):
    lig_name = lig_sdf
    lig_mol2 = lig_sdf[:-3]+'mol2'
    now_cwd = os.getcwd()
    lig_sdf = osp.join(work_dir, lig_sdf)
    cwd_mol2 = osp.join(now_cwd, lig_mol2)
    work_mol2 = osp.join(work_dir, lig_mol2)
    command = '''obabel {lig} -O {lig_mol2}'''.format(lig=lig_sdf,
                                                        lig_mol2 = work_mol2)
    obabel_result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if verbose:
        if obabel_result.returncode == 0:
            if verbose:
                print('obabel successfully')
        else:
            print('command: ', command)
            print(obabel_result.stderr)
        
    shutil.copy(work_mol2, now_cwd)

    command = '''prepare_ligand -l {lig_mol2} -A hydrogens'''.format(lig_mol2=cwd_mol2)
    prep_lig_result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if prep_lig_result.returncode == 0:
            print('prepare_ligand successfully')
        else:
            print('command: ', command)
            print(prep_lig_result.stderr)

    lig_pdbqt = lig_name[:-3]+'pdbqt'
    cwd_pdbqt = osp.join(now_cwd, lig_pdbqt)
    work_pdbqt = osp.join(work_dir, lig_pdbqt)
    os.remove(cwd_mol2)
    os.remove(work_mol2)
    if osp.exists(work_pdbqt):
        os.remove(work_pdbqt)
    shutil.move(cwd_pdbqt, work_dir)
    if os.path.exists(lig_pdbqt):
        if verbose:
            print('prepare successfully !')
        else:
            print('generation failed!')
    return work_pdbqt

def prepare_ligand_obabel(work_dir, ligand, out_lig_name=None, mode='ph'):
    ligand = osp.join(work_dir, ligand)

    if not osp.exists(ligand):
        raise ValueError('ligand file does not exist')
    
    lig_pdbqt = ligand.replace('sdf','pdbqt')
    if mode == 'ph':
        command = 'obabel {lig_sdf} -O {lig_pdbqt} -p 7.4'.format(lig_sdf=ligand, lig_pdbqt=lig_pdbqt)
        #command = 'prepare_ligand -l {lig_sdf} -o {lig_pdbqt} -A hydrogens'.format(lig_sdf=lig_sdf, lig_pdbqt=lig_pdbqt)
    elif mode == 'noh':
        command = 'obabel {lig_sdf} -O {lig_pdbqt}'.format(lig_sdf=ligand, lig_pdbqt=lig_pdbqt)
    elif mode == 'allH':
        mol = Chem.SDMolSupplier(ligand)[0]
        ligand_allH = ligand[:-4]+'_allH.sdf'
        write_sdf(mol, ligand_allH)
        ligand = ligand_allH
        command = 'obabel {lig_sdf} -O {lig_pdbqt}'.format(lig_sdf=ligand, lig_pdbqt=lig_pdbqt)
    
    pre_lig_result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if pre_lig_result.returncode == 0:
        pass
    else:
        print('command: ', command)
        print(pre_lig_result.stderr)

    return out_lig_name

def docking_with_smina(work_dir, protein_pdbqt, lig_pdbqt, centroid, verbose=1, out_lig_sdf=None, save_pdbqt=False):
    '''
    work_dir: is same as the prepare_target
    protein_pdbqt: .pdbqt file
    lig_sdf: ligand .sdf format file
    '''
    # prepare target
    lig_pdbqt = osp.join(work_dir, lig_pdbqt)
    protein_pdbqt = osp.join(work_dir, protein_pdbqt)
    cx, cy, cz = centroid
    out_lig_sdf_dirname = osp.dirname(lig_pdbqt)
    out_lig_pdbqt_filename = osp.basename(lig_pdbqt).split('.')[0]+'_out.pdbqt'
    out_lig_pdbqt = osp.join(out_lig_sdf_dirname, out_lig_pdbqt_filename) 
    if out_lig_sdf is None:
        out_lig_sdf_filename = osp.basename(lig_pdbqt).split('.')[0]+'_out.sdf'
        out_lig_sdf = osp.join(out_lig_sdf_dirname, out_lig_sdf_filename) 
    else:
        out_lig_sdf = osp.join(work_dir, out_lig_sdf)

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
    
    if dock_result.returncode == 0:
        if verbose:
            print('docking successfully')
        pass
    else:
        print('command: ', command)
        print(dock_result.stderr)
        raise ValueError('docking failed')


    os.remove(lig_pdbqt)
    if not save_pdbqt:
        os.remove(out_lig_pdbqt)
    
    return out_lig_sdf

def docking_with_sminaflex(work_dir, protein_pdbqt, lig_pdbqt, centroid, verbose=1, out_lig_sdf=None, save_pdbqt=False):
    '''
    work_dir: you need to put the protein_pdbqt and lig_sdf in the same directory, which is work_dir
    protein_pdbqt: .pdbqt file name
    lig_sdf: ligand .sdf file name
    '''
    # prepare target
    lig_pdbqt = osp.join(work_dir, lig_pdbqt)
    protein_pdbqt = osp.join(work_dir, protein_pdbqt)
    cx, cy, cz = centroid
    out_lig_sdf_dirname = osp.dirname(lig_pdbqt)
    out_lig_pdbqt_filename = osp.basename(lig_pdbqt).split('.')[0]+'_out.pdbqt'
    out_lig_pdbqt = osp.join(out_lig_sdf_dirname, out_lig_pdbqt_filename) 
    if out_lig_sdf is None:
        out_lig_sdf_filename = osp.basename(lig_pdbqt).split('.')[0]+'_out.sdf'
        out_lig_sdf = osp.join(out_lig_sdf_dirname, out_lig_sdf_filename) 
    else:
        out_lig_sdf = osp.join(work_dir, out_lig_sdf)

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
    
    if dock_result.returncode == 0:
        if verbose:
            print('docking successfully')
        pass
    else:
        print('command: ', command)
        print(dock_result.stderr)
        raise ValueError('docking failed')


    os.remove(lig_pdbqt)
    if not save_pdbqt:
        os.remove(out_lig_pdbqt)
    
    return out_lig_sdf