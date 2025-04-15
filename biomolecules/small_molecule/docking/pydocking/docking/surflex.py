import subprocess
import re
from .chem import read_sdf, write_sdf, set_mol_position, pdb2mol2, sdf2mol2
import os
import shutil
import os.path as osp

os.environ['TA_LICENSE'] = '/home/haotian/Molecule_Generation/Docking/classical/software/sybylx2.1.1/AdminTools11.6'
surflex_path = '/home/haotian/Molecule_Generation/Docking/classical/software/sybylx2.1.1/sybylx2.1.1/bin/linuxem64/surflex-dock.exe'

def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()

def get_protomol(protein_mol2, ligand_mol2, out_file=None, use_exitprotomol=True, verbose=True):
    '''
    Protomol is the special defination used in surflex-dock, it is the active site determined by the 
    orignial crystal ligand or by the pocket searching method (e.g. fpocket)
    '''
    if out_file is None:
        out_file = protein_mol2.replace('.mol2', '')
    if os.path.exists(out_file+'-protomol.mol2') and use_exitprotomol:
        print(f'{out_file}'+'-protomol.mol2: Exists!')
        return out_file+'-protomol.mol2'
    
    redundant_files = out_file + '-*.pdb'
    command = f'{surflex_path} proto {protein_mol2} {ligand_mol2} {out_file} && rm -rf {redundant_files}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(f'{protein_mol2}  Have been converted to protomol file!')
    else:
        print(result.stderr)
    return out_file+'-protomol.mol2'

def surflex_dock(protein_mol2, ligand_mol2, ori_mol2, log_file=None, n_confs=5, verbose=True):
    if log_file is None:
        log_file = ligand_mol2.replace('.mol2', '')
    # create a protomol file
    protomol_file = get_protomol(protein_mol2, ori_mol2, verbose=verbose)
    # create a list that stores the ligand mol2 file name 
    list_file = osp.join(osp.dirname(log_file),'list')
    write_file(outline=ligand_mol2, output_file=list_file)
    # run the surflex-dock command
    command = f'{surflex_path} -ndock_final {n_confs} +misc_premin +fastsearch -multistart 10 dock_list {list_file} {protomol_file} {protein_mol2} {log_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print(f'{ligand_mol2}  Have been docked!')
        else:
            print(result.stderr)
    result_file = log_file + '-results.mol2'
    # process the result file to the normal sdf file
    result_sdf = process_docked_mol2(result_file, log_file)
    result_tab = log_file + '-results_tab.log'
    
    os.remove(ligand_mol2)
    os.remove(result_file)
    os.remove(list_file)
    os.remove(log_file)
    os.remove(result_tab)
    return result_sdf

def surflex_score(protein_mol2, ligand_mol2, ori_mol2, log_file=None, verbose=True, save=True):
    if log_file is None:
        log_file = ligand_mol2.replace('.mol2', '')
    # create a protomol file
    protomol_file = get_protomol(protein_mol2, ori_mol2, verbose=verbose)
    # create a list that stores the ligand mol2 file name 
    list_file = osp.join(osp.dirname(log_file),'list')
    write_file(outline=ligand_mol2, output_file=list_file)
    # run the surflex-dock command
    command = f'{surflex_path} +misc_premin +misc_premin score_list {list_file} {protomol_file} {protein_mol2} {log_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print(f'{ligand_mol2}  Have been docked!')
        else:
            print(result.stderr)
    affin = open(log_file,'r').readline().split('Non-opt: ')[1].split()[0]
    result_file = log_file + '-results.mol2'

    os.remove(result_file)
    os.remove(ligand_mol2)
    os.remove(list_file)
    os.remove(log_file)
    if save:
        write_file(ligand_mol2.replace('.mol2','_surflex_score'), affin)
    return affin


def process_docked_mol2(docked_mol2file, docking_log, verbose=True):
    '''
    You will get the log file and corresponding results file after docking
    process_docked_mol2 will save the mol2 file in the sdf format, with marking its scored affinity. (all the docked results are saved as this format)
    '''
    pattern = re.compile(r'obj01_\d{3}:\s([-\d.]+)')
    affins_list = []
    with open(docking_log, 'r') as file:
        for line in file:
            matches = pattern.findall(line)
            affins_list.extend(matches)
    # print(affins_list)
    sdf_file = docked_mol2file.replace('-results.mol2', '_surflex.sdf') 
    command = f'obabel {docked_mol2file} -O {sdf_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose==True:
        if result.returncode == 0:
            print(f'{docked_mol2file}  Have been converted to sdf file!')
        else:
            print(result.stderr)

    mols = read_sdf(sdf_file)
    if len(mols) != len(affins_list):
        raise ValueError('The number of successful molecule and docking results does not match!')
    
    for affin, mol in zip(affins_list, mols):
        surflex_remark = f' SurfLex RESULT:      {affin}      0.000      0.000'
        mol.SetProp('REMARK', surflex_remark)
    write_sdf(mols, sdf_file)
    
    return sdf_file


import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, default='./test/4fny_protein.pdb')
    parser.add_argument('--sdf_file', type=str, default='./test/4fny_ori.sdf')
    parser.add_argument('--ori_file', type=str, default='./test/4fny_ori.sdf')
    parser.add_argument('--n_confs', type=int, default=5)
    args = parser.parse_args()
    
    # docking 
    protein_mol2 = pdb2mol2(args.pdb_file)
    ligand_mol2 = sdf2mol2(args.sdf_file)
    ori_mol2 = sdf2mol2(args.ori_file)
    result_sdf = surflex_dock(protein_mol2, ligand_mol2, ori_mol2, args.n_confs)

    # score 
    protein_mol2 = pdb2mol2(args.pdb_file)
    ligand_mol2 = sdf2mol2(args.sdf_file)
    affin = surflex_score(protein_mol2, ligand_mol2) # it also saved as ligand_mol2.replace('.mol2','_surflex_score')
