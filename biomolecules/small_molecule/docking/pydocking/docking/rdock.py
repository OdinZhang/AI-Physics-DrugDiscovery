# mamba install -c bioconda rdock
from .chem import read_sdf, write_sdf, set_mol_position
import tempfile
import shutil
import os.path as osp
import subprocess
import os 
import timeout_decorator

def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()

def process_rdock(rdock_sdf):
    new_mols = []
    mols = read_sdf(rdock_sdf)
    for mol in mols:
        affin = mol.GetProp('SCORE')
        rdock_remark = f' rDock RESULT:      {affin}      0.000      0.000'
        mol.SetProp('REMARK', rdock_remark)
        new_mols.append(mol)
    write_sdf(new_mols, rdock_sdf)

def rdock_virtualscreen(protein_mol2, dock_ligand_list, ref_ligand, out_file=None, n_conf=10, verbose=True, rm_tmp=True):
    '''
    rdock is not intelligent enough to handle the path, so we need to copy the files to a temporary directory, and then copy the results back
    '''
    run = False
    for dock_ligand in dock_ligand_list:
        docked_basefile = osp.basename(dock_ligand.replace('.sdf', '_rdock'))
        # if the out_file is None, then the docked file will be saved in the same directory as the dock_ligand
        if out_file is None:
            if osp.isabs(dock_ligand):
                output_path = osp.dirname(dock_ligand)
            else:
                output_path = osp.join(os.getcwd(), osp.dirname(dock_ligand))
            docked_target = osp.join(output_path, docked_basefile+'.sdf')
        else:
            docked_target = out_file
        
        if not osp.exists(docked_target):
            run = True
            break
    if not run:
        return 'All the docked files exist'

    tmpdirname = tempfile.mkdtemp()
    shutil.copy(protein_mol2, tmpdirname)
    shutil.copy(ref_ligand, tmpdirname)
    for dock_ligand in dock_ligand_list:
        shutil.copy(dock_ligand, tmpdirname)
    
    cavity_params = f'''  
    RBT_PARAMETER_FILE_V1.00
    TITLE {osp.basename(protein_mol2)}_dock

    RECEPTOR_FILE {osp.basename(protein_mol2)}
    ### RECEPTOR_FLEX 3.0

    ##################################################################
    ### CAVITY DEFINITION: REFERENCE LIGAND METHOD
    ##################################################################
    SECTION MAPPER
        SITE_MAPPER RbtLigandSiteMapper
        REF_MOL {osp.basename(ref_ligand)}
            RADIUS 6.0
            SMALL_SPHERE 1.5
    ##      LARGE_SPHERE 4.0
            MAX_CAVITIES 1
            MIN_VOLUME 100
            VOL_INCR 0.0
            GRIDSTEP 0.5
    END_SECTION

    #################################
    #CAVITY RESTRAINT PENALTY
    #################################
    SECTION CAVITY
        SCORING_FUNCTION RbtCavityGridSF
        WEIGHT 1.0
    END_SECTION
    '''
    # remove the first line and remove the space at the beginning of each line 
    processed_params = '\n'.join([line[4:] for line in cavity_params.split('\n')][1:]) 
    prm_file = osp.join(tmpdirname, 'cavity.prm')
    write_file(prm_file, processed_params)
    
    command = f'cd {tmpdirname} && rbcavity -was -d -r {prm_file}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print('successfully prepare the cavity')
        else:
            print('failed to prepare the cavity: ', result.stderr)
    return_docked_list = []
    for dock_ligand in dock_ligand_list:
        docked_basefile = osp.basename(dock_ligand.replace('.sdf', '_rdock'))

        # if the out_file is None, then the docked file will be saved in the same directory as the dock_ligand
        if out_file is None:
            if osp.isabs(dock_ligand):
                output_path = osp.dirname(dock_ligand)
            else:
                output_path = osp.join(os.getcwd(), osp.dirname(dock_ligand))
            docked_target = osp.join(output_path, docked_basefile+'.sdf')
        else:
            docked_target = out_file
        
        if osp.exists(docked_target):
            return_docked_list.append(docked_target)
            continue

        dock_ligand_basefile = osp.basename(dock_ligand)

        command = f'cd {tmpdirname} && rbdock -i {dock_ligand_basefile} -o {docked_basefile} -r {prm_file} -p dock.prm -n {n_conf} && mv {docked_basefile}.sd {docked_target}'
        # print(command)
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if verbose:
            if result.returncode == 0:
                print('successfully docked', docked_target)
            else:
                print('failed to dock: ', result.stderr)
        process_rdock(docked_target)
        return_docked_list.append(docked_target)
    shutil.rmtree(tmpdirname)
    
    return return_docked_list

import argparse
if __name__ == '__main__':
    from .chem import pdb2mol
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, default='./test/4fny_protein.pdb')
    parser.add_argument('--sdf_file', type=str, default='./test/4fny_ori.sdf')
    parser.add_argument('--ori_lig', type=str, default='./test/4fny_ori.sdf')
    args = parser.parse_args()
    
    protein_mol = pdb2mol(args.pdb_file)
    docked_file = rdock_virtualscreen(args.pdb_file, [args.sdf_file], args.ori_lig, n_conf=10, verbose=True, rm_tmp=True)
    