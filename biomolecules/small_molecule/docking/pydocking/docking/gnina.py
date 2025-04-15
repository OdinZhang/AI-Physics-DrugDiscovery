from .chem import read_sdf, write_sdf, set_mol_position, sdf2centroid
import os
import os.path as osp
import subprocess

gnina_bin = '/home/haotian/Molecule_Generation/Docking/classical/software/gnina'

def process_gnina_sdf(sdf_file):
    mols = read_sdf(sdf_file, sanitize=False)
    for mol in mols:
        affin = mol.GetProp('CNNscore')
        gnina_remark = f' GNINA RESULT:      {affin}      0.000      0.000'
        mol.SetProp('REMARK', gnina_remark)
    write_sdf(mols, sdf_file)


def gnina_dock(protein_pdb, ligand_sdf, centroid, out_file=None, verbose=True):
    if out_file is None:
        out_file = ligand_sdf.replace('.sdf', '_gnina.sdf')
    cx, cy, cz = centroid
    # check whether the docked sdf exists
    if osp.exists(out_file):
        return out_file

    command = '''{gnina_bin}\
            -r {protein_pdb}\
            -l {ligand_sdf}\
            --center_x {center_x: .4f}\
            --center_y {center_y: .4f}\
            --center_z {center_z: .4f}\
            --size_x 20 --size_y 20 --size_z 20\
            -o {out_file}'''.format(gnina_bin=gnina_bin,
                                    protein_pdb=protein_pdb,
                                    ligand_sdf=ligand_sdf,
                                    center_x=cx,
                                    center_y=cy,
                                    center_z=cz,
                                    out_file=out_file)
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print('Gnina Docking finished')
        else:
            print('Gnina Docking failed: {}'.format(result.stderr))
    process_gnina_sdf(out_file)
    return out_file

def gnina_flexdock(protein_pdb, ligand_sdf, centroid, out_file=None, verbose=True):
    if out_file is None:
        out_file = ligand_sdf.replace('.sdf', '_flexgnina.sdf')
    cx, cy, cz = centroid
    # check whether the docked sdf exists
    if osp.exists(out_file):
        return out_file
    
    command = '''{gnina_bin}\
            -r {protein_pdb}\
            -l {ligand_sdf}\
            --center_x {center_x: .4f}\
            --center_y {center_y: .4f}\
            --center_z {center_z: .4f}\
            --size_x 20 --size_y 20 --size_z 20\
            --flexdist_ligand {ligand_sdf}\
            --flexdist 4.0\
            -o {out_file}'''.format(gnina_bin=gnina_bin,
                                    protein_pdb=protein_pdb,
                                    ligand_sdf=ligand_sdf,
                                    center_x=cx,
                                    center_y=cy,
                                    center_z=cz,
                                    out_file=out_file)
    
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        if result.returncode == 0:
            print('Gnina Flex Docking finished: ', out_file)
        else:
            print('Gnina Flex Docking failed: {}'.format(result.stderr))
    process_gnina_sdf(out_file)
    return out_file

import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, default='./test/4fny_protein.pdb')
    parser.add_argument('--sdf_file', type=str, default='./test/4fny_ori.sdf')
    parser.add_argument('--ori_lig', type=str, default='./test/4fny_ori.sdf')
    args = parser.parse_args()
    
    centeroid = sdf2centroid(args.ori_lig)
    # docking 
    sdf_file = gnina_dock(args.pdb_file, args.sdf_file, centeroid) # 0.5min

    # flex docking 
    sdf_file = gnina_flexdock(args.pdb_file, args.sdf_file, centeroid) # 2.5min


# def gnina_dock(protein_pdb, ligand_sdf, out_file=None, verbose=True):
#     if out_file is None:
#         pdbqt_out = ligand_sdf.replace('.sdf', '_gnina.pdbqt')
#     else:
#         pdbqt_out = out_file.replace('.sdf', '.pdbqt')
#     sdf_out = pdbqt_out.replace('.pdbqt', '.sdf')
        
#     # check whether the docked sdf exists
#     if osp.exists(sdf_out):
#         return sdf_out

#     command = f'{gnina_bin} -r {protein_pdb} -l {ligand_sdf} --autobox_ligand {ligand_sdf} -o {pdbqt_out}'
#     print(command)
#     result = subprocess.run(command, shell=True, capture_output=True, text=True)
#     if verbose:
#         if result.returncode == 0:
#             command = f'obabel {pdbqt_out} -O {sdf_out}'
#             result = subprocess.run(command, shell=True, capture_output=True, text=True)
#             print('Gnina Docking finished')
#         else:
#             print('Gnina Docking failed: {}'.format(result.stderr))
#     return result

# def gnina_flexdock(protein_pdb, ligand_sdf, out_file=None, verbose=True):
#     if out_file is None:
#         pdbqt_out = ligand_sdf.replace('.sdf', '_gnina.pdbqt')
#     else:
#         pdbqt_out = out_file.replace('.sdf', '.pdbqt')
#     sdf_out = pdbqt_out.replace('.pdbqt', '.sdf')
        
#     # check whether the docked sdf exists
#     if osp.exists(sdf_out):
#         return sdf_out
    
#     command = f'{gnina_bin} -r {protein_pdb} -l {ligand_sdf} --autobox_ligand {ligand_sdf} --flexdist_ligand {ligand_sdf} --flexdist {4.0} -o {pdbqt_out}'
#     result = subprocess.run(command, shell=True, capture_output=True, text=True)
#     if verbose:
#         if result.returncode == 0:
#             command = f'obabel {pdbqt_out} -O {sdf_out}'
#             result = subprocess.run(command, shell=True, capture_output=True, text=True)
#             print('Gnina Flex Docking finished')
#         else:
#             print('Gnina Flex Docking failed: {}'.format(result.stderr))
    
#     return sdf_out