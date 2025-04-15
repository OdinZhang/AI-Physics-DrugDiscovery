import numpy as np
import os.path as osp


# chainID
chainID_idx = 22 - 1
# resName
resName_sidx = 18 - 1
resName_eidx = 20
# resSeq
resSeq_sidx = 23 - 1
resSeq_eidx = 26
# element
ele_sidx = 77 - 1
ele_eidx = 78
# Number of heavy atoms of the real ligand
num_valid_ats = 3

def is_nmr(pdb):
    '''
    Determine if pdb is an NMR structure
    '''
    for i, line in enumerate(pdb):
        if line.startswith("MODEL"):
            break
    if i == len(pdb) - 1:
        return False
    else:
        return True

def sel_last_model(pdb):
    '''
    Select the last MODEL of the NMR structure, including the last CONECT record.
    '''
    pdb = pdb[::-1]  
    last_model = []
    for line in pdb:
        if line.startswith("MODEL"):
            break
        else:
            last_model.append(line)
    return last_model[::-1]

def sel_real_heattm(pdb):
    '''
     Select the real HETATM record, including the last CONECT record, the judgement criterion is that there is a TER record in the line before the HETATM record.
    '''
    st_num = None
    for idx, line in enumerate(pdb):
        if line.startswith("TER"):
            if pdb[idx+1].startswith("HETATM"):
                st_num = idx+1
                break
    if st_num:
        return pdb[st_num:], st_num  # Returns the real HETATM record and the subscript where it started.
    else:
        return None, None
    

def sel_conect(pdb):
    '''
    Selecting CONECT records
    '''
    st_num = None
    for idx, line in enumerate(pdb):
        if line.startswith("CONECT"):
            st_num = idx
            break
    if st_num:
        return pdb[st_num:], st_num   # Returns the CONECT record and the subscript at its beginning.
    else:
        return None, None

def get_num_heavat(records): # records is list
    count = 0
    for line in records:
        if line[ele_sidx:ele_eidx].strip() != 'H':
            count = count + 1
    return count

def is_hoh(records): # records is list
    is_hoh = False
    for line in records:
        if line[resName_sidx:resName_eidx]=='HOH':
            is_hoh = True
            break
    return is_hoh

def get_ligs(heattm):
    '''
    Get the unique identifier and its corresponding subscript for each HETATM start record
    '''
    # Get the chain ID, residue name, residue number, and subscript of the HETATM record for each record in the real HETATM record that starts with HETATM.
    chainIDs, resNames, resSeqs, corr_idxs = [], [], [], []
    for i, ht in enumerate(heattm):
        if ht.startswith('HETATM'):  # There may be non-HETATM start fields, e.g. 'ANISOU', 'CONECT', 'MASTER', 'END'
            chainIDs.append(ht[chainID_idx])
            resNames.append(ht[resName_sidx:resName_eidx])
            resSeqs.append(ht[resSeq_sidx:resSeq_eidx])
            corr_idxs.append(i)  

    resNames = [_.strip() for _ in resNames]  # Remove spaces
    resSeqs = [_.strip() for _ in resSeqs]  # Remove spaces

    # Get a unique identifier for each HETATM start record
    ligs = []
    for chainID, resName, resSeq in zip(chainIDs,resNames,resSeqs):
        ligs.append('%s_%s_%s' % (chainID, resName, resSeq))
    return ligs, corr_idxs


def get_num_ter_heatam(pdb):
    '''
    Get the number of TER records followed by HETATM records.
    If the number is greater than 1, the ligand will not be separated. For example, PDBID: 6C8N
    '''
    num = 0
    for idx, line in enumerate(pdb):
        if line.startswith("TER"):
            if pdb[idx+1].startswith("HETATM"):
                num = num + 1
    return num

def exsit_atom_record(pdb):
    '''
    Check if there is an ATOM record in the pdb, some pdb ids, such as 2GQ7, do not exist
    '''
    flag = False
    for idx, line in enumerate(pdb):
        if line.startswith("ATOM"):
            flag = True
            break
    return flag

def pdb_lig_split(pdb_file, pdb_id=None, out_dir='./'):
    """
    Splits PDB file into separate files for protein and each unique ligand.
    
    Args:
        pdb_file (str): Path to the PDB file.
        pdb_id (str, optional): Identifier for the PDB file. Defaults to the filename without extension.
        out_dir (str): Output directory to save split files.
        
    Returns:
        tuple: Path to the saved protein file and list of paths to the saved ligand files.
    """
    if pdb_id is None:
        pdb_id = pdb_file.split('/')[-1].split('.')[0]

    with open(pdb_file, 'r') as f:
        pdb = f.readlines()

    # determine if the pdb is an NMR structure
    if is_nmr(pdb):
        pdb = sel_last_model(pdb)

    heattm, st_num = sel_real_heattm(pdb)
    num_ter_heatam = get_num_ter_heatam(pdb)
    atom_record_flag = exsit_atom_record(pdb)

    if (not heattm) or (num_ter_heatam > 1) or (not atom_record_flag):
        print('skipping...%s' % pdb_id)
        return None, None
    

    conect, st_num_ = sel_conect(pdb)
    if not conect:
        print('There is no real CONECT records in %s' % pdb_id)
    ligs, corr_idxs = get_ligs(heattm)

    # Resolve residues in HETATM records
    uni_ligs =  list(np.unique(ligs))
    num_uni_ligs = len(uni_ligs)
    num_ligs = len(ligs)
    uni_ligs_idxs = []
    for lig in uni_ligs:
        _ = [corr_idxs[i] for i in range(num_ligs) if ligs[i] == lig]  # Get record subscripts for all atoms of each residue
        uni_ligs_idxs.append(_)

    # Count heavy atoms per residue
    num_heavat = [get_num_heavat(list(np.array(heattm)[uni_ligs_idxs[i]])) for i in range(num_uni_ligs)]

    flags = np.array(['ORIGIN_ORIGIN' for _ in heattm])

    # Determine 'VALID' and 'INVALID' records
    for i, num in enumerate(num_heavat):
        if num >= num_valid_ats:
            flags[uni_ligs_idxs[i]] = uni_ligs[i]  
        if num<num_valid_ats:
            flags[uni_ligs_idxs[i]] = "INVAL"
    
    # write ligand files
    lig_num = 0
    saved_ligs = []
    for ligs in uni_ligs:
        if np.any(flags == ligs):
            name = ligs
            saved_lig = osp.join(out_dir, f'{pdb_id}_{name}.pdb')
            saved_ligs.append(saved_lig)
            with open(saved_lig, 'w') as f:
                if conect:
                    f.write(''.join(list(np.array(heattm)[flags==ligs])+conect))
                    lig_num += 1
                else:
                    f.write(''.join(list(np.array(heattm)[flags==ligs])))
                    lig_num += 1
    
    # write protein file
    saved_protein = osp.join(out_dir, f'{pdb_id}_protein.pdb')
    # Write remaining 'VALID' and 'ORIGINAL' records to the protein file
    if np.any(flags=='ORIGIN_ORIGIN'):
        # _ = np.logical_or(flags == 'VAL', flags == 'ORIGIN')
        with open(saved_protein, 'w') as f:
            f.write(''.join(list(pdb[:st_num])+list(np.array(heattm)[flags == 'ORIGIN_ORIGIN'])))

    return saved_protein, saved_ligs