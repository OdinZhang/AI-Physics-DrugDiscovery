import requests
import os
from six.moves.urllib.request import urlopen
import json
from .pdb_utils import align_and_rmsd


def pdb2uniprot(pdb):
    '''
    Transform the pdbid to the uniprot id
    '''
    pdb = pdb.lower()
    try:
        content = urlopen('https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/' + pdb).read()
    except:
        print(pdb,"PDB Not Found (HTTP Error 404). Skipped.")
        return None
    content = json.loads(content.decode('utf-8'))
    uniprotid = list(content[pdb]['UniProt'].keys())[0]
    return uniprotid

def download_alphafold_prediction(uniprot_id, output_dir, out_name=None):
    '''
    Use the uniprot id to download AlphaFold-generated conformations
    '''
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id.upper()}-F1-model_v4.pdb"
    response = requests.get(url)
    if response.status_code == 404:
        print('Network Error')
    if response.status_code == 200:
        if out_name:
            output_path = os.path.join(output_dir, out_name)
        else:
            output_path = os.path.join(output_dir, f"AF-{uniprot_id.upper()}-F1-model_v4.pdb")
        with open(output_path, "w") as output_file:
            output_file.write(response.text)
        print(f"Downloaded AlphaFold prediction for {uniprot_id} to {output_path}")
    else:
        print(f"AlphaFold prediction not found for {uniprot_id}")



if __name__ == '__main': 

    import torch
    from tqdm import tqdm
    import os.path as osp
    import shutil
    cross_base = '/home/haotian/Molecule_Generation/Graph-BP/data/crossdock2020'
    new_base = '/home/haotian/Molecule_Generation/Res2Mol/data/alphafold/alphafold_pro'

    split = '../split_by_name.pt'
    split = torch.load(split)
    test_split = split['test']

    for pair in tqdm(test_split):
        try:
            pdb_file = pair[0]
            data_dir = pdb_file.split('/')[0]
            if not osp.exists(osp.join(new_base,data_dir)):
                os.makedirs(osp.join(new_base,data_dir), exist_ok=True)
                sdf_file = pair[1]
                shutil.copy(osp.join(cross_base,sdf_file), osp.join(new_base,data_dir))
                pdbid = pdb_file.split('/')[-1][:4]
                uniprotid = pdb2uniprot(pdbid)
                download_alphafold_prediction(uniprotid, osp.join(new_base,data_dir), out_name=pdbid+'_af.pdb')
                af_pdb = osp.join(new_base,data_dir, pdbid+'_af.pdb')

                ori_pdbpath = osp.join(cross_base,data_dir,pdb_file.split('/')[-1][:10]+'.pdb')
                align_and_rmsd(af_pdb, ori_pdbpath)
        except:
            print('failed',pdb_file)
    