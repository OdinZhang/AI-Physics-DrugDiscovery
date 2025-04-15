import multiprocessing
from rdkit.Chem import AllChem
from rdkit import Chem
from tqdm import tqdm
import multiprocessing
import signal
from rdkit.Chem import AllChem
from rdkit import Chem
from tqdm import tqdm
def read_sdf(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    return [mol for mol in suppl if mol is not None]

def suppress_rdkit_warnings():
    from rdkit import RDLogger
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.CRITICAL)

# Call this function before running the rest of your script
suppress_rdkit_warnings()

def write_sdf(mol_list, file):
    writer = Chem.SDWriter(file)
    for i in mol_list:
        writer.write(i)
    writer.close()
    
def sanitize_filename(name):
    return name.replace('/', '_').replace('\\', '_')

# Timeout decorator
class TimeoutException(Exception):
    pass

def timeout(seconds):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutException()

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result
        return wrapper
    return decorator

def generate_and_save_conformer(mol):
    mol = Chem.AddHs(mol)
    if mol.HasProp('IDNUMBER'):
        name = sanitize_filename(mol.GetProp('IDNUMBER'))
    else:
        print("Warning: 'IDNUMBER' not found for a molecule, skipping.")
        return None  # Skip this molecule or you can set a default name
    ps = AllChem.ETKDGv2()
    id = AllChem.EmbedMolecule(mol, ps)
    if id == -1:
        print('rdkit coords could not be generated without using random coords. Using random coords now.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol, ps)
    AllChem.MMFFOptimizeMolecule(mol, confId=0)

    write_sdf([mol], f'SPECS/{name}.sdf')
    return name

def generate_and_save_conformer(mol, name=None):
    if name is None:
        name = sanitize_filename(mol.GetProp('IDNUMBER'))
    else:
        print('No name found')
        return None
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv2()
    id = AllChem.EmbedMolecule(mol, ps)
    if id == -1:
        print('Using random coords for embedding.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol, ps)
    AllChem.MMFFOptimizeMolecule(mol, confId=0)
    # os.makedirs('CHEMDIV', exist_ok=True)
    write_sdf([mol], f'CHEMDIV/{name}.sdf')
    return name


import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fast similarity search using OpenEye")
    parser.add_argument('--mols_sdf', type=str)
    args = parser.parse_args()
    
    print('read mols....')
    mols = read_sdf(args.mols_sdf)

    print('Reading Done')
    for mol in mols:
        mol.RemoveAllConformers()
    total = len(mols)
    for i, mol in enumerate(mols):
        ratio = i/total
        if i % 1000 == 0:
            print(f'{i} has done, ratio: {ratio}')
        generate_and_save_conformer(mol)
        
