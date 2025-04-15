from rdkit import Chem
import numpy as np


# cub file default unit is BOHR, we need to transform the BOHR to Angstrom
BOHR = 1.8897259886

def _putline(*args):
    """
    Generate a line to be written to a cube file where
    the first field is an int and the remaining fields are floats.

    params:
        *args: first arg is formatted as int and remaining as floats

    returns: formatted string to be written to file with trailing newline
    """
    s = "{0:^ 8d}".format(args[0])
    s += "".join("{0:< 12.6f}".format(arg) for arg in args[1:])
    return s + "\n"

def write_cube(fname, meta, data):
    with open(fname, "w") as cube:
        # first two lines are comments
        cube.write(" Cubefile created by cubetools.py\n  source: none\n")
        natm = len(meta['atoms'])
        nx, ny, nz = data.shape
        cube.write(_putline(natm, *meta['org'])) # 3rd line #atoms and origin
        cube.write(_putline(nx, *meta['xvec']))
        cube.write(_putline(ny, *meta['yvec']))
        cube.write(_putline(nz, *meta['zvec']))
        for atom_mass, atom_pos in meta['atoms']:
            cube.write(_putline(atom_mass, *atom_pos)) #skip the newline
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if (i or j or k) and k%6==0:
                        cube.write("\n")
                    cube.write(" {0: .5E}".format(data[i,j,k]))

def write_dualcube(pkt_mol, lig_mol, pkt_data, lig_data):

    lig_conf = lig_mol.GetConformer(0).GetPositions() *BOHR 
    lig_atoms = [(int(atom.GetMass()), lig_conf[atom.GetIdx()]) for atom in lig_mol.GetAtoms()]
    
    pkt_conf = pkt_mol.GetConformer(0).GetPositions() *BOHR
    pkt_atoms = [(int(atom.GetMass()), pkt_conf[atom.GetIdx()]) for atom in pkt_mol.GetAtoms()]
    pkt_center = pkt_conf.mean(axis=0)

    # in the grid computation method, we use the pkt_center as the system center  
    pkt_meta = {'org':np.array([-12,-12,-12])*BOHR + pkt_center, \
        'xvec':[1.0*BOHR, 0.0, 0.0],'yvec':[0.0, 1.0*BOHR, 0.0],'zvec':[0.0, 0.0, 1.0*BOHR],\
            'nx':np.array(24),'ny':np.array(24),'nz':np.array(24),\
                'atoms':pkt_atoms}

    lig_meta = {'org':np.array([-12,-12,-12])*BOHR + pkt_center, \
        'xvec':[1.0*BOHR, 0.0, 0.0],'yvec':[0.0, 1.0*BOHR, 0.0],'zvec':[0.0, 0.0, 1.0*BOHR],\
            'nx':np.array(24),'ny':np.array(24),'nz':np.array(24),\
                'atoms':lig_atoms}
    
    write_cube('./pkt.cub',pkt_meta,pkt_data)
    write_cube('./lig.cub',lig_meta,lig_data)

if __name__ == '__main__':
    pkt_data = np.load('./pkt_occup.npy')
    lig_data = np.load('./lig_occup.npy')
    lig_mol = Chem.MolFromMolFile('./lig.sdf')
    pkt_mol = Chem.MolFromMolFile('./pkt.sdf')

    # write cub file
    write_dualcube(pkt_mol, lig_mol, pkt_data, lig_data)

    # save the PDB file for vmd visualization, because the SDF format could not be readed directly. 
    Chem.MolToPDBFile(pkt_mol, './pkt.pdb')
    Chem.MolToPDBFile(lig_mol, './lig.pdb')


