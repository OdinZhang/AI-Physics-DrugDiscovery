from utils.xtb_density import CDCalculator, interplot_ecloud
from utils.grid import *
from utils.chem import *
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.tools.voxeldescriptors import _getOccupancyC

calculater = CDCalculator(xtb_command='xtb')


pkt_mol = Chem.MolFromPDBFile('./pkt.pdb')
lig_mol = read_sdf('./lig.sdf')[0]

# radius is the radius of the grid sphere
radius = 10

boundary = 2 * radius 
resolution = 0.5 # boundy / size
size = int(boundary/resolution)
assert size % 1 == 0, print('Error: size must be an integer')
size += 1

N = [size, size, size]
llc = (np.zeros(3) - float(size * resolution / 2)) + resolution / 2
# Now, the box is 24×24×24 A^3
expanded_pcenters = BuildGridCenters(llc, N, resolution)

rotation = True

pkt_smallmol = SmallMol(pkt_mol)


lig_coords = lig_mol.GetConformer().GetPositions()
lig_center = lig_coords.mean(axis=0)

# define the pkt channel
pkt_sigmas, pkt_coords, pkt_center = generate_sigmas(pkt_smallmol)

# use the pkt_center as the whole center
pkt_grids = expanded_pcenters + pkt_center
lig_grids = expanded_pcenters + pkt_center


# Do the rotation
if rotation:
    rrot = uniformRandomRotation()  # Rotation
    lig_coords = rotate(lig_coords, rrot, center=pkt_center)
    pkt_coords_ = rotate(pkt_mol.GetConformer().GetPositions(), rrot, center=pkt_center)
    pkt_coords = rotate(pkt_coords, rrot, center=pkt_center)



pkt_channel = _getOccupancyC(pkt_coords.astype(np.float32),
                            pkt_grids.reshape(-1, 3),
                            pkt_sigmas).reshape(size, size, size, 8)


rotated_lig_mol = set_mol_position(lig_mol, lig_coords)
lig_ecloud = calculater.calculate(rotated_lig_mol)
lig_density = interplot_ecloud(lig_ecloud, lig_grids.transpose(3,0,1,2)).reshape(N)

rotated_pkt_mol = set_mol_position(pkt_mol, pkt_coords_)

from utils.xtb_density import write_new_cube
from utils.cubtools import write_cube
BOHR = 1.8897259886

new_X, new_Y, new_Z = lig_grids.transpose(3,0,1,2)

vec_x = new_X[1, 0, 0] - new_X[0, 0, 0]
vec_y = new_Y[0, 1, 0] - new_Y[0, 0, 0]
vec_z = new_Z[0, 0, 1] - new_Z[0, 0, 0]

lig_meta = {'org': llc*BOHR, #llc
        'xvec': np.array([vec_x, 0, 0])*BOHR,
        'yvec': np.array([0, vec_y, 0])*BOHR,
        'zvec': np.array([0, 0, vec_z])*BOHR,
        'nx': size,
        'ny': size,
        'nz': size,
        'atoms': lig_ecloud['meta']['atoms'],
}
write_cube(lig_density/(BOHR**3), lig_meta, './rotated_lig.cube')
pkt_conf = pkt_mol.GetConformer(0).GetPositions() *BOHR
pkt_atoms = [(int(atom.GetMass()), pkt_conf[atom.GetIdx()]) for atom in pkt_mol.GetAtoms()]
pkt_meta = {'org': llc*BOHR, #llc
        'xvec': np.array([vec_x, 0, 0])*BOHR,
        'yvec': np.array([0, vec_y, 0])*BOHR,
        'zvec': np.array([0, 0, vec_z])*BOHR,
        'nx': size,
        'ny': size,
        'nz': size,
        'atoms': pkt_atoms,
}
write_cube(pkt_channel[:,:,:,4],pkt_meta,'./rotated_pkt.cub')
Chem.MolToPDBFile(rotated_pkt_mol, 'rotated_pkt.pdb')
Chem.MolToPDBFile(rotated_lig_mol, 'rotated_lig.pdb')