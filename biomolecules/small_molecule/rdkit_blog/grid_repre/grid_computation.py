import numpy as np
from moleculekit.tools.voxeldescriptors import _getOccupancyC
from moleculekit.util import uniformRandomRotation
from htmd2moleculekit import _getChannelRadii

def _getGridCenters(llc, N, resolution):
    xrange = [llc[0] + resolution * x for x in range(0, N[0])]
    yrange = [llc[1] + resolution * x for x in range(0, N[1])]
    zrange = [llc[2] + resolution * x for x in range(0, N[2])]
    centers = np.zeros((N[0], N[1], N[2], 3))
    for i, x in enumerate(xrange):
        for j, y in enumerate(yrange):
            for k, z in enumerate(zrange):
                centers[i, j, k, :] = np.array([x, y, z])
    return centers

resolution = 1.
size = 24
N = [size, size, size]
bbm = (np.zeros(3) - float(size * 1. / 2))
global_centers = _getGridCenters(bbm, N, resolution)

def rotate(coords, rotMat, center=(0,0,0)):
    """
    Rotate a selection of atoms by a given rotation around a center
    """

    newcoords = coords - center
    return np.dot(newcoords, np.transpose(rotMat)) + center


def get_aromatic_groups(in_mol):
    """
    Obtain groups of aromatic rings
    """
    groups = []
    ring_atoms = in_mol.GetRingInfo().AtomRings()
    for ring_group in ring_atoms:
        if all([in_mol.GetAtomWithIdx(x).GetIsAromatic() for x in ring_group]):
            groups.append(ring_group)
    return groups


def generate_sigmas(mol):
    """
    Calculates sigmas for elements as well as pharmacophores.
    Returns sigmas, coordinates and center of ligand.
    """
    coords = mol._coords[: , : , 0]
    n_atoms = len(coords)
    lig_center = mol.getCenter()

    # Calculate all the channels
    multisigmas = _getChannelRadii(mol)[:, [0, 1, 2, 3, 7]] #pick the sepecif feature we need

    aromatic_groups = get_aromatic_groups(mol._mol)
    aromatics = [coords[np.array(a_group)].mean(axis=0) for a_group in aromatic_groups]
    aromatics = np.array(aromatics) #(1, 3)
    if len(aromatics) == 0:  # Make sure the shape is correct
        aromatics = aromatics.reshape(aromatics.shape[0], 3)

    # Generate the pharmacophores
    aromatic_loc = aromatics + (np.random.rand(*aromatics.shape) - 0.5)
 
    acceptor_ph = (multisigmas[:, 2] > 0.01)
    donor_ph = (multisigmas[:, 3] > 0.01)

    # Generate locations
    acc_loc = coords[acceptor_ph]
    acc_loc = acc_loc + (np.random.rand(*acc_loc.shape) - 0.5)
    donor_loc = coords[donor_ph]

    donor_loc = donor_loc + (np.random.rand(*donor_loc.shape) - 0.5)
    coords = np.vstack([coords, aromatic_loc, acc_loc, donor_loc])

    final_sigmas = np.zeros((coords.shape[0], 8))
    final_sigmas[:n_atoms, :5] = multisigmas
    pos1 = n_atoms + len(aromatic_loc)  # aromatics end

    final_sigmas[n_atoms:(pos1), 5] = 2.
    pos2 = pos1 + len(acc_loc)
    final_sigmas[pos1:pos2, 6] = 2.
    final_sigmas[pos2:, 7] = 2.

    return final_sigmas, coords, lig_center
    

def voxelize_pkt_lig(dual_sigmas, dual_coords, dual_center, displacement=2., rotation=True):
    """
    Generates molecule representation.
    """
    pkt_sigmas, lig_sigmas = dual_sigmas
    pkt_coords, lig_coords = dual_coords
    pkt_center, lig_center = dual_center
    # Do the rotation
    if rotation:
        rrot = uniformRandomRotation()  # Rotation
        lig_coords = rotate(lig_coords, rrot, center=lig_center)
        pkt_coords = rotate(pkt_coords, rrot, center=pkt_center)

    # Do the translation
    # center = center + (np.random.rand(3) - 0.5) * 2 * displacemen

    lig_centers2D = global_centers + pkt_center
    pkt_centers2D = global_centers + pkt_center 

    pkt_occupancy = _getOccupancyC(pkt_coords.astype(np.float32),
                               pkt_centers2D.reshape(-1, 3),
                               pkt_sigmas).reshape(size, size, size, 8)

    lig_occupancy = _getOccupancyC(lig_coords.astype(np.float32),
                               lig_centers2D.reshape(-1, 3),
                               lig_sigmas).reshape(size, size, size, 8)
    return pkt_occupancy.astype(np.float32).transpose(3, 0, 1, 2,), lig_occupancy.astype(np.float32).transpose(3, 0, 1, 2,)

def vox_from_mol(pkt_mol, lig_mol):
    pkt_sigmas, pkt_coords, pkt_center = generate_sigmas(pkt_mol)
    lig_sigmas, lig_coords, lig_center = generate_sigmas(lig_mol)

    pkt_vox, lig_vox = voxelize_pkt_lig((pkt_sigmas, lig_sigmas), (pkt_coords, lig_coords), (pkt_center,lig_center), rotation=False)
    return pkt_vox, lig_vox


if __name__ == '__main__':
    sdf_file = './lig.sdf'
    pkt_file = './pkt.pdb'
    lig_mol = Chem.MolFromMolFile(sdf_file)
    pkt_mol = Chem.MolFromPDBFile(pkt_file)
    lig_mol, pkt_mol = align_pkt_lig_to_zero(lig_mol, pkt_mol)
    lig_mol = SmallMol(lig_mol)
    pkt_mol = SmallMol(pkt_mol)
    pkt_vox, lig_vox = vox_from_mol(pkt_mol, lig_mol)
    # 4 denotes the occup 
    np.save('pkt_occup.npy', pkt_vox[4])
    np.save('lig_occup.npy', lig_vox[4])