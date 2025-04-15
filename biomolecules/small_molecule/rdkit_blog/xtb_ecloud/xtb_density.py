"""
Calculate charge density using xTB, only works for GNU/Linux.
"""

import os
import uuid
import numpy
import shutil
import rdkit as rd
import rdkit.Chem
try:
    import cubtools
    from cubtols import write_cube
except: 
    import utils.cubtools as cubtools
    from utils.cubtools import write_cube

import numpy as np


BOHR = 1.8897259886
XTBINP = r"""
$cube
    step={:f}
$end
$write
    density=true
    spin density=false
    fod=false
    charges=false
    mulliken=false
    geosum=false
    inertia=false
    mos=false
    wiberg=false
$end
""".format(0.5 * BOHR)


class CDCalculator(object):
    '''
    simple wrapper for xTB to calculate charge density
    Constructor:
        xtb_command: the command to run xTB, default is 'xtb', or you can use the binary installation path
    Functions:
        calculate(mol): calculate the charge density of a molecule
        clean(): clean the temporary files
    '''

    def __init__(self, xtb_command: str = 'xtb') -> None:
        self.clean()
        self.rootdir = os.getcwd()
        self.workdir = os.getcwd() + '/' + uuid.uuid4().hex
        os.mkdir(self.workdir)
        self.xtb_command = xtb_command
        os.environ['KMP_STACKSIZE'] = '1G'
        os.environ['OMP_NUM_THREADS'] = '1'
        os.environ['MKL_NUM_THREADS'] = '1'
        with open(self.workdir + '/xtb.inp', 'w') as fp:
            print(XTBINP, file=fp, flush=True)
        self._errno = 0

    def __del__(self):
        self.clean()

    def clean(self):
        try:
            shutil.rmtree(self.workdir)
        except:
            pass

    def calculate(self, mol: rd.Chem.rdchem.Mol) -> dict:
        os.chdir(self.workdir)
        rd.Chem.rdmolfiles.MolToXYZFile(mol, 'tmp.xyz')
        command = self.xtb_command + ' --norestart --input xtb.inp tmp.xyz' + \
                                     ' 2> ' + os.devnull + ' > log'
        self._errno = os.system(command)
        if (self._errno == 0):
            try:
                density, meta = cubtools.read_cube('density.cub')
                meta['org'] = numpy.array(meta['org']) / BOHR
                meta['len'] = (numpy.array(meta['xvec']) * meta['nx'] +
                               numpy.array(meta['yvec']) * meta['ny'] +
                               numpy.array(meta['zvec']) * meta['nz']) / BOHR
                # meta.pop('atoms')
                result = {'density': density * (BOHR ** 3), 'meta': meta}
            except Exception as e:
                print(e)
                result = {}
        else:
            result = {}
        os.chdir(self.rootdir)
        return result

    @property
    def err_msg(self):
        if (self._errno == 0):
            return ''
        else:
            return open(self.workdir + 'log', 'r').read()

from scipy.interpolate import RegularGridInterpolator

def ecloud2grid(cub):
    """
    compute the grid coordinates of the ecloud dict
    Parameters:
        cub: ecloud dict
    """

    x_len, y_len, z_len = cub['meta']['len'] # three dimensional length of the box
    x_llc, y_llc, z_llc = cub['meta']['org']
    x_cell, y_cell, z_cell = cub['meta']['nx'], cub['meta']['ny'], cub['meta']['nz']


    X, Y, Z = np.meshgrid(np.arange(float(x_cell)), np.arange(float(y_cell)), np.arange(float(z_cell)))

    X *= x_len / x_cell 
    Y *= y_len / y_cell  
    Z *= z_len / z_cell 
    # grid_coords = np.stack([X, Y, Z], axis=-1)
    return X+x_llc, Y+y_llc, Z+z_llc

def interplot_ecloud(ecloud, new_gridcoord):
    """ 
    Assign the density of the ecloud to the new coordinates
    Parameters:
        ecoud: ecloud dict
        new_gridcoord: new coordinates of the ecloud
    """
    X, Y, Z = ecloud2grid(ecloud)
    X_1d = X[0, :, 0]
    Y_1d = Y[:, 0, 0]
    Z_1d = Z[0, 0, :]
    
    interpolator = RegularGridInterpolator((X_1d, Y_1d, Z_1d), ecloud['density'], 
                                        bounds_error=False, fill_value=0.0)
    new_X,new_Y, new_Z = new_gridcoord
    new_points = np.array([new_X.flatten(), new_Y.flatten(), new_Z.flatten()]).T
    new_density = interpolator(new_points)

    return new_density

def write_new_cube(new_gridcoord, new_density, meta, fname):
    """
    Write the new cube file
    Parameters:
        new_gridcoord: new coordinates of the ecloud (Unit: A) shape = (3, nx, ny, nz)
        new_density: new density of the ecloud
        meta: meta data of the ecloud
        fname: name of the new cube file
    """
    new_X,new_Y, new_Z = new_gridcoord
    vec_x = new_X[1, 0, 0] - new_X[0, 0, 0]
    vec_y = new_Y[0, 1, 0] - new_Y[0, 0, 0]
    vec_z = new_Z[0, 0, 1] - new_Z[0, 0, 0]

    new_meta = {'org':np.array([np.min(new_X), np.min(new_Y), np.min(new_X)])*BOHR, #llc
            'xvec': np.array([vec_x, 0, 0])*BOHR,
            'yvec': np.array([0, vec_y, 0])*BOHR,
            'zvec': np.array([0, 0, vec_z])*BOHR,
            'nx': new_gridcoord[0].shape[0],
            'ny': new_gridcoord[0].shape[1],
            'nz': new_gridcoord[0].shape[2],
            'atoms': meta['atoms'],
    }

    new_density = new_density.reshape(new_X.shape) 
    write_cube(new_density, new_meta, fname)

def BuildGridCenters(llc, N, step):
    """
    llc: lower left corner
    N: number of cells in each direction
    step: step size
    """

    if type(step) == float:
        xrange = [llc[0] + step * x for x in range(0, N[0])]
        yrange = [llc[1] + step * x for x in range(0, N[1])]
        zrange = [llc[2] + step * x for x in range(0, N[2])]
    elif type(step) == list or type(step) == tuple:
        xrange = [llc[0] + step[0] * x for x in range(0, N[0])]
        yrange = [llc[1] + step[1] * x for x in range(0, N[1])]
        zrange = [llc[2] + step[2] * x for x in range(0, N[2])]

    centers = np.zeros((N[0], N[1], N[2], 3))
    for i, x in enumerate(xrange):
        for j, y in enumerate(yrange):
            for k, z in enumerate(zrange):
                centers[i, j, k, :] = np.array([x, y, z])
    return centers

step = 0.5
size = 24
cell_number = [size, size, size]
llc = (np.zeros(3) - float(size * step / 2))
# Now, the box is 24×24×24 A^3

expanded_pcenters = BuildGridCenters(llc, cell_number, step)
# Now, the expanded_pcenters is the coordinates of the grid points


if __name__ == '__main__':
    # reconstruct cube file
    from rdkit import Chem
    from cubtools import write_cube
    import numpy as np

    claculator = CDCalculator(xtb_command='xtb')
    mol = Chem.MolFromMolFile('./90.sdf')
    ecloud = claculator.calculate(mol)
    recon_meta = {'org': ecloud['meta']['org']*BOHR,
                'xvec': np.array(ecloud['meta']['xvec']),
                'yvec': np.array(ecloud['meta']['yvec']),
                'zvec': np.array(ecloud['meta']['zvec']),
                'nx': ecloud['meta']['nx'],
                'ny': ecloud['meta']['ny'],
                'nz': ecloud['meta']['nz'],
                'len': ecloud['meta']['len']*BOHR,
                'atoms': ecloud['meta']['atoms']}

    write_cube(ecloud['density']/(BOHR**3), recon_meta, './recon.cub')
    # recon_x, recon_y, recon_z = cub2grid(result)
    # np.sum(ecloud['density'] != result['density'])
    # TEST: Interpolate the ecloud to the new coordinates
    new_density = interplot_ecloud(ecloud, expanded_pcenters.transpose(3, 0, 1, 2))
    write_new_cube(expanded_pcenters.transpose(3,0,1,2), new_density/ (BOHR**3), ecloud['meta'], './new.cub')



# Provide a example of the 90.sdf file
"""


     RDKit          3D

 51 54  0  0  0  0  0  0  0  0999 V2000
   -5.1035   -1.7655   -0.4640 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9370   -1.6108    0.2183 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1679   -0.5627   -0.4487 C   0  0  1  0  0  0  0  0  0  0  0  0
   -4.8805    0.7129   -0.9008 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9268    1.9018   -0.8241 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4458    2.8259   -1.0836 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1118    1.7545   -1.5353 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3495    2.0255    0.5892 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6973    0.7948    0.9705 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3820    0.5679    1.1444 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9598   -0.5599    1.3603 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3986    1.7388    1.0121 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.1404    1.8917   -0.5009 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1808    2.5847   -0.8246 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2841    2.0299   -0.0859 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.0734    0.6945   -0.6722 S   0  0  2  0  0  0  0  0  0  0  0  0
    4.3763    0.7110   -0.0717 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9715    0.7271   -2.0983 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3085   -0.8326   -0.0589 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9620   -1.5131    0.9247 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4805   -2.7319    1.4664 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3288   -3.3095    1.0239 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6331   -2.6459   -0.0043 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4871   -2.9383   -0.6336 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6879   -1.9636   -1.5065 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2611   -1.0295   -1.4869 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.1234   -1.4007   -0.5605 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9415   -4.2299    1.4276 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0532   -3.1995    2.2526 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8834   -1.0914    1.2958 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1155    2.1152    1.3417 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8635    1.3788    1.8125 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0212    0.3032    1.7243 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6882    1.5964    2.8681 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0608    3.1767    1.6087 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0099    1.6883    1.7977 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1313    3.6530   -0.5845 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4148    2.4661   -1.8826 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1179    0.9053   -0.9658 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9545    2.4570   -0.9566 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8160    2.6646    1.4155 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6042   -0.3397    0.9574 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0465   -1.2154    1.2910 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4210   -0.1442    1.6610 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6620    2.8668    0.6453 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1687    2.2061    1.2969 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2458    0.5953   -1.9222 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7435    0.8927   -0.2544 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3284   -0.7581   -1.1242 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5705   -2.6630   -0.1614 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5017   -1.9220   -1.4633 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1 50  1  0
  1 51  1  0
  3  4  1  0
  3 42  1  0
  3 49  1  6
  4  5  1  0
  4 47  1  0
  4 48  1  0
  5  6  1  0
  5  7  1  0
  5  8  1  0
  8  9  1  0
  8 45  1  0
  8 46  1  0
  9 10  1  0
  9 42  1  0
 10 11  2  0
 10 12  1  0
 12 13  1  0
 12 32  1  0
 12 41  1  1
 13 14  1  0
 13 39  1  0
 13 40  1  0
 14 15  1  0
 14 37  1  0
 14 38  1  0
 16 15  1  1
 15 31  1  0
 16 17  2  0
 16 18  2  0
 16 19  1  0
 19 20  2  0
 19 27  1  0
 20 21  1  0
 20 30  1  0
 21 22  2  0
 21 29  1  0
 22 23  1  0
 22 28  1  0
 23 24  2  0
 23 27  1  0
 24 25  1  0
 25 26  1  0
 26 27  2  0
 31 32  1  0
 31 35  1  0
 31 36  1  0
 32 33  1  0
 32 34  1  0
 42 43  1  0
 42 44  1  0
M  END
$$$$
"""