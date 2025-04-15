rdkit rotate molecules https://iwatobipen.wordpress.com/2019/02/27/rotate-molecule-and-visualize-it-rdkit-chemoinformatics/
rdkit rtz https://github.com/wutobias/r2z
https://chemistry.stackexchange.com/questions/136930/computational-chemistry-software-to-generate-files-of-cartesian-coordinates-and
rdkit切分fragment
https://zhuanlan.zhihu.com/p/390678258
rdkit显示键的编号
https://zhuanlan.zhihu.com/p/163827512
二面角扫描
https://pschmidtke.github.io/blog/torsion/dihedral/oss/opensource/rdkit/xtb/energy/2021/02/16/torsion-angle-scans-xtb.html
https://github.com/vermasrijan/redial-2020/blob/5cfd09d5645301654373cb9008ded1b2552bd920/mayachemtools/bin/RDKitPerformTorsionScan.py#L588
fragment_utils
https://github.com/deepchem/deepchem/blob/9664aeab16940e1ee46030cac08d90e12aa24751/deepchem/utils/test/test_fragment_utils.py

不同rdkit版本的容错率还不一样，如果直接从mol文件当中读取分子，比如
mol = Chem.MolFromMolFile(lig_path) 的话，很多时候会报valence的错误
这个时候mol = Chem.MolFromMolFile(lig_path,sanitize=False) 就可以继续操作了

读取sdf文件
suppl = Chem.SDMolSupplier(lig_path, sanitize=False)
mols = [mol for mol in suppl if mol]	
mol = mols[0]

读取mol2文件
mol = Chem.MolFromMol2File(mol2_path, sanitize=False)

屏蔽Rdkit的Error输出
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
smi = 'CO(C)C'
mol = Chem.MolFromSmiles(smi)
try:
    Chem.SanitizeMol(mol)
except:
    print('OK')
   
# rdkit内部会有拓扑匹配的事情，所以要非常小心地处理拓扑问题。拿到新的分子最好的方式就是先去H然后再加氢
def cannaical_mol(mol, 3d=True, addhs=True):
    mol_ = copy.deepcopy(mol)
    mol_ = Chem.RemoveAllHs(mol_)
    smi = Chem.MolToSmiles(mol_)
    if 3D is True and addhs is True:
        mol_ = Chem.AddHs(mol_, addCoords=True)
    if 3D is not True and addhs is True:
        mol_ = Chem.AddHs(mol_)
    return mol_, smi
    return mol_, smi

def rm_radical(mol):
    mol = copy.deepcopy(mol)
    for atom in mol.GetAtoms():
        atom.SetNumRadicalElectrons(0)
    return mol

#rdkit的各种blog
#https://greglandrum.github.io/rdkit-blog/prototypes/technical/2020/01/25/trying-the-tautomer-canonicalization-code.html

#rdkit显示原子序号
#https://www.jianshu.com/p/ec4f3b9e57f9
# method1
mol = Chem.MolFromSmiles('c1ccccc(C(N)=O)1')
for atom in mol.GetAtoms():
    atom.SetProp("atomNote", str(atom.GetIdx()))
mol
# method2
from rdkit import Chem

def showAtomNum(smi):
    try:
        mol=Chem.MolFromSmiles(smi)
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetProp('molAtomMapNumber',str(i))
        return Chem.MolToSmiles(mol)
    except Exception as e:
        print (e)
        return None

#https://molvs.readthedocs.io/en/latest/
#Google Summer Code Project

# 重新规划rdkit的顺序
from rdkit.Chem import rdmolfiles
new_order1 = rdmolfiles.CanonicalRankAtoms(m1)
mol1 = rdmolops.RenumberAtoms(m1, new_order1)

# MolToGrid, and save the svg file
def plot_rdkit_svg_grid(mols, mols_per_row=5, filename=None, **kwargs):
    """
    Plots a grid of RDKit molecules in SVG.
    :param mols: a list of RDKit molecules
    :param mols_per_row: size of the grid
    :param filename: save an image with the given filename
    :param kwargs: additional arguments for `RDKit.Chem.Draw.MolsToGridImage`
    :return: the SVG as a string
    """
    svg = Draw.MolsToGridImage(mols, molsPerRow=mols_per_row, useSVG=True, **kwargs)
    if filename is not None:
        if not filename.endswith('.svg'):
            filename += '.svg'
        with open(filename, 'w') as f:
            f.write(svg)
    return svg 

#rdkit进阶教程
#http://rdkit.chenzhaoqiang.com/mediaManual.html#id8

#rdkit 计算shape and electrostatic similarity 
#https://iwatobipen.wordpress.com/2021/12/10/compare-shape-and-electrostatic-similarity-of-molecules-rdkit-espsim-python/

#单个分子highlight
Draw.MolsToGridImage([mol],highlightAtomLists=[anchor_idx])
#存储highlight
Draw.MolToFile(mol, 'test.png',highlightAtoms=anchor_idx )
