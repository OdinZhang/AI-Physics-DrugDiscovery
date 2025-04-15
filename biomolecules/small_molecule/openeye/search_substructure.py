# 2023.12.09: Odin Zhang 
from openeye import oechem
import argparse
import os.path as osp
oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Error)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_sdf', type=str, default='./query.sdf')
    parser.add_argument('--library', type=str, default='./library/chemdiv.sdf')
    parser.add_argument('--saved_file', type=str, default=None)
    parser.add_argument('--aliphatic', action='store')
    parser.add_argument('--bondtop', action='store')
    args = parser.parse_args()

    qfile = oechem.oemolistream(args.query_sdf)
    tfile = oechem.oemolistream(args.library)
    query_name = osp.basename(args.query_sdf).split('.')[0]
    library_name = osp.basename(args.library).split('.')[0]

    if args.saved_file is None:
        saved_file = f'{query_name}_{library_name}_searched.sdf'
    else:
        saved_file = args.saved_file
    
    # set the same aromaticity model for the query and the target file
    aromodel = oechem.OEIFlavor_Generic_OEAroModelMDL
    qflavor = qfile.GetFlavor(qfile.GetFormat())
    qfile.SetFlavor(qfile.GetFormat(), (qflavor | aromodel))
    tflavor = tfile.GetFlavor(tfile.GetFormat())
    tfile.SetFlavor(tfile.GetFormat(), (tflavor | aromodel))

    # read MDL query and initialize the substructure search
    opts = oechem.OEMDLQueryOpts_Default | oechem.OEMDLQueryOpts_SuppressExplicitH
    if args.aliphatic:
        opts = oechem.OEMDLQueryOpts_Default | oechem.OEMDLQueryOpts_SuppressExplicitH | oechem.OEMDLQueryOpts_AddBondAliphaticConstraint
    if args.bondtop:
        opts = oechem.OEMDLQueryOpts_Default | oechem.OEMDLQueryOpts_SuppressExplicitH |oechem.OEMDLQueryOpts_AddBondTopologyConstraint
    qmol = oechem.OEQMol()

    oechem.OEReadMDLQueryFile(qfile, qmol, opts)
    ss = oechem.OESubSearch(qmol)

    # loop over target structures
    ofs = oechem.oemolostream(saved_file)
    tindex = 1
    for tmol in tfile.GetOEGraphMols():
        oechem.OEPrepareSearch(tmol, ss)
        if ss.SingleMatch(tmol):
            print("hit target =", tindex)
            oechem.OEWriteMolecule(ofs, tmol)
        tindex += 1
        if tindex%10000 == 0:
            print(f'Finished Searching {tindex} Files')

    # Close the file streams
    qfile.close()
    tfile.close()
    ofs.close()