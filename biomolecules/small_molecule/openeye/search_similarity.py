# Odin Zhang 20240410 :(
# python search_similarity.py --query_sdf ./frag2.sdf --molfname ./library/Specs_ExAcD_Aug_2020.sdf --fpdbfname ./library/Specs_ExAcD_Aug_2020.fpbin --saved_sdf ./case_examples/sim_searched.sdf --num_return 50
import sys
from openeye import oechem
from openeye import oegraphsim
import os.path as osp
import argparse
oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Error)

if __name__ == '__main__':
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Fast similarity search using OpenEye")
    parser.add_argument('--query_sdf', type=str, default='./query.sdf', help='Query SDF file')
    parser.add_argument('--molfname', type=str, default='./library/chemdiv.sdf', help='Molecule library SDF file')
    parser.add_argument('--fpdbfname', type=str, default='./library/chemdiv.fpbin', help='Fingerprint database file')
    parser.add_argument('--saved_sdf', type=str, default=None, help='Output SDF file for saving results')
    parser.add_argument('--num_return', type=int, default=1, help='Number of molecules to return')
    args = parser.parse_args()

    # Initialize databases
    timer = oechem.OEWallTimer()
    timer.Start()

    ifs = oechem.oemolistream()
    if not ifs.open(args.query_sdf):
        oechem.OEThrow.Fatal("Cannot open input file!")

    ofs = oechem.oemolostream()
    if not ofs.open(args.saved_sdf):
        oechem.OEThrow.Fatal("Cannot open output file!")

    query = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, query):
        oechem.OEThrow.Fatal("Cannot read query molecule!")

    moldb = oechem.OEMolDatabase()
    if not moldb.Open(args.molfname):
        oechem.OEThrow.Fatal("Cannot open molecule database!")

    fpdb = oegraphsim.OEFastFPDatabase(args.fpdbfname)
    if not fpdb.IsValid():
        oechem.OEThrow.Fatal("Cannot open fingerprint database!")
    nrfps = fpdb.NumFingerPrints()

    if not oegraphsim.OEAreCompatibleDatabases(moldb, fpdb):
        oechem.OEThrow.Fatal("Databases are not compatible!")

    print("%5.2f sec to initialize databases" % timer.Elapsed())

    fptype = fpdb.GetFPTypeBase()
    print("Using fingerprint type %s" % fptype.GetFPTypeString())

    # Configure search options
    opts = oegraphsim.OEFPDatabaseOptions(args.num_return, oegraphsim.OESimMeasure_Tanimoto)

    # Search fingerprint database
    timer.Start()
    scores = fpdb.GetSortedScores(query, opts)
    print("%5.2f sec to search %d fingerprints" % (timer.Elapsed(), nrfps))

    # Write hits
    timer.Start()
    cnt = 0
    for si in scores:
        hit = oechem.OEGraphMol()
        if moldb.GetMolecule(hit, si.GetIdx()):
            oechem.OESetSDData(hit, "Similarity score", "%.2f" % si.GetScore())
            oechem.OEWriteMolecule(ofs, hit)
            cnt += 1
    print("%5.2f sec to write %d hits" % (timer.Elapsed(), cnt))
