#!/usr/bin/env python3
#python searchfastfp.py -fpdb DC01_400000.fpbin -molfname DC01_400000.sdf -query KRD4.sdf -out hits.sdf -memorytype in-memory
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

#############################################################################
# Searches fast a fingerprint database
#############################################################################

import sys
from openeye import oechem
from openeye import oegraphsim


def main(argv=[__name__]):

    itf = oechem.OEInterface()
    oechem.OEConfigure(itf, InterfaceData)

    defopts = oegraphsim.OEFPDatabaseOptions(1, oegraphsim.OESimMeasure_Tanimoto)
    oegraphsim.OEConfigureFPDatabaseOptions(itf, defopts)
    oegraphsim.OEConfigureFPDatabaseMemoryType(itf)

    if not oechem.OEParseCommandLine(itf, argv):
        return 0

    qfname = itf.GetString("-query")
    mfname = itf.GetString("-molfname")
    ffname = itf.GetString("-fpdbfname")
    ofname = itf.GetString("-out")

    # initialize databases

    timer = oechem.OEWallTimer()
    timer.Start()

    ifs = oechem.oemolistream()
    if not ifs.open(qfname):
        oechem.OEThrow.Fatal("Cannot open input file!")

    ofs = oechem.oemolostream()
    if not ofs.open(ofname):
        oechem.OEThrow.Fatal("Cannot open output file!")

    query = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, query):
        oechem.OEThrow.Fatal("Cannot read query molecule!")

    moldb = oechem.OEMolDatabase()
    if not moldb.Open(mfname):
        oechem.OEThrow.Fatal("Cannot open molecule database!")

    memtype = oegraphsim.OEGetFPDatabaseMemoryType(itf)

    fpdb = oegraphsim.OEFastFPDatabase(ffname, memtype)
    if not fpdb.IsValid():
        oechem.OEThrow.Fatal("Cannot open fingerprint database!")
    nrfps = fpdb.NumFingerPrints()
    memtypestr = fpdb.GetMemoryTypeString()

    if not oegraphsim.OEAreCompatibleDatabases(moldb, fpdb):
        oechem.OEThrow.Fatal("Databases are not compatible!")

    print("%5.2f sec to initialize databases" % timer.Elapsed())

    fptype = fpdb.GetFPTypeBase()
    print("Using fingerprint type %s" % fptype.GetFPTypeString())

    opts = oegraphsim.OEFPDatabaseOptions()
    oegraphsim.OESetupFPDatabaseOptions(opts, itf)

    # search fingerprint database

    timer.Start()
    scores = fpdb.GetSortedScores(query, opts)
    print("%5.2f sec to search %d fingerprints %s" % (timer.Elapsed(), nrfps, memtypestr))

    timer.Start()
    hit = oechem.OEGraphMol()
    for si in scores:
        if moldb.GetMolecule(hit, si.GetIdx()):
            oechem.OESetSDData(hit, "Similarity score", "%.2f" % si.GetScore())
            oechem.OEWriteMolecule(ofs, hit)
    print("%5.2f sec to write %d hits" % (timer.Elapsed(), opts.GetLimit()))

    return 0


InterfaceData = """
!BRIEF [-query] <molfile> [-molfname] <molfile> [-fpdbfname] <fpfile>  [-out] <molfile>

!CATEGORY "input/output options"

  !PARAMETER -query
    !ALIAS -q
    !TYPE string
    !REQUIRED true
    !KEYLESS 1
    !VISIBILITY simple
    !BRIEF Input query filename
  !END

  !PARAMETER -molfname
    !ALIAS -mol
    !TYPE string
    !REQUIRED true
    !KEYLESS 2
    !VISIBILITY simple
    !BRIEF Input molecule filename
  !END

  !PARAMETER -fpdbfname
    !ALIAS -fpdb
    !TYPE string
    !REQUIRED true
    !KEYLESS 3
    !VISIBILITY simple
    !BRIEF Input fast fingerprint database filename
  !END

  !PARAMETER -out
    !ALIAS -o
    !TYPE string
    !REQUIRED true
    !KEYLESS 4
    !VISIBILITY simple
    !BRIEF Output molecule filename
  !END

!END
"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))