#!/usr/bin/env python3
#python makefastfp.py -in DC03_400000.sdf -fpdbfname DC03_400000.fpbin
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
# Generates a binary fingerprint file for fast fingerprint search
#############################################################################

import sys
import os
from openeye import oechem
from openeye import oegraphsim


def main(argv=[__name__]):

    itf = oechem.OEInterface()
    oechem.OEConfigure(itf, InterfaceData)
    oegraphsim.OEConfigureFingerPrint(itf, oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree))

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    ifname = itf.GetString("-in")
    ffname = itf.GetString("-fpdb")

    if oechem.OEGetFileExtension(ffname) != "fpbin":
        oechem.OEThrow.Fatal("Fingerprint database file should have '.fpbin' file extension!")

    idxfname = oechem.OEGetMolDatabaseIdxFileName(ifname)

    if not os.path.exists(idxfname):
        if not oechem.OECreateMolDatabaseIdx(ifname):
            oechem.OEThrow.Warning("Unable to create %s molecule index file", idxfname)

    print("Using %s index molecule file" % idxfname)

    moldb = oechem.OEMolDatabase()
    if not moldb.Open(ifname):
        oechem.OEThrow.Fatal("Cannot open molecule database file!")

    fptype = oegraphsim.OESetupFingerPrint(itf)
    print("Using fingerprint type %s" % fptype.GetFPTypeString())

    opts = oegraphsim.OECreateFastFPDatabaseOptions(fptype)
    opts.SetTracer(oechem.OEDots(100000, 1000, "fingerprints"))
    nrprocs = itf.GetInt("-numprocessors")
    opts.SetNumProcessors(nrprocs)

    print("Generating fingerprints with %d threads" % opts.GetNumProcessors())

    timer = oechem.OEWallTimer()
    if not oegraphsim.OECreateFastFPDatabaseFile(ffname, ifname, opts):
        oechem.OEThrow.Fatal("Cannot create fingerprint database file!")

    nrmols = moldb.GetMaxMolIdx()
    print("%5.2f secs to generate %d fingerprints" % (timer.Elapsed(), nrmols))

    return 0


InterfaceData = """
!BRIEF [-in] <input> [-fpdbfname] <output>

!CATEGORY "input/output options" 1

  !PARAMETER -in
    !ALIAS -i
    !TYPE string
    !REQUIRED true
    !KEYLESS 1
    !VISIBILITY simple
    !BRIEF Input molecule filename
  !END

  !PARAMETER -fpdbfname
    !ALIAS -fpdb
    !TYPE string
    !REQUIRED true
    !KEYLESS 2
    !VISIBILITY simple
    !BRIEF Output fingerprint database filename
  !END

!END

!CATEGORY "fingerprint generation options" 2

  !PARAMETER -numprocessors
    !ALIAS -numprocs
    !TYPE int
    !REQUIRED false
    !DEFAULT 0
    !VISIBILITY simple
    !BRIEF Number of processors used for generate fingerprints
    !DETAIL
       0 means that the number of processors available in the machines will be
       used to generate fingerprints
  !END

!END


"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))