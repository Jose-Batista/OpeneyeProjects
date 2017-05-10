#!/usr/bin/env python
#############################################################################
# Copyright (C) 2017 OpenEye Scientific Software, Inc.
#############################################################################
# Generates binary fingerprint file for fast fingerprint search
#############################################################################

import sys
import os
from openeye.oechem import *
from openeye.oegraphsim import *
from operator import itemgetter

def searchfastfp( query, mfname, molffname, itf ):

    timer = OEWallTimer()
    timer.Start()

    memtype = OEFastFPDatabaseMemoryType_MemoryMapped
    memtypestr = "memory-mapped"
    fpdb = OEFastFPDatabase(molffname, memtype)
    if not fpdb.IsValid():
        OEThrow.Fatal("Cannot open fingerprint database!")
    nrfps = fpdb.NumFingerPrints()

    moldb = OEMolDatabase()
    if not moldb.Open(mfname):
        OEThrow.Fatal("Cannot open molecule database!")

    if not OEAreCompatibleDatabases(moldb, fpdb):
        OEThrow.Fatal("Databases are not compatible!")

    OEThrow.Info("%5.2f sec to initialize databases" % timer.Elapsed())

    fptype = fpdb.GetFPTypeBase()
    OEThrow.Info("Using fingerprint type %s" % fptype.GetFPTypeString())

    opts = OEFPDatabaseOptions()
    OESetupFPDatabaseOptions(opts, itf)

    # search fingerprint database

    timer.Start()
    scores = fpdb.GetSortedScores(query, opts)
    OEThrow.Info("%5.2f sec to search %d fingerprints %s" % (timer.Elapsed(), nrfps, memtypestr))

    timer.Start()
    nrhits = 0
    hit = OEGraphMol()
    hitlist = list()
    for si in scores:
        if moldb.GetMolecule(hit, si.GetIdx()):
            nrhits += 1
            OESetSDData(hit, "Similarity score", "%.2f" % si.GetScore())
            hitlist.append( [hit.GetTitle(),si.GetScore()] )
    OEThrow.Info("%5.2f sec to write %d hits" % (timer.Elapsed(), nrhits))
    return hitlist

def createFastFP( ifname, ffname, itf):

    if OEGetFileExtension(ffname) != "fpbin":
        OEThrow.Fatal("Fingerprint database file should have '.fpbin' file extension!")

    idxfname = OEGetMolDatabaseIdxFileName(ifname)

    if not os.path.exists(idxfname):
        if not OECreateMolDatabaseIdx(ifname):
            OEThrow.Warning("Unable to create %s molecule index file", idxfname)

    OEThrow.Info("Using %s index molecule file" % idxfname)

    moldb = OEMolDatabase()
    if not moldb.Open(ifname):
        OEThrow.Fatal("Cannot open molecule database file!")

    nrmols = moldb.GetMaxMolIdx()

    fptype = OESetupFingerPrint(itf)
    OEThrow.Info("Using fingerprint type %s" % fptype.GetFPTypeString())

    opts = OECreateFastFPDatabaseOptions(fptype)
    opts.SetTracer(OEDots(100000, 1000, "fingerprints"))
    OEThrow.Info("Generating fingerprints with %d threads" % opts.GetNumProcessors())

    timer = OEWallTimer()
    if not OECreateFastFPDatabaseFile(ffname, ifname, opts):
        OEThrow.Fatal("Cannot create fingerprint database file!")

    OEThrow.Info("%5.2f secs to generate %d fingerprints" % (timer.Elapsed(), nrmols))


def GetQuery( ifn ):
    ifs = oemolistream()
    if not ifs.open( ifn ):
        OEThrow.Fatal( "Unable to open file %s for reading" % ifn )
    
    mol = OEGraphMol()
    OEReadMolecule( ifs, mol )

    ifs.close()
    return mol

def searchslowfp(qmol, ifname):

    ifs = oemolistream()
    if not ifs.open( ifname ):
        OEThrow.Fatal( "Unable to open file %s for reading" % ifname )

    qfp = OEFingerPrint()    
    OEMakeFP(qfp, qmol, OEFPType_Tree) 
    OEThrow.Info("Slow sarch is using fingerprint type %s" % qfp.GetFPTypeBase().GetFPTypeString() )

    hitlist = list()
    for mol in ifs.GetOEGraphMols():
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, OEFPType_Tree)
        hitlist.append([mol.GetTitle(), OETanimoto(qfp, fp)])
        
    tmp = sorted(hitlist, key=itemgetter(1), reverse=True)
    return tmp

def main(argv=[__name__]):

    itf = OEInterface()
    OEConfigure(itf, InterfaceData)
    OEConfigureFingerPrint(itf, OEGetFPType(OEFPType_Tree))

    if not OEParseCommandLine(itf, argv):
        return 1

    ifname = itf.GetString("-in")
    fpname = itf.GetString("-fpdb")
    qfname = itf.GetString("-query")

    qmol = GetQuery(qfname)
    createFastFP( ifname, fpname, itf )
    fastfphitlist = searchfastfp(qmol,ifname,fpname,itf)  
    
    for i in range(10):
        print(fastfphitlist[ i ])

    fphitlist = searchslowfp(qmol,ifname)

    for i in range(10):
        print(fphitlist[ i ])

    return 0


InterfaceData = """
!BRIEF [-in] <input> [-fpdbfname] <output>

!CATEGORY "input/output options"

  !PARAMETER -in
    !ALIAS -i
    !TYPE string
    !REQUIRED true
    !KEYLESS 1
    !VISIBILITY simple
    !BRIEF Input molecule filename
  !END

  !PARAMETER -query
    !ALIAS -q
    !TYPE string
    !REQUIRED true
    !KEYLESS 2
    !VISIBILITY simple
    !BRIEF query molecule filename
  !END

  !PARAMETER -fpdbfname
    !ALIAS -fpdb
    !TYPE string
    !REQUIRED true
    !KEYLESS 3
    !VISIBILITY simple
    !BRIEF Output fingerprint database filename
  !END

!END
"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))
