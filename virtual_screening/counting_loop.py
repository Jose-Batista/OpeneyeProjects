#!/usr/bin/env python
from __future__ import print_function
from openeye.oechem import *

import sys
import os


def count_database(database):
    ifs = oemolistream()

    if not ifs.open(database):
        OEThrow.Fatal("Unable to open inputfile" )

    count = 0
    for mol in ifs.GetOEMols():
        count += 1
    return count

def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    ind = itf.GetString("-in_database")

    count = count_database(ind)
    print(count)

InterfaceData = """
!PARAMETER -in_database
    !ALIAS -id
    !TYPE string
    !BRIEF Input Database of Molecules
    !REQUIRED true
    !KEYLESS 1
!END

"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
