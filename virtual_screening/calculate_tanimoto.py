
#!/usr/bin/env python
from __future__ import print_function
from openeye.oechem import *
from openeye.oegraphsim import *
import sys

def read_db(input_query, fptype):
    ifs = oemolistream()

    if not ifs.open(input_query):
        OEThrow.Fatal("Unable to open inputfile" )

    mol_list = []
    for mol in ifs.GetOEGraphMols():
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, fptype)
        mol.SetData(str(fptype), fp)
        mol_list.append(mol.CreateCopy())
    return mol_list


def calculate_tanimoto(db_list, mol, fptype):
    fp = mol[0].GetData(str(fptype))
    print("Molecule : %s   Smiles : %s" %(mol[0].GetTitle(), OEMolToSmiles(mol[0]))) 
    print("Results")
    tanimoto_list = list()
    for dbmol in db_list:
        dbfp = dbmol.GetData(str(fptype))
        tanimoto = OETanimoto(fp, dbfp)
        print("%s, %s  Tanimoto value : %4f" %(dbmol.GetTitle(), OEMolToSmiles(dbmol), tanimoto))
        tanimoto_list.append((dbmol.GetTitle(), OEMolToSmiles(dbmol), tanimoto))
        tanimoto_list = sorted(tanimoto_list, reverse = True, key = lambda mol:mol[2])
    return tanimoto_list

def write_output(mol, tanimoto_list, out):
    tanimoto_results = open(out, "w")
    tanimoto_results.write("Molecule : %s   Smiles : %s\n" %(mol[0].GetTitle(), OEMolToSmiles(mol[0])))
    for tanimoto in tanimoto_list:
        tanimoto_results.write("%s, %s Tanimoto value : %.4f\n" %(tanimoto[0], tanimoto[1], tanimoto[2]))
    tanimoto_results:close()

def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    ind = itf.GetString("-in_database")
    inm = itf.GetString("-in_molecule")
    out = itf.GetString("-output")
    fptype = itf.GetInt("-fprint")

    db_list = read_db(ind, fptype)
    mol = read_db(inm, fptype)
    tanimoto_list = calculate_tanimoto(db_list, mol, fptype)
    write_output(mol, tanimoto_list, out)

InterfaceData = """

!PARAMETER -in_database
  !ALIAS -id
  !TYPE string
  !BRIEF Input Database of Molecule file
  !REQUIRED true
  !KEYLESS 1
!END

!PARAMETER -in_molecule
  !ALIAS -im
  !TYPE string
  !BRIEF Input of Molecule file
  !REQUIRED true
  !KEYLESS 2
!END

!PARAMETER -output
  !ALIAS -o
  !TYPE string
  !BRIEF Output File
  !REQUIRED false
  !KEYLESS 3
!END

!PARAMETER -fprint
  !ALIAS -fp
  !TYPE int
  !BRIEF Fingerprint Type (101 for MACCS, 102 for Path, 103 for Lingo, 104 for Circular, 105 for Tree)
  !REQUIRED true
  !KEYLESS 4
!END
"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))
