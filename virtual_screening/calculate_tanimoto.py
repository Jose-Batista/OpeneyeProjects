
#!/usr/bin/env python
from __future__ import print_function
from openeye.oechem import *
from openeye.oegraphsim import *
import sys

import json,requests
import urllib.parse as parse

def read_db(input_query, fptype):
    ifs = oemolistream()

    if not ifs.open(input_query):
        OEThrow.Fatal("Unable to open inputfile" )

    fp_list = list()
    mol_list = []
    for mol in ifs.GetOEGraphMols():
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, fptype)
        mol_list.append(mol.CreateCopy())
        fp_list.append(fp)
    return mol_list, fp_list


def calculate_tanimoto(dbfp_list, db_list, fp, mol, fptype):
    fp = fp[0]
    print("Molecule : %s   Smiles : %s" %(mol[0].GetTitle(), OEMolToSmiles(mol[0]))) 
    print("Results")
    tanimoto_list = list()
    for i, dbmol in enumerate(db_list):
        dbfp = dbfp_list[i]
        tanimoto = OETanimoto(fp, dbfp)
        print("%s, %s  Tanimoto value : %4f" %(dbmol.GetTitle(), OEMolToSmiles(dbmol), tanimoto))
        tanimoto_list.append((dbmol.GetTitle(), OEMolToSmiles(dbmol), tanimoto))
        tanimoto_list = sorted(tanimoto_list, reverse = True, key = lambda mol:mol[2])
        if len(tanimoto_list) > 100:
            tanimoto.pop()
    return tanimoto_list

def write_output(mol, tanimoto_list, out):
    tanimoto_results = open(out, "a")
    tanimoto_results.write("Molecule : %s   Smiles : %s\n" %(mol[0].GetTitle(), OEMolToSmiles(mol[0])))
    for tanimoto in tanimoto_list:
        tanimoto_results.write("%s, %s Tanimoto value : %.4f\n" %(tanimoto[0], tanimoto[1], tanimoto[2]))
    tanimoto_results.write("\n")
    tanimoto_results.close()

def request_tanimoto(mol, baseurl, data):
    ranking = list()
    smiles = OEMolToSmiles(mol[0])
    safe_smiles = parse.quote(smiles)
    url = "%s/%s/hitlist?smiles=%s&oformat=csv" %(baseurl, data['databases'][0], safe_smiles) 
    response = requests.get( url )
    hitlist = response.content.decode().split('\n')
    hitlist.pop(0)
    hitlist.pop()
    for mol in hitlist:
        cur_mol = mol.split(',')
        ranking.append((cur_mol[1], cur_mol[0], float(cur_mol[4])))
        print((cur_mol[1], cur_mol[0], float(cur_mol[4])))
    return ranking


def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    ind = itf.GetString("-in_database")
    inm = itf.GetString("-in_molecule")
    out = itf.GetString("-output")
    fptype = itf.GetInt("-fprint")

    (db_list, dbfp_list) = read_db(ind, fptype)
    (mol, fp) = read_db(inm, fptype)
    tanimoto_list = calculate_tanimoto(dbfp_list, db_list, fp, mol, fptype)
    write_output(mol, tanimoto_list, out)

    baseurl = "http://10.0.1.22:8069"
    response = requests.get( baseurl )
    data = response.json()

    tanimoto_list = request_tanimoto(mol, baseurl, data)
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
  !REQUIRED true
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
