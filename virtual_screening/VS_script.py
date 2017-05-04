#!/usr/bin/env python
from __future__ import print_function
from openeye.oechem import *
from openeye.oegraphsim import *

import sys
import os
import random
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_database(database):
	ifs = oemolistream()

	if not ifs.open(database):
		OEThrow.Fatal("Unable to open inputfile" )

	mol_list = []
	for mol in ifs.GetOEGraphMols()
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, fptype)
        mol.SetData(str(fptype),fp)
		mol_list.append(mol.CreateCopy())
	return mol_list

def ReadIndex(index_input):
	index_log = open(index_input, "r")
    index = index_log.read()
    index_list = index.split('set N°')
    index_list = index_list[1:]
    for set_id, random_set in enumerate(index_list):
        random_set = random_set.split(' ')
        random_set = random_set[1:-1]
        for i, idx in enumerate(random_set):
            random_set[i] = int(idx)
        index_list[set] = random_set

    index_log.close()
	return index_lists

def RankActives(act_list, index_list, fptype, topn, nb_ka):
    ranking_list = list()

    for i, baitset in enumerate(index_list):
        ranking = pd.DataFrame(columns=["Molecule", "idx", "Tanimoto", "Rank", "KA"])
        c = 0
        for idx in baitset:
            while c < idx:
                fp = act_list[idx].GetData(str(fptype))
                simval = GetSimValAgainstAC(fp, act_list, baitset, fptype)
                OESetSDData(act_list[idx], "Similarity Value (Tanimoto) :", str(simval))
                OESetSDData(act_list[idx], "Trial Set :", str(i))
                OESetSDData(act_list[idx], "Known Active :", "1")
                ranking = UpdateRanking(act_list[idx], idx, simval, True, ranking, topn)
                c += 1
            c += 1
        while c < len(act_list):
            fp = act_list[c].GetData(str(fptype))
            simval = GetSimValAgainstAC(fp, act_list, baitset, fptype)
            ranking = UpdateRanking(act_list[c], idx, simval, True, ranking, topn)
            c += 1

        ranking_list.append(ranking)
    return ranking_list

def GetSimValAgainstAC(fp, act_list, baitset, fptype):
	maxval = 0
	for idx in baitset:
		fp_act = act_list[idx].GetData(str(fptype))
		tanimoto = OETanimoto(fp, fp_act)
		if tanimoto > maxval:
			maxval = tanimoto
	return maxval

def UpdateRanking(mol, idx, tanimoto, KA, ranking, topn):
	ranking.loc[len(ranking)] = [mol, idx, tanimoto, np.NaN, KA]
	ranking["Rank"] = ranking["Tanimoto"].rank(method = "min", ascending = 0)
	ranking = ranking[ranking["Rank"] < topn + 1]
	ranking = ranking.sort_values('Rank')
	ranking = ranking.reset_index()
	ranking = ranking.drop(['index'], axis = 1)

	return ranking

def RankingAnalysis(ranking, nb_ka):
    results = pd.DataFrame()
	ranking["Nb_KA"] = ranking.KA.cumsum()
	ranking["Count"] = 1
	ranking["RR"] = 100 * ranking.Nb_KA/nb_ka
	ranking["HR"] = 100 * ranking.Nb_KA/ranking.Count.cumsum()

	return ranking

def PlotResults(ranking, plot_output):

    results.plot(y = 'RR', label = "RR Set ")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("RR Rates")
    path = plot_output + "RR_plot.svg"
    plt.savefig(path)
    
    results.plot(y = 'HR', label = "HR Set ")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("HR Rates")
    path = plot_output + "HR_plot.svg"
    plt.savefig(path)
    
    results.plot(y = 'Average RR', label = "Average RR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average RR Rates")
    path = plot_output + "Average_RR_plot.svg"
    plt.savefig(path)

    results.plot(y = 'Average HR', label = "Average HR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average HR Rates")
    path = plot_output + "Average_HR_plot.svg"
    plt.savefig(path)
    
    #plt.show()

def write_output(ranking, act_list, dec_list, out, output_dir):
	ofs = oemolostream()
	output_path = out

	if not ofs.open(output_path):
		OEThrow.Warning( "Unable to create output file")

	for index, mol in ranking.iterrows():
		if mol["KA"] == 1 :
			top_mol = act_list[int(mol["idx"])]
		else:
			top_mol = dec_list[int(mol["idx"])]
		OESetSDData(top_mol, "Similarity Value (Tanimoto) :", str(mol["Tanimoto"]))
		OESetSDData(top_mol, "Trial Set :", str(mol["Set"]))
		print('%s has a similarity of %.3f' % (OEMolToSmiles(top_mol), mol["Tanimoto"]))
		OEWriteMolecule(ofs, top_mol)

	path = output_dir + "ranking.csv"
	ranking[["Set", "Molecule ID", "Tanimoto"]].to_csv(path)
	print(ranking)
	PlotResults(ranking, output_dir)

def main(argv=[__name__]):
	itf = OEInterface(InterfaceData, argv)

    ina = itf.GetString("-in_act_database")
    ind = itf.GetString("-in_decoys")
    ini = itf.GetString("-in_index_set")
    out = itf.GetString("-output")
    od = itf.GetString("-output_directory")
    topn = itf.GetInt("-topN")
    fptype = itf.GetInt("-fprint")

    print("Reading inputs")
    index_list = ReadIndex(ini)
    act_list = read_database(ina, fptype)
    
    nb_ka = len(act_list) - len(index_list[0])
    iteration = len(index_list)

    print("Ranking the Known Actives")
    ranking_list = RankActives(act_list, index_list, fptype, topn)

    print("Ranking the decoys")
    ifs = oemolistream()
    if not ifs.open(ind):
        OEThrow.Fatal("Unable to open inputfile" )

    dbfp = OEFingerPrint()
    for mol in ifs.GetOEMols():
        OEMakeFP(dbfp, mol, fptype)
        mol.SetData(str(fptype), dbfp)

        for i in range(iteration):
            simval = GetSimValAgainstAC(dbfp, act_list, index_list[i], fptype)

            OESetSDData(mol, "Similarity Value (Tanimoto) :", str(simval))
            OESetSDData(mol, "Trial Set :", str(i))
            OESetSDData(mol, "Known Active :",'0' )
            ranking_list[i] = (UpdateRanking(mol, simval, False, ranking_list[i], topn))
        
    print("Analysing")
    results = RankingAnalysis(ranking_list, nb_ka, iteration)
    print("Printing output")
    write_output(ranking_list, results, iteration, out, od)

	#ranking = pd.DataFrame(columns=["Set", "Molecule ID", "idx", "Tanimoto", "Rank", "KA", "Nb_KA", "Count", "RR", "HR"])

	#ranking = pd.concat([ranking, RankDatabase(act_list, dec_list, index_set, i, fptype, topn, nb_ka, start_time)])
        
    for idx in range(len(dec_list)):
        dbfp = dec_list[idx].GetData(str(fptype))
        mol_id = dec_list[idx].GetTitle()
        KA = 0
        simval = GetSimValAgainstAC(dbfp, act_list, index_set, fptype)
        ranking = UpdateRanking(set_nb, mol_id, idx, simval, KA, ranking, topn)

InterfaceData = """
!PARAMETER -in_act_database
	!ALIAS -ia
	!TYPE string
	!BRIEF Input Database of Active Molecules
	!REQUIRED true
	!KEYLESS 1
!END

!PARAMETER -in_decoys
	!ALIAS -id
	!TYPE string
	!BRIEF Input Database of Decoy Molecules
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


!PARAMETER -out_index_set
	!ALIAS -os
	!TYPE string
	!BRIEF Output Random Index Set Log
	!REQUIRED true
	!KEYLESS 4
!END

!PARAMETER -output_directory
	!ALIAS -od
	!TYPE string
	!BRIEF Output Directory for the plots and the ranking list
	!REQUIRED true
	!KEYLESS 5
!END

!PARAMETER -topN
	!ALIAS -t
	!TYPE int
	!BRIEF Number of top Molecules
	!REQUIRED true
	!KEYLESS 6
!END

!PARAMETER -fprint
	!ALIAS -fp
	!TYPE int
	!BRIEF Fingerprint Type (101 for MACCS, 102 for Path, 103 for Lingo, 104 for Circular, 105 for Tree)
	!REQUIRED true
	!KEYLESS 7
!END

!PARAMETER -iteration
	!ALIAS -iter
	!TYPE int
	!BRIEF Number of Iterations of the Test
	!REQUIRED true
	!KEYLESS 8
!END

!PARAMETER -ratio
	!ALIAS -r
	!TYPE int
	!BRIEF Ratio between Baits and Known Actives
	!REQUIRED false
	!DEFAULT 2
	!KEYLESS 9
!END
"""

if __name__ == "__main__":
		sys.exit(main(sys.argv))
