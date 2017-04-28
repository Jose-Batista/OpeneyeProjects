#!/usr/bin/env python
from __future__ import print_function
from openeye.oechem import *
from openeye.oegraphsim import *

import sys
import os
import random

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_database(database):
	ifs = oemolistream()

	if not ifs.open(database):
		OEThrow.Fatal("Unable to open inputfile" )

	mol_list = []
	for mol in ifs.GetOEGraphMols():
		mol_list.append(mol.CreateCopy())
	return mol_list

def RandomIndex(total, nb_index, index_log, iteration):
	i = 0
	index_set = set()
	while i < nb_index :
		index = random.randint(0, total - 1)
		if index in index_set:
			continue
		else : 
			index_set.add(index)
			i += 1

	index_log = open(index_log, "a")
	index_log.write("set N°%d: " % iteration)
	for index in index_set:
		index_log.write(str(index) + " ")
	index_log.write("\n")
	index_log.close()

	return index_set

def CalculateFP(mol_list, fptype):
	for idx in range(len(mol_list)):
		fp = OEFingerPrint()
		OEMakeFP(fp, mol_list[idx], fptype)
		mol_list[idx].SetData(str(fptype), fp)

def RankDatabase(act_list, dec_list, index_set, set_nb, fptype, topn, nb_ka):
		ranking = []
		for idx in range(len(act_list)):
			if idx not in index_set:
				dbfp = act_list[idx].GetData(str(fptype))
				mol_id = act_list[idx].GetTitle()
				KA = 1
				simval = GetSimValAgainstAC(dbfp, act_list, index_set, fptype)
				ranking = UpdateRanking(set_nb, mol_id, idx, simval, KA, ranking, topn)
		for idx in range(len(dec_list)):
			dbfp = dec_list[idx].GetData(str(fptype))
			mol_id = dec_list[idx].GetTitle()
			KA = 0
			simval = GetSimValAgainstAC(dbfp, act_list, index_set, fptype)
			ranking = UpdateRanking(set_nb, mol_id, idx, simval, KA, ranking, topn)

		ranking = RankingAnalysis(ranking, nb_ka)
		return ranking

def GetSimValAgainstAC(dbfp, act_list, index_set, fptype):
	maxval = 0
	for idx in index_set:
		fp_act = act_list[idx].GetData(str(fptype))
		if not fp_act.IsValid():
			print("fp_act uninitialized fingerprint")
		tanimoto = OETanimoto(dbfp, fp_act)
		if tanimoto > maxval:
			maxval = tanimoto
	return maxval

def UpdateRanking(set_nb, mol_id, idx, tanimoto, KA, ranking, topn):
	index = 0
	for top_mol in ranking:
		if tanimoto < top_mol[3]:
			index = ranking.index(top_mol) + 1
		else:
			break

	upper = ranking[:index]
	lower = ranking[index:]
	ranking = upper + [(set_nb, mol_id, idx, tanimoto, KA)] + lower

	i = topn - 1
	while i < len(ranking) - 1:
		if ranking[i][3] != ranking[i + 1][3]:
			ranking = ranking[:i + 1]

			break
		else:
			i += 1

	return ranking


def RankingAnalysis(ranking, nb_ka):
	data = []
	count = 0
	count_ka = 0
	for mol in ranking:
		count += 1
		if mol[4] == 1:
			count_ka += 1
		rr = 100 * count_ka/nb_ka
		hr = 100 * count_ka/count
		data.append((rr, hr))

	return data

def PlotResults(ranking, plot_output):
	ranking_by_set = ranking.groupby("Set")
	plt.figure(1)
	for Set, group in ranking_by_set:
		plt.plot(group['RR'], label = "RR Set " + str(int(Set)))
	plt.xlabel('Top Molecules')
	plt.ylabel('Rate (%)')
	plt.legend( loc='best')
	plt.title("RR Rates")
	path = plot_output + "RR_plot.svg"
	plt.savefig(path)
	
	plt.figure(2)
	for Set, group in ranking_by_set:
		plt.plot(group['HR'], label = "HR Set " + str(int(Set)))
	plt.xlabel('Top Molecules')
	plt.ylabel('Rate (%)')
	plt.legend( loc='best')
	plt.title("HR Rates")
	path = plot_output + "HR_plot.svg"
	plt.savefig(path)
	
	plt.show()
		

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
	out = itf.GetString("-output")
	outs = itf.GetString("-out_index_set")
	od = itf.GetString("-output_directory")
	topn = itf.GetInt("-topN")
	fptype = itf.GetInt("-fprint")
	iteration = itf.GetInt("-iteration")
	ratio = itf.GetInt("-ratio")

	act_list = read_database(ina)
	dec_list = read_database(ind)
	nb_act = len(act_list)
	nb_baits = nb_act//(ratio + 1)
	nb_ka = nb_act - nb_baits

	CalculateFP(act_list, fptype)
	CalculateFP(dec_list, fptype)

	i = 0
	with open(outs, 'w'): pass
	ranking = pd.DataFrame(columns=["Set", "Molecule ID", "idx", "Tanimoto", "Rank", "KA", "Nb_KA", "Count", "RR", "HR"])

	for i in range(iteration):
		print("Calculating iteration %d..." % i)
		index_set = RandomIndex(nb_act, nb_baits, outs, i)
		ranking = pd.concat([ranking, RankDatabase(act_list, dec_list, index_set, i, fptype, topn, nb_ka)])
        
       #average_result = pd.DataFrame(ranking.reset_index().groupby("index")["RR"].mean()
	write_output(ranking, act_list, dec_list, out, od)


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