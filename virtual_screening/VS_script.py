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

def read_database(database, fptype):
	ifs = oemolistream()

	if not ifs.open(database):
		OEThrow.Fatal("Unable to open inputfile" )

	mol_list = []
	for mol in ifs.GetOEMols():
		fp = OEFingerPrint()
		OEMakeFP(fp, mol, fptype)
		mol.SetData(str(fptype), fp)
		mol_list.append(mol.CreateCopy())
	return mol_list

def ReadIndex(index_input):
	index_log = open(index_input, 'r')
	index = index_log.read()
	index_list = index.split('set N°')
	index_list = index_list[1:]
	for set, random_set in enumerate(index_list):
		random_set = random_set.split(' ')
		random_set = random_set[1:-1]
		for i, idx in enumerate(random_set):
			random_set[i] = int(idx)
		index_list[set] = random_set

	index_log.close()
	return index_list


def RankDatabase(act_list, dec_database, index_list, set_nb, fptype, topn, nb_ka):
		ranking = []

		print("start")

		for idx in range(len(act_list)):
			if idx not in index_list[set_nb]:
				dbfp = act_list[idx].GetData(str(fptype))

				simval = GetSimValAgainstAC(dbfp, act_list, index_list, set_nb, fptype)

				#OESetSDData(act_list[idx], "Similarity Value (Tanimoto) :", str(simval))
				#OESetSDData(act_list[idx], "Trial Set :", str(set_nb))
				mol_id = act_list[idx].GetTitle()
				KA = 1
				#OESetSDData(act_list[idx], "Known Active :", str(KA))
				
				ranking = UpdateRanking(set_nb, mol_id, simval, KA, ranking, topn)
		
		print("start decoys")
		ifs = oemolistream()

		if not ifs.open(dec_database):
			OEThrow.Fatal("Unable to open inputfile" )

		for mol in ifs.GetOEMols():
			dbfp = OEFingerPrint()
			OEMakeFP(dbfp, mol, fptype)
			#mol.SetData(str(fptype), dbfp)
			
			simval = GetSimValAgainstAC(dbfp, act_list, index_list, set_nb, fptype)

			#OESetSDData(mol, "Similarity Value (Tanimoto) :", str(simval))
			#OESetSDData(mol, "Trial Set :", str(set_nb))
			mol_id = mol.GetTitle()
			KA = 0
			#OESetSDData(mol, "Known Active :", str(KA))

			ranking = UpdateRanking(set_nb, mol_id, simval, KA, ranking, topn)

		print("start analysis")
		data = RankingAnalysis(ranking, nb_ka)
		print ("end")
		return ranking, data

def GetSimValAgainstAC(dbfp, act_list, index_list, set_nb, fptype):
	maxval = 0
	for idx in index_list[set_nb]:
		fp_act = act_list[idx].GetData(str(fptype))
		if not fp_act.IsValid():
			print("fp_act uninitialized fingerprint")
		tanimoto = OETanimoto(dbfp, fp_act)
		if tanimoto > maxval:
			maxval = tanimoto
	return maxval

def UpdateRanking(set_nb, mol_id, tanimoto, KA, ranking, topn):
	index = len(ranking)
	for top_mol in reversed(ranking):
		if tanimoto > top_mol[2]:
			index = ranking.index(top_mol) 
		else:
			break

	upper = ranking[:index]
	lower = ranking[index:]
	ranking = upper + [(set_nb, mol_id, tanimoto, KA)] + lower

	i = topn - 1
	while i < len(ranking) - 1:
		if ranking[i][2] != ranking[i + 1][2]:
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
		if mol[3] == 1:
			count_ka += 1
		rr = 100 * count_ka/nb_ka
		hr = 100 * count_ka/count
		data.append((rr, hr))

	return data

def PlotResults(data, iteration, plot_output):

	plt.figure(1)
	for data_set in data:
		rr_set = [result[0] for result in data_set[1]]
		plt.plot(rr_set, label = "RR Set " + str(data_set[0]))
	plt.xlabel('Top Molecules')
	plt.ylabel('Rate (%)')
	plt.legend( loc='best')
	plt.title("RR Rates")
	path = plot_output + "RR_plot.svg"
	plt.savefig(path)
	
	plt.figure(2)
	for data_set in data:
		hr_set = [result[1] for result in data_set[1]]
		plt.plot(hr_set, label = "HR Set " + str(data_set[0]))
	plt.xlabel('Top Molecules')
	plt.ylabel('Rate (%)')
	plt.legend( loc='best')
	plt.title("HR Rates")
	path = plot_output + "HR_plot.svg"
	plt.savefig(path)
	
	plt.show()
		

def write_output(ranking, data, iteration, out, output_dir):
	#ofs = oemolostream()
	#output_path = out

	#if not ofs.open(output_path):
		#OEThrow.Warning( "Unable to create output file")

	#for mol in ranking:
	#	OEWriteMolecule(ofs, mol[1])

	path = output_dir + "ranking.txt"
	ranking_save = open(path, "w")
	for mol in ranking:
		mol_data = str(mol[0]) + " " + mol[1] + " " + str(mol[2])
		ranking_save.write(mol_data)
	ranking_save.close()

	PlotResults(data, iteration, output_dir)

def main(argv=[__name__]):
	itf = OEInterface(InterfaceData, argv)

	ina = itf.GetString("-in_act_database")
	ind = itf.GetString("-in_decoys")
	ini = itf.GetString("-in_index_set")
	out = itf.GetString("-output")
	od = itf.GetString("-output_directory")
	topn = itf.GetInt("-topN")
	fptype = itf.GetInt("-fprint")
	iteration = itf.GetInt("-iteration")
	ratio = itf.GetInt("-ratio")

	index_list = ReadIndex(ini)
	act_list = read_database(ina, fptype)
	nb_act = len(act_list)
	nb_baits = nb_act//(ratio + 1)
	nb_ka = nb_act - nb_baits

	ranking = []
	data = []

	for i in range(iteration):
		print("Calculating iteration %d..." % i)
	
		(new_ranking, new_data) = RankDatabase(act_list, ind, index_list, i, fptype, topn, nb_ka)
		ranking = ranking + new_ranking
		data.append((i, new_data))
        
       #Average
	write_output(ranking, data, iteration, out, od)


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

!PARAMETER -in_index_set
	!ALIAS -ii
	!TYPE string
	!BRIEF Input Random Index Set Log
	!REQUIRED true
	!KEYLESS 3
!END

!PARAMETER -output
	!ALIAS -o
	!TYPE string
	!BRIEF Output File
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