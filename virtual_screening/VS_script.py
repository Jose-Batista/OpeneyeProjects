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
    for set_id, random_set in enumerate(index_list):
        random_set = random_set.split(' ')
        random_set = random_set[1:-1]
        for i, idx in enumerate(random_set):
            random_set[i] = int(idx)
        index_list[set_id] = random_set

    index_log.close()
    return index_list


def RankActives(act_list, index_list, fptype, topn):

    ranking_list = list()

    for i, baitset in enumerate(index_list):
        ranking = list()
        c = 0
        for idx in baitset:
            while c < idx:
                fp = act_list[c].GetData(str(fptype))
                simval = GetSimValAgainstAC(fp, act_list, baitset, fptype)
                OESetSDData(act_list[idx], "Similarity Value (Tanimoto) :", str(simval))
                OESetSDData(act_list[idx], "Trial Set :", str(i))
                OESetSDData(act_list[idx], "Known Active :", "1")
                ranking = UpdateRanking(act_list[c], simval, True, ranking, topn)
                c += 1
            c += 1
        while c < len(act_list):
            fp = act_list[c].GetData(str(fptype))
            simval = GetSimValAgainstAC(fp, act_list, baitset, fptype)
            ranking = UpdateRanking(act_list[c], simval, True, ranking, topn)
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

def UpdateRanking(mol, tanimoto, KA, ranking, topn):
    index = len(ranking)
    if (index >= topn and tanimoto < ranking[index-1][1]):
        return ranking
    else:    
        for top_mol in reversed(ranking):
            if tanimoto > top_mol[1]:
                index = ranking.index(top_mol) 
            else:
                break

        upper = ranking[:index]
        lower = ranking[index:]
        ranking = upper + [(mol.CreateCopy(), tanimoto, KA)] + lower

        i = topn - 1
        while i < len(ranking) - 1:
            if ranking[i][1] != ranking[i + 1][1]:
                ranking = ranking[:i + 1]

                break
            else:
                i += 1

        return ranking

def RankingAnalysis(ranking_list, nb_ka):
    results = pd.DataFrame()
    for ranking in ranking_list:
        set_results = pd.DataFrame(columns = ['RR', 'HR', 'Set'])
        count = 0
        count_ka = 0
        for row, mol in enumerate(ranking):
            count += 1
            if mol[2] == 1:
                count_ka += 1
            rr = 100 * count_ka/nb_ka
            hr = 100 * count_ka/count
            set_results.loc[row] = [rr, hr, i]
        results = pd.concat([results, set_results])

    results['Average RR'] = results['RR'].mean(axis=1)
    results['Average HR'] = results['HR'].mean(axis=1)
    #print(results)
    return results

def PlotResults(results, plot_output):

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


def write_output(ranking_list, results, iteration, out, output_dir):
    ofs = oemolostream()
    output_path = out

    if not ofs.open(output_path):
        OEThrow.Warning( "Unable to create output file")

    for ranking in ranking_list:
        for mol in ranking:
            OEWriteMolecule(ofs, mol[0])

    path = output_dir + "ranking.txt"
    ranking_save = open(path, "w")
    for i, ranking in enumerate(ranking_list):
        for mol in ranking:
            mol_data = str(i) + " " + mol[0].GetTitle() + " " + str(mol[1])
            ranking_save.write(mol_data)
    ranking_save.close()

    PlotResults(results, output_dir)

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

"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
