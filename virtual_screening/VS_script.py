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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_database(database, fptype):
    ifs = oemolistream()

    if not ifs.open(database):
        OEThrow.Fatal("Unable to open inputfile" )

    mol_list = list()
    fp_list = list()
    for mol in ifs.GetOEMols():
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, fptype)
        fp_list.append(fp)
        mol.SetData(str(fptype), fp)
        mol_list.append(mol.CreateCopy())
    return mol_list, fp_list

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


def RankActives(act_list, index_list, fp_list, fptype, topn):

    ranking_list = list()

    for i, baitset in enumerate(index_list):
        ranking = list()
        c = 0
        for idx in baitset:
            while c < idx:
                fp = act_list[c].GetData(str(fptype))
                simval = GetSimValAgainstAC(fp, fp_list, baitset, fptype)
                OESetSDData(act_list[idx], "Similarity Value (Tanimoto) :", str(simval))
                OESetSDData(act_list[idx], "Trial Set :", str(i))
                OESetSDData(act_list[idx], "Known Active :", "1")
                ranking = UpdateRanking(act_list[c], simval, True, ranking, topn)
                c += 1
            c += 1
        while c < len(act_list):
            fp = act_list[c].GetData(str(fptype))
            simval = GetSimValAgainstAC(fp, fp_list, baitset, fptype)
            ranking = UpdateRanking(act_list[c], simval, True, ranking, topn)
            c += 1

        ranking_list.append(ranking)


    return ranking_list

def GetSimValAgainstAC(fp, fp_list, baitset, fptype):
    maxval = 0
    for idx in baitset:
        tanimoto = OETanimoto(fp, fp_list[idx])
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

def RankingAnalysis(ranking_list, nb_ka, topn, FPType):
    results = pd.DataFrame()
    for i, ranking in enumerate(ranking_list):
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
    
    results_avg = pd.DataFrame()
    results_avg['Average RR ' + FPType] = results.groupby(results.index)['RR'].mean()
    results_avg['Average HR ' + FPType] = results.groupby(results.index)['HR'].mean()
    results_avg = results_avg.head(topn)

    return results_avg

def PlotResults(results_avg, plot_output, FPType):

    results_avg.plot(y = 'Average RR', label = "Average RR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average RR Rates " + FPType)
    path = plot_output + "Average_RR_plot_" + FPType + ".svg"
    plt.savefig(path)

    results_avg.plot(y = 'Average HR', label = "Average HR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average HR Rates " + FPType)
    path = plot_output + "Average_HR_plot_" + FPType + ".svg"
    plt.savefig(path)
    
    #plt.show()


def write_output(ranking_list, results_avg, FPType, output_dir):
    ofs = oemolostream()
    path = output_dir + "ranking_" + FPType + ".oeb"

    if not ofs.open(path):
        OEThrow.Warning( "Unable to create output file")

    for ranking in ranking_list:
        for mol in ranking:
            OEWriteMolecule(ofs, mol[0])

    path = output_dir + "ranking_" + FPType + ".txt"
    ranking_save = open(path, "w")
    for i, ranking in enumerate(ranking_list):
        for mol in ranking:
            mol_data = str(i) + " " + mol[0].GetTitle() + " " + str(mol[1]) + " " + str(mol[2])
            ranking_save.write(mol_data)
    ranking_save.close()

    path = output_dir + "results_" + FPType + ".csv"
    results_avg.to_csv(path)
    
    PlotResults(results_avg, output_dir, FPType)

def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    ina = itf.GetString("-in_act_database")
    ind = itf.GetString("-in_decoys")
    ini = itf.GetString("-in_index_set")
    od = itf.GetString("-output_directory")
    topn = itf.GetInt("-topN")
    fptype = itf.GetInt("-fprint")

    start_time = time.time()

    print("Reading inputs")
    index_list = ReadIndex(ini)
    (act_list, fp_list) = read_database(ina, fptype)
    
    nb_ka = len(act_list) - len(index_list[0])
    iteration = len(index_list)

    print("Ranking the Known Actives")
    ranking_list = RankActives(act_list, index_list, fp_list, fptype, topn)

    print("Ranking the decoys")
    ifs = oemolistream()
    if not ifs.open(ind):
        OEThrow.Fatal("Unable to open inputfile" )

    dbfp = OEFingerPrint()
    count = 0
    for mol in ifs.GetOEMols():
        count += 1
        #if count < 10:
            #print("Before FP : ", time.time() - start_time)
        OEMakeFP(dbfp, mol, fptype)
        mol.SetData(str(fptype), dbfp)
        #if count < 10:
            #print("After FP : ", time.time() - start_time)

        for i in range(iteration):
            #if count < 10:
                #print("Before SimVal : ", time.time() - start_time)
            simval = GetSimValAgainstAC(dbfp, fp_list, index_list[i], fptype)
            #if count < 10:
                #print("After SimVal : ", time.time() - start_time)

            OESetSDData(mol, "Similarity Value (Tanimoto) :", str(simval))
            OESetSDData(mol, "Trial Set :", str(i))
            OESetSDData(mol, "Known Active :",'0' )
            #if count < 10:
                #print("Before Ranking : ", time.time() - start_time)
            ranking_list[i] = (UpdateRanking(mol, simval, False, ranking_list[i], topn))
        
    print("Analysing")
    fptypes = {101 : 'MACCS', 102 : 'path', 104 : 'circular', 105 : 'tree'}
    FPType = fptypes[fptype]
    results_avg = RankingAnalysis(ranking_list, nb_ka, topn, FPType)
    print("Printing output")
    write_output(ranking_list, results_avg, FPType, od)


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

!PARAMETER -output_directory
    !ALIAS -od
    !TYPE string
    !BRIEF Output Directory for the plots and the ranking list
    !REQUIRED true
    !KEYLESS 4
!END

!PARAMETER -topN
    !ALIAS -t
    !TYPE int
    !BRIEF Number of top Molecules
    !REQUIRED true
    !KEYLESS 5
!END

!PARAMETER -fprint
    !ALIAS -fp
    !TYPE int
    !BRIEF Fingerprint Type (101 for MACCS, 102 for Path, 103 for Lingo, 104 for Circular, 105 for Tree)
    !REQUIRED true
    !KEYLESS 6
!END

"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
