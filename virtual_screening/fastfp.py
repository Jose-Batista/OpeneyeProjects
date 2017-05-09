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

import json,requests

def read_database(database, fptype):
    ifs = oemolistream()

    if not ifs.open(database):
        OEThrow.Fatal("Unable to open inputfile" )

    mol_list = list()
    for mol in ifs.GetOEMols():
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


def CreateRankings(act_list, index_list, baseurl, data):
    ranking_list = list()
    for baitset in index_list:
        ranking = list()
        for idx in baitset:
            print(idx)
            url = "%s/%s/hitlist?smiles=%s&oformat=csv" %(baseurl, data['databases'][0], OEMolToSmiles(act_list[idx]))
            response = requests.get( url )
            print(response.content)
            hitlist = response.content.decode().split('\n')
            hitlist = hitlist[1:-1]
            cur_rank = list()
            for mol in hitlist:
                cur_mol = mol.split(',')
                print(cur_mol[0], cur_mol[1], cur_mol[2])
                cur_rank.append((cur_mol[0], cur_mol[1], float(cur_mol[2]), False))
            if len(ranking) == 0:
                ranking = cur_rank
            else:
                ranking = MergeRankings(ranking, cur_rank)
        ranking_list.append(ranking)
    return ranking_list

def MergeRankings(ranking_1, ranking_2):
    merged_list = list()
    i = 0
    j = 0
    while i < len(ranking_1):
        while j < len(ranking_2) and ranking_2[j][2] > ranking_1[i][2]:
            merged_list.append(ranking_2[j])
            j += 1
        merged_list.append(ranking_1[i])  
        i += 1

    while j < len(ranking_2):
        merged_list.append(ranking_2[j])
        j += 1

    return merged_list

def InsertKnownActives(ranking_list, act_list, index_list, baseurl, data):
    for i, baitset in enumerate(index_list):
        c = 0
        for idx in baitset:
            while c < idx:
                url = "%s/%s/neighbor?smiles=%s&oformat=csv" %(baseurl, data['databases'][0], OEMolToSmiles(act_list[c]))
                response = requests.get( url )
                print(response.content)
                neighbor = response.content.decode().split('\n')
                neighbor = neighbor[1:-1]
                neighbor = neighbor[0].split(',')
                known_act = list((OEMolToSmiles(act_list[c]), OEMol.GetTitle(), float(neighbor[2]), True )
                ranking_list[i] = MergeRankings(ranking_list[i], known_act)
                c += 1
            c += 1
        while c < len(act_list):
            url = "%s/%s/neighbor?smiles=%s&oformat=csv" %(baseurl, data['databases'][0], OEMolToSmiles(act_list[c]))
            response = requests.get( url )
            print(response.content)
            neighbor = response.content.decode().split('\n')
            neighbor = neighbor[1:-1]
            neighbor = neighbor[0].split(',')
            known_act = list((OEMolToSmiles(act_list[c]), OEMol.GetTitle(), float(neighbor[2]), True )
            ranking_list[i] = MergeRankings(ranking_list[i], known_act)
            c += 1

    return ranking_list

def RankingAnalysis(ranking_list, nb_ka):
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
    results_avg['Average RR'] = results.groupby(results.index)['RR'].mean()
    results_avg['Average HR'] = results.groupby(results.index)['HR'].mean()

    return results_avg

def PlotResults(results_avg, plot_output, fptype):

    results_avg.plot(y = 'Average RR', label = "Average RR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average RR Rates FP" + str(fptype))
    path = plot_output + "Average_RR_plot_FP" + str(fptype) + ".svg"
    plt.savefig(path)

    results_avg.plot(y = 'Average HR', label = "Average HR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average HR Rates FP" + str(fptype))
    path = plot_output + "Average_HR_plot_FP" + str(fptype) + ".svg"
    plt.savefig(path)
    
    #plt.show()


def write_output(ranking_list, results_avg, fptype, output_dir):
    ofs = oemolostream()
    path = output_dir + "ranking_FP" + str(fptype) + ".oeb"

    if not ofs.open(path):
        OEThrow.Warning( "Unable to create output file")

    for ranking in ranking_list:
        for mol in ranking:
            OEWriteMolecule(ofs, mol[0])

    path = output_dir + "ranking_FP" + str(fptype) + ".txt"
    ranking_save = open(path, "w")
    for i, ranking in enumerate(ranking_list):
        for mol in ranking:
            mol_data = str(i) + " " + mol[0].GetTitle() + " " + str(mol[1])
            ranking_save.write(mol_data)
    ranking_save.close()

    path = output_dir + "results_FP" + str(fptype) + ".csv"
    results_avg.to_csv(path)
    
    PlotResults(results_avg, output_dir, fptype)
    
def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    ina = itf.GetString("-in_act_database")
    ini = itf.GetString("-in_index_set")
    od = itf.GetString("-output_directory")
    topn = itf.GetInt("-topN")
    fptype = itf.GetInt("-fprint")

    baseurl = "http://130.180.63.34:8089"
    response = requests.get( baseurl )
    data = response.json()

    print("Reading inputs")
    index_list = ReadIndex(ini)
    act_list = read_database(ina)
    
    nb_ka = len(act_list) - len(index_list[0])

    ranking_list = CreateRanking(act_list, index_list, baseurl, data)
    ranking_list = InsertKnownActives(ranking_list, act_list, index_list, baseurl)

    print("Analysing")
    results_avg = RankingAnalysis(ranking_list, nb_ka)
    print("Printing output")
    write_output(ranking_list, results_avg, fptype, od)


InterfaceData = """
!PARAMETER -in_act_database
    !ALIAS -ia
    !TYPE string
    !BRIEF Input Database of Active Molecules
    !REQUIRED true
    !KEYLESS 1
!END

!PARAMETER -in_index_set
    !ALIAS -ii
    !TYPE string
    !BRIEF Input Random Index Set Log
    !REQUIRED true
    !KEYLESS 2
!END

!PARAMETER -output_directory
    !ALIAS -od
    !TYPE string
    !BRIEF Output Directory for the plots and the ranking list
    !REQUIRED true
    !KEYLESS 3
!END

!PARAMETER -topN
    !ALIAS -t
    !TYPE int
    !BRIEF Number of top Molecules
    !REQUIRED true
    !KEYLESS 4
!END

!PARAMETER -fprint
    !ALIAS -fp
    !TYPE int
    !BRIEF Fingerprint Type (101 for MACCS, 102 for Path, 103 for Lingo, 104 for Circular, 105 for Tree)
    !REQUIRED true
    !KEYLESS 5
!END

"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
