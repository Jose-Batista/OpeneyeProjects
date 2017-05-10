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

def CreateRankings(act_list, index_list, baseurl, data, topn):
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
                ranking = MergeRankings(ranking, cur_rank, topn)
        ranking_list.append(ranking)
    return ranking_list

def MergeRankings(ranking_1, ranking_2, topn):
    merged_list = list()
    i = 0
    j = 0
    count = 0
    while i < len(ranking_1):
        while j < len(ranking_2) and ranking_2[j][2] > ranking_1[i][2]:
            if count < topn or ranking_2[j][2] == merged_list[count-1][2]:
                merged_list.append(ranking_2[j])
                count += 1
                j += 1
            else:
                break
        if count < topn or ranking_2[j][2] == merged_list[count-1][2]:
            merged_list.append(ranking_1[i])  
            count += 1
            i += 1
        else:
            break

    while j < len(ranking_2):
        if count < topn or ranking_2[j][2] == merged_list[count-1][2]:
            merged_list.append(ranking_2[j])
            count += 1
            j += 1
        else:
            break

    return merged_list

baseurl = "http://130.180.63.34:8089"
response = requests.get( baseurl )
data = response.json()

topn = 100

index_list = [[0,1,2], [0,3,5], [0,2,4]]

act_list = list()
mol = OEMol()
OESmilesToMol(mol, "c1cccc1c")
act_list.append(mol.CreateCopy())

mol = OEMol()
OESmilesToMol(mol, "C1CCCC1C")
act_list.append(mol.CreateCopy())

mol = OEMol()
OESmilesToMol(mol, "O=Cc1ccc(O)c(OC)c1")
act_list.append(mol.CreateCopy())

mol = OEMol()
OESmilesToMol(mol, "NC(Cl)(Br)C(=O)O")
act_list.append(mol.CreateCopy())

mol = OEMol()
OESmilesToMol(mol, "CN1CCC[C@H]1c2cccnc2")
act_list.append(mol.CreateCopy())

mol = OEMol()
OESmilesToMol(mol, "N[C@@H](C)C(=O)O")
act_list.append(mol.CreateCopy())

for mol in act_list:
    print('Mol : ', OEMolToSmiles(mol))

ranking_list = CreateRankings(act_list, index_list, baseurl, data, topn)

for ranking in ranking_list:
    print('Ranking')
    for mol in ranking:
        print('Molecule : %-20s Simval : %.3f   KA : %d' % (mol[1], mol[2], mol[3]))

