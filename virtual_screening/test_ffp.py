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
import urllib.parse as parse

def CreateRankings(act_list, index_list, baseurl, data, topn):
    ranking_list = list()
    for baitset in index_list:
        ranking = list()
        for idx in baitset:
            print(idx)
            smiles = OEMolToSmiles(act_list[idx])
            #safe_smiles = parse.quote(smiles, safe='~@#$&()*!+=:;,.?/\'')
            safe_smiles = parse.quote(smiles)
            url = "%s/%s/hitlist?smiles=%s&oformat=csv" %(baseurl, data['databases'][0], safe_smiles) 
            print(url)
            response = requests.get( url )
            hitlist = response.content.decode().split('\n')
            print(type(hitlist))
            hitlist.pop(0)
            hitlist.pop()
            cur_rank = list()
            for mol in hitlist:
                cur_mol = mol.split(',')
                cur_rank.append((cur_mol[0], cur_mol[1], float(cur_mol[5]), False))
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
    id_set = set()
    while i < len(ranking_1):
        while j < len(ranking_2) and ranking_2[j][2] > ranking_1[i][2]:
            if ranking_2[j][1] not in id_set: 
                if count < topn or ranking_2[j][2] == merged_list[count-1][2]:
                    merged_list.append(ranking_2[j])
                    count += 1
                    id_set.add(ranking_2[j][1])
                    j += 1
                else:
                    break
            else:
                j += 1

        if ranking_1[i][1] not in id_set: 
            if ranking_1[i] not in id_set and (count < topn or ranking_1[i][2] == merged_list[count-1][2]):
                merged_list.append(ranking_1[i])  
                count += 1
                id_set.add(ranking_1[i][1])
                i += 1
            else:
                break
        else:
            i += 1

    while j < len(ranking_2):
        if ranking_2[j][1] not in id_set: 
            if ranking_2[j] not in id_set and (count < topn or ranking_2[j][2] == merged_list[count-1][2]):
                merged_list.append(ranking_2[j])
                count += 1
                id_set.add(ranking_2[j][1])
                j += 1
            else:
                break
        else:
            i += 1

    return merged_list

baseurl = "http://10.0.1.22:8069"
response = requests.get( baseurl )
data = response.json()

topn = 100

index_list = [[0,1,2], [0,3,5], [0,2,4]]

act_list = list()
mol = OEMol()
OESmilesToMol(mol, "N#N")
act_list.append(mol.CreateCopy())

mol = OEMol()
OESmilesToMol(mol, "[OH3+]")
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

