#!/usr/bin/env python
from __future__ import print_function
from openeye.oechem import *

import sys
import os

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import ttest_rel
from scipy.stats import bayes_mvs

def TTest(results, subdirs, folders, output):
    sample_100 = pd.DataFrame(columns=folders)
    sample_500 = pd.DataFrame(columns=folders)
    for i, AC in enumerate(subdirs):
        tree_1024 = [col for col in results.columns if ('RR' in col and AC in col and 'tree' in col and '1024' in col)]
        tree_2048 = [col for col in results.columns if ('RR' in col and AC in col and 'tree' in col and '2048' in col)]
        tree_4096 = [col for col in results.columns if ('RR' in col and AC in col and 'tree' in col and '4096' in col)]
        tree_8192 = [col for col in results.columns if ('RR' in col and AC in col and 'tree' in col and '8192' in col)]
        sample_100.loc[i]=[results[tree_1024[0]].iloc[99], results[tree_2048[0]].iloc[99], results[tree_4096[0]].iloc[99], results[tree_8192[0]].iloc[99]]#, results[rocs[0]].iloc[99]]
        sample_500.loc[i]=[results[tree_1024[0]].iloc[499], results[tree_2048[0]].iloc[499], results[tree_4096[0]].iloc[499], results[tree_8192[0]].iloc[499]]#, results[rocs[0]].iloc[499]]

    sample_100['1024_2048']=sample_100['1024']-sample_100['2048']
    sample_100['1024_4096']=sample_100['1024']-sample_100['4096']
    sample_100['1024_8192']=sample_100['1024']-sample_100['8192']
    sample_100['2048_4096']=sample_100['2048']-sample_100['4096']
    sample_100['2048_8192']=sample_100['2048']-sample_100['8192']
    sample_100['8192_4096']=sample_100['8192']-sample_100['4096']

    sample_500['1024_2048']=sample_500['1024']-sample_500['2048']
    sample_500['1024_4096']=sample_500['1024']-sample_500['4096']
    sample_500['1024_8192']=sample_500['1024']-sample_500['8192']
    sample_500['2048_4096']=sample_500['2048']-sample_500['4096']
    sample_500['2048_8192']=sample_500['2048']-sample_500['8192']
    sample_500['8192_4096']=sample_500['8192']-sample_500['4096']

    print(sample_100)
    print(sample_500)

    path = output + "FPsize_Ttest.txt"
    stream = open(path, 'w')

    print("Recovery Rate at TOP100")
    tree_1024_tree_2048_100 = ttest_rel(sample_100['1024'], sample_100['2048'])
    print(' T-Test results on RR for tree_1024 against tree_2048 for top100 : ', tree_1024_tree_2048_100)
    stream.write(' T-Test results on RR for tree_1024 against tree_2048 for top100 : ' + str(tree_1024_tree_2048_100) + '\n')
    mean, var, std = bayes_mvs(sample_100['1024_2048'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_1024_tree_4096_100 = ttest_rel(sample_100['1024'], sample_100['4096'])
    print(' T-Test results on RR for tree_1024 against tree_4096 for top100 : ', tree_1024_tree_4096_100)
    stream.write(' T-Test results on RR for tree_1024 against tree_4096 for top100 : ' + str(tree_1024_tree_4096_100) + '\n')
    mean, var, std = bayes_mvs(sample_100['1024_4096'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_1024_tree_8192_100 = ttest_rel(sample_100['1024'], sample_100['8192'])
    print(' T-Test results on RR for tree_1024 against tree_8192 for top100 : ', tree_1024_tree_8192_100)
    stream.write(' T-Test results on RR for tree_1024 against tree_8192 for top100 : ' + str(tree_1024_tree_8192_100) + '\n')
    mean, var, std = bayes_mvs(sample_100['1024_8192'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_2048_tree_4096_100 = ttest_rel(sample_100['2048'], sample_100['4096'])
    print(' T-Test results on RR for tree_2048 against tree_4096 for top100 : ', tree_2048_tree_4096_100)
    stream.write(' T-Test results on RR for tree_2048 against tree_4096 for top100 : ' + str(tree_2048_tree_4096_100) + '\n')
    mean, var, std = bayes_mvs(sample_100['2048_4096'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_2048_tree_8192_100 = ttest_rel(sample_100['2048'], sample_100['8192'])
    print(' T-Test results on RR for tree_2048 against tree_8192 for top100 : ', tree_2048_tree_8192_100)
    stream.write(' T-Test results on RR for tree_2048 against tree_8192 for top100 : ' + str(tree_2048_tree_8192_100) + '\n')
    mean, var, std = bayes_mvs(sample_100['2048_8192'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_8192_tree_4096_100 = ttest_rel(sample_100['8192'], sample_100['4096'])
    print(' T-Test results on RR for tree_8192 against tree_4096 for top100 : ', tree_8192_tree_4096_100)
    stream.write(' T-Test results on RR for tree_8192 against tree_4096 for top100 : ' + str(tree_8192_tree_4096_100) + '\n')
    mean, var, std = bayes_mvs(sample_100['8192_4096'])
    print(mean)
    stream.write(str(mean) + '\n')

    print("\nRecovery Rate at TOP500")
    tree_1024_tree_2048_500 = ttest_rel(sample_500['1024'], sample_500['2048'])
    print(' \n T-Test results on RR for tree_1024 against tree_2048 for top500 : ', tree_1024_tree_2048_500)
    stream.write(' T-Test results on RR for tree_1024 against tree_2048 for top500 : ' + str(tree_1024_tree_2048_500) + '\n')
    mean, var, std = bayes_mvs(sample_500['1024_2048'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_1024_tree_4096_500 = ttest_rel(sample_500['1024'], sample_500['4096'])
    print(' T-Test results on RR for tree_1024 against tree_4096 for top500 : ', tree_1024_tree_4096_500)
    stream.write(' T-Test results on RR for tree_1024 against tree_4096 for top500 : ' + str(tree_1024_tree_4096_500) + '\n')
    mean, var, std = bayes_mvs(sample_500['1024_4096'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_1024_tree_8192_500 = ttest_rel(sample_500['1024'], sample_500['8192'])
    print(' T-Test results on RR for tree_1024 against tree_8192 for top500 : ', tree_1024_tree_8192_500)
    stream.write(' T-Test results on RR for tree_1024 against tree_8192 for top500 : ' + str(tree_1024_tree_8192_500) + '\n')
    mean, var, std = bayes_mvs(sample_500['1024_8192'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_2048_tree_4096_500 = ttest_rel(sample_500['2048'], sample_500['4096'])
    print(' T-Test results on RR for tree_2048 against tree_4096 for top500 : ', tree_2048_tree_4096_500)
    stream.write(' T-Test results on RR for tree_2048 against tree_4096 for top500 : ' + str(tree_2048_tree_4096_500) + '\n')
    mean, var, std = bayes_mvs(sample_500['2048_4096'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_2048_tree_8192_500 = ttest_rel(sample_500['2048'], sample_500['8192'])
    print(' T-Test results on RR for tree_2048 against tree_8192 for top500 : ', tree_2048_tree_8192_500)
    stream.write(' T-Test results on RR for tree_2048 against tree_8192 for top500 : ' + str(tree_2048_tree_8192_500) + '\n')
    mean, var, std = bayes_mvs(sample_500['2048_8192'])
    print(mean)
    stream.write(str(mean) + '\n')

    tree_8192_tree_4096_500 = ttest_rel(sample_500['8192'], sample_500['4096'])
    print(' T-Test results on RR for tree_8192 against tree_4096 for top500 : ', tree_8192_tree_4096_500)
    stream.write(' T-Test results on RR for tree_8192 against tree_4096 for top500 : ' + str(tree_8192_tree_4096_500) + '\n')
    mean, var, std = bayes_mvs(sample_500['8192_4096'])
    print(mean)
    stream.write(str(mean) + '\n')

    stream.close()


def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    results = pd.DataFrame()

    path = itf.GetString("-in_dir")
    folders = ['1024', '2048', '4096', '8192']

    for folder in folders:
        subdirs = [x for x in os.listdir( path ) if (os.path.isdir(path + x) and os.path.isdir(path + x + '/' + folder + '/'))]

        for AC in subdirs:
            result_files = [x for x in os.listdir( path + AC + '/' + folder + '/') if x.endswith('.csv')]
            for result_file in result_files:
                result = pd.read_csv(path + AC + '/' + folder + '/' + result_file)
                result.rename(columns=lambda x: x + " " + AC + " " + folder, inplace=True)
                results = pd.concat([results, result], axis = 1)

    output = itf.GetString("-output")
    #PlotResults(results, subdirs, folder, output)
    #CompareFPResults(results, subdirs, folder, output)
    TTest(results, subdirs, folders, output)

InterfaceData = """

!PARAMETER -output
    !ALIAS -o
    !TYPE string
    !BRIEF Output Directory
    !REQUIRED true
    !KEYLESS 1
!END

!PARAMETER -in_dir
    !ALIAS -i
    !TYPE string
    !BRIEF Input Directory of Results
    !REQUIRED true
    !KEYLESS 2
!END

"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
