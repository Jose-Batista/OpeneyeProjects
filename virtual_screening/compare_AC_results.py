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

def PlotResults(results, subdirs, folder, plot_output):
    AC_count = 0
    AC_left = len(subdirs)

    path = plot_output + folder + "_AC_comparison.pdf"
    pdf_file = PdfPages(path)

    fig, axes = plt.subplots(nrows=min(3, AC_left), ncols=2, figsize=(7,9))
    plt.subplots_adjust(hspace=0.4)

    for AC in subdirs:
        cols = [col for col in results.columns if ('RR' in col and AC in col)]
        for col in cols:
            if 'tree' in col:
                tree = col
            if 'path' in col:
                path = col
            if 'circular' in col:
                circular = col

        results.plot(ax=axes[AC_count, 0], y = [col for col in results.columns if ('RR' in col and AC in col)], sharex = True, sharey = True) 
        axes[AC_count, 0].set_xlabel('Top Ranked Molecules')
        axes[AC_count, 0].set_ylabel('Rate (%)')
        axes[AC_count, 0].set_title("Avg Recovery Rates " + AC) 
        for line in axes[AC_count, 0].get_lines():
            if 'tree' in line.get_label():
                line.set_color("C0")
            if 'path' in line.get_label():
                line.set_color("C1")
            if 'circular' in line.get_label():
                line.set_color("C2")
            if 'FR' in line.get_label():
                line.set_color("red")
        axes[AC_count, 0].legend([col.split(" ")[2] for col in results.columns if ('RR' in col and AC in col)])

        results.plot(ax=axes[AC_count, 1], y = [col for col in results.columns if ('HR' in col and AC in col)], sharex = True, sharey = True)
        axes[AC_count, 1].set_xlabel('Top Ranked Molecules')
        axes[AC_count, 1].set_title("Avg Hit Rates " + AC) 
        for line in axes[AC_count, 1].get_lines():
            if 'tree' in line.get_label():
                line.set_color("C0")
            if 'path' in line.get_label():
                line.set_color("C1")
            if 'circular' in line.get_label():
                line.set_color("C2")
            if 'FR' in line.get_label():
                line.set_color("red")
        axes[AC_count, 1].legend([col.split(" ")[2] for col in results.columns if ('HR' in col and AC in col)])

        AC_count += 1
        if AC_count == 3 and AC_left !=3: 
            pdf_file.savefig(fig)
            AC_count = 0
            AC_left -= 3
            fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(7, 9))
            plt.subplots_adjust(hspace=0.4)

    while AC_count < 3:
        axes[AC_count, 0].axis('off')
        axes[AC_count, 1].axis('off')
        AC_count += 1

    pdf_file.savefig(fig)
    pdf_file.close()

def CompareFPResults(results, subdirs, folder, plot_output):
    AC_count = 0
    AC_left = len(subdirs)

    path = plot_output + folder + "_FP_comparison.pdf"
    pdf_file = PdfPages(path)

    fig, axes = plt.subplots(nrows=min(3, AC_left), ncols=2, figsize=(7,9))
    plt.subplots_adjust(hspace=0.4)

    compare = pd.DataFrame()
    compare = results.iloc[[49, 99, 149, 499]]
    compare['Index'] = [50, 100, 150, 500]

    for AC in subdirs:
        cols = [col for col in results.columns if ('RR' in col and AC in col)]
        bars = compare.plot(kind='bar', ax=axes[AC_count, 0], x = compare['Index'], y = [col for col in results.columns if ('RR' in col and AC in col)], sharex = True, sharey = True)

        #for col in cols:
            #bar = compare.plot(kind='bar', ax=axes[AC_count, 0], y = col, label = 'tree', sharex = True, sharey = True) 
            #bar.set_label(col.split(" ")[2])
        axes[AC_count, 0].set_xlabel('Top Ranked Molecules')
        axes[AC_count, 0].set_ylabel('Rate (%)')
        axes[AC_count, 0].set_title("Avg Recovery Rates " + AC) 
        axes[AC_count, 0].legend([col.split(" ")[2] for col in results.columns if ('RR' in col and AC in col)])
        barlist=filter(lambda x: isinstance(x, matplotlib.patches.Rectangle), axes[AC_count, 0].get_children())
        for bar in barlist:
            print(bar.get_label())
            if 'tree' in bar.get_label():
                bar.set_color("cyan")
            if 'path' in bar.get_label():
                bar.set_color("orange")
            if 'circular' in bar.get_label():
                bar.set_color("green")
            if 'FR' in bar.get_label():
                bar.set_color("red")
        axes[AC_count, 0].legend([col.split(" ")[2] for col in results.columns if ('RR' in col and AC in col)])

        compare.plot(kind='bar', ax=axes[AC_count, 1], y = [col for col in results.columns if ('HR' in col and AC in col)], sharex = True, sharey = True)
        axes[AC_count, 1].set_xlabel('Top Ranked Molecules')
        axes[AC_count, 1].set_title("Avg Hit Rates " + AC) 
        barlist=filter(lambda x: isinstance(x, matplotlib.patches.Rectangle), axes[AC_count, 1].get_children())
        for bar in barlist:
            print(bar)
            if 'tree' in bar.get_label():
                line.set_color("cyan")
            if 'path' in bar.get_label():
                line.set_color("orange")
            if 'circular' in bar.get_label():
                line.set_color("green")
            if 'FR' in bar.get_label():
                line.set_color("red")
        axes[AC_count, 1].legend([col.split(" ")[2] for col in results.columns if ('HR' in col and AC in col)])

        AC_count += 1
        if AC_count == 3 and AC_left !=3: 
            pdf_file.savefig(fig)
            AC_count = 0
            AC_left -= 3
            fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(7, 9))
            plt.subplots_adjust(hspace=0.4)

    while AC_count < 3:
        axes[AC_count, 0].axis('off')
        axes[AC_count, 1].axis('off')
        AC_count += 1

    pdf_file.savefig(fig)
    pdf_file.close()

def TTest(results, subdirs, folder, output):
    columns=['tree', 'path', 'circular']
    #columns=['tree', 'path', 'circular', 'rocs']
    sample_100 = pd.DataFrame(columns=columns)
    sample_500 = pd.DataFrame(columns=columns)
    for i, AC in enumerate(subdirs):
        tree = [col for col in results.columns if ('RR' in col and AC in col and 'tree' in col)]
        path = [col for col in results.columns if ('RR' in col and AC in col and 'path' in col)]
        circular = [col for col in results.columns if ('RR' in col and AC in col and 'circular' in col)]
        rocs = [col for col in results.columns if ('RR' in col and AC in col and 'FR' in col)]
        sample_100.loc[i]=[results[tree[0]].iloc[99], results[path[0]].iloc[99], results[circular[0]].iloc[99]]#, results[rocs[0]].iloc[99]]
        sample_500.loc[i]=[results[tree[0]].iloc[499], results[path[0]].iloc[499], results[circular[0]].iloc[499]]#, results[rocs[0]].iloc[499]]

    #sample_500['tree_rocs']=sample_500['tree']-sample_500['rocs']

    path = output + folder + "Ttest.txt"
    stream = open(path, 'w')

    print(sample_100)
    print(sample_500)
    tree_path_100 = ttest_rel(sample_100['tree'], sample_100['path'])
    tree_circular_100 = ttest_rel(sample_100['tree'], sample_100['circular'])
    circular_path_100 = ttest_rel(sample_100['circular'], sample_100['path'])
    #tree_rocs_100 = ttest_rel(sample_100['tree'], sample_100['rocs'])
    #path_rocs_100 = ttest_rel(sample_100['path'], sample_100['rocs'])
    #circular_rocs_100 = ttest_rel(sample_100['circular'], sample_100['rocs'])
    print(' T-Test results on RR for tree against path for top100 : ', tree_path_100)
    print(' T-Test results on RR for tree against circular for top100 : ', tree_circular_100)
    print(' T-Test results on RR for circular against path for top100 : ', circular_path_100)
    #print(' T-Test results on RR for tree against rocs for top100 : ', tree_rocs_100)
    #print(' T-Test results on RR for path against rocs for top100 : ', path_rocs_100)
    #print(' T-Test results on RR for circular against rocs for top100 : ', circular_rocs_100)
    stream.write(' T-Test results on RR for tree against path for top100 : ' + str(tree_path_100) + '\n')
    stream.write(' T-Test results on RR for tree against circular for top100 : ' + str(tree_circular_100) + '\n')
    stream.write(' T-Test results on RR for circular against path for top100 : ' + str(circular_path_100) + '\n')
    #stream.write(' T-Test results on RR for tree against rocs for top100 : ' + str(tree_rocs_100) + '\n')
    #stream.write(' T-Test results on RR for path against rocs for top100 : ' + str(path_rocs_100) + '\n')
    #stream.write(' T-Test results on RR for circular against rocs for top100 : ' + str(circular_rocs_100) + '\n')

    tree_path_500 = ttest_rel(sample_500['tree'], sample_500['path'])
    tree_circular_500 = ttest_rel(sample_500['tree'], sample_500['circular'])
    circular_path_500 = ttest_rel(sample_500['circular'], sample_500['path'])
    #tree_rocs_500 = ttest_rel(sample_500['tree'], sample_500['rocs'])
    #path_rocs_500 = ttest_rel(sample_500['path'], sample_500['rocs'])
    #circular_rocs_500 = ttest_rel(sample_500['circular'], sample_500['rocs'])
    print(' T-Test results on RR for tree against path for top500 : ', tree_path_500)
    print(' T-Test results on RR for tree against circular for top500 : ', tree_circular_500)
    print(' T-Test results on RR for circular against path for top500 : ', circular_path_500)
    #print(' T-Test results on RR for tree against rocs for top500 : ', tree_rocs_500)
    #print(' T-Test results on RR for path against rocs for top500 : ', path_rocs_500)
    #print(' T-Test results on RR for circular against rocs for top500 : ', circular_rocs_500)
    stream.write(' \n T-Test results on RR for tree against path for top500 : ' + str(tree_path_500) + '\n')
    stream.write(' T-Test results on RR for tree against circular for top500 : ' + str(tree_circular_500) + '\n')
    stream.write(' T-Test results on RR for circular against path for top500 : ' + str(circular_path_500) + '\n')
    #stream.write(' T-Test results on RR for tree against rocs for top500 : ' + str(tree_rocs_500) + '\n')
    #stream.write(' T-Test results on RR for path against rocs for top500 : ' + str(path_rocs_500) + '\n')
    #stream.write(' T-Test results on RR for circular against rocs for top500 : ' + str(circular_rocs_500) + '\n')

    #mean, var, std = bayes_mvs(sample_500['tree_rocs'])
    #print(mean)
    #stream.write(str(mean) + '\n')
    stream.close()

def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    results = pd.DataFrame()

    path = itf.GetString("-in_dir")
    folder = itf.GetString("-method")

    subdirs = [x for x in os.listdir( path ) if (os.path.isdir(path + x) and os.path.isdir(path + x + '/' + folder + '/'))]

    for AC in subdirs:
        result_files = [x for x in os.listdir( path + AC + '/' + folder + '/') if x.endswith('.csv')]
        for result_file in result_files:
            result = pd.read_csv(path + AC + '/' + folder + '/' + result_file)
            result.rename(columns=lambda x: x + " " + AC, inplace=True)
            results = pd.concat([results, result], axis = 1)

    output = itf.GetString("-output")
    PlotResults(results, subdirs, folder, output)
    CompareFPResults(results, subdirs, folder, output)
    TTest(results, subdirs, folder, output)

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

!PARAMETER -method
    !ALIAS -m
    !TYPE string
    !BRIEF Input Calculation method
    !REQUIRED true
    !KEYLESS 3
!END
"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
