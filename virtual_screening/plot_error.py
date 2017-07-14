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

def ReadResultsFile(path, nb_ka):
    datas = pd.DataFrame()
    raw_datas = pd.DataFrame()

    results = open(path, 'rb')
    chunk = results.read()
    chunk = chunk.decode('utf-8')
    results.close()

    Sets = chunk.split("Set n°")
    Sets.pop(0)
    for i, Set in enumerate(Sets):
        Set = Set.split("\n")
        Set.pop(0)
        Set.pop()

        count = 0
        count_ka = 0
        for row, line in enumerate(Set):
            set_results = pd.DataFrame(columns = ['RR', 'HR'])
            count += 1
            if 'True' in line: 
                count_ka += 1
            rr = 100 * count_ka/nb_ka
            hr = 100 * count_ka/count
            set_results.loc[row] = [rr, hr]
            raw_datas = pd.concat([raw_datas, set_results])
    print(raw_datas.loc[50])
    print(raw_datas.loc[499])

    datas['Average RR'] = raw_datas.groupby(raw_datas.index)['RR'].mean()
    datas['Average HR'] = raw_datas.groupby(raw_datas.index)['HR'].mean()
    datas['Minimum RR'] = raw_datas.groupby(raw_datas.index)['RR'].min()
    datas['Minimum HR'] = raw_datas.groupby(raw_datas.index)['HR'].min()
    datas['Maximum RR'] = raw_datas.groupby(raw_datas.index)['RR'].max()
    datas['Maximum HR'] = raw_datas.groupby(raw_datas.index)['HR'].max()
    datas = datas.head(500)
    print('Values at top 100 : Average = ', datas['Average RR'][99], ' ; Min = ', datas['Minimum RR'][99], ' ; Max = ', datas['Maximum RR'][99])
    print('Values at top 500 : Average = ', datas['Average RR'][499], ' ; Min = ', datas['Minimum RR'][499], ' ; Max = ', datas['Maximum RR'][499])

    return datas

def PlotDatas(output, datas):
    path = output
    idx = [i for i in range(len(datas)) if i%50==0]
    data_points=datas.iloc[idx]
    print(data_points)
    pdf_file = PdfPages(path)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7,3))

    #min_error = data_points['Average RR']-data_points['Minimum RR']
    #max_error = data_points['Maximum RR']-data_points['Average RR']
    min_error = datas['Average RR']-datas['Minimum RR']
    max_error = datas['Maximum RR']-datas['Average RR']
    datas.plot(ax=axes[0], y = [col for col in datas.columns if ('RR' in col)])
    #data_points.plot(kind='box', ax=axes[0], y = 'Average RR', yerr=[2, 1], sharex = True, sharey = True) 
    plt.errorbar(np.arange(500), datas['Average RR'], yerr=[min_error, max_error], errorevery=50, lw=1)
    axes[1].legend(['Average RR with error bars'])
    pdf_file.savefig(fig)
    pdf_file.close()
        #result.rename(columns=lambda x: x + " " + AC, inplace=True)
        #results = pd.concat([results, result], axis = 1)

def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    datas = pd.DataFrame()

    input_results = itf.GetString("-input")
    nb_ka = itf.GetInt("-num_ka")

    datas = ReadResultsFile(input_results, nb_ka)

    output = itf.GetString("-output")
    PlotDatas(output, datas)
    #PlotResults(results, subdirs, folder, output)
    #CompareFPResults(results, subdirs, folder, output)

InterfaceData = """

!PARAMETER -input
    !ALIAS -i
    !TYPE string
    !BRIEF Input Results file
    !REQUIRED true
    !KEYLESS 1
!END

!PARAMETER -output
    !ALIAS -o
    !TYPE string
    !BRIEF Output file
    !REQUIRED true
    !KEYLESS 2
!END

!PARAMETER -num_ka
    !ALIAS -ka
    !TYPE int
    !BRIEF Number of Known actives in the set
    !REQUIRED true
    !KEYLESS 3
!END
"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
