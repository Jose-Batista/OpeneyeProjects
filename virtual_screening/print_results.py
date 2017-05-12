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

def PlotResults(results, plot_output):
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9,8))
    plt.subplots_adjust(hspace=0.4)

    results.plot(ax=axes[0], y = [col for col in results.columns if 'RR' in col])
    axes[0].set_xlabel('Top Rank Molecules')
    axes[0].set_ylabel('Rate (%)')
    axes[0].set_title("Average RR Rates") 
    #path = plot_output + "Average_RR_plot.svg"
    #plt.savefig(path)
    
    results.plot(ax=axes[1], y = [col for col in results.columns if 'HR' in col])
    axes[1].set_xlabel('Top Rank Molecules')
    axes[1].set_ylabel('Rate (%)')
    axes[1].set_title("Average HR Rates") 
    
    path = plot_output + "Average_HR_plot.svg"
    plt.savefig(path)
        
        #plt.show()

def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    results = pd.DataFrame()
    for result_file in itf.GetStringList("-in_data"):
        result = pd.read_csv(result_file)
        results = pd.concat([results, result], axis = 1)

    print(results)

    output = itf.GetString("-output")
    
    PlotResults(results, output)

InterfaceData = """

!PARAMETER -output
    !ALIAS -o
    !TYPE string
    !BRIEF Output file
    !REQUIRED true
    !KEYLESS 1
!END

!PARAMETER -in_data
    !ALIAS -i
    !TYPE string
    !BRIEF Input Database of Active Molecules
    !REQUIRED true
    !LIST true
    !KEYLESS 2
!END
"""

if __name__ == "__main__":
        sys.exit(main(sys.argv))
