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

    results.plot(y = [col for col in results.columns if 'RR' in col], label = "Average RR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average RR Rates FP") 
    path = plot_output + "Average_RR_plot.svg"
    plt.savefig(path)
    
    results.plot(y =  [col for col in results.columns if 'HR' in col], label = "Average HR")
    plt.xlabel('Top Rank Molecules')
    plt.ylabel('Rate (%)')
    plt.legend( loc='best')
    plt.title("Average HR Rates FP")
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
