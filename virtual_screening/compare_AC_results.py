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

def PlotResults(results, subdirs, plot_output):
    AC_count = 0
    AC_left = len(subdirs)

    path = plot_output + "AC_comparison.pdf"
    pdf_file = PdfPages(path)

    fig, axes = plt.subplots(nrows=min(3, AC_left), ncols=2, figsize=(9,8))
    plt.subplots_adjust(hspace=0.4)

    for AC in subdirs:
        if AC_left == 1 :
            results.plot(ax=axes[0], y = [col for col in results.columns if ('RR' in col and AC in col)], sharex = True, sharey = True)
            axes[0].set_xlabel('Top Rank Molecules')
            axes[0].set_ylabel('Rate (%)')
            axes[0].set_title("Avg Recovery Rates " + AC) 
            axes[0].legend([col.split(" ")[2] for col in results.columns if ('RR' in col and AC in col)])
            
            results.plot(ax=axes[1], y = [col for col in results.columns if ('HR' in col and AC in col)], sharex = True, sharey = True)
            axes[1].set_xlabel('Top Rank Molecules')
            axes[1].set_title("Avg Hit Rates " + AC) 
            axes[1].legend([col.split(" ")[2] for col in results.columns if ('HR' in col and AC in col)])

        else :
            results.plot(ax=axes[AC_count, 0], y = [col for col in results.columns if ('RR' in col and AC in col)], sharex = True, sharey = True) 
            axes[AC_count, 0].set_xlabel('Top Rank Molecules')
            axes[AC_count, 0].set_ylabel('Rate (%)')
            axes[AC_count, 0].set_title("Avg Recovery Rates " + AC) 
            axes[AC_count, 0].legend([col.split(" ")[2] for col in results.columns if ('RR' in col and AC in col)])
            
            results.plot(ax=axes[AC_count, 1], y = [col for col in results.columns if ('HR' in col and AC in col)], sharex = True, sharey = True)
            axes[AC_count, 1].set_xlabel('Top Rank Molecules')
            axes[AC_count, 1].set_title("Avg Hit Rates " + AC) 
            axes[AC_count, 1].legend([col.split(" ")[2] for col in results.columns if ('HR' in col and AC in col)])
        
        AC_count += 1
        if AC_count == 3: 
            pdf_file.savefig(fig)
            AC_count = 0
            AC_left -= 3
            fig, axes = plt.subplots(nrows=min(3, AC_left), ncols=2, figsize=(9,8))
            plt.subplots_adjust(hspace=0.4)
        
    pdf_file.savefig(fig)
    pdf_file.close()
    #plt.savefig(path)
        
        #plt.show()

def main(argv=[__name__]):
    itf = OEInterface(InterfaceData, argv)

    results = pd.DataFrame()

    path = itf.GetString("-in_dir")

    subdirs = [x for x in os.listdir( path ) if (os.path.isdir(path + x) and os.path.isdir(path + x + '/FFP/'))]

    for AC in subdirs:
        result_files = [x for x in os.listdir( path + AC + '/FFP/') if x.endswith('.csv')]
        for result_file in result_files:
            result = pd.read_csv(path + AC + '/FFP/' + result_file)
            result.rename(columns=lambda x: x + " " + AC, inplace=True)
            results = pd.concat([results, result], axis = 1)

    output = itf.GetString("-output")
    
    PlotResults(results, subdirs, output)

InterfaceData = """

!PARAMETER -output
    !ALIAS -o
    !TYPE string
    !BRIEF Output file
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
