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

def PlotResults(results, subdirs, folder, plot_output):
    AC_count = 0
    AC_left = len(subdirs)

    path = plot_output + folder + "_AC_comparison.pdf"
    pdf_file = PdfPages(path)

    fig, axes = plt.subplots(nrows=min(3, AC_left), ncols=2, figsize=(7,9))
    plt.subplots_adjust(hspace=0.4)

    for AC in subdirs:
        results.plot(ax=axes[AC_count, 0], y = [col for col in results.columns if ('RR' in col and AC in col)], sharex = True, sharey = True) 
        axes[AC_count, 0].set_xlabel('Top Rank Molecules')
        axes[AC_count, 0].set_ylabel('Rate (%)')
        axes[AC_count, 0].set_title("Avg Recovery Rates " + AC) 
        axes[AC_count, 0].legend([col.split(" ")[2] for col in results.columns if ('RR' in col and AC in col)])
        for line in axes[AC_count, 0].get_lines():
            if 'tree' in line.get_label():
                line.set_color("C0")
            if 'path' in line.get_label():
                line.set_color("C1")
            if 'circular' in line.get_label():
                line.set_color("C2")
            if 'FR' in line.get_label():
                line.set_color("red")

        results.plot(ax=axes[AC_count, 1], y = [col for col in results.columns if ('HR' in col and AC in col)], sharex = True, sharey = True)
        axes[AC_count, 1].set_xlabel('Top Rank Molecules')
        axes[AC_count, 1].set_title("Avg Hit Rates " + AC) 
        axes[AC_count, 1].legend([col.split(" ")[2] for col in results.columns if ('HR' in col and AC in col)])
        for line in axes[AC_count, 1].get_lines():
            if 'tree' in line.get_label():
                line.set_color("C0")
            if 'path' in line.get_label():
                line.set_color("C1")
            if 'circular' in line.get_label():
                line.set_color("C2")
            if 'FR' in line.get_label():
                line.set_color("red")

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
