#!/usr/bin/env python
from __future__ import print_function
from openeye.oechem import *

import sys
import os
import random

def CountDatabase(database):
	ifs = oemolistream()

	if not ifs.open(database):
		OEThrow.Fatal("Unable to open inputfile" )

	count = 0

	for mol in ifs.GetOEMols():
		count += 1
	return count

def RandomIndex(total, ratio, iteration, index_out):
	nb_index = total//(ratio + 1)

	with open(index_out, "w"): pass
	
	for iter in range (iteration):
		i = 0

		index_set = set()
		while i < nb_index :
			index = random.randint(0, total - 1)
			if index in index_set:
				continue
			else : 
				index_set.add(index)
				i += 1

		index_log = open(index_out, 'a')
		index_log.write("set N°%d: " % iter)
		for index in index_set:
			index_log.write(str(index) + " ")
		index_log.write("\n")
		index_log.close()

	return index_set

total = 2040
nb_index = total//(3)
index_out = 'clean_databases/FA10/index_log_FA10.txt'
iteration = 5

def main(argv=[__name__]):
	itf = OEInterface(InterfaceData, argv)

	ina = itf.GetString("-act_database")
	ratio = itf.GetInt("-ratio")
	iteration = itf.GetInt("-iteration")
	index_out = itf.GetString("-output")

	tot_range = CountDatabase(ina)

	RandomIndex(tot_range, ratio, iteration, index_out)

InterfaceData = """
!PARAMETER -act_database
	!ALIAS -ad
	!TYPE string
	!BRIEF Database to be Counted
	!REQUIRED true
	!KEYLESS 1
!END

!PARAMETER -ratio
	!ALIAS -r
	!TYPE int
	!BRIEF Ratio of index generated
	!REQUIRED true
	!KEYLESS 2
!END

!PARAMETER -iteration
	!ALIAS -it
	!TYPE int
	!BRIEF Number of Iterations wanted for the Test
	!REQUIRED true
	!KEYLESS 3
!END

!PARAMETER -output
	!ALIAS -o
	!TYPE string
	!BRIEF Output File
	!REQUIRED true
	!KEYLESS 4
!END

"""

if __name__ == "__main__":
		sys.exit(main(sys.argv))