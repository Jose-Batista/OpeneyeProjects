#!/usr/bin/env python

import sys
import os
import random

def RandomIndex(total, nb_index, iteration, index_out):
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

	tot_range = itf.GetInt("-range")
	ratio = itf.GetInt("-ratio")
	iteration = itf.GetInt("-iteration")
	index_out = itf.GetString("-output")

	nb_index = tot_range//(ratio + 1)

	RandomIndex(tot_range, nb_index, iteration, index_out)

InterfaceData = """
!PARAMETER -range
	!ALIAS -ran
	!TYPE int
	!BRIEF Range of the Randomised set
	!REQUIRED true
	!KEYLESS 1
!END

!PARAMETER -ratio
	!ALIAS -rat
	!TYPE int
	!BRIEF Ratio of index generated
	!REQUIRED true
	!KEYLESS 2
!END

!PARAMETER -iteration
	!ALIAS -iter
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