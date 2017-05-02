#!/usr/bin/env python

import sys
import os
import random

def RandomIndex(total, nb_index, index_out, iteration):
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

RandomIndex(total, nb_index, index_out, iteration)