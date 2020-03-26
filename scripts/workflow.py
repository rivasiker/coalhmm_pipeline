from Bio import AlignIO
from Bio.AlignIO import MafIO
import os
import pandas as pd
import sys
from gwf import Workflow


###################### FILTERING OF THE MAF FILE ######################


###################### SLICING OF THE FILTERED MAF FILE ######################


##	Required parameters:
##		target_seqname: string with reference sequence and chromosome in the
##						form of 'species.chromosome'.
##		window_size:	integer corresponding to the preferred window size
##						for slicing the maf file for the coalHMM program.
##
##	Required files:
##		'../tmp/filtered.maf':	a filtered maf file as outputted from maffilter.

from start_end import start_end
import pickle

if not os.path.exists('../tmp/slice_dct.txt'):
	# Load the alignment
	alignment = AlignIO.parse('../tmp/filtered.maf', 'maf')
	# Save slice list
	slice_lst = start_end(alignment, target_seqname, window_size)
	# Save slice list as temporary file
	pickle.dump(slice_lst, open('../tmp/slice_dct.txt', 'wb'))
else:
	# Load slice dictionary
	slice_lst = pickle.load(open('../tmp/slice_dct.txt', 'rb'))


###################### CREATE INFO TABLES AND RUN COALHMM ######################


##	Required parameters:
##		target_seqname: string with reference sequence and chromosome in the
##						form of 'species.chromosome'.
##
##	Required files:
##		'../tmp/filtered.maf':	a filtered maf file as outputted from maffilter.
##      '../tmp/slice_dct.txt': a text file containing the dictionary of slices.




# It needs to be created above, in the first target
########## os.mkdir('../tmp/info_tables')


# Create maf indexing
idx = MafIO.MafIndex('../tmp/filtered.mafindex', '../tmp/filtered.maf', target_seqname)
# For each slice of the maffilter
for run in range(len(slice_lst)):
	# Create temporary directory with coalHMM run index
	os.mkdir('../tmp/run_{}'.format(run))
	# Create input and output lists for the gwf run
	inputs = ['../tmp/filtered.mafindex', '../tmp/filtered.maf']
	outputs = ['../tmp/run_{}/fasta_{}.fa'.format(run, j) for j in range(slice_lst[2])]
	outputs.append('../tmp/info_tables/run_{}.csv'.format(run))
	# Save individual fasta files and info df
	gwf.target('run_{}'.format(run), inputs=inputs, outputs=outputs) << """
	python create_fasta_and_info_table.py {} {} {} {}
	""".format(run, target_seqname, slice_lst[run][0], slice_lst[run][1])
	# Run coalHMM





gwf = Workflow()

gwf.target('MyTarget', inputs=[], outputs=[]) << """
echo hello world
"""


















