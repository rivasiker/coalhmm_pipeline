from Bio import AlignIO
from Bio.AlignIO import MafIO
import os
import pandas as pd
import sys
from gwf import Workflow

gwf = Workflow()

if not os.path.isdir('../tmp/'):
	os.mkdir('../tmp/')

###################### FILTERING OF THE MAF FILE ######################

## 	Required files: Maffilter stand alone program required (can't be installed via conda)

##	Required parameters: strings of : big_maf_file, species1, species2, species3, species4

big_maf_file = '../../chr3_subset/subset_chr3.maf.gz'
species1 = 'Homo_sapiens'
species2 = 'Pan_troglodytes'
species3 = 'Gorilla_gorilla_gorilla'
species4 = 'Pongo_abelii'

gwf.target('Mafffilter_control_file', 
           inputs=[], 
		   outputs=['../tmp/control_file'],
		   cores=1,
    	   memory='2g',
		   walltime= '00:10:00') << """
./maffilter_controlfile_generation.sh {} {} {} {} {}
""".format(big_maf_file, species1, species2, species3, species4)


gwf.target('Maffilter', 
		   inputs=['../tmp/control_file'], 
		   outputs=['../tmp/filtered.maf', '../tmp/maf_filtering.log'],
		   cores=1,
    	   memory='10g',
		   walltime= '02:00:00') << """
./maffilter param=../tmp/control_file
"""

gwf.target('coalhmm_params', 
		   inputs=[], 
		   outputs=['../tmp/params.file'],
		   cores=1,
    	   memory='10g',
		   walltime= '02:00:00') << """
./coalhmm_paramfile_generation.sh {} {} {} {}
""".format(species1, species2, species3, species4)

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

# Define the reference species and chromosome
target_seqname = 'Homo_sapiens.chr3'
# Define the window size
window_size = 1000000

# TRICK FOR GWF TO WORK
# If the splice file exists
if not os.path.isfile('../tmp/slice_dct.txt'): 
	# If the filtered maf exists
	if os.path.isfile('../tmp/filtered.maf'):
		# Load the alignment
		alignment = AlignIO.parse('../tmp/filtered.maf', 'maf')
		# Save slice list
		slice_lst = start_end(alignment, target_seqname, window_size)
		# Save slice list as temporary file
		pickle.dump(slice_lst, open('../tmp/slice_dct.txt', 'wb'))

		# Create/load maf indexing
		idx = MafIO.MafIndex('../tmp/filtered.mafindex', '../tmp/filtered.maf', target_seqname)
	
	# If the filtered maf does not exist
	else:
		# Save an empty list
		slice_lst = []
# If the splice list exists
else:
	# Load slice dictionary
	slice_lst = pickle.load(open('../tmp/slice_dct.txt', 'rb'))
	# Create maf indexing
	idx = MafIO.MafIndex('../tmp/filtered.mafindex', '../tmp/filtered.maf', target_seqname)



###################### CREATE INFO TABLES AND RUN COALHMM ######################


##	Required parameters:
##		target_seqname: string with reference sequence and chromosome in the
##						form of 'species.chromosome'.
##
##	Required files:
##		'../tmp/filtered.maf':	a filtered maf file as outputted from maffilter.
##      '../tmp/slice_dct.txt': a text file containing the dictionary of slices.



if not os.path.isdir('../tmp/info_tables'):
	# Create temporary directory with coalHMM run index
	os.mkdir('../tmp/info_tables')


# If the filtered maf file exists, then run coalHMM
if os.path.isfile('../tmp/filtered.maf'):
	# For each slice of the maffilter
	for run in range(len(slice_lst)):
		if not os.path.isdir('../tmp/run_{}'.format(run)):
			# Create temporary directory with coalHMM run index
			os.mkdir('../tmp/run_{}'.format(run))
		# Create input and output lists for the gwf run
		inputs = ['../tmp/filtered.mafindex', '../tmp/filtered.maf']
		outputs = ['../tmp/run_{}/fasta_{}.fa'.format(run, j) for j in range(slice_lst[run][2])]
		outputs.append('../tmp/info_tables/run_{}.csv'.format(run))
		# Save individual fasta files and info df
		gwf.target('run_{}'.format(run), 
				inputs=inputs, outputs=outputs,
				cores=4,
				memory='4g',
				walltime= '01:00:00') << """
		python create_fasta_and_info_table.py {} {} {} {}
		""".format(run, target_seqname, slice_lst[run][0], slice_lst[run][1])
	#run coalhmm
	gwf.target('coalhmm_run_{}'.format(run), input=, output=cores=4, memory='4g', walltime= '02:00:00') << """
	coalhmm --noninteractive=yes param=../tmp/params.file species1={} species2={} species3={} outgroup={} input.sequence.multiparts=yes input.sequence.format=Fasta input.sequence.list={} input.sequence.multiparts.prefix=#XXX?(A path to be added before all paths in the file list) input.sequence.multiparts.reset=yes optimization.profiler={}.profiler optimization.message_handler={}.messages output.posterior.states=none output.hidden_states={}.states output.hidden_states.divergences={}.divergences output.posterior.values={}.posteriors output.estimated.parameters={}.params output.userfriendly.parameters={}.estimates input.sequence.sites_to_use=all input.sequence.max_gap_allowed=50%
""".format(species1, species2, species3, species4, INPUT_LIST, PATH??, BASENAME)
















