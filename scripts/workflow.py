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

gwf.target('Maffilter_control_file', 
           inputs=[], 
		   outputs=['../tmp/control_file', '../tmp/params.file'],
		   cores=1,
    	   memory='2g',
		   walltime= '00:10:00') << """
./maffilter_controlfile_generation.sh {} {} {} {} {}
./coalhmm_paramfile_generation.sh {} {} {} {}
""".format(big_maf_file, species1, species2, species3, species4, species1, species2, species3, species4)


gwf.target('Maffilter', 
		   inputs=['../tmp/control_file'], 
		   outputs=['../tmp/filtered.maf', '../tmp/maf_filtering.log'],
		   cores=1,
    	   memory='10g',
		   walltime= '02:00:00') << """
./maffilter param=../tmp/control_file
"""

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
	os.mkdir('../tmp/info_tables')
if not os.path.isdir('../tmp/fasta_names'):
	os.mkdir('../tmp/fasta_names')
if not os.path.isdir('../tmp/inputs'):
	os.mkdir('../tmp/inputs')
if not os.path.isdir('../tmp/outputs'):
	os.mkdir('../tmp/outputs')
if not os.path.isdir('../tmp/results/'):
	os.mkdir('../tmp/results/')


# If the filtered maf file exists, then run coalHMM
if os.path.isfile('../tmp/filtered.maf'):
	if os.path.isfile('../tmp/slice_dct.txt'): 
		# For each slice of the maffilter
		for run in range(len(slice_lst)):
			
			if not os.path.isdir('../tmp/inputs/run_{}'.format(run)):
				# Create temporary directory with coalHMM run index
				os.mkdir('../tmp/inputs/run_{}'.format(run))
			if not os.path.isdir('../tmp/outputs/run_{}'.format(run)):
				# Create temporary directory with coalHMM run index
				os.mkdir('../tmp/outputs/run_{}'.format(run))


			# Create input and output lists for the gwf run
			inputs = ['../tmp/filtered.mafindex', '../tmp/filtered.maf']
			outputs = ['../tmp/inputs/run_{}/fasta_{}.fa'.format(run, j) for j in range(slice_lst[run][2])]
			outputs = outputs + ['../tmp/info_tables/run_{}.csv'.format(run), '../tmp/fasta_names/run_{}.txt'.format(run)]
			# Save individual fasta files and info df
			gwf.target('run_{}'.format(run), 
					   inputs=inputs, 
					   outputs=outputs,
					   cores=4,
					   memory='4g',
					   walltime= '01:00:00') << """
			python create_fasta_and_info_table.py {} {} {} {}
			""".format(run, target_seqname, slice_lst[run][0], slice_lst[run][1])


			# Define format of output paths
			out = '../tmp/outputs/run_{}/'.format(run)
			results = ['../tmp/outputs/run_{}/estimates'.format(run),
			 		   '../tmp/outputs/run_{}/posterior_states'.format(run),
			           '../tmp/outputs/run_{}/posteriors'.format(run), 
					   '../tmp/outputs/run_{}/hidden_states'.format(run),
					   '../tmp/outputs/run_{}/divergences'.format(run)]
			# Run coalhmm
			gwf.target('coalhmm_run_{}'.format(run), 
					   inputs=outputs, 
					   outputs=results,
					   cores=4, 
					   memory='4g', 
					   walltime= '01:00:00') << """
			./coalhmm --noninteractive=yes param=../tmp/params.file species1={} species2={} species3={} outgroup={} \
			input.sequence.multiparts=yes input.sequence.format=Fasta input.sequence.list={} input.sequence.multiparts.prefix={} \
			input.sequence.multiparts.reset=yes optimization.profiler={}profiler optimization.message_handler={}messages \
			output.posterior.states={}posterior_states output.hidden_states={}hidden_states output.hidden_states.divergences={}divergences \
			output.posterior.values={}posteriors output.estimated.parameters={}params output.userfriendly.parameters={}estimates \
			input.sequence.sites_to_use=all input.sequence.max_gap_allowed=50%
			""".format(species1, species2, species3, species4, 
			           '../tmp/fasta_names/run_{}.txt'.format(run), 
					   '../tmp/inputs/run_{}/'.format(run), 
					   out, out, out, out, out, out, out, out)

			# Collect results from each of the runs and combine them with the coordinates
			gwf.target('collect_run_{}'.format(run), 
					   inputs=['../tmp/info_tables/run_{}.csv'.format(run), '../tmp/outputs/run_{}/posteriors'.format(run)], 
					   outputs=['../tmp/results/run_{}.csv'.format(run)],
					   cores=1,
					   memory='4g',
					   walltime= '01:00:00') << """
			python collect_posteriors.py {}
			""".format(run)


		gwf.target('final_table', 
					inputs=['../tmp/results/run_{}.csv'.format(i) for i in range(len(slice_lst))], 
					outputs=['../tmp/final_table.csv'],
					cores=1,
					memory='4g',
					walltime= '01:00:00') << """
		OutFileName="../tmp/final_table.csv"      			 # Fix the output name
		i=0                                       			 # Reset a counter
		for filename in ../tmp/results/*.csv; do 
			if [ "$filename"  != "$OutFileName" ] ;      	 # Avoid recursion 
			then 
				if [[ $i -eq 0 ]] ; then 
					head -1  "$filename" >   "$OutFileName"  # Copy header if it is the first file
				fi
				tail -n +2  "$filename" >>  "$OutFileName"   # Append from the 2nd line each file
				i=$(( $i + 1 ))                              # Increase the counter
			fi
		done
		"""














