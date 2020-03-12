from Bio import AlignIO
from Bio.AlignIO import MafIO
import os
import pandas as pd
import sys
import getopt

# Path to the input uncompressed maf file
input_file = sys.argv[1]
# Path to an existing output directory
output_dir = sys.argv[2].strip('/')+'/'
# Name of the output file 
output_name = sys.argv[3]
# Target sequence in the form of species.chromosome
target_seqname = sys.argv[4]
# Size of the window for slicing the maf file in bp
window_size = sys.argv[5]



# Parse the maf alignment from the first argument of the script
alignment = AlignIO.parse(input_file, 'maf')
# Create or load a mafindex
idx = MafIO.MafIndex(output_dir+output_name+".mafindex", input_file, target_seqname)


len(alignment)