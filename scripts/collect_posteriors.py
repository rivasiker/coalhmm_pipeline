import pandas as pd
import sys

# Save the run identifier
run = int(sys.argv[1])

# Load the info table with the coordinates and the gap information
info_table = pd.read_csv('../tmp/info_tables/run_{}.csv'.format(run))
# Load the posterior probabilities
posteriors = pd.read_csv('../tmp/outputs/run_{}/posteriors'.format(run), delim_whitespace=True)

# Create empty dictionaries for the position and the chromosome
pos = {}
chrom = {}
# Add species as keys to the dictionaries
for species in list(set(info_table['species'])):
    pos[species] = []
    chrom[species+'.chr'] = []
    
# For each individual fasta file in the run
for n in range(max(info_table['file'])+1):
    # Filter the sub-table for the fasta file
    tab = info_table[info_table['file'] == n]
    # Re-set binary vector for the columns for which the posteriors have been computed
    a = ''
    # Define list of gap information from the subtable
    lst = list(tab['gaps'])
    # For each index in the list of gap info
    for i in range(len(lst[0])):
        # If at least 50% of the rows are non-gaps (from the coalHMM param file)
        if int(lst[0][i])+int(lst[1][i])+int(lst[2][i])+int(lst[3][i]) >= 2:
            # Add a 1 to the bianry vector
            a += '1'
        else:
            # Add a 0 to the binary vector
            a += '0'

    # For each of the species
    for count, species in enumerate(list(tab['species'])):
        # Save the start coordinate
        start = int(tab[tab['species'] == species]['start'])
        # Save the chromosome name
        chromosome = tab[tab['species'] == species]['chr'].to_string(index = False).strip()
        # For each element in the binary vector
        for i in range(len(a)):
            # If at least 50% of the rows are non-gaps
            if a[i] == '1':
                # If the site is not a gap
                if lst[count][i] == '1':
                    # Add the coordinate
                    pos[species].append(start)
                    # Add the chromosome
                    chrom[species+'.chr'].append(chromosome)
                    # Update coordinates
                    start += 1
                # If the site is a gap
                else:
                    # Add None to both dictionaries
                    pos[species].append(None)
                    chrom[species+'.chr'].append(None)
            # If less of 50% of the rows are non-gaps
            elif a[i] == '0':
                # If the site is not a gap
                if lst[count][i] == '1':
                    # Update the coordinate 
                    start += 1

# Save the concatenated csv file
pd.concat([pd.DataFrame.from_dict(pos, dtype='int64'),
           pd.DataFrame.from_dict(chrom),
           posteriors.reset_index(drop=True).drop(['Chunk'], axis=1)], 
          axis = 1).to_csv('../tmp/results/run_{}.csv'.format(run), index = False)