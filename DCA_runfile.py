from DCA_Module import *
import numpy as np
import numpy.linalg as npl
import random as rand
import sys
import csv

# Import Data
[concat_list, concat_matrix, concat_list_raw] = load_txt_list('sequence_list')

# Instantiate Constants
AA_out='0'
AA_matrix=['+','-']
x=0.7
q=3
[M,L] = np.shape(concat_matrix)

# IN FUTURE INSERT Km IMPLEMENTATION
km_matrix = np.array([1]*M)
# IN FUTURE INSERT Km IMPLEMENTATION 

# Calculate single frequencies (fiA)
[fiA_matrix_matrix, fiAB, fiA_vector1, fiA_out_list] = single_residue_freq(concat_matrix, AA_matrix, q, AA_out, km_matrix)

# Calculate residue pair frequencies (fijAB)
fijAB_matrix = residue_pair_freq(concat_matrix, km_matrix)

# Generate covariance matrix
Cij_matrix = generate_covariance_matrix(fijAB_matrix, fiAB, fiA_matrix_matrix)

# Generate coupling energies
eij_matrix = generate_coupling_energies(Cij_matrix)

# Save relevant output files
np.savetxt('stock_km_matrix', km_matrix)
np.savetxt('stocks_fiA', fiA_vector1)
np.savetxt('stocks_fiA_out', fiA_out_list)
np.savetxt('stocks_eij', eij_matrix)
np.savetxt('Covariant Matrix', Cij_matrix)

