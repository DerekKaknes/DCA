from DCA_Module import *
import numpy as np
import time

global_start = time.time()

# Generate covariance matrix?  
generate_Cij = "Yes"

if generate_Cij == "Yes":

  # Import Data
  ts = time.time()
  [concat_list, concat_matrix, concat_list_raw] = load_txt_list('test_sequence.txt')
  runtime = time.time() - ts
  print "Executed in ", runtime, "seconds"

  # Instantiate Constants
  AA_out = '0'
  AA_matrix = ['+','-']
  x = 0.7
  q = len(AA_matrix) + len(AA_out)
  [M, L] = np.shape(concat_matrix)

  # IN FUTURE INSERT Km IMPLEMENTATION
  km_matrix = np.array([1]*M)
  # IN FUTURE INSERT Km IMPLEMENTATION 

  # Calculate single frequencies (fiA)
  ts = time.time()
  [fiA_matrix_matrix, fiAB, fiA_vector1, fiA_out_list] = single_residue_freq(concat_matrix, AA_matrix, AA_out, km_matrix)
  runtime = time.time() - ts
  print "Executed in ", runtime, "seconds"

  # Calculate residue pair frequencies (fijAB)
  ts = time.time()
  fijAB_matrix = residue_pair_freq(concat_matrix, km_matrix, AA_matrix, AA_out)
  runtime = time.time() - ts
  print "Executed in ", runtime, "seconds"

  # Generate covariance matrix
  ts = time.time()
  Cij_matrix = generate_covariance_matrix(fijAB_matrix, fiAB, fiA_matrix_matrix)
  runtime = time.time() - ts
  print "Executed in ", runtime, "seconds"

  # Generate coupling energies
  ts = time.time()
  eij_matrix = generate_coupling_energies(Cij_matrix)
  runtime = time.time() - ts
  print "Executed in ", runtime, "seconds"

  # Save relevant output files
  np.savetxt('stock_km_matrix', km_matrix)
  np.savetxt('stocks_fiA', fiA_vector1)
  np.savetxt('stocks_fiA_out', fiA_out_list)
  np.savetxt('stocks_eij', eij_matrix)
  np.savetxt('Covariant Matrix', Cij_matrix)

else:
  fiA_vector1 = np.loadtxt('stocks_fiA')
  eij_matrix = np.loadtxt('stocks_eij')
  fiA_out_list = np.loadtxt('stocks_fiA_out')
  L = len(fiA_vector1)/2


# Reorganize fiA and fiA_out_list
ts = time.time()
fiA_matrix = join_fiA(fiA_vector1, fiA_out_list)
runtime = time.time() - ts
print "Executed in ", runtime, "seconds"

# Calc DI for each residue
ts = time.time()
DI_matrix = direct_information(L,q,eij_matrix,fiA_matrix)
runtime = time.time() - ts
print "Executed in ", runtime, "seconds"

np.savetxt('fiA_list_list', fiA_matrix)
np.savetxt('stocks_DI_scores', DI_matrix)

print "End of program - executed in ", time.time() - global_start, "seconds"
print "DI Scores: ", DI_matrix
