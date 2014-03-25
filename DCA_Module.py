import numpy as np
import numpy.linalg as npl
import random as rand
import sys
import csv


# Instantiate Constants
AA_out='0'
AA_matrix=['+','-']
x=0.7
q=3



#######  FUNCTION IMPLEMENTATIONS ##########

def load_txt_list(filename, dtype = "str"):
    print "Executing load_txt_list"
    raw_text_input = np.loadtxt(filename, dtype)
    line_parsed_list = []
    for seq in raw_text_input:
      seq = list(seq)
      line_parsed_list.append(seq)
    data_matrix = np.array(line_parsed_list)
    return [line_parsed_list, data_matrix, raw_text_input]


def single_residue_freq(concat_matrix, AA_matrix, q, AA_out, km_matrix):
  print "Executing single_residue_freq"

  [M,L] = np.shape(concat_matrix)
  Meff = sum(km_matrix)
  lam = Meff
  fiA_matrix_matrix=[]
  fiA_out_list=[]
  insert_state_list = []

  # iterate over length of data sequence (L is number of timeseries points)
  for i in range(L):
      count_matrix=[0]*2
      out_count=0

      # iterate over length of actors (M is number of companies in data set)
      for j in range(M):
          km = km_matrix[j]
          if concat_matrix[j][i] != AA_out:
              index=AA_matrix.index(concat_matrix[j][i])
              count_matrix[index] += km
          else:
              out_count += km

      fiA_matrix = []
      for k in range(2):
          fiA = (1 / float(lam + Meff)) * ((lam / float(q)) + count_matrix[k])
          fiA_matrix.append(fiA)


      #fiA_out is (1/(2*M)(M/q)) + out_count
      fiA_out = (1 / float(lam + Meff)) * ((lam / float(q)) + out_count)
      fiA_matrix_matrix += fiA_matrix
      fiA_out_list.append(fiA_out)

  #returns from while loops: fiA_matrix_matrix and fiAB (only values that persist in later calcs)
  fiA_vector1 = np.array(fiA_matrix_matrix)
  fiA_vector2 = np.reshape(fiA_vector1, (2*L,1))
  fiAB=fiA_vector1 * fiA_vector2
  return [fiA_matrix_matrix, fiAB, fiA_vector1, fiA_out_list]

def residue_pair_freq(concat_matrix, km_matrix):
  print "Executing residue_pair_freq"
  [M,L] = np.shape(concat_matrix)
  Meff = sum(km_matrix)
  lam = Meff

  fijAB_matrix=np.zeros((2*L,2*L))

  for i in range(L):
    for j in range(L):
      for k in range(M):
        i_val = concat_matrix[k][i]
        j_val = concat_matrix[k][j]
        if i_val != AA_out and j_val != AA_out:
          i_index = AA_matrix.index(i_val)
          j_index = AA_matrix.index(j_val)
          pair_address_x = (i) * 2 + i_index
          pair_address_y = (j) * 2 + j_index
          fijAB_matrix[pair_address_y][pair_address_x] += km_matrix[k]

  fijAB_matrix = fijAB_matrix + (lam / float(q**2))
  fijAB_matrix = fijAB_matrix / float(lam + Meff)
  return fijAB_matrix

def generate_covariance_matrix(fijAB_matrix, fiAB, fiA_matrix_matrix):
  print "Executing generate_covariance_matrix"
  L = len(fiAB)/2
  Cij_matrix = fijAB_matrix - fiAB

  for i in range(L):
    for j in range(2):
      for k in range(2):
        Cij_matrix[2*i+j][2*j+k] = -(fiA_matrix_matrix[2*i+j] * fiA_matrix_matrix[2*j+k])
  for i in range(2*L):
    Cij_matrix[i][i] = fiA_matrix_matrix[i] * (1- fiA_matrix_matrix[i])
  return Cij_matrix

def generate_coupling_energies(Cij_matrix):
  print "Executing generate_coupling_energies"
  eij_matrix1 = npl.inv(Cij_matrix)
  eij_matrix = -eij_matrix1
  return eij_matrix