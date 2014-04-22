import numpy as np
import numpy.linalg as npl
import random as rand
import sys
import csv
from math import *
import time


# Instantiate Constants

#######  FUNCTION IMPLEMENTATIONS ##########

def load_txt_list(filename, dtype = "str"):
    print "Executing load_txt_list..."
    raw_text_input = np.loadtxt(filename, dtype)
    line_parsed_list = []
    for sequence in raw_text_input:
      line_parsed_list.append(list(sequence))
      print "Seq: ", sequence
    data_matrix = np.array(line_parsed_list)
    return [line_parsed_list, data_matrix, raw_text_input]


def single_residue_freq(concat_matrix, AA_matrix, AA_out, km_matrix):
  print "Executing single_residue_freq..."

  [M,L] = np.shape(concat_matrix)
  Meff = sum(km_matrix)
  lam = Meff
  q = len(AA_matrix) + len(AA_out)
  fiA_matrix_matrix = np.zeros((L,2), int)
  fiA_out_list = np.zeros(L, int)

  # iterate over each amino acid position in the protein sequence (L is the length of the sequence)
  for i in range(L):
      count_matrix = np.zeros(2, int)
      out_count = 0

      # iterate over each sample (M is number of experimentally determined sequences)
      for j in range(M):
          km = km_matrix[j] # km would normally down-weight sequences that are over-represented in the available experimental data
          if concat_matrix[j][i] != AA_out:
              index = AA_matrix.index(concat_matrix[j][i])
              count_matrix[index] += km
          else:
              out_count += km

      position_frequency_count = np.zeros(len(count_matrix), int)
      for k in range(len(count_matrix)):
          fiA = (1 / float(lam + Meff)) * ((lam / float(q)) + count_matrix[k])
          position_frequency_count[k] = fiA

      #fiA_out is (1/(2*M)(M/q)) + out_count
      fiA_out = (1 / float(lam + Meff)) * ((lam / float(q)) + out_count)
      fiA_out = out_count
      fiA_matrix_matrix[i] = position_frequency_count
      fiA_out_list[i] = fiA_out


  #returns from while loops: fiA_matrix_matrix and fiAB (only values that persist in later calcs)
  fiA_vector = np.reshape(fiA_matrix_matrix, 2 * L)
  fiA_trans = np.reshape(fiA_vector, (2*L, 1))
  fiAB = fiA_trans * fiA_vector 
  return [fiA_matrix_matrix, fiAB, fiA_vector, fiA_out_list]

def residue_pair_freq(concat_matrix, km_matrix, AA_matrix, AA_out):
  print "Executing residue_pair_freq..."
  [M,L] = np.shape(concat_matrix)
  Meff = sum(km_matrix)
  lam = Meff
  q = len(AA_matrix) + len(AA_out)

  fijAB_matrix = np.zeros((2*L,2*L))

  time1 = time.time()
  for i in range(L):
    time2 = time.time()
    diff = time2 - time1
    time1 = time2
    #print "Executed residue ",i," in ",diff," seconds"
    for j in range(L):
      for day in range(M): # M is very long, 3801 items, reps the number of trading day data points
        i_return = concat_matrix[day][i] # the return (+,-,0) of company i on day k
        j_return = concat_matrix[day][j] # the return (+,-,0) of company j on day k
        if i_return != AA_out and j_return != AA_out: # if the return of CompA and CompB is non-zero, then:
          i_index = AA_matrix.index(i_return)  # i_index is the return value (+,-) - could later include larger range of returns
          j_index = AA_matrix.index(j_return)
          pair_address_x = (i) * 2 + i_index # paX equals double the company index (i) plus the amino acid index
          pair_address_y = (j) * 2 + j_index
          fijAB_matrix[pair_address_y][pair_address_x] += km_matrix[day]
          #print "Companies (j,i, day) ",[j,i,day], "Pair Address: ",[pair_address_y, pair_address_x]," Returns: ",[j_return, i_return], " fijAB value: ",fijAB_matrix[pair_address_y][pair_address_x]
          #print "Sum of frequency matrix = ", sum(sum(fijAB_matrix))

  fijAB_matrix = fijAB_matrix + (lam / float(q**2))
  fijAB_matrix = fijAB_matrix / float(lam + Meff)
  return fijAB_matrix

def generate_covariance_matrix(fijAB_matrix, fiAB, fiA_matrix_matrix):
  print "Executing generate_covariance_matrix..."
  L = len(fiAB)/2
  Cij_matrix = fijAB_matrix - fiAB
  fiA_matrix_matrix = np.reshape(fiA_matrix_matrix, (np.size(fiA_matrix_matrix)))

  for i in range(L):
    for j in range(2):
      for k in range(2):
        Cij_matrix[2*i+j][2*j+k] = -(fiA_matrix_matrix[2*i+j] * fiA_matrix_matrix[2*j+k])
  for i in range(2*L):
    Cij_matrix[i][i] = fiA_matrix_matrix[i] * (1- fiA_matrix_matrix[i])
  return Cij_matrix

def generate_coupling_energies(Cij_matrix):
  print "Executing generate_coupling_energies..."
  eij_matrix1 = npl.inv(Cij_matrix)
  eij_matrix = -eij_matrix1
  return eij_matrix

def join_fiA(fiA, fiA_out):
  # This function takes in fiA and fiA_out and returns a new array, fiA_list_list,
  # which combines the input arrays in the pattern of [fiA[0], fiA[1], fiA_out[0]]
  print "Executing join_fiA..."
  result = []
  for i in range(len(fiA_out)):
      for j in range(2):
          result.append(fiA[2*i+j])
      result.append(fiA_out[i])
  return result


def eij_matrix(eij, i, j, q):
  # This function takes the (2L x 2L) matrix eij and maps to a qxq matrix equal to: 
  # [[e00,e01,1],
  #  [e10,e11,1],
  #  [1,1,1]]

  result = []
  for k in range(2):
      eij_list=[]
      for m in range(2):
          eij1 = exp(eij[2*i+k][2*j+m])
          #print "eij[",2*i+k,"][",2*j+m,"] = ",eij1,": i,j = ",i,",",j
          eij_list.append(eij1)
      # Set last column to 1
      eij_list.append(1)
      result.append(eij_list)
  # Set last row to 1's
  eij_out = [1]*q
  result.append(eij_out)
  return np.array(result)

def local_biases(pi,pj,W, epsilon = 0.0001):
  #This function takes pi, pj and W and returns the local_biases hiA and hiB
  # THIS FUNCTION IS THE PRIMARY COMPUTATIONAL EXPENSE
  q = len(pi)
  mu1 = np.ones(q)/q
  mu2 = np.ones(q)/q
  diff = 1

  while diff > epsilon:
      
      scra1 = np.dot(mu2,np.transpose(W))
      scra2 = np.dot(mu1,W)

      new1 = pi/scra1
      new1 = new1/sum(new1)

      new2 = pj/scra2
      new2 = new2/sum(new2)

      sub_diff1 = new1 - mu1
      sub_diff2 = new2 - mu2

      diff = max(max(np.absolute(sub_diff1)), max(np.absolute(sub_diff2)))
      
      mu1 = new1
      mu2 = new2

  hiA = mu1
  hiB = mu2
  return [hiA, hiB]

def individual_DI(fiA_list_list, pdir,i,j,q):
  DI = 0
  for k in range(q):
      for m in range(q):
          pdir_sub = pdir[k][m]
          DI += pdir_sub*log(pdir_sub/float(fiA_list_list[(i)*q + k]*fiA_list_list[(j)*q + m]))

  return DI

def direct_information(L,q,eij, fiA_list_list):
  print "Executing direct_information..."
  DI_matrix = []

  for i in range(L):
      DI_list = []
      for j in range(L):

          W = eij_matrix(eij,i,j,q)

          pi = np.array(fiA_list_list[(i)*q:(i)*q + q])
          pj = np.array(fiA_list_list[(j)*q:(j)*q + q])

          [hiA, hiB] = local_biases(pi,pj,W)

          pdir = W*(np.reshape(hiA, (q, 1))*hiB)
          pdir = pdir / pdir.sum()

          DI = individual_DI(fiA_list_list, pdir,i,j,q)
          DI_list.append(DI)

      DI_matrix.append(DI_list)

  return DI_matrix

