import numpy as np
import numpy.linalg as npl
import random as rand

#load concatenated sequences

concat_list_raw=np.loadtxt('sequence_list_abbrev', dtype = 'str')


concat_list = []
for seq in concat_list_raw:

    seq = list(seq[:50])
    concat_list.append(seq)

concat_matrix=np.array(concat_list)

print np.shape(concat_matrix)

#compile sequence matrix

#constants

AA_out='0'
x=0.7
L=len(concat_list[0])
M=len(concat_list)
q=3
AA_matrix=['+','-']

L = len(concat_matrix[0])

# calculate km

print 'calculating km...'

#t=0
#km_matrix=[]
#while t<M:
#    t=t+1
#    t1=0
#    km=0
#    while t1<M:
#        t1=t1+1
#        t2=0
#        Kronecker=0
#        while t2<L:
#            t2=t2+1
#            A=concat_matrix[t-1][t2-1]
#            B=concat_matrix[t1-1][t2-1]
#            if A==B:
#                Kronecker=Kronecker+1

#        theta=Kronecker-x*L
#        if theta>0:
#            km=km+1
#    km = 1/float(km)
#    km_matrix.append(km)
#    print t, km

km_matrix = [1]*M

km_matrix = np.array(km_matrix)

np.savetxt('stock_km_matrix', km_matrix)

# calculate km

# calculate number of effective sequences (Meff)

Meff = sum(km_matrix)
#print 'Meff = ', Meff

lam=Meff

# calculate numer of effective sequences (Meff)

# calculate single amino acid frequencies (fiA)

t=0
fiA_matrix_matrix=[]
fiA_out_list=[]
insert_state_list = []
while t<L:
    t=t+1
    t1=0
    count_matrix=[0]*2
    out_count=0
    while t1<M:
        t1=t1+1
        count=km_matrix[t1-1]
        if concat_matrix[t1-1][t-1] != AA_out:
            index=AA_matrix.index(concat_matrix[t1-1][t-1])
            count_matrix[index]=count_matrix[index]+count
        else:
            out_count=out_count+count

    t2=0
    fiA_matrix=[]
    while t2<2:
        t2=t2+1
        fiA=(1/float(lam+Meff))*((lam/float(q))+count_matrix[t2-1])
        fiA_matrix.append(fiA)
    fiA_out=(1/float(lam+Meff))*((lam/float(q))+out_count)
    fiA_matrix_matrix = fiA_matrix_matrix+fiA_matrix
    fiA_out_list.append(fiA_out)

fiA_vector1=np.array(fiA_matrix_matrix)
fiA_vector2=np.reshape(fiA_vector1, (2*L,1))
fiAB=fiA_vector1*fiA_vector2

#calculate single amino acid frequencies (fiA)

# calculate pair frequencies (fijAB)

print 'calculating Cij matrix...'

dimension=2*L
fijAB_matrix=np.zeros((dimension,dimension))

i=0
while i<L:
    i=i+1
    #print '     ',i
    j=0
    while j<L:
        j=j+1
        t=0
        while t<M:
            t=t+1
            i_val=concat_matrix[t-1][i-1]
            j_val=concat_matrix[t-1][j-1]
            if i_val != AA_out and j_val != AA_out:
                i_index=AA_matrix.index(i_val)
                j_index=AA_matrix.index(j_val)
                pair_address_x=(i-1)*2+i_index
                pair_address_y=(j-1)*2+j_index
                fijAB_matrix[pair_address_y][pair_address_x]=fijAB_matrix[pair_address_y][pair_address_x]+km_matrix[t-1]

fijAB_matrix=fijAB_matrix+(lam/float(q**2))
fijAB_matrix=fijAB_matrix/float(lam+Meff)

# calculate pair frequencies (fijAB)

# calculate empirical correlation matrix

Cij_matrix=fijAB_matrix-fiAB

t=0
while t<L:
    t=t+1

    t1=0
    while t1 < 2:
        t1 = t1+1

        t2 = 0
        while t2 < 2:
            t2 = t2+1

            Cij_matrix[2*(t-1) + t1-1][2*(t-1) + t2-1] = -fiA_matrix_matrix[2*(t-1)+t1-1]*fiA_matrix_matrix[2*(t-1)+t2-1]

t=0
while t<dimension:
    t=t+1
    Cij_matrix[t-1][t-1]=fiA_matrix_matrix[t-1]*(1-fiA_matrix_matrix[t-1])

# calculate couplings (eij)

print 'calculating couplings...'

eij_matrix1=npl.inv(Cij_matrix)
eij_matrix=-eij_matrix1

# calculate couplings (eij)

# save

np.savetxt('stocks_fiA', fiA_vector1)
np.savetxt('stocks_fiA_out', fiA_out_list)
np.savetxt('stocks_eij', eij_matrix)
np.savetxt('Covariant Matrix', Cij_matrix)

# save

print 'end'

