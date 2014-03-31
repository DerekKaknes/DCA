import numpy as np
#import scipy.optimize as so
from math import *

q = 3

print 'loading...'

fiA=np.loadtxt('stocks_fiA')
eij=np.loadtxt('stocks_eij')
fiA_out=np.loadtxt('stocks_fiA_out')


print 'calculating fields...'

L=len(fiA)/2
L = 5

t=0
fiA_list_list=[]
while t<L:
    t=t+1
    t1=0
    fiA_list=[]
    while t1<2:
        t1=t1+1
        fiA_list.append(fiA[(t-1)*2+t1-1])
    fiA_list.append(fiA_out[t-1])
    fiA_list_list=fiA_list_list+fiA_list

print len(fiA_list_list)
        
t=0
DI_list_list = []
while t<L:
    t=t+1
    print "i = ",t
    t1=0
    hiA_matrix=[]
    DI_list = []
    while t1<L:
        t1=t1+1
        t2=0
        eij_list_list=[]
        while t2<2:
            t2=t2+1
            t3=0
            eij_list=[]
            while t3<2:
                t3=t3+1
                eij1=exp(eij[(t-1)*2+t2-1][(t1-1)*2+t3-1])
                eij_list.append(eij1)
            eij_list.append(1)
            eij_list_list.append(eij_list)
        eij_out=[1]*q
        eij_list_list.append(eij_out)

        W = np.array(eij_list_list)

        epsilon = 0.0001
        diff = 1
        mu1 = [1]*q
        mu1 = np.array(mu1)
        mu1 = mu1/float(q)
        mu2 = [1]*q
        mu2 = np.array(mu2)
        mu2 = mu2/float(q)
        pi = fiA_list_list[(t-1)*q:(t-1)*q + q]
        pi = np.array(pi)
        pj = fiA_list_list[(t1-1)*q:(t1-1)*q + q]
        pj = np.array(pj)
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

        pdir = W*(np.reshape(mu1, (q, 1))*mu2)
        pdir = pdir / pdir.sum()

        t2 = 0
        DI_sub_list = []
        while t2 < q:
            t2 = t2 + 1

            t3 = 0
            while t3 < q:
                t3 = t3+1

                pdir_sub = pdir[t2-1][t3-1]

                DI_sub = pdir_sub*log(pdir_sub/float(fiA_list_list[(t-1)*q + t2-1]*fiA_list_list[(t1-1)*q + t3 - 1]))
                print "index: ",[t,t1,t2,t3]," DI_sub: ", DI_sub
                DI_sub_list.append(DI_sub)

        DI = sum(DI_sub_list)
        DI_list.append(DI)
        print "DI List: ", np.shape(DI_list)

    DI_list_list.append(DI_list)

print "DI List List: ", DI_list_list
print "DI Shape: ", np.shape(DI_list_list)
np.savetxt('stocks_DI_scores', DI_list_list)
            

                
                    
                
                
