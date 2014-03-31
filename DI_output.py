import numpy as np
#import get_sector
from pylab import *

symbol_list_all = np.loadtxt('symbol_list', "str")
sector_list_all = np.loadtxt('sector_list', "str")
industry_list_all = np.loadtxt('industry_list', "str", delimiter="!!!")


q=3

DI_matrix = np.loadtxt('stocks_DI_scores')
eij_matrix = np.loadtxt('stocks_eij')

DI_list = []
for i in range(len(DI_matrix)):
    for j in range(len(DI_matrix[0])):

        symbol_list_all_index_t = symbol_list_all.where(symbol_list[i].capitalize())
        symbol_list_all_index_t1= symbol_list_all.where(symbol_list[j].capitalize())

        sort_list = [DI_matrix[i][j], 
            symbol_list[i], 
            sector_list_all[symbol_list_all_index_t],
            industry_list_all[symbol_list_all_index_t],
            symbol_list[j], 
            sector_list_all[symbol_list_all_index_t1], 
            industry_list_all[symbol_list_all_index_t1]]

        sort_list = [DI_matrix[i][j], symbol_list[i], symbol_list[j]]

        DI_list.append(sort_list)

DI_list_sorted = sorted(DI_list, reverse=True)

for i in range(len(DI_list_sorted)):
    print DI_list_sorted[i]



        
