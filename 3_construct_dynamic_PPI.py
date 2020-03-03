# -*-  coding: utf-8 -*-
# @Time: 2020/2/15  15:47
# @Author: LiuJ
# @File: 3_construct_dynamic_PPI.py

import numpy as np
import pandas as pd

str1 = './yeast_PPI/species_yeast1.csv'
str2 = './yeast_PPI/translate.csv'
str3 = './GSE3431/GSE3431_gene_filtered_binary.csv'
str4 = './GSE3431/GSE3431_geneSyb_filter.csv'

df1 = pd.read_csv(str1)
df2 = pd.read_csv(str2)
df4 = pd.read_csv(str4)
list1 = []
list2 = []
dic = {}
for row in df2.itertuples():
    dic[getattr(row, 'From')] = getattr(row, 'To')
for row in df1.itertuples():
    if dic.get(getattr(row, 'Identify_A'), 0) != 0:
        list1.append(dic[getattr(row, 'Identify_A')])
    else:
        list1.append('NaN')
    if dic.get(getattr(row, 'Identify_B'), 0) != 0:
        list2.append(dic[getattr(row, 'Identify_B')])
    else:
        list2.append('NaN')
df1.loc[:, 'Identify_A'] = list1
df1.loc[:, 'Identify_B'] = list2
print('the number of PPIs before filtered: ' + str(len(list1)))
list1 = []
for row in df1.itertuples():
    if getattr(row, 'Identify_A') is 'NaN' or getattr(row, 'Identify_B') is 'NaN':
        list1.append(getattr(row, 'Index'))
df1 = df1.drop(list1)
print('the number of PPIs after filtering NaN: ' + str(len(df1.index.to_list())))

df3 = pd.read_csv(str3)
sum=0
for i in range(36):
    list1 = []
    dic = {}
    for row in df3.itertuples():
        dic[getattr(row, 'Gene_Symbol')] = getattr(row, 'GSM77' + str(298 + i))
    for row1 in df1.itertuples():
        if dic.get(getattr(row1, 'Identify_A'), -1) == 1 and dic.get(getattr(row1, 'Identify_B'), -1) == 1:
            list1.append((getattr(row1, 'Identify_A'), getattr(row1, 'Identify_B')))
    print('number of PPIs in time step %d: %d' % (i + 1, len(list1)))
    sum+=len(list1)
print('average number of edges are: %d'%(sum/36))

    # for ll in list1:
    #     print(ll)
    # print('\n')
dic = {}
cnt = 0
for row in df2.itertuples():
    dic[getattr(row, 'To')] = getattr(row, 'Index')
for row in df4.itertuples():
    if dic.get(getattr(row, 'Gene_Symbol'), -1) != -1:
        cnt += 1

print('The number of proteins covered: %d   the percent of cover: %.4f%%' % (cnt, cnt / 5053 * 100))
