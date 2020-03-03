# -*-  coding: utf-8 -*-
# @Time: 2020/1/9  16:07
# @Author: LiuJ
# @File: 1_read_file.py


import pandas as pd


str1 = "./GSE3431/GSE3431_gene_symbol.csv"
str2 = "./GSE3431/GSE3431_series_matrix.csv"
str3 = "./GSE3431/GSE3431_geneSyb_filter.csv"
str4 = "./GSE3431/GSE3431_serMtx_filter.csv"

df1 = pd.read_csv(str1, index_col=False)
df2 = pd.read_csv(str2, index_col=False)

df1.drop_duplicates()
df2.drop_duplicates()
cnt = 0
list = []
sch = []
for row in df1.itertuples():
    if getattr(row, 'Gene_Symbol') == '---':
        list.append(getattr(row, 'Index'))
    elif getattr(row, 'Gene_Symbol').find('///') != -1:
        fstSyb = getattr(row, 'Gene_Symbol').split('///')[0]
        sch.append([getattr(row, 'Index'), fstSyb])

# print('the sum number is: %d' % cnt)
# print(len(list))
df1 = df1.drop(index=list)  # Delete those Prob_Set whose geneSymbols are '---'

nn = len(sch)
for i in range(nn):
    df1.loc[sch[i][0], 'Gene_Symbol'] = sch[i][1]  # Select the 1st geneSymbol of multiple geneSymbols
print('the number of gene expressions after drop ---: %d'%(len(df1.index.to_list())))

col_nam = df2.columns.tolist()
col_nam.insert(1, 'Gene_Symbol')
df2 = df2.reindex(columns=col_nam)
df2['Gene_Symbol'] = 'NaN'  # A new column for geneSymbol

list2 = []
dic={}
for row1 in df1.itertuples():
    dic[getattr(row1,'Probe_Set_ID')]=getattr(row1,'Gene_Symbol')

for row2 in df2.itertuples():
    if dic.get(getattr(row2, 'ID_REF'),0)!=0:
        a = [getattr(row2, 'Index'), dic[getattr(row2, 'ID_REF')]]
        list2.append(a)
print(len(list2))
for i in range(len(list2)):
    df2.loc[list2[i][0], 'Gene_Symbol'] = list2[i][1]  # Find and insert geneSymbol in df2


cnt = 0
list3=[]
for row in df2.itertuples():
    if getattr(row, 'Gene_Symbol') == 'NaN':
        cnt += 1
        list3.append(getattr(row,'Index'))
print('the number of NaN is: %d'%cnt)
df2 = df2.drop(index=list3)

df1.to_csv(path_or_buf=str3, sep=',', header=True, index=False)
df2.to_csv(path_or_buf=str4, sep=',', header=True, index=False)
