# -*-  coding: utf-8 -*-
# @Time: 2020/2/7  17:07
# @Author: LiuJ
# @File: 2_filter_gene.py

import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib
import matplotlib.pyplot as plt
from plotfile import myplot, myplot_2,myplot_3

str1 = './GSE3431/GSE3431_serMtx_filter.csv'
str2 = './GSE3431/GSE3431_gene_filtered.csv'
str3 = './GSE3431/GSE3431_gene_filtered_binary.csv'

p = 6  # the len of step
M = 36  # number of time stamps of GSE3431
p_thshd = 0.01
mean_thshd = 0.5
kk = 3


df1 = pd.read_csv(filepath_or_buffer=str1, index_col=['ID_REF', 'Gene_Symbol'])


i = 0
j = 0
sum_list=[]
for i in range(len(df1.index.to_list())):
    sum_list.append(True)

for p in range(6,7):
    print('#Round %d start'%(p-5))
    list = []
    V=M - p
    Y = np.zeros(shape=(V, 1), dtype=np.float)
    X = np.zeros(shape=(V, p + 1), dtype=np.float)
    beta = np.zeros(shape=(p + 1, 1), dtype=np.float)
    X[:, 0] = 1
    for row in df1.itertuples():  # traverse all the genes
        for i in range(p, M):
            for j in range(i - (p - 1), i + 1):
                X[i - p, j - (i - p)] = getattr(row, 'GSM77' + str(j + 297))
            Y[i - p] = getattr(row, 'GSM77' + str(j + 297 + 1))  # for one gene, compose X and Y

        beta = np.dot(np.dot(np.linalg.inv(np.dot(X.T, X)), X.T), Y)
        sigma = np.power(np.linalg.norm(Y - np.dot(X, beta), axis=0, keepdims=True), 2) / V

        beta_c = np.sum(Y, axis=0) / V
        sigma_c = np.dot((Y - beta_c).T, (Y - beta_c)) / V

        F = (sigma_c / sigma - 1) * (M - 2 * p - 1) / p
        p_value = st.f.sf(abs(F), p, M - 2 * p - 1)
        # print('the p_value of'+str(getattr(row, 'Index'))+'is: '+str(p_value.tolist()))
        if p_value < p_thshd:
            list.append('time_dependent')
        else:
            list.append('time_independent')
    for i in range(len(list)):
        if list[i] is 'time_dependent':
            sum_list[i] = sum_list[i] and False
        else:
            sum_list[i] = sum_list[i] and True
    print('#Round %d complete, p=%d'%(p-5,p))
list=[]
for ll in sum_list:
    if ll is False:
        list.append( 'time_dependent')
    else:
        list.append( 'time_independent')

cnt = 0
for k in list:
    if k is 'time_dependent':
        cnt += 1

print('The Thresholds Set:\np: [6:%d]\tp_threshold: %.3f\t mean_threshold: %.3f\tk: %.3f' % (p, p_thshd, mean_thshd, kk))
print('number_time_dependent: ' + str(cnt) + '  number_time_independent: ' + str(len(list) - cnt))
print('percent_time_dependent: %.3f%%     percent_time_independent: %.3f%%' % (
    cnt / len(list) * 100, (1 - cnt / len(list)) * 100))

df1.insert(0, column='Independence', value=list)
df2 = pd.DataFrame(data=df1)
list2 = []
list3 = []
for row in df1.itertuples():
    if getattr(row, 'Independence') is 'time_dependent':
        list2.append(getattr(row, 'Index'))
    else:
        list3.append(getattr(row, 'Index'))
df1 = df1.drop(list2, axis=0)  # time_independent genes left
df2 = df2.drop(list3, axis=0)  # time_dependent genes left

list4 = df1.mean(axis=1).to_list()
df1.insert(0, column='Mean', value=list4)
df1 = df1.sort_values(by='Mean', ascending=True)
list4 = df2.mean(axis=1).to_list()
df2.insert(0, column='Mean', value=list4)
df2 = df2.sort_values(by='Mean', ascending=True)

list4 = []
for row in df1.itertuples():
    list4.append(getattr(row, 'Index'))
    if getattr(row, 'Mean') > mean_thshd:
        break
print('the number of filter genes: %d  the number of genes left: %d' % (len(list4), len(list) - len(list4)))
print('percent_low_mean_time_independent_genes/time_independent_genes: %.4f%%' % (len(list4) / (len(list) - cnt) * 100))
print('percent_low_mean_time_independent_genes/all_genes: %.4f%%' % (len(list4) / (len(list)) * 100))

df1 = df1.drop(list4, axis=0)
df1 = pd.concat([df1, df2])

one_cycle = np.zeros(shape=(len(list) - len(list4), 12), dtype=np.float)
# one_cycle = np.zeros(shape=(len(list) , 12), dtype=np.float)
temp = np.zeros(shape=(12, 1), dtype=np.float)
cnt = 0
for row in df1.itertuples():
    temp = np.zeros(shape=(12, 1), dtype=np.float)
    for i in range(36):
        if 12 <= i < 24:
            j = i % 12
        elif i >= 24:
            j = i % 24
        else:
            j = i
        temp[j] += getattr(row, 'GSM77' + str(298 + i))
    temp /= 3
    one_cycle[cnt, :] = temp.reshape(12)
    cnt += 1

mu = one_cycle.mean(axis=1)
std = one_cycle.std(axis=1,ddof=1)
for i in range(len(std)):
    if i>5000:
        std[i]/=mu[i]

active_thshd = mu + kk * std * (1 - 1/(1+std*std))
active_thshd = active_thshd.tolist()


df1.insert(0, column='Active_threshold', value=active_thshd)
df1.insert(0,column='Sigma',value=std)
df1.to_csv(path_or_buf=str2, sep=',', header=True, index=True)
bina = np.zeros(shape=(len(df1.index.to_list()), 36), dtype=np.float)

df1=df1.sort_values(by='Mean',ascending=True)
cnt = 0
for row in df1.itertuples():
    for i in range(36):
        if getattr(row, 'Active_threshold') > getattr(row, 'GSM77' + str(i + 298)):
            bina[cnt, i] = 0
        else:
            bina[cnt, i] = 1
    cnt += 1


count=0
liy=[]
for i in range(bina.shape[0]):
    count=0
    for j in range(bina.shape[1]):
        if bina[i,j]==1:
            count+=1
    liy.append(count)
lix=df1.loc[:,'Mean'].to_list()
myplot_2(lix,liy)

num1=0
sum=0
list4=[]
for j in range(bina.shape[1]):
    num1=0
    for i in range(bina.shape[0]):
        if bina[i,j] == 1:
            num1+=1
    print('active genes in time step %d is: %d/4031'%(j+1,num1))
    list4.append(num1)
    sum+=num1
print('average number of nodes is: %d'%(sum/36))

df1 = df1.drop(['Mean', 'Active_threshold', 'Independence','Sigma'], axis=1)
col_list = df1.columns.to_list()
df1.loc[:, col_list] = bina
df1.to_csv(path_or_buf=str3, sep=',', header=True, index=True)


