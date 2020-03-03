# -*-  coding: utf-8 -*-
# @Time: 2020/3/3  11:40
# @Author: LiuJ
# @File: 4_plotfile.py

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from scipy.stats import norm
import numpy as np
# 从pyplot导入MultipleLocator类，这个类用于设置刻度间隔

def myplot(liy6):
    fig=plt.figure()
    ax2=fig.add_subplot(111)
    x_major_locator = MultipleLocator(2)
    _,bins,_=ax2.hist(liy6, bins=100, density=1, facecolor="maroon", edgecolor="black", alpha=0.7,label='mean=[16,etc]')
    mu = np.mean(np.array(liy6))
    std = np.std(np.array(liy6))
    y = norm.pdf(bins, mu, std)  # 拟合一条最佳正态分布曲线
    ax2.plot(bins, y, 'b--')  # 绘制y的曲线
    ax2.xaxis.set_major_locator(x_major_locator)
    ax2.set_xlim(6, 30)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Frequency', fontsize=20)
    ax2.set_xlabel('Expression Value',fontsize=20)
    ax2.legend(loc='upper center',fontsize=15)
    plt.show()

def myplot_2(lix,liy):
    fig = plt.figure()

    ax1 = fig.add_subplot(111)
    ax1.bar(lix,liy,edgecolor='cyan',facecolor='black',width=0.5)
    ax1.set_ylabel('Active numbers', fontsize=20)
    ax1.set_title('Mean-active_number Curve', fontsize=24)
    ax1.set_xlabel('Mean values', fontsize=20)
    ax1.set_ylim(0,36)


    x_major_locator = MultipleLocator(2)
    # ax2 = ax1.twinx()  # this is the important function
    # ax1.hist(lix, bins=100, normed=0, facecolor="deepskyblue", edgecolor="black", alpha=0.7)
    # ax1.xaxis.set_major_locator(x_major_locator)
    ax1.set_xlim(0,60)
    # ax1.set_ylim(0,2000)
    # ax1.set_ylabel('Frequency',fontsize=20)
    plt.show()

def myplot_3(liy):
    fig =plt.figure()
    ax=fig.add_subplot(111)
    lix=[x+1 for x in range(36)]
    ax.plot(lix,liy,':',linewidth=2,color='deeppink',marker='o')
    ax.set_xlabel('Time',fontsize=20)
    ax.set_ylabel('Number of nodes',fontsize=20)
    ax.set_xlim(0,37)
    ax.legend(loc='upper left',fontsize=15)

    plt.show()