import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import numexpr as npr
from itertools import chain
from matplotlib_venn import venn2
from math import pi

# 读取数据,tsv转为csv格式
STRING = pd.read_csv('E:\Python training\data\downstream\string_node_degrees.tsv', sep='\t')
STRING.to_csv('example.csv', index=False)
#print(STRING)

#提取数据
table1 = STRING.iloc[:,0:1]
#print(table1)
tab1 = np.asarray(table1)  # 转成数组形式
t1 = STRING['#node'].tolist()   #list形式
print(t1)

call = ['GeneSymbol','foldchange','pvalue']

#读取转录组数据
MC1_DE = pd.read_excel('E:\Python training\data\Con-VS-MC1_DE.xlsx',engine='openpyxl')
MC1_gene = pd.DataFrame(MC1_DE,columns =['GeneSymbol','foldchange','pvalue'])
#print(MC1_gene)
MC03_DE = pd.read_excel('E:\Python training\data\Con-VS-MC0.3_DE.xlsx',engine='openpyxl')
MC03_gene = pd.DataFrame(MC03_DE,columns =['GeneSymbol','foldchange','pvalue'])
#print(MC03_gene)
MC0_DE = pd.read_excel('E:\Python training\data\Con-VS-MC0_DE.xlsx',engine='openpyxl')
MC0_gene = pd.DataFrame(MC0_DE,columns =['GeneSymbol','foldchange','pvalue'])
#print(MC0_gene)
MCM1_DE = pd.read_excel('E:\Python training\data\MG-MG+C_M1_DE_anno.xlsx',engine='openpyxl')
MCM1_gene = pd.DataFrame(MCM1_DE,columns =['GeneSymbol','foldchange','pvalue'])
#print(MC1_gene)
MCM03_DE = pd.read_excel('E:\Python training\data\MG-MG+C_M0.3_DE_anno.xlsx',engine='openpyxl')
MCM03_gene = pd.DataFrame(MCM03_DE,columns =['GeneSymbol','foldchange','pvalue'])
#print(MC1_gene)
MCM1_C_DE = pd.read_excel('E:\Python training\data\Con-MG+C_M1_DE_anno.xlsx',engine='openpyxl')
MCM1_C_gene = pd.DataFrame(MCM1_C_DE,columns =['GeneSymbol','foldchange','pvalue'])
#print(MC1_gene)
MCM03_C_DE = pd.read_excel('E:\Python training\data\Con-MG+C_M0.3_DE_anno.xlsx',engine='openpyxl')
MCM03_C_gene = pd.DataFrame(MCM03_C_DE,columns =['GeneSymbol','foldchange','pvalue'])
#print(MC1_gene)

#MC1差异基因
MC1_r=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MC1_gene.loc[(MC1_gene['GeneSymbol'] == Name)&(MC1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC1_r.extend(FFF)
MC1_r=pd.DataFrame(columns=call,data=MC1_r)
#MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('STRING_MC1：'.format(MC1_r),MC1_r)

#MC0.3差异基因
MC03_r=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MC03_gene.loc[(MC03_gene['GeneSymbol'] == Name)&(MC03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC03_r.extend(FFF)
MC03_r=pd.DataFrame(columns=call,data=MC03_r)
#MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('STRING_MC03：'.format(MC03_r),MC03_r)

#MC0差异基因
MC0_r=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MC0_gene.loc[(MC0_gene['GeneSymbol'] == Name)&(MC0_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC0_r.extend(FFF)
MC0_r=pd.DataFrame(columns=call,data=MC0_r)
#MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('STRING_MC0：'.format(MC0_r),MC0_r)

#MCM1差异基因
MCM1_r=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MCM1_gene.loc[(MCM1_gene['GeneSymbol'] == Name)&(MCM1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM1_r.extend(FFF)
MCM1_r=pd.DataFrame(columns=call,data=MCM1_r)
#MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('STRING_MCM1：'.format(MCM1_r),MCM1_r)

#MCM03差异基因
MCM03_r=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MCM03_gene.loc[(MCM03_gene['GeneSymbol'] == Name)&(MCM03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM03_r.extend(FFF)
MCM03_r=pd.DataFrame(columns=call,data=MCM03_r)
#MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('STRING_MCM03：'.format(MCM03_r),MCM03_r)

#MCM1_C差异基因
MCM1_C_r=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MCM1_C_gene.loc[(MCM1_C_gene['GeneSymbol'] == Name)&(MCM1_C_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM1_C_r.extend(FFF)
MCM1_C_r=pd.DataFrame(columns=call,data=MCM1_C_r)
#MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('STRING_MCM1_C：'.format(MCM1_C_r),MCM1_C_r)

#MCM03_C差异基因
MCM03_C_r=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MCM03_C_gene.loc[(MCM03_C_gene['GeneSymbol'] == Name)&(MCM03_C_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM03_C_r.extend(FFF)
MCM03_C_r=pd.DataFrame(columns=call,data=MCM03_C_r)
#MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('STRING_MCM0.3_C：'.format(MCM03_C_r),MCM03_C_r)

print('____________________________________')

# 设置数据
df = pd.DataFrame({
    'group': [         'EPGN',     'EREG',  'ICAM1','ERBB2','ERBB3'],
    'MG0-(MG+C)/MG1': [  0,        0,         0.709850, 0,        0],
    'MG0-(MG+C)/MG0.3': [0,        0.704500,  0.647004, 0,       0],
    'Con-MG+C1': [15.446455,        2.156673,   0,        0,       0],
    'Con-(MG+C)/MG1': [13.839042,  1.639992,   0.798557, 0.772943, 0.504478],
    'Con-(MG+C)/MG0.3': [7.390744, 0,         0.729308, 0.825874, 0.554262],
    'Con-MG+C0': [       11.319941, 0,         0,        0,        0],

})

# 目标数量
categories = list(df)[1:]
N = len(categories)

# 角度
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

# 初始化
ax = plt.subplot(111, polar=True)

# 设置第一处
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)

# 添加背景信息
plt.xticks(angles[:-1], categories)
ax.set_rlabel_position(0)
plt.yticks([0, 1.25, 2.5], ["0", "1.25", "2.5"], color="grey", size=7)
plt.ylim(0, 2.5)

# 添加数据图

# 第一个
values1 = df.loc[0].drop('group').values.flatten().tolist()
values1 += values1[:1]
ax.plot(angles, values1, 'blue',linewidth=1.5, linestyle='solid', label="EPGN")
ax.fill(angles, values1, 'blue', alpha=0.1)

# 第二个
values2 = df.loc[1].drop('group').values.flatten().tolist()
values2 += values2[:1]
ax.plot(angles, values2, 'green',linewidth=1.5, linestyle='solid', label="EREG")
ax.fill(angles, values2, 'green', alpha=0.1)

# 第三个
values3 = df.loc[2].drop('group').values.flatten().tolist()
values3 += values3[:1]
ax.plot(angles, values3,  'red',linewidth=1.5, linestyle='solid', label="ICAM1")
ax.fill(angles, values3, 'red', alpha=0.1)

# 第四个
values4 = df.loc[3].drop('group').values.flatten().tolist()
values4 += values4[:1]
ax.plot(angles, values4,  'yellow',linewidth=1.5, linestyle='solid', label="ERBB2")
ax.fill(angles, values4, 'yellow', alpha=0.1)

# 第五个
values5 = df.loc[4].drop('group').values.flatten().tolist()
values5 += values5[:1]
ax.plot(angles, values5,'black', linewidth=1.5, linestyle='solid', label="ERBB3")
ax.fill(angles, values5, 'black', alpha=0.1)

# 添加图例
plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

# 显示
plt.show()
