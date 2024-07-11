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
    'group': [         'EPGN',     'EREG', 'ALDOC', 'ICAM1','EGFR','ERBB2','ERBB3'],
   'MG0-(MG+C)/MG0.3': [0,         0.704500,   0,     0.647004,   0,       0,       0],
   'Con-(MG+C)/MG1': [  13.839042, 1.639992,   0,     0.798557, 1.097102,0.772943,0.504478],
   'Con-(MG+C)/MG0.3': [7.390744,  0,      0.812304,  0.729308 ,  0,    0.825874, 0.554262],
   'MG0-(MG+C)/MG1': [   0,        0,          0,      0.709850,   0,       0,       0],
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
plt.yticks([0, 1, 2], ["0", "1", "2"], color="grey", size=7)
plt.ylim(0, 2)

# 添加数据图

# 第一个
values = df.loc[0].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label="EPGN")
ax.fill(angles, values, 'blue', alpha=0.1)

# 第二个
values = df.loc[1].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label="EREG")
ax.fill(angles, values, 'r', alpha=0.1)

# 第三个
values = df.loc[2].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label="ALDOC")
ax.fill(angles, values, 'red', alpha=0.1)

# 第四个
values = df.loc[3].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label="ICAM1")
ax.fill(angles, values, 'yellow', alpha=0.1)

# 第五个
values = df.loc[4].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label="EGFR")
ax.fill(angles, values, 'black', alpha=0.1)

# 第六个
values = df.loc[5].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label="ERBB2")
ax.fill(angles, values, 'purple', alpha=0.1)

# 第七个
values = df.loc[6].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid', label="ERBB3")
ax.fill(angles, values, 'green', alpha=0.1)


# 添加图例
plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

# 显示
plt.show()
