import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import numexpr as npr
from itertools import chain
from matplotlib_venn import venn2
from math import pi

# 读取数据
Jas_U = pd.read_csv('E:\Python training\AREG\Jas-UCSC.csv')
Pro_T = pd.read_csv('E:\Python training\AREG\Pro-Transfac.csv')
hTF_TJ = pd.read_csv('E:\Python training\AREG\hTFDB-Trans-Jas.csv')

#提取数据
table1 = Jas_U.iloc[:,3:4]
table2 = Pro_T.iloc[:,3:4]
table3 = hTF_TJ.iloc[:,3:4]
#table4 = GTRD.iloc[:,3:4]
# #print(table2) #print(table3)
tab1 = np.asarray(table1)  # 转成数组形式
tab2 = np.asarray(table2)  # 转成数组形式
tab3 = np.asarray(table3)  # 转成数组形式
#tab4 = np.asarray(table4)  # 转成数组形式
t1 = Jas_U['to_gene'].tolist()   #list形式
t2 = Pro_T['to_gene'].tolist()   #list形式
t3 = hTF_TJ['to_gene'].tolist()  #list形式

#求交集
a = np.intersect1d(tab1, tab3, assume_unique=True)
a_list = a.tolist()
#print('JASPAR-UCSC和hTFDB-TRANSFAC-JASPAR交集:'.format(a_list),a_list)
b = np.intersect1d(tab1, tab2, assume_unique=True)
b_list = b.tolist()
#print('JASPAR-UCSC和PROMO-TRANSFAC交集:'.format(b_list),b_list)
c = np.intersect1d(tab2, tab3, assume_unique=True)
c_list = c.tolist()
#print('PROMO-TRANSFAC和hTFDB-TRANSFAC-JASPAR交集:'.format(c_list),c_list)

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

call = ['GeneSymbol','foldchange','pvalue']
#交集的差异基因检索
#MC1差异基因筛选：a和b交集没有筛选出结果，c交集结果‘STAT4’
MC1_re=[]
for i in range(len(c_list)):
    Name = c_list[i]
    FF = MC1_gene.loc[(MC1_gene['GeneSymbol'] == Name)&(MC1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC1_re.extend(FFF)
MC1_re=pd.DataFrame(columns=call,data=MC1_re)
print('PROMO∩hTFDB-MC1_DE：'.format(MC1_re),MC1_re)
#MC0差异基因筛选：c交集结果‘STAT4’
MC0_re=[]
for i in range(len(c_list)):
    Name = c_list[i]
    FF = MC0_gene.loc[(MC0_gene['GeneSymbol'] == Name)&(MC0_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC0_re.extend(FFF)
MC0_re=pd.DataFrame(columns=call,data=MC0_re)
print('PROMO∩hTFDB-MC0_DE：'.format(MC0_re),MC0_re)
#MC0.3差异基因筛选均没有结果
MC03_re=[]
for i in range(len(c_list)):
    Name = c_list[i]
    FF = MC03_gene.loc[(MC03_gene['GeneSymbol'] == Name)&(MC03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC03_re.extend(FFF)
MC03_re=pd.DataFrame(columns=call,data=MC03_re)
print('PROMO∩hTFDB-MC0.3_DE：'.format(MC03_re),MC03_re)

#MCM0.3差异基因筛选：a和b交集没有筛选出结果，c交集结果‘MAZ,STAT4’
MCM03c_re=[]
for i in range(len(c_list)):
    Name = c_list[i]
    FF = MCM03_gene.loc[(MCM03_gene['GeneSymbol'] == Name)&(MCM03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:,0:3].values
    #print(FFF)
    MCM03c_re.extend(FFF)
MCM03c_re=pd.DataFrame(columns=call,data=MCM03c_re)
print('PROMO∩hTFDB-MCM0.3c_DE：'.format(MCM03c_re),MCM03c_re)

#MCM1差异基因筛选：a和b交集没有筛选出结果，c交集结果‘ATF,MAZ,YY1’
MCM1c_re=[]
for i in range(len(c_list)):
    Name = c_list[i]
    FF = MCM1_gene.loc[(MCM1_gene['GeneSymbol'] == Name)&(MCM1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:,0:3].values
   # print(FFF)
    MCM1c_re.extend(FFF)
MCM1c_re=pd.DataFrame(columns=call,data=MCM1c_re)
print('PROMO∩hTFDB-MCM1c_DE：'.format(MCM1c_re),MCM1c_re)

#MCM0.3-Control差异基因筛选：a和b交集没有筛选出结果，c交集结果‘MAZ,STAT4’
MCM03_Cc_re=[]
for i in range(len(c_list)):
    Name = c_list[i]
    FF = MCM03_C_gene.loc[(MCM03_C_gene['GeneSymbol'] == Name)&(MCM03_C_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:,0:3].values
    #print(FFF)
    MCM03_Cc_re.extend(FFF)
MCM03_Cc_re=pd.DataFrame(columns=call,data=MCM03_Cc_re)
print('PROMO∩hTFDB-MCM0.3_Cc_DE：'.format(MCM03_Cc_re),MCM03_Cc_re)

#MCM1-Control差异基因筛选：a和b交集没有筛选出结果，c交集结果‘ATF,MAZ,YY1’
MCM1_Cc_re=[]
for i in range(len(c_list)):
    Name = c_list[i]
    FF = MCM1_C_gene.loc[(MCM1_C_gene['GeneSymbol'] == Name)&(MCM1_C_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:,0:3].values
   # print(FFF)
    MCM1_Cc_re.extend(FFF)
MCM1_Cc_re=pd.DataFrame(columns=call,data=MCM1_Cc_re)
print('PROMO∩hTFDB-MCM1_Cc_DE：'.format(MCM1_Cc_re),MCM1_Cc_re)

print('____________________________________')

# 设置数据
#df = pd.DataFrame({
#    'group': ['A', 'B', 'C', 'D'],
#    'var1': [38, 1.5, 30, 4],
#    'var2': [29, 10, 9, 34],
#    'var3': [8, 39, 23, 24],
#    'var4': [7, 31, 33, 14],
#    'var5': [28, 15, 32, 14]
#})
df = pd.DataFrame({
    'group': [         'STAT4',   'VDR',   'MAZ', 'ATF3',  'YY1',  'STAT3','USF2'],
   'Con-MG+C0': [        2.15008,   0,         0,      0,        0,        0,       0],
   'Con-MG+C1': [        3.361371,2.725614,    0,      0,        0,        0,       0],
   'MG0-(MG+C)/MG0.3': [0.535794,0.657309,1.184437,   0,        0,        0,       0],
   'Con-(MG+C)/MG1': [  1.946288, 1.25395, 0.748229, 0.748569, 0.844321, 0.790886, 0.762300],
    'Con-(MG+C)/MG0.3': [0,       1.257652, 0.812304, 0,       0.894249, 0.814439, 0.806312],
    'MG0-(MG+C)/MG1': [  0,        0.656085, 1.092083,0.819491, 0.904841,    0,        0],
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
plt.yticks([0, 1.7, 3.4], ["0", "1.7", "3.4"], color="grey", size=7)
plt.ylim(0, 3.4)

# 添加数据图

# 第一个
values = df.loc[0].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, 'blue',linewidth=1.5, linestyle='solid', label="STAT4")
ax.fill(angles, values, 'blue', alpha=0.1)

# 第二个
values = df.loc[1].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, 'orange',linewidth=1.5, linestyle='solid', label="VDR")
ax.fill(angles, values, 'orange', alpha=0.1)

# 第三个
values = df.loc[2].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values,'red',  linewidth=1.5, linestyle='solid', label="MAZ")
ax.fill(angles, values, 'red', alpha=0.1)

# 第四个
values = df.loc[3].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, 'yellow',linewidth=1.5, linestyle='solid', label="ATF3")
ax.fill(angles, values, 'yellow', alpha=0.1)

# 第五个
values = df.loc[4].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, 'black',linewidth=1.5, linestyle='solid', label="YY1")
ax.fill(angles, values, 'black', alpha=0.1)

# 第六个
values = df.loc[5].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, 'purple',linewidth=1.5, linestyle='solid', label="STAT3")
ax.fill(angles, values, 'purple', alpha=0.1)

# 第七个
values = df.loc[6].drop('group').values.flatten().tolist()
values += values[:1]
ax.plot(angles, values, 'green',linewidth=1.5, linestyle='solid', label="USF2")
ax.fill(angles, values, 'green', alpha=0.1)


# 添加图例
plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

# 显示
plt.show()