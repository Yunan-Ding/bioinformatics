import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import numexpr as npr
from itertools import chain
from matplotlib_venn import venn2


# 读取数据
Jas_U = pd.read_csv('E:\Python training\AREG\Jas-UCSC.csv')
Pro_T = pd.read_csv('E:\Python training\AREG\Pro-Transfac.csv')
hTF_TJ = pd.read_csv('E:\Python training\AREG\hTFDB-Trans-Jas.csv')
#GTRD = pd.read_csv('E:\Python training\AREG\GTRD1.csv')
#print(Jas_U) #print(Pro_T) #print(hTF_TJ)

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
print('JASPAR-UCSC和hTFDB-TRANSFAC-JASPAR交集:'.format(a_list),a_list)
b = np.intersect1d(tab1, tab2, assume_unique=True)
b_list = b.tolist()
print('JASPAR-UCSC和PROMO-TRANSFAC交集:'.format(b_list),b_list)
c = np.intersect1d(tab2, tab3, assume_unique=True)
c_list = c.tolist()
print('PROMO-TRANSFAC和hTFDB-TRANSFAC-JASPAR交集:'.format(c_list),c_list)
#d = np.intersect1d(tab1, tab4, assume_unique=True)
#d_list = d.tolist()
#print('JASPAR-UCSC和GTRD交集:'.format(d),d)
#e = np.intersect1d(tab2, tab4, assume_unique=True)
#e_list = e.tolist()
#print('PROMO-TRANSFAC和GTRD交集:'.format(e),e)
#f = np.intersect1d(tab3, tab4, assume_unique=True)
#f_list = f.tolist()
#countf = pd.value_counts(f_list)
#print(countf)
#print('hTFDB-TRANSFAC-JASPAR和GTRD交集:'.format(f),f)

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

print('-------------------------------------------------------')
#原始数据库的差异基因检索

#MC1差异基因
MC1_r2=[]
for i in range(len(t2)):
    Name = t2[i]
    FF = MC1_gene.loc[(MC1_gene['GeneSymbol'] == Name)&(MC1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC1_r2.extend(FFF)
MC1_r2=pd.DataFrame(columns=call,data=MC1_r2)
MC1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC1_DE.csv',index=False)
print('PRO-MC1_DE：'.format(MC1_r2),MC1_r2)

MC1_r3=[]
for i in range(len(t3)):
    Name = t3[i]
    FF = MC1_gene.loc[(MC1_gene['GeneSymbol'] == Name)&(MC1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC1_r3.extend(FFF)

#MC1_all = list(chain(MC1_r1, MC1_r2,MC1_r3))
MC1_r3=pd.DataFrame(columns=call,data=MC1_r3)
#MC1_all=pd.DataFrame(columns=call,data=MC1_all)
MC1_r3.to_csv('E:\Python training\AREG\everydatabase\hTFDB-MC1_DE.csv',index=False)
#MC1_all.to_csv('E:\Python training\AREG\MC1_ALL.csv',index=False)
print('hTFDB-MC1_DE：'.format(MC1_r3),MC1_r3)

#MC0.3差异基因
MC03_r3=[]
for i in range(len(t3)):
    Name = t3[i]
    FF = MC03_gene.loc[(MC03_gene['GeneSymbol'] == Name)&(MC03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC03_r3.extend(FFF)
MC03_r3=pd.DataFrame(columns=call,data=MC03_r3)
#MC03_all = list(chain(MC03_r1, MC03_r2,MC03_r3))
#MC03_all=pd.DataFrame(columns=call,data=MC03_all)
MC03_r3.to_csv('E:\Python training\AREG\everydatabase\hTFDB-MC0.3_DE.csv',index=False)
#MC03_all.to_csv('E:\Python training\AREG\MC03_ALL.csv',index=False)
print('hTFDB-MC0.3_DE：'.format(MC03_r3),MC03_r3)

#MC0差异基因
MC0_r2=[]
for i in range(len(t2)):
    Name = t2[i]
    FF = MC0_gene.loc[(MC0_gene['GeneSymbol'] == Name)&(MC0_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC0_r2.extend(FFF)
MC0_r2=pd.DataFrame(columns=call,data=MC0_r2)
MC0_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC0_DE.csv',index=False)
print('PRO-MC0_DE：'.format(MC0_r2),MC0_r2)

MC0_r3=[]
for i in range(len(t3)):
    Name = t3[i]
    FF = MC0_gene.loc[(MC0_gene['GeneSymbol'] == Name)&(MC0_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MC0_r3.extend(FFF)
MC0_r3=pd.DataFrame(columns=call,data=MC0_r3)
#MC0_all = list(chain(MC0_r1, MC0_r2,MC0_r3))
#MC0_all=pd.DataFrame(columns=call,data=MC0_all)
MC0_r3.to_csv('E:\Python training\AREG\everydatabase\hTFDB-MC0_DE.csv',index=False)
#MC0_all.to_csv('E:\Python training\AREG\MC0_ALL.csv',index=False)
print('hTFDB-MC0_DE：'.format(MC0_r3),MC0_r3)

#MCM0.3差异基因
MCM03_r1=[]
for i in range(len(t1)):
    Name = t1[i]
    FF = MCM03_gene.loc[(MCM03_gene['GeneSymbol'] == Name)&(MCM03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM03_r1.extend(FFF)
MCM03_r1=pd.DataFrame(columns=call,data=MCM03_r1)
MCM03_r1.to_csv('E:\Python training\AREG\everydatabase\JAS-MC0.3_DE.csv',index=False)
print('JAS-MC0.3_DE：'.format(MCM03_r1),MCM03_r1)

MCM03_r2=[]
for i in range(len(t2)):
    Name = t2[i]
    FF = MCM03_gene.loc[(MCM03_gene['GeneSymbol'] == Name)&(MCM03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM03_r2.extend(FFF)
MCM03_r2=pd.DataFrame(columns=call,data=MCM03_r2)
MCM03_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MC0.3_DE.csv',index=False)
print('PRO-MC0.3_DE：'.format(MCM03_r2),MCM03_r2)

MCM03_r3=[]
for i in range(len(t3)):
    Name = t3[i]
    FF = MCM03_gene.loc[(MCM03_gene['GeneSymbol'] == Name)&(MCM03_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM03_r3.extend(FFF)
MCM03_r3=pd.DataFrame(columns=call,data=MCM03_r3)
MCM03_r3.to_csv('E:\Python training\AREG\everydatabase\hTFDB-MCM0.3_DE.csv',index=False)
print('hTFDB-MCM0.3_DE：'.format(MCM03_r3),MCM03_r3)
#MCM03_all = list(chain(MCM03_r1, MCM03_r2,MCM03_r3))
#MCM03_all=pd.DataFrame(columns=call,data=MCM03_all)
#MCM03_all.to_csv('E:\Python training\AREG\MCM03_ALL.csv',index=False)

#MCM1差异基因
MCM1_r2=[]
for i in range(len(t2)):
    Name = t2[i]
    FF = MCM1_gene.loc[(MCM1_gene['GeneSymbol'] == Name)&(MCM1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM1_r2.extend(FFF)
MCM1_r2=pd.DataFrame(columns=call,data=MCM1_r2)
MCM1_r2.to_csv('E:\Python training\AREG\everydatabase\PRO-MCM1_DE.csv',index=False)
print('PRO-MCM1_DE：'.format(MCM1_r2),MCM1_r2)

MCM1_r3=[]
for i in range(len(t3)):
    Name = t3[i]
    FF = MCM1_gene.loc[(MCM1_gene['GeneSymbol'] == Name)&(MCM1_gene['pvalue'] <0.05)]
    FFF = FF.iloc[:, 0:3].values
    MCM1_r3.extend(FFF)
MCM1_r3=pd.DataFrame(columns=call,data=MCM1_r3)
MCM1_r3.to_csv('E:\Python training\AREG\everydatabase\hTFDB-MCM1_DE.csv',index=False)
print('hTFDB-MCM1_DE：'.format(MCM1_r3),MCM1_r3)
#MCM1_all = list(chain(MCM1_r2,MCM1_r3))
#MCM1_all=pd.DataFrame(columns=call,data=MCM1_all)
#MCM1_all.to_csv('E:\Python training\AREG\MCM1_ALL.csv',index=False)

