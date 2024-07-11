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

#绘制venn
v=venn2([set(t2), set(t3)],set_labels=('PROMO','hTFDB'))
# 单独改变A的标签
#v.get_label_by_id('t2').set_text('PROMO')
#v.get_label_by_id('t3').set_text('hTFDB')
# 画圆，linestyle线条类型，linewith线宽，color线条颜色
#c = venn3_circles(subsets = (10, 8, 22, 6,9,4,2), linestyle='dashed', linewidth=1, color="grey")
plt.show()