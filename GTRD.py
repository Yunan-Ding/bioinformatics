import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import numexpr as npr

# 读取数据
df = pd.read_excel('E:\Python training\AREG\GTRD.xlsx',engine='openpyxl')
#print(df)

# 按照转录因子名字出现次数进行计数
name = df.groupby('Type').size()
name_num = df.groupby('Type').count()
# print(name_num)
count = pd.Series(name.index, index=name.values)
count_num = np.asarray(count)  # 转成数组形式
list_count_num = count_num.tolist()  # 数组转化为list
#print('TFDB出现计数结果'.format(list_count_num),list_count_num)

# 取指定列
##建一个表格
result = name_num.iloc[:, :1]
# print(result)
# 表格新建一列
result['weight'] = name_num.iloc[:, :1]
print(result)
result['from_gene'] = 'AREG'
print(result)
result['to_gene'] = count_num
print(result)
print('整理后的表格'.format(result),result)

# 创建一个空的无向图/创建一个空的加权图
G = nx.Graph()
# 添加节点
nodes = result['to_gene']
edge_data = result
for node in nodes:
    G.add_node(node)

# 添加带权边
for edge in edge_data.itertuples():
    src = edge.from_gene
    dst = edge.to_gene
    weight = edge.weight
    G.add_edge(src, dst, weight=weight)

# 绘制网络图

nx.draw(G, with_labels=True,edge_color='grey',node_color='chartreuse',font_size =4,
         node_size=100,
        width=[float(d['weight'] * 0.1) for (u, v, d) in G.edges(data=True)])

plt.show()

result.to_csv('E:\Python training\AREG\GTRD1.csv',index=False)