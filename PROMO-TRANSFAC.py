import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import numexpr as npr

# 读取数据
df = pd.read_excel('E:\Python training\AREG\PROMO-TRANSFAC.xlsx',engine='openpyxl')
#print(df)

# 按照转录因子名字出现次数进行计数
name = df.groupby('Name').size()
name_num = df.groupby('Name').count()
# print(name_num)
count = pd.Series(name.index, index=name.values)
count_num = np.asarray(count)  # 转成数组形式
list_count_num = count_num.tolist()  # 数组转化为list
#print('TFDB出现计数结果'.format(list_count_num),list_count_num)

re=[]
for i in range(len(list_count_num)):
    a = list_count_num[i]
    FF = df.loc[(df['Name'] == a)]
    #print(FF)
    b = FF.iloc[:,1:3]
    score_max = b.iloc[:,1:2].max()
    c = [a,score_max]
    cc = score_max.values
    re.extend(cc)
print(re)

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
result['Score'] =re
print('整理后的表格'.format(result),result)

# 创建一个空的无向图/创建一个空的加权图
G = nx.Graph()
# 添加节点
nodes = result['to_gene']
result1= result.iloc[:3,:]
print(result1)
result2= result.iloc[4:,:]
result3 = pd.concat([result1, result2],axis=0)
print('调整权重后'.format(result3),result3)
edge_data = result3
for node in nodes:
    G.add_node(node)

# 添加带权边
for edge in edge_data.itertuples():
    src = edge.from_gene
    dst = edge.to_gene
    weight = edge.weight
    G.add_edge(src, dst, weight=weight)


# 绘制网络图
dd = dict(result.Score)
print(dd)

nx.draw(G, with_labels=True,edge_color='grey',node_color='skyblue',font_size =5,
        nodelist=dd, node_size=[v ** 2.1 for v in dd.values()],
        width=[float(d['weight'] * 0.1) for (u, v, d) in G.edges(data=True)])

plt.show()

result3.to_csv('E:\Python training\AREG\Pro-Transfac.csv',index=False)