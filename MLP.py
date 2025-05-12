# -*- coding: utf-8 -*-
# @Time    : 2022/2/27 16:46
# @Author  : gepeng
# @File    : MLP.py
#用MLP逼近 XOR 函数

from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from collections import Counter

# import neurolab as nl
# import matplotlib.pyplot as plt
# y = [0, 1, 1, 0]
# X = [[0, 0], [0, 1], [1, 0], [1,1]]
# classifier = MLPClassifier(solver='lbfgs', activation='logistic',
#                            hidden_layer_sizes=(2, ), random_state=20)
# classifier.fit(X, y)
# y_hat = classifier.predict(X)
# print(classifier.score(X, y))
# for i, p in enumerate(y_hat):
#     print('真实值: %s, 预测值: %s' % (y[i], p))

# #MLP模型的应用案例
# plt.rcParams['font.sans-serif'] = ['SimHei']  ## 用来正常显示中文标签
# plt.rcParams['axes.unicode_minus'] = False  ## 用来正常显示负号
# ## 虚构练数据
# min = -10 ##最小值
# max = 10 ##最大值
# N = 100 ##观测值个数
# np.random.seed(1)  ##设定伪随机数种子
# ## 生成自变量数据
# x = np.linspace(min, max, N)
# ## 生成结果变量数据
# y = 2*np.square(x) + 10
# y = y / np.linalg.norm(y)
# X = x.reshape(100, 1) ##训练的特征数据
# labels = y.reshape(100, 1) ##训练的标签数据
# ## 将训练数据可视化
# plt.figure()
# plt.scatter(X, labels)
# plt.xlabel('X轴')
# plt.ylabel('Y轴')
# plt.title('训练数据')
# plt.show()
# mlp = nl.net.newff([[min, max]], [10, 10, 10, 10, 10, 1])## 定义一个带有5个隐藏层、每个隐藏层包括10个神经元的MLP,最后输出层包括一个神经元
# mlp.trainf = nl.train.train_gd## 训练的算法为梯度下降法
# error = mlp.train(X, labels, epochs=1000, show=100, goal=0.01)## 训练mlp的误差
# yhat=mlp.sim(X)## 基于训练好的mlp进行预测
# ## 可视化训练误差过程
# plt.figure()
# plt.plot(error)
# plt.xlabel('Epoches时期数')
# plt.ylabel('Error训练误差')
# plt.title('训练的误差过程')
# plt.show()
# #可视化预测结果
# x2=np.linspace(min, max, N*2)
# y2=mlp.sim(x2.reshape(x2.size, 1)).reshape(x2.size)
# y3=yhat.reshape(N)
# plt.figure()
# plt.plot(x2, y2,'o')
# plt.plot(x, y,'.')
# plt.plot(x, y3,'p')
# plt.title('真实值 vs 预测值')
# plt.show()
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import matplotlib.colors as colour
d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene2kinds.csv',sep=',',nrows=1,index_col=0)
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgenecankinds.csv',sep=',',nrows=1,index_col=0)
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene3kinds.csv',sep=',',nrows=1,index_col=0)
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene.csv',sep=',',header=None,nrows=1,index_col=0)
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgenecankinds.csv',sep=',',nrows=1,index_col=0)
a=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维**样本
c=np.array(d)#1*48
c=c.ravel()#转置成（47，）格式,不可以用np.transpose()
f = Counter(c)#统计分类个数，最小分类不能小于2
print(f)
b=np.array(a)
b=np.transpose(b)#转置，得是np.array类型
print(b.shape)#(16383, 2540)
print(c.shape)#(16383,)
h=set(c)
h=list(h)
h.sort(reverse=False)
print(h)
X_train, X_test, y_train, y_test = train_test_split(b, c, random_state=1,stratify=c)#random_state为随机种子数，test_size默认为0.25,抽样分布参考y的分布
# mlp = MLPClassifier(solver='adam', activation='relu',alpha=1e-4,hidden_layer_sizes=(50,50), random_state=1,max_iter=100,verbose=10,learning_rate_init=.1)
# # 使用solver='lbfgs',准确率为79%，比较适合小(少于几千)数据集来说，且使用的是全训练集训练，比较消耗内存
# # 使用solver='adam'，准确率只有67%,但是默认solver ‘adam’在相对较大的数据集上效果比较好（几千个样本或者更多）
# # 使用solver='sgd'，准确率为98%，且每次训练都会分batch，消耗更小的内存
# mlp.fit(X_train,y_train)
# print(mlp.score(X_test,y_test))
# print(mlp.n_layers_)
# print(mlp.n_iter_)
# print(mlp.loss)
# print(mlp.out_activation_)

#加网格搜索
parameters = {"hidden_layer_sizes": [ (3000, )],
                             "solver": ['sgd'],
                             "max_iter": [200],
                             "verbose": [True],
                             }
mlp = MLPClassifier(random_state=0)
estimator = GridSearchCV(mlp,parameters, n_jobs=-1)
estimator.fit(X_train, y_train)
print(estimator.get_params().keys())
print(estimator.best_params_)
print('ACC:',estimator.best_score_)
y_predprob = estimator.predict_proba(X_test)
print('AUC:',roc_auc_score(y_test,y_predprob[:,1]))

# roc_ovr=roc_auc_score(y_test, y_predprob , multi_class='ovr')
# print('--AUC-ovr:',roc_ovr)
# roc_ovo=roc_auc_score(y_test, y_predprob, multi_class='ovo')
# print('--AUC-ovo:',roc_ovo)



def Find_Optimal_Cutoff(tper,fper,threshold):
    y=tper-fper
    Youden_index=np.argmax(y)
    optimal_threshold=threshold[Youden_index]
    point=[fper[Youden_index],tper[Youden_index]]
    return optimal_threshold,point
prob = y_predprob[:, 1]
fper, tper, thresholds = roc_curve(y_test, prob,pos_label='Tumors')
auc =auc(fper, tper)
optimal_th,optimal_point=Find_Optimal_Cutoff(tper,fper,thresholds)
plt.plot(fper, tper, color='blue', label='ROC curve (area = %0.2f)' % auc)
plt.plot([0, 1], [0, 1], color='green', linestyle='--')
plt.plot(optimal_point[0], optimal_point[1], marker='o', color='r')
plt.text(optimal_point[0],optimal_point[1], f'Threshold:{optimal_th:.2f}')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('MLP ROC Curve')
plt.legend(loc="lower right")
plt.savefig('./MLP.pdf',dpi=600)


#多分类
# n_classes = len(h)
# fpr = dict()#创建字典
# tpr = dict()
# roc_auc = dict()
# for i in range(n_classes):
#     fpr[i], tpr[i], _= roc_curve(y_test, y_predprob[:,i],pos_label=h[i])
#     roc_auc[i] = auc(fpr[i], tpr[i])
# # Compute macro-average ROC curve and ROC area（方法一）
# # First aggregate all false positive rates
# all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
#
# # Then interpolate all ROC curves at this points
# mean_tpr = np.zeros_like(all_fpr)
# for i in range(n_classes):
#     mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])
#
# # Finally average it and compute AUC
# mean_tpr /= n_classes
# fpr["macro"] = all_fpr
# tpr["macro"] = mean_tpr
# roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
#
# # Plot all ROC curves
# lw=2
# plt.figure()
# plt.plot(fpr["macro"], tpr["macro"],
#          label='macro-average ROC curve (area = {0:0.2f})'
#                ''.format(roc_auc["macro"]),
#          color='navy', linestyle=':', linewidth=4)
# from itertools import cycle
# # colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
# colors=colour.cnames
# for i, color in zip(range(n_classes), colors):
#     plt.plot(fpr[i], tpr[i], color=color, lw=lw,
#              label='{0} (area = {1:0.2f})'
#              ''.format(h[i], roc_auc[i]))
#
# plt.plot([0, 1], [0, 1], 'k--', lw=lw)
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('MLP ROC Curve to 3-class')
# plt.legend(loc='lower right',fontsize=7)
# plt.savefig('./MLP3kinds.pdf',dpi=600)











# #加网格搜索+两次分类
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene.csv',sep=',',header=None,nrows=1,index_col=0)
# a=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维47样本
# c=np.array(d)#1*48
# c=c.ravel()#转置成（47，）格式,不可以用np.transpose()
# # f = Counter(c)#统计分类个数，最小分类不能小于2
# # print(f)
# b=np.array(a)
# b=np.transpose(b)#转置，得是np.array类型
# # print(b.shape)#(23484, 2540)
# # print(c.shape)#(23484,)
# indices = np.arange(23484)
# X_train, X_test, y_train, y_test,idX_train,idX_test = train_test_split(b, c, indices,random_state=1,stratify=c)#random_state为随机种子数，test_size默认为0.25,抽样分布参考y的分布
# parameters = {"hidden_layer_sizes": [(3000,),(2500,)],
#                              "solver": ['sgd'],
#                              "max_iter": [300],
#                              "verbose": [True],
#                              }
# mlp = MLPClassifier(random_state=0)
# estimator = GridSearchCV(mlp,parameters, n_jobs=-1)
# estimator.fit(X_train, y_train)
# print(estimator.get_params().keys())
# print(estimator.best_params_)
# print('ACC:',estimator.best_score_)
# y_predprob = estimator.predict_proba(X_test)
# roc_ovr=roc_auc_score(y_test, y_predprob , multi_class='ovr')
# print('--AUC-ovr:',roc_ovr)
# roc_ovo=roc_auc_score(y_test, y_predprob, multi_class='ovo')
# print('--AUC-ovo:',roc_ovo)

# c=np.array(d)#1*48
# c=c.ravel()#转置成（47，）格式,不可以用np.transpose()
# # f = Counter(c)#统计分类个数，最小分类不能小于2
# # print(f)
# b=np.array(a)
# b=np.transpose(b)#转置，得是np.array类型
# # print(b.shape)#(16383, 2540)
# # print(c.shape)#(16383,)
# X_train2, X_test2, y_train2, y_test2 = train_test_split(b, c, random_state=1,stratify=c)#random_state为随机种子数，test_size默认为0.25,抽样分布参考y的分布,shuffle参数默认是True，洗牌，会在拆分前重组数据顺序。
# parameters = {"hidden_layer_sizes": [(100,50,30), (100, 30),(100,)],
#                              "solver": ['adam', 'sgd', 'lbfgs'],
#                              "max_iter": [200],
#                              "verbose": [True],
#                              }
# mlp = MLPClassifier(random_state=0)
# estimator = GridSearchCV(mlp,parameters, n_jobs=-1)
# estimator.fit(X_train, y_train)
# print(estimator.get_params().keys())
# print(estimator.best_params_)
# print('ACC:',estimator.best_score_)
