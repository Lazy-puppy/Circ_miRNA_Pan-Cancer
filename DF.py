# -*- coding: utf-8 -*-
# @Time    : 2022/2/15 21:15
# @Author  : gepeng
# @File    : DF.py
# from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from deepforest import CascadeForestClassifier #安装deep-forest
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from collections import Counter
# X, y = load_digits(return_X_y=True)#（手写数字识别）X为1797*64（1797个样本64维），y为(1797,),X,y是numpy.ndarray类型

from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import matplotlib.colors as colour
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene2kinds.csv',sep=',',nrows=1,index_col=0)
# a=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维**样本
# c=np.array(d)#1*48
# c=c.ravel()#转置成（47，）格式,不可以用np.transpose()
# # f = Counter(c)#统计分类个数，最小分类不能小于2
# # print(f)
# b=np.array(a)
# b=np.transpose(b)#转置，得是np.array类型
# print(f)
# print(b.shape)#(16383, 2540)
# print(c.shape)#(16383,)
# X_train, X_test, y_train, y_test = train_test_split(b, c, random_state=1,stratify=c)#random_state为随机种子数，test_size默认为0.25,抽样分布参考y的分布
# model = CascadeForestClassifier(random_state=0,n_estimators=10,n_trees=236,max_depth=22, n_jobs = -1,use_predictor=True, predictor='forest', backend="sklearn")
# model.fit(X_train, y_train)
# y_pred = model.predict(X_test)
# y_predprob = model.predict_proba(X_test)
# acc = accuracy_score(y_test, y_pred) * 100
# print("\nTesting Accuracy: {:.3f} %".format(acc))
# print('AUC:',roc_auc_score(y_test, y_predprob[:,1]))#正标签会放在第二列所以应[:,1]
# def Find_Optimal_Cutoff(tper,fper,threshold):
#     y=tper-fper
#     Youden_index=np.argmax(y)
#     optimal_threshold=threshold[Youden_index]
#     point=[fper[Youden_index],tper[Youden_index]]
#     return optimal_threshold,point
# prob = y_predprob[:, 1]
# fper, tper, thresholds = roc_curve(y_test, prob,pos_label='Tumors')
# auc =auc(fper, tper)
# optimal_th,optimal_point=Find_Optimal_Cutoff(tper,fper,thresholds)
# plt.plot(fper, tper, color='blue', label='ROC curve (area = %0.2f)' % auc)
# plt.plot([0, 1], [0, 1], color='green', linestyle='--')
# plt.plot(optimal_point[0], optimal_point[1], marker='o', color='r')
# plt.text(optimal_point[0],optimal_point[1], f'Threshold:{optimal_th:.2f}')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('DF ROC Curve')
# plt.legend(loc="lower right")
# plt.savefig('./DF.jpg')

# #三分类
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene3kinds.csv',sep=',',nrows=1,index_col=0)
# a=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维47样本
# c=np.array(d)#1*48
# c=c.ravel()#转置成（47，）格式,不可以用np.transpose()
# # f = Counter(c)#统计分类个数，最小分类不能小于2
# # print(f)
# b=np.array(a)
# b=np.transpose(b)#转置，得是np.array类型
# # print(b.shape)#(16383, 2540)
# # print(c.shape)#(16383,)
# X_train, X_test, y_train, y_test = train_test_split(b, c, random_state=1,stratify=c)#random_state为随机种子数，test_size默认为0.25,抽样分布参考y的分布
# model = CascadeForestClassifier(random_state=0,n_estimators=10,n_trees=236,max_depth=22)
# model.fit(X_train, y_train)
# y_pred = model.predict(X_test)
# y_predprob = model.predict_proba(X_test)
# print('y_predprob:',y_predprob[0:3,:])
# acc = accuracy_score(y_test, y_pred) * 100
# print("\nTesting Accuracy: {:.3f} %".format(acc))
# roc_ovr=roc_auc_score(y_test, y_predprob , multi_class='ovr')#正标签会放在第二列所以应[:,1]
# print('--roc-ovr:',roc_ovr)
# roc_ovo=roc_auc_score(y_test, y_predprob, multi_class='ovo')
# print('--roc-ovo:',roc_ovo)

# #多分类
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgenecankinds.csv',sep=',',nrows=1,index_col=0)
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene2kinds.csv',sep=',',nrows=1,index_col=0)
# d=pd.read_csv(r'/home/MF21300005/3Dgene/g1.csv',sep=',',nrows=1,index_col=0,header=None)
d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene3kinds.csv',sep=',',nrows=1,index_col=0)
a=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维47样本
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
model = CascadeForestClassifier(random_state=0,n_estimators=10,n_trees=236,max_depth=22)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
y_predprob = model.predict_proba(X_test)
acc = accuracy_score(y_test, y_pred) * 100
print("\nTesting Accuracy: {:.3f} %".format(acc))
# roc_ovr=roc_auc_score(y_test, y_predprob , multi_class='ovr')#正标签会放在第二列所以应[:,1]
# print('--roc-ovr:',roc_ovr)
# roc_ovo=roc_auc_score(y_test, y_predprob, multi_class='ovo')
# print('--roc-ovo:',roc_ovo)
n_classes = len(h)
fpr = dict()#创建字典
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _= roc_curve(y_test, y_predprob[:,i],pos_label=h[i])
    roc_auc[i] = auc(fpr[i], tpr[i])
# Compute macro-average ROC curve and ROC area（方法一）
# First aggregate all false positive rates
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])

# Finally average it and compute AUC
mean_tpr /= n_classes
fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# Plot all ROC curves
lw=2
plt.figure()
plt.plot(fpr["macro"], tpr["macro"],
         label='macro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle=':', linewidth=4)

# colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
colors=colour.cnames
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='{0} (area = {1:0.2f})'
             ''.format(h[i], roc_auc[i]))

plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('DF ROC Curve to 3-class')
plt.legend(loc='lower right',fontsize=8)
plt.savefig('./qcb3Dgene_DF_ROC3kinds.pdf',dpi=600)

