#coding=utf-8

from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from sklearn.model_selection import GridSearchCV
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from collections import Counter
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import matplotlib.colors as colour
from sklearn.preprocessing import label_binarize
from itertools import cycle
import joblib
# d=pd.read_csv(r'/home/MF21300005/pandora-seq/1_41/潘多拉大样本汇总/bigsample2kinds.csv',sep=',',nrows=1,index_col=0,header=None)
# a=pd.read_csv(r'/home/MF21300005/pandora-seq/1_41/潘多拉大样本汇总/bigsample2kinds.csv',sep=',',skiprows=0,index_col=0)#2540维**样本
# # d=pd.read_csv(r'/home/MF21300005/3Dgene/agilent3d2.csv',sep=',',nrows=1,index_col=0,header=None)
# a=pd.read_csv(r'/home/MF21300005/3Dgene/agilent3d3.csv',sep=',',skiprows=0,index_col=0)#2540维**样本
# d=pd.read_csv(r'/home/MF21300005/otherE-MTAB-8026.csv',sep=',',nrows=1,index_col=0,header=None)
# a=pd.read_csv(r'/home/MF21300005/otherE-MTAB-8026.csv',sep=',',skiprows=0,index_col=0)#2540维**样本

#2分类
# d=pd.read_csv(r'/home/MF21300005/newqcb3Dgene2kindsn.csv',sep=',',nrows=1,index_col=0)
# a=pd.read_csv(r'/home/MF21300005/newqcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维**样本

#3分类
# d=pd.read_csv(r'/home/MF21300005/newqcb3Dgene3kindsn.csv',sep=',',nrows=1,index_col=0)
# a=pd.read_csv(r'/home/MF21300005/newqcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维**样本


#25分类
d=pd.read_csv(r'/home/MF21300005/3Dgene/2023.9.283Dgene/qcb3Dgenecankinds.csv',sep=',',nrows=1,index_col=0)
a=pd.read_csv(r'/home/MF21300005/3Dgene/2023.9.283Dgene/qcb3Dgene.csv',sep=',',skiprows=0,index_col=0)#2540维**样本

c=np.array(d)#1*48
c=c.ravel()#转置成（47，）格式,不可以用np.transpose()


f = Counter(c)#统计分类个数，最小分类不能小于2
print(f)
h=set(c)
h=list(h)
h.sort(reverse=False)
print("h:",h)
b=np.array(a)
from sklearn.preprocessing import MinMaxScaler

scaler = MinMaxScaler( )
b=scaler.fit_transform(b)
b=np.transpose(b)#转置，得是np.array类型
print(b.shape)#(23484, 2540)
print(c.shape)#(23484,)

x_train, x_test, y_train, y_test = train_test_split(b, c, random_state=1,stratify=c)#random_state为随机种子数，test_size默认为0.25,抽样分布参考y的分布
print(x_train.shape)#(12287, 2540)
print(y_train.shape)#(12287,)
# print(y_train.isnull().any())#判断哪些”列”存在缺失值

# params = {'n_estimators' : [50,60,70,80,90,100,110,120,130,140,150,200,250,300,400],
#            'max_depth' : list(range(2,46,5)),
#            'min_samples_leaf' : list(range(1, 22, 2)),
#            'min_samples_split' : list(range(2, 9, 2))
#             }

params = {'n_estimators' : [150,300],
           'max_depth' : [42],
           'min_samples_leaf' : [1],
           'min_samples_split' : [2]
            }
rf_clf = RandomForestClassifier(random_state = 0, n_jobs = -1,class_weight='balanced')
grid_cv = GridSearchCV(rf_clf, param_grid = params, cv = 5, n_jobs = -1,refit=True)
grid_cv.fit(x_train, y_train)
print('Optimal Hyper Parameter, RF: ', grid_cv.best_params_)
print('Maximum Accuracy, RF: {:.4f}'.format(grid_cv.best_score_))
rf2 = grid_cv.best_estimator_
rf2score=rf2.score(x_test,y_test)
print('Test Accuracy, RF: {:.4f}'.format(rf2score))

e=pd.read_csv(r'/home/MF21300005/3Dgene/2023.9.283Dgene/qcb3Dgene.csv',sep=',',skiprows=0)
print(e.shape)
e=np.array(e)
e=np.transpose(e)
feat_labels=e[0]
y_pred = rf2.predict(x_test)
y_predprob = rf2.predict_proba(x_test)
# print("y_test20:",y_test[1:20])
# print("y_predprob",y_predprob[0:10,:])
# print(y_predprob.shape)
# print("y_pred20:",y_pred[1:20])

# AUC=roc_auc_score(y_test, y_predprob[:,1])
# print('AUC:',AUC)

importance =rf2.feature_importances_
imp_result = np.argsort(importance)[::-1][:20]
for i in range(len(imp_result)):
    print("%2d. %-*s %f" % (i + 1, 30, feat_labels[imp_result[i]], importance[imp_result[i]]))




#
joblib.dump(rf2, "qcb3Dgenetrain_modelcankinds.m") #存储
clf = joblib.load("qcb3Dgenetrain_modelcankinds.m") #调用
# clf = joblib.load("agilent3d3train_model2kinds.m") #调用
#

#测试GSE113486
d=pd.read_csv(r'/home/MF21300005/3Dgene/2023.9.283Dgene/g1GSE113486.csv',sep=',',nrows=1,index_col=0)
# d=pd.read_csv(r'/home/MF21300005/3Dgene/qcb3Dgenecankinds.csv',sep=',',nrows=1,index_col=0)
a=pd.read_csv(r'/home/MF21300005/3Dgene/2023.9.283Dgene/g1GSE113486.csv',sep=',',skiprows=0,index_col=0)#2540维**样本
c=np.array(d)#1*48
c=c.ravel()#转置成（47，）格式,不可以用np.transpose()
f = Counter(c)#统计分类个数，最小分类不能小于2
print(f)
h=set(c)
h=list(h)
h.sort(reverse=False)
print(h)
b=np.array(a)
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler( )
b=scaler.fit_transform(b)
b=np.transpose(b)#转置，得是np.array类型
print(b.shape)#(16383, 2540)
print(c.shape)#(16383,)
clfscore=clf.score(b,c)
print('Test Accuracy, RF: {:.4f}'.format(clfscore))
# e=pd.read_csv(r'/home/MF21300005/3Dgene/qcbagilent32kinds.csv',sep=',',skiprows=0)
# print(e.shape)
# e=np.array(e)
# e=np.transpose(e)
# feat_labels=e[0]
y_pred = clf.predict(b)

y_predprob = clf.predict_proba(b)

print(y_predprob[0:10,:])

y_trainpred=rf2.predict(x_train)
from sklearn.metrics import confusion_matrix
C1 = confusion_matrix(y_train, y_trainpred, labels=[0,1])  # 可将'1'等替换成自己的类别，如'cat'。
plt.matshow(C1, cmap=plt.cm.Reds)  # 根据最下面的图按自己需求更改颜色
# plt.colorbar()
for i in range(len(C1)):
    for j in range(len(C1)):
        plt.annotate(C1[j, i], xy=(i, j), horizontalalignment='center', verticalalignment='center')

# plt.tick_params(labelsize=15) # 设置左边和上面的label类别如0,1,2,3,4的字体大小。

plt.ylabel('True label')
plt.xlabel('Predicted label')
# plt.ylabel('True label', fontdict={'family': 'Times New Roman', 'size': 20}) # 设置字体大小。
# plt.xlabel('Predicted label', fontdict={'family': 'Times New Roman', 'size': 20})
plt.xticks(range(0,2), labels=['0','1']) # 将x轴或y轴坐标，刻度 替换为文字/字符
plt.yticks(range(0,2), labels=['0','1'])
#_____________________*******************____________________________________
plt.savefig('3DgeneRFcankinds_train0confusematrix.pdf',bbox_inches ='tight',dpi=600)
plt.close()

C2 = confusion_matrix(y_test, y_pred, labels=[0,1])  # 可将'1'等替换成自己的类别，如'cat'。
plt.matshow(C2, cmap=plt.cm.Reds)  # 根据最下面的图按自己需求更改颜色
# plt.colorbar()
for i in range(len(C2)):
    for j in range(len(C2)):
        plt.annotate(C2[j, i], xy=(i, j), horizontalalignment='center', verticalalignment='center')

# plt.tick_params(labelsize=15) # 设置左边和上面的label类别如0,1,2,3,4的字体大小。

plt.ylabel('True label')
plt.xlabel('Predicted label')
# plt.ylabel('True label', fontdict={'family': 'Times New Roman', 'size': 20}) # 设置字体大小。
# plt.xlabel('Predicted label', fontdict={'family': 'Times New Roman', 'size': 20})
plt.xticks(range(0,2), labels=['0','1']) # 将x轴或y轴坐标，刻度 替换为文字/字符
plt.yticks(range(0,2), labels=['0','1'])
#_____________________*******************____________________________________
plt.savefig('3DgeneRFcankinds_test0confusematrix.pdf',bbox_inches ='tight',dpi=600)
plt.close()

# print("c:",c)
# print(y_predprob.shape)
# print("y_pred:",y_pred)
# AUC=roc_auc_score(c, y_predprob[:,1])
# print('AUC:',AUC)
# importance = clf.feature_importances_
# imp_result = np.argsort(importance)[::-1][:10]
# for i in range(len(imp_result)):
#     print("%2d. %-*s %f" % (i + 1, 30, feat_labels[imp_result[i]], importance[imp_result[i]]))

# roc_ovr=roc_auc_score(c, y_predprob, multi_class='ovr')
# print('--AUC-ovr:',roc_ovr)
# roc_ovo=roc_auc_score(c, y_predprob, multi_class='ovo')
# print('--AUC-ovo:',roc_ovo)

#2分类ROC
# def Find_Optimal_Cutoff(tper,fper,threshold):
#     y=tper-fper
#     Youden_index=np.argmax(y)
#     optimal_threshold=threshold[Youden_index]
#     point=[fper[Youden_index],tper[Youden_index]]
#     return optimal_threshold,point
# prob = y_predprob[:, 1]
# fper, tper, thresholds = roc_curve(y_test, prob,pos_label=1)
# auc =auc(fper, tper)
# optimal_th,optimal_point=Find_Optimal_Cutoff(tper,fper,thresholds)
# plt.plot(fper, tper, color='blue', label='ROC curve (area = %0.2f)' % auc)
# plt.plot([0, 1], [0, 1], color='green', linestyle='--')
# plt.plot(optimal_point[0], optimal_point[1], marker='o', color='r')
# plt.text(optimal_point[0],optimal_point[1], f'Threshold:{optimal_th:.2f}')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('RF ROC Curve')
# plt.legend(loc="lower right")
# plt.savefig('./new3DgeneROC2kindsRF.pdf',dpi=300)
#
#
#
# from sklearn.metrics import precision_recall_curve,average_precision_score
# # precision, recall值的计算
# precision, recall, _ = precision_recall_curve(y_test,y_pred)
# # average_precision值的计算
# PRC = average_precision_score(y_test,y_pred)
# print("PRC:",PRC)
# # PRC曲线绘制
# plt.figure(figsize=(6,6))
# plt.title('PR curves')
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.xlim([-0.05,1.05])
# plt.ylim([-0.05,1.05])
# plt.step(recall, precision, color='b', label=' (PRC={:.4f})'.format(PRC))
# plt.plot([0, 1], [1, 0], color='m', linestyle='--')
# plt.legend(loc='lower right')
# # 保存图片(常用格式如下)
# plt.savefig('new3DgenePRcurves2kinds.pdf',dpi=300)




#多分类ROC

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
plt.title('RF ROC Curve to multi-class')
plt.legend(loc='lower right',fontsize=8)
plt.savefig('./qcb3Dgene_GSE113486_ROC25kinds.pdf',dpi=600)