#!/usr/bin/python
import sys
import csv
import numpy
from sklearn.linear_model import Lasso
from numpy import linalg as LA
import os
from scipy import sparse

# this script reads the data, constructs the training set, creates the feature representation and performs z-score regression.
# finally, the model parameters are saved as the coefficient and the intercept.
# the arguments for this python script, sample is an index of a random subsample of training data, train_size is the number of viruses in the data, regularization is the parameter for regression.
# the current folder is sample_and_train_sample_train_size which has a file train_v_set.txt which contains the list of unique virus sequences in this random subsample.
# this script matches these unique sequences with the training data: antibodies.txt, viruses.txt and z_scores.txt (these files are in the parent folder of the current folder).
# the files antibodies.txt, viruses.txt and z_scores.txt contain the ab and v binding site strings ad the corresponding z_scores respectively.

sample=int(sys.argv[1]) 
train_size=int(sys.argv[2])
regularization=float(sys.argv[3])


# mapping the 20 amino acids to indices
def obtain_index(acid):
    acid_vector=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    for i in range(0, len(acid_vector)):
        if acid==acid_vector[i]:
           return i

# the following function performs logistic regression, scikit-learn Lasso in particular
def logistic_regression_with_transformed_features(x_train,x_test,y_train,y_test,reg_weight):
    logistic=Lasso(alpha=reg_weight).fit(x_train,y_train) 
    corr_coef= numpy.corrcoef(y_test,logistic.predict(x_test))
    print corr_coef[0,1]
   
    intercept=[] 
    intercept.append(logistic.intercept_)
    print intercept
    corr_coef_var=[]
    corr_coef_var.append(corr_coef[0,1])
    print corr_coef_var
    return logistic.coef_.ravel(),corr_coef_var,intercept
    

# this function shuffles the training data randomly, important in case of cross-validation
def predict_logistic(train_to_test_ratio,M):
    features=open("./sample_and_train_%d_%d/features_train_%d_%d.txt" %(sample,train_size,sample,train_size))
    f=features.readlines()
    
    K=int(train_to_test_ratio*M)
    overall_cum_coef=[]
    reg_weights=[regularization]
    #reg_weights=[0.01]
    for param in range(0,len(reg_weights)):
        C=reg_weights[param]
        cum_coef=0
        for cv in range(0,1):
            arr = numpy.arange(M)
            numpy.random.shuffle(arr)
            c=numpy.zeros((M,num_f+1))
            for j in range(0,M):
                #print j
                index=arr[j]
                #index=j
                line=list(f[index])
                one_line=[]
                for jj in range(0,len(line)-1):
                    one_line.append(int(line[jj]))
                one_line=numpy.asarray(one_line)
                current=numpy.hstack((one_line,x[index]))
                c[j,:]=current
            ff=c[:,0:num_f]
            y=c[:,num_f:num_f+1]
            y=numpy.reshape(y,(M,))
            x_train=ff[0:K]
            x_test=ff[K:M]
            y_train=y[0:K]
            y_test=y[K:M]
            x_train=sparse.csr_matrix(x_train)
            coef1,coef,intercept= logistic_regression_with_transformed_features(x_train,x_train,y_train,y_train,C)
            
        overall_cum_coef.append( float(cum_coef/1))
    numpy.savetxt('./sample_and_train_%d_%d/coefficients.txt' %(sample,train_size),coef1)
    numpy.savetxt('./sample_and_train_%d_%d/intercept.txt' %(sample,train_size),numpy.array(intercept))
    numpy.savetxt('./sample_and_train_%d_%d/corr_coef.txt' %(sample,train_size),coef)




# this section reads the unique virus sequences (equal to the training_size) in the training set and matches those to ab-v-z-score data points in the training data to create the training and the test set for the particular random subsample; this section could have been written much more optimally 

import random 
train_ab_file=open("./sample_and_train_%d_%d/train_set_ab_%d_%d.txt" %(sample,train_size,sample,train_size),"w")
train_v_file=open("./sample_and_train_%d_%d/train_set_v_%d_%d.txt"%(sample,train_size,sample,train_size),"w")
train_scores=open("./sample_and_train_%d_%d/train_set_scores_%d_%d.txt"%(sample,train_size,sample,train_size),"w")


test_ab_file=open("./sample_and_train_%d_%d/test_set_ab_%d_%d.txt" %(sample,train_size,sample,train_size),"w")
test_v_file=open("./sample_and_train_%d_%d/test_set_v_%d_%d.txt"%(sample,train_size,sample,train_size),"w")
test_scores=open("./sample_and_train_%d_%d/test_set_scores_%d_%d.txt"%(sample,train_size,sample,train_size),"w")
native_file=open("./sample_and_train_%d_%d/native_list.txt" %(sample,train_size),"w")

all_ab=open("../antibodies.txt")
all_ab=all_ab.readlines()
all_v=open("../viruses.txt")
all_v=all_v.readlines()
all_score=open("../z_scores.txt")
all_score=all_score.readlines()


train_v_set=open("./sample_and_train_%d_%d/train_v_set.txt"%(sample,train_size),"r")
train_v_set=train_v_set.readlines()


for j in range(0,train_size):
    match=0
    for native_index in range(0,180):
        v=open("../all_virus/%d.txt" %native_index)
        v=v.readlines()
        virus=v[0]
        #temp=list(virus)
       
        if virus == train_v_set[j]:
           #print virus
           match=1
           match_index=native_index
           break
    
    
    v=open("../all_virus/%d.txt" %match_index)
    v=v.readlines() 
    for jj in range(0,4):
        virus=v[jj]
        temp=list(virus)
        temp=temp[0:len(temp)-1]
        virus="".join(temp)
        if jj==0:
           native_file.write("%s" %v[jj])
       
        found=0
        for  k in range(0,len(all_v)):
             temp= list(all_v[k])
             temp=temp[0:len(temp)-2]
             current_v="".join(temp)
            
             if current_v==virus:
                   found=found+1
                   train_ab_file.write("%s" %all_ab[k])
                   train_v_file.write("%s" %all_v[k])
                   train_scores.write("%s" %all_score[k]) 


test_v_set=open("./sample_and_train_%d_%d/test_v_set.txt"%(sample,train_size),"r")
test_v_set=test_v_set.readlines()

for j in range(0,len(test_v_set)):
    train=0
   
    for native_index in range(0,180):
        v=open("../all_virus/%d.txt" %native_index)
        v=v.readlines()
        virus=v[0]
       
        if virus == test_v_set[j]:
           train=1
           match_index=native_index
           break

    if train==1:
      
       v=open("../all_virus/%d.txt" %match_index)
       v=v.readlines()
       for jj in range(0,4):
           virus=v[jj]
           temp=list(virus)
           temp=temp[0:len(temp)-1]
           virus="".join(temp)
           
           found=0
           for  k in range(0,len(all_v)):
              temp= list(all_v[k])
              temp=temp[0:len(temp)-2]
              current_v="".join(temp)
              
              if current_v==virus:
                  
                   found=found+1
                   test_ab_file.write("%s" %all_ab[k])
                   test_v_file.write("%s" %all_v[k])
                   test_scores.write("%s" %all_score[k])
       

train_ab_file.close()
train_v_file.close()
train_scores.close()

native_file.close()
test_ab_file.close()
test_v_file.close()
test_scores.close()
    
# following section reads the training data as antibody and virus strings and constructs one-hot encoding features
a_file=open("./sample_and_train_%d_%d/train_set_ab_%d_%d.txt" %(sample,train_size,sample,train_size))
v_file=open("./sample_and_train_%d_%d/train_set_v_%d_%d.txt" %(sample,train_size,sample,train_size))

feature_file=open("./sample_and_train_%d_%d/features_train_%d_%d.txt" %(sample,train_size,sample,train_size),"w")

lines_a=a_file.readlines()
lines_v=v_file.readlines()

L=len(lines_a)
print L

for data_point in range(0,len(lines_a)):
    features_a=lines_a[data_point]
    features_v=lines_v[data_point]
   

    count=0

    for i in range (0,27):
        #the antibody side
        acid=features_a[i]
        acid_index=obtain_index(acid)
        for j in range(0,20):
           if acid_index==j:
              feature_file.write("%d" %1)
           else:
              feature_file.write("%d" %0)
           count=count+1
     
    for i in range (0,32):
        #the antibody side
        acid=features_v[i]
        acid_index=obtain_index(acid)
        for j in range(0,20):
            if acid_index==j:
               feature_file.write("%d" %1)
            else:
               feature_file.write("%d" %0)
            count=count+1

    
    for i in range (0,27):
        acid_a=features_a[i]
        acid_index_a=obtain_index(acid_a)
        for m in range(0,32):
            acid_v=features_v[m]
            acid_index_v=obtain_index(acid_v)
            for j in range(0,20):
                for k in range(0,20):
                    if acid_index_a==j and acid_index_v==k:
                       feature_file.write("%d" %1)
                    else:
                       feature_file.write("%d" %0)
                    count=count+1
    
    feature_file.write("\n")
feature_file.close()


# loads the z_scores corresponding to the training set and calls the shuffle data function which calls the regression function
x=numpy.loadtxt("./sample_and_train_%d_%d/train_set_scores_%d_%d.txt" %(sample,train_size,sample,train_size))    
num_f=27*20+32*20+27*32*20*20
predict_logistic(1,L)
