#!/usr/bin/python
import numpy
import sys
import pickle
import os
import cplex
from cplex.exceptions import CplexError
import pylab
import matplotlib.pyplot as plt


# this script requires the following arguments: the random subsample index as sample, the corresponding unique virus sequences as train_size, and the number of mutations (\alpha in the paper) as alpha. 
# the current folder is sample_and_train_sample_train_size which has the training data, the z_scores and the parameters saved from regression.
# it loads the coefficient and the intercept parameters from the z-score regression model trained previously on this subsample. 
# finally, it constructs the MILP 6 in the paper for robust ab design and uses cplex python interface to compute the solution.


sample=int( sys.argv[1])
train_size=int( sys.argv[2])
alpha=int( sys.argv[3])

acid_vector=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","X"]
M=20

def obtain_index(acid):
    acid_vector=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    for i in range(0, len(acid_vector)):
        if acid==acid_vector[i]:
           return i
N_ab=27
N_v=32
Z=10

# count matrix is the probability matrix  computed from mutation history
count_matrix=numpy.loadtxt("mutation_counts.txt")
count_matrix=count_matrix>0

# native list is the list of unique virus sequences in the training set, equal to train_size
virus_file=open("sample_and_train_%d_%d/native_list.txt" %(sample,train_size))

virus=virus_file.readlines()

# native antibody is vrc23
antibody=open("vrc23.txt")
antibody=antibody.readlines()
native_antibody=antibody[0]
native_antibody=native_antibody[0:N_ab]
native_ab_feature_vector=[]
for i in range (0,N_ab):
    acid=native_antibody[i]
    acid_index=obtain_index(acid)
    for j in range(0,20):
        if acid_index==j:
           native_ab_feature_vector.append(1)
        else:
           native_ab_feature_vector.append(0)


# loading the saved coefficient and intercept from the regression model
weights=numpy.loadtxt("sample_and_train_%d_%d/coefficients.txt" %(sample,train_size))

intercept=numpy.loadtxt("sample_and_train_%d_%d/intercept.txt" %(sample,train_size))


# ab design MILP from paper (Equation 6) in cplex python interface
ab_weights=weights[0:N_ab*20]
v_weights=weights[N_ab*20:N_v*20+N_ab*20]



train_size=len(virus)
variables=[]
objective=[]
upper_bounds=[]
lower_bounds=[]
v_types=[]

for i in range(0,N_ab):
    for j in range(0,M):
        variables.append("a_"+str(i)+"_"+str(j))
        objective.append(train_size*ab_weights[i*20+j])
        upper_bounds.append(1)
        lower_bounds.append(0)
        v_types.append("I")

for k in range(0,train_size):
    for i in range(0,N_v):
        
            variables.append("d_" +str(k)+"_"+str(i))
            objective.append(1)
            upper_bounds.append(cplex.infinity)
            lower_bounds.append(0)
            v_types.append("C")
for k in range(0,train_size):
    variables.append("d_" +str(k)+"_"+str(N_v))
    objective.append(alpha-N_v)
    upper_bounds.append(cplex.infinity)
    lower_bounds.append(0)
    v_types.append("C")

for k in range(0,train_size):
    for i in range(0,N_v):
        for j in range(0,M):
            variables.append("s2_" +str(k)+"_"+str(i)+"_"+str(j))
            objective.append(1)
            upper_bounds.append(cplex.infinity)
            lower_bounds.append(0)
            v_types.append("C")
        
for k in range(0,train_size):
    for i in range(0,N_v):
        for j in range(0,M):
            variables.append("s3_" +str(k)+"_"+str(i)+"_"+str(j))
            objective.append(Z*(count_matrix[i][j]))
            upper_bounds.append(cplex.infinity)
            lower_bounds.append(0)
            v_types.append("C")


coef_weights=numpy.zeros((N_v,M))
virus_weights=weights[N_ab*M:N_ab*M+N_v*M]
count=0
for i in range(0,N_v):
    for j in range(0,20):
        index=i*20+j
        antibody_coef_sum=0
        for k in range(0,N_ab):
            for u in range(0,M):
                matrix=weights[((N_ab*20+N_v*20)+N_v*20*20*k+20*20*i):((N_ab*20+N_v*20)+N_v*20*20*k+20*20*(i+1))]
                Q_ki=numpy.reshape(matrix,(20,20))
                antibody_coef_sum=antibody_coef_sum+native_ab_feature_vector[20*k+u]*Q_ki[u][j]
                #coef_weights[i][j]=antibody_coef_sum


prob = cplex.Cplex()
prob.objective.set_sense(prob.objective.sense.minimize)
prob.variables.add(obj = objective, ub = upper_bounds,lb=lower_bounds, names = variables,types=v_types)


for t in range(0,train_size):
    
    native_virus=virus[t]
    native_virus=native_virus[0:32]
    native_feature_vector=[]
    for i in range (0,N_v):
        acid=native_virus[i]
        acid_index=obtain_index(acid)
        for j in range(0,20):
            if acid_index==j:
               native_feature_vector.append(1)
            else:
               native_feature_vector.append(0)
    for i in range(0,N_v):
        for j in range(0,M):
             
            var=[]
            coef=[]
            right_hand_side=v_weights[i*20+j]
            var.append("s2_" +str(t)+"_"+str(i)+"_"+str(j))
            coef.append(1)
            var.append("s3_" +str(t)+"_"+str(i)+"_"+str(j))
            coef.append(1)
            var.append("d_"+str(t)+"_"+str(i))
            coef.append(1)
            var.append("d_" +str(t)+"_"+str(N_v))
            coef.append(-1*native_feature_vector[i*20+j])
            for k in range(0,N_ab):
                for u in range(0,M):
                    #print "inside ab variables loop"
                    matrix=weights[((N_ab*20+N_v*20)+N_v*20*20*k+20*20*i):((N_ab*20+N_v*20)+N_v*20*20*k+20*20*(i+1))]
                    Q_ki=numpy.reshape(matrix,(20,20))
                    var.append("a_"+str(k)+"_"+str(u))
                    coef.append(-1*Q_ki[u][j])
            prob.linear_constraints.add(lin_expr =[cplex.SparsePair(ind = var, val = coef)],senses = ["G"],rhs = [right_hand_side])

for k in range(0,N_ab):
    var=[]
    coef=[]
    for u in range(0,M):
        var.append("a_"+str(k)+"_"+str(u))
        coef.append(1)
    prob.linear_constraints.add(lin_expr =[cplex.SparsePair(ind = var, val = coef)],senses = ["E"],rhs = [1])
    
try:
   prob.write("bilevel_%d_%d_%d.lp" %(alpha,sample,train_size))
   print "written"
   prob.solve()
except CplexError, exc:
   print exc
    
print ()
# solution.get_status() returns an integer code
print "Solution status = " , prob.solution.get_status(), ":",
# the following line prints the corresponding string
print(prob.solution.status[prob.solution.get_status()])
print("Solution value  = ", prob.solution.get_objective_value())

x  = prob.solution.get_values()
    
numcols = prob.variables.get_num()
#for j in range(numcols):
#    print("Column %d:  Value = %10f %s " % (j, x[j],variables[j]))            
x=x[0:N_ab*M]


# write solution antibody to file
write_file=open("optimized_ab_%d_%d_%d.txt" %(alpha,sample,train_size),"w")
best_ab=[]
for i in range(0,N_ab):
    for j in range(0,M):
        if x[i*M+j]==1:
           acid=acid_vector[j]
           best_ab.append(acid) 
           write_file.write("%s" %acid)

write_file.close()

print best_ab

