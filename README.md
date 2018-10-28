The scripts in this repository correspond to the paper Game-theoretic Antibody Design.

ILP_virus_escape.py loads the z-score regression model parameters and solves the ILP 4 in the paper to compute the escaping virus sequence 
starting from a native virus sequence, given an antibody, within \alpha mutations.

--This script requires the following arguments: the random subsample index as sample, the corresponding unique virus sequences as 
--train_size, and the number of mutations (\alpha in the paper) as alpha. 
--The current folder is sample_and_train_sample_train_size which has the training data, the z_scores and the parameters saved from   
--regression.
--It loads the coefficient and the intercept parameters from the z-score regression model trained previously on this subsample. 
--Finally, it constructs the ILP 4 in the paper for solving virus escape and uses cplex python interface to compute the solution.

MILP_ab_design.py loads the z-score regression model parameters and solves the MILP 6 in the paper to compute the robust antibody  within \alpha mutations.

--This script requires the following arguments: the random subsample index as sample, the corresponding unique virus sequences as 
--train_size, and the number of mutations (\alpha in the paper) as alpha. 
--The current folder is sample_and_train_sample_train_size which has the training data, the z_scores and the parameters saved from
--regression.
--It loads the coefficient and the intercept parameters from the z-score regression model trained previously on this subsample. 
--Finally, it constructs the MILP 6 in the paper for robust ab design and uses cplex python interface to compute the solution.

z_score_regression_training.py performs z_score regression on a given random subsample of the data.

--This script reads the data, constructs the training set, creates the feature representation and performs z-score regression.
--Finally, the model parameters are saved as the coefficient and the intercept.
--The arguments for this python script, sample is an index of a random subsample of training data, train_size is the number of viruses in
--the data, regularization is the parameter for regression.
--The current folder is sample_and_train_sample_train_size which has a file train_v_set.txt which contains the list of unique virus
--sequences in this random subsample.
--This script matches these unique sequences with the training data: antibodies.txt, viruses.txt and z_scores.txt (these files are in the
--parent folder of the current folder).
--The files antibodies.txt, viruses.txt and z_scores.txt contain the ab and v binding site strings ad the corresponding z_scores
--respectively.
--Sparse matrixes are used to speed up training
