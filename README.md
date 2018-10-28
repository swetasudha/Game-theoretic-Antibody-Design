The scripts in this repository correspond to the paper Game-theoretic Antibody Design.

z_score_regression_training.py performs z_score regression on a given random subsample of the data.

ILP_virus_escape.py loads the z-score regression model parameters and solves the ILP 4 in the paper to compute the escaping virus sequence 
starting from a native virus sequence, given an antibody, within \alpha mutations.

MILP_ab_design.py loads the z-score regression model parameters and solves the MILP 6 in the paper to compute the robust antibody  within \alpha mutations.
