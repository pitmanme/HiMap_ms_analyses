import lomap

"""
Driver code for calculation of Frobenius norms:

A = X - Y
d(X, Y) = sqrt(Tr{A'A})

Refs:
[1] https://inst.eecs.berkeley.edu/~ee127/sp21/livebook/l_mats_norms.html
[2] https://search.r-project.org/CRAN/refmans/SMFilter/html/FDist2.html
"""
# *****************************************************************************
# Written by Dr. Mary Pitman. 2022
# You must have HiMap installed to run, and include mod from Ref 2 above
# This is an 'academic script' toggled via commenting.
# *****************************************************************************

#-------------------------------------------------------#
# Generate similarity scores.
#-------------------------------------------------------#
# Generate random weights for 20 ligands
n_arr = lomap.rand_sim_scores(20)

#-------------------------------------------------------#
# Optimization.
#-------------------------------------------------------#
n_trials = [1, 2, 3]

# n = 10
k_inputs = ['min', '1n', 15, '2n', 'nlnn', 30, 35, 40, 'max']

# n = 20
# k_inputs = ['min', 20, 30, 40, 50, 'nlnn', 70, 80, 90, 100, 120, 140, 160, 180, 'max']

# n = 40
# k_inputs = [42, 45, 50, 60, '2n', 100, 'nlnn', 200, 300, 500, 600, 700, 'max']

# n = 80
# k_inputs = ['min', '1n', 15, '2n', 'nlnn', 30, 35, 40, 'max']

# Send the user selected clusters for optimization.
for k in k_inputs:
    # Repeat operation 20 times for averaging
    for i in n_trials:
        lomap.Optimize(n_arr, num_edges = k, optim_types = ['random', 'random'])


