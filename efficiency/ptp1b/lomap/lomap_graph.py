# Import lomap
import lomap
import sys
import networkx as nx
import numpy as np

# Create the molecule database by using .mol2 files
# The DBMolecule class must be created with a valid
# directory name
        
db_mol = lomap.DBMolecules('/Users/mpitman/work/gits/packages/Lomap-main/examples/PLB_optim_data/lomap_similarities/ptp1b', output=True, radial=False )
#use the radial option with hub ligand set as 
#db_mol = lomap.DBMolecules('test/radial/', output=True, radial=True, hub="ejm_46.mol2")
    
# Generate the strict and loose symmetric similarity
# score matrices.
strict, loose = db_mol.build_matrices()
    
# Convert the matrices in standard numpy matrices
strict_numpy = strict.to_numpy_2D_array()
loose_numpy = loose.to_numpy_2D_array()

# Networkx graph generation based on the similarity 
# score matrices
          
nx_graph = db_mol.build_graph() 


# Calculate the Maximum Common Subgraph (MCS) between 
# the first two molecules in the molecule database 
# ignoring hydrogens and depicting the mapping in a file
    
MC = lomap.MCS.getMapping(db_mol[0].getMolecule(), db_mol[1].getMolecule(), hydrogens=False, fname='mcs.png')

