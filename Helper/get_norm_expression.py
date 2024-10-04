from scipy import sparse
from sklearn.utils import sparsefuncs
import pandas as pd
import numpy as np

def normalization(dataframe):
    X = sparse.csr_matrix(dataframe.values, dtype='float64')
    counts_per_cell = X.sum(axis = 1)  # original counts per cell
    counts_per_cell = np.ravel(counts_per_cell)
    counts_greater_than_zero = counts_per_cell[counts_per_cell > 0]
    # Use 1 per 10k 
    after = 1e4
    counts_per_cell += counts_per_cell == 0
    counts_per_cell = counts_per_cell / after
    # Scale the data
    sparsefuncs.inplace_row_scale(X, 1 / counts_per_cell)
    # Apply natural log transformation
    result = pd.DataFrame.sparse.from_spmatrix(X)
    result.columns = dataframe.columns
    result = result.sparse.to_dense()
    return result