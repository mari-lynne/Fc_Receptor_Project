#Linear dependencies check ###

# Linear depenent cols can't be added to make zero
# Can be made as the product of another column in the matrix

# 1) Single observation problem:
# If a column only has observations in x rows that are not present in any other column, they are linear dependent.

# 2) Null vector problem
# If a column has no observations then it can cause linear dependencies

# 3) Inverse observation problem:
# If one column is the inverse of another column, the sum would = 0, and therefore dependent

# 4) Opposite observation problem:
# If all each observations in a column results in the opposite observation in another column, they are also linearly dependent. For example, imagine a time point with both TD and nTD samples, e.g baseline, if comparing to observations in the diagnosis time point, TD time will allways be negative for the nTD samples, and positive for TD samples.

# PROOF
a = c(1, 0, 0, 1)
b = c(0, 1, 1, 0)
c = c(1, 1, 0, 1)
d = c(0, 0, 1, 1)
dmat <- matrix(c(a, b, c, d),
              nrow = 4,
              ncol = 4,
              dimnames = list(c("a", "b", "c", "d")))
print(mat)
print(mat[, 1] + mat[, 2]) #dependent cols
print(mat[, 1] + mat[, 3]) #independent cols

# If [,i] is always the opposite of [,j] when added there will only be one observation
mat <-
  apply(mat, MARGIN = 2, function(x)
    (x + 2)) #adding 2 to matrix to change example
print(mat[, 1] + mat[, 2]) # Only one observation

# Function 4

#filter - poss Option to filter ld columns from matrix
# recode non-numeric vars to numeric so we can add :)

ld_check <- function(mat) {
  #/ 1) Checks:
  # If mat is matrix; else convert
  
  #/ 2) Matrix Addition
  ind <-
    t(combn(ncol(mat), 2)) # All combos of col indexes to form pairwise-col comparisons
  c_names <-
    paste0(ind[, 1], "_v_", ind[, 2]) #Name comparisons cols
  out <-
    apply(ind, 1, function(x)
      (mat[, x[1]] + mat[, x[2]])) #Add all pairwise vectors
  colnames(out) <- c_names # Update sum table colnames
  
  #/ 3) Validate unique observations 
  for (i in seq_along(c_names)) {
    if (length(unique(out[, i])) == 1) {
      print(paste0("Linear dependent matrix columns =", c_names[i]))
    }
  }
}

dmat = mat
ld_check(mat=dmat)
debug(ld_check)
  
 
