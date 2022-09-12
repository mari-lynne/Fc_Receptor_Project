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
mat <- matrix(c(a, b, c, d),
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

# LD check Function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ld_check <- function(mat,
                     filter=c("1")) { #Default filter first col out of LDs
                                      #To switch to remove 2nd col type 2
  message("Default:  Filter = 1
          If LDs are present, remove the first column
          Filter = 2
          If LDs are present, remove the second column
          Filter = 3
          No removal of LD or null columns")
  
  #/ 1) Checks:
  # If mat is matrix; else convert
  
  #/ 1) Rank check
  print(paste0("ncols = ",ncol(mat)))
  print(paste0("rank = ",qr(mat)$rank))
  
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
  for (i in 1:ncol(out)) {
    if (filter == 1 && length(unique(out[, i])) == 1) {
      print(paste0("Linear dependent cols =", c_names[i]))
      x <- str_extract(c_names[1],"\\d")
    } else if ((filter == 2||3) && length(unique(out[, i])) == 1) {
      print(paste0("Linear dependent cols =", c_names[i]))
      x <- str_extract(c_names[1],"\\d$")
    }
  }
  
  #/ 4) Highlight zero vars
  for (i in 1:ncol(mat)){
    if(sum(mat[,i]) == 0){
      print(paste0("Null vectors present - col ",i))
      y <- i
    }
  }
  #/ 5) Highlight vars with only one observations
  for (i in 1:ncol(mat)){
    if(sum(mat[,i]) == 1){
      print(paste0("Single observation - col ",i))
       z <- i
    }
  }
  
  #/ 6) Column removal check
  t1 <- c(x,y,z)
  print(t1)
  message(paste("removing problem columns: ", t1," from matrix"))
  ind_mat <- mat[,t1]
  print(paste0("Updated ncols = ",ncol(mat)))
  print(paste0("Updated rank = ",qr(mat)$rank))
  return(ind_mat)
}


ld_check(mat=dmat)

#debug(ld_check)

mat[,c(x,y,z)]

t <- c(x,y,z)

mat[,t]

paste0(qt, t, c(qt,//,))
rm(z)


addQuotes_<-function(str,qt="'"){
  return(strp<-paste(qt,str,qt,sep=''));
}



#Updates 

#filter - poss Option to filter ld columns from matrix
# recode non-numeric vars to numeric so we can add :)
# If there are lots of LD columns - just show the most common denominator
# If col x or + by other col (pairwise) = value of column in matrix remove
# Null check
# Rank check
colSums(mat) == 0

#Testing ####
a = c(1, 0, 0, 1)
b = c(0, 1, 1, 0)
c = c(1, 1, 0, 1)
d = c(0, 0, 1, 1)
mat <- matrix(c(a, b, c, d),
               nrow = 4,
               ncol = 4,
               dimnames = list(c("a", "b", "c", "d")))
mat <- cbind(mat,(c(0,0,0,0)))
mat <- cbind(mat,(c(0,0,0,1)))

str_extract(c_names[1],"\\d") #just take first mismatch by default
#will therefore remove first col
#change to second d as option
str_extract(c_names[1],"\\d$") #works :)

