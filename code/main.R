# Library for KDE
library(np)
# Library for Non-Parametric Single Level Regular Vines
source("nprv.R")
# Library for Gaussian Single Level Regular Vines
source("grv.R")

# Do all files (they are already stored in copula form)
for(i in dir(pattern="*txt"))
{
  # Read file

  x <- read.table(i)

  # 700 train, 300 test

  s <- x[sample(1:nrow(x), 1000, replace=T),]
  s_train <- sample(1:1000, 300)
  s_test  <- setdiff(1:1000, s_train)
  s_train <- s[s_train,]
  s_test  <- s[s_test,]

  # Non-Parametric Single Level Regular Vine
  nprv_ken     <- nprv(s_train)
  res_nprv_ken <- mean(log(nprv_eval(nprv_ken,s_test)))

  # Gaussian Single Level Regular Vine
  grv_ken     <- grv(s_train)
  res_grv_ken <- mean(log(grv_eval(grv_ken,s_test)))

  # Regular Kernel Density Estimator (KDE)
  res_kde <- mean(log(npudens(tdat=s_train, edat=s_test, bwmethod="normal-reference")$dens))

  # Output results dto densities.res
  print(c(i, res_kde, res_grv_ken, res_nprv_ken))
}

