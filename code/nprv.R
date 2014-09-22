#
# Non-Parametric Single Level Regular Vine
#

library(np)
library(kernlab)

#
# Maximum Mean Discrepancy
#

mmd <- function(sampleX, sampleY, a = 0.05)
{
  sampleX <- as.matrix(sampleX)
  sampleY <- as.matrix(sampleY)
  kmmd(sampleX, sampleY, alpha = a)@H0
}

#
# Prim's Algorithm
#

prim <- function(G)
{
  V     <- 1:nrow(G) 
  V_new <- sample(1:length(V), 1)
  E_new <- NULL

  best_idx <- 1

  while(setequal(V, V_new) == FALSE)
  {
    best <- -Inf

    for(i in intersect(V, V_new))
    {
      for(j in setdiff(V, V_new))
      {
        if(G[i,j] > best)
        {
          best <- G[i,j]
          best_idx <- c(i,j)
        }
      }
    }

    V_new <- c(V_new, best_idx[2])
    E_new <- rbind(E_new, best_idx)
  }

  E_new
}

#
# Kernel Density Estimation routine
#

ke <- function(x)
{
  bwdens <- npudensbw(as.matrix(x), bwmethod="nor")
  bwdist <- npudistbw(as.matrix(x), bwmethod="nor")

  list(
    pdf = function(x) npudens(bwdens, edat=x)$dens,
    cdf = function(x) npudist(bwdist, edat=x)$dist,
    x   = x
  )
}

#
# Kernel Estimation of Marginal PDF 
#

NPPDF_fit <- function(d)
{
  ke(as.matrix(d))$pdf
}

#
# Kernel Estimation of Marginal CDF 
#

NPCDF_fit <- function(d)
{
  ke(as.matrix(d))$cdf
}

#
# Fit NPCopula object
#

NPCopula_fit <- function(sample)
{
  sampleR <- qnorm(as.matrix(sample))
  kest    <- ke(sampleR)

  pdf     <- function(u,v)
  {
    kest$pdf(cbind(qnorm(u), qnorm(v))) / dnorm(qnorm(u)) / dnorm(qnorm(v))
  }

  d       <- sampleR

  list(pdf = pdf, d = sampleR)
}

#
# Build Maximum-Spanning-Trees
#

NPVine_select_edges <- function(d)
{
  cparams <- NULL
  dd      <- ncol(d)
  
  H <- cor(d,method="kendall")

  T <- prim(H)

  for(e in 1:nrow(T))
  {
    cparams <- c(cparams, list(sort(T[e,])))
  }

  cparams
}

#
# Non-Parametric R-Vine evaluation at sites d
#

nprv_eval <- function(npv, d)
{
  res <- 1

  for(i in 1:npv$n)
  {
    res   <- res * npv$pdfs[[i]](d[,i])
    d[,i] <- npv$cdfs[[i]](d[,i])
  }
  
  d[d >= 0.99988] <- 0.99988
  d[d <= 0.00001] <- 0.00001

  for(cp in npv$cparams)
  {
    res <- res * npv$copulas[cp[1],cp[2]][[1]]$pdf(d[,cp[1]],d[,cp[2]])
  }

  res
}

#
# Fit R-Vine with kernel estimated bivariate Copula functions. 
#

nprv <- function(d, unlabeled=NULL)
{
  n    <- ncol(d)
  pdfs <- NULL
  cdfs <- NULL
  miny <- min(d[,n])-2*sd(d[,n])
  maxy <- max(d[,n])+2*sd(d[,n])

  X <- rbind(d[,1:(n-1)],unlabeled)

  for(i in 1:(n-1))
  {
    pdfs <- c(pdfs, NPPDF_fit(X[,i]))
    cdfs <- c(cdfs, NPCDF_fit(X[,i]))
    X[,i] <- cdfs[[i]](X[,i])
    d[,i] <- cdfs[[i]](d[,i])
  }

  pdfs  <- c(pdfs, NPPDF_fit(d[,n]))
  cdfs  <- c(cdfs, NPCDF_fit(d[,n]))
  d[,n] <- cdfs[[n]](d[,n])

  cparams <- NPVine_select_edges(d)
  copulas <- array(list(), c(n,n))
  
  for(i in cparams)
  {
    if(n == i[2])
    {
      copulas[i[1],i[2]] <- list(NPCopula_fit(cbind(d[,i[1]], d[,n])))
    }
    
    else
    {
      copulas[i[1],i[2]] <- list(NPCopula_fit(cbind(X[,i[1]], X[,i[2]])))
    }
  }

  list(pdfs=pdfs,cdfs=cdfs,copulas=copulas,cparams=cparams,n=n,miny=miny,maxy=maxy)
}

