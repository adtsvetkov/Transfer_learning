args = commandArgs(trailingOnly=TRUE)

require(pracma)

#вектор B мы убрали

#радиально-базисная функция - функция Гаусса

#работает
RBF <- function(eps, x1_vec, x2_vec) 
{
  return(exp(-(eps*Norm(x1_vec-x2_vec))^2))
}

#работает

get_Kij <- function(xs_matr, eps)
{
  r <- ncol(xs_matr)-1
  K <- matrix(nrow = r, ncol = r)
  for (i in 1:r) 
  {
    for (j in 1:r) 
    {
      K[i, j] <- RBF(eps, xs_matr[, i], xs_matr[, j])
    }
  }
  return(K)
}

#работает

get_ki <- function(xs_matr, xt_matr, eps)
{
  k <- vector()
  s <- ncol(xs_matr)-1
  t <- ncol(xt_matr)-1
  for (i in 1:s)
  {
    vec <- vector()
    for(j in 1:t) 
    {
      vec <- append(vec, RBF(eps, xs_matr[, i], xt_matr[, j]))
    }
    k <- append(k, (s/t)*sum(vec))
  }
  return(k)
}

#V - матрица n x m
#n_s = n = n_t*t

KMM <- function(x_t, x_s, eps)
{
  m <- ncol(x_t)-1
  V <- matrix()
  n_t = nrow(x_t)
  n_s = nrow(x_s)
  
  K <- get_Kij(x_s, eps)
  print(K)
  k <- get_ki(x_s, x_t, eps)
  print(k)
  
  ####
    
  #grid <- 10^seq(10, -2, length = 100)
  #x <- model.matrix(K)
  #lasso <- glmnet(x, K, alpha = 1, lambda = grid)
  #plot(lasso, xvar = "lambda", label = T, lwd = 2)
    
  ###
}

if (length(args) == 0) 
{
  path1 <- "C:\\source_domain.csv"
  path2 <- "C:\\target_domain.csv"
} else if (length(args)==2) {
  path1 <- args[1]
  path2 <- args[2]
} else if (length(args)!=2)
{
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

source <- as.matrix(read.csv(path1, row.names = 1))
target <- as.matrix(read.csv(path2, row.names = 1))

V <- KMM(target, source, 0.1)
