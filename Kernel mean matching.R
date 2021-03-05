#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

require(pracma)
require(quadprog)
require(glmnet)
require(xlsx)

#радиально-базисная функция - функция Гаусса

RBF <- function(eps, x1_vec, x2_vec) 
{
  return(exp(-(eps*Norm(x1_vec-x2_vec))^2))
}

minimize <- function(beta, K, k, j, n_s){
  return((1/2)*(t(beta)%*%K%*%beta) - t(k)%*%beta + j*abs(sum(beta)-n_s))
}

get_Kij <- function(xs_matr, eps)
{
  r <- nrow(xs_matr)
  K <- matrix(nrow = r, ncol = r)
  for (i in 1:r) 
  {
    for (j in 1:r) 
    {
      K[i, j] <- RBF(eps, xs_matr[i, ], xs_matr[j, ])
    }
  }
  return(K)
}

get_ki <- function(xs_matr, xt_matr, eps)
{
  k <- vector()
  s <- nrow(xs_matr)
  t <- nrow(xt_matr)
  for (i in 1:s)
  {
    vec <- vector()
    for(j in 1:t) 
    {
      vec <- append(vec, RBF(eps, xs_matr[i, ], xt_matr[j, ]))
    }
    k <- append(k, (s/t)*sum(vec))
  }
  return(k)
}

#V - матрица n x m
#n_s = n = n_t*t

answer <- data.frame()
answer_indexes <- data.frame()

KMM <- function(x_t, x_s, eps)
{
  #задаем размерности
  m <- ncol(x_t) - 1
  V <- matrix()
  n_t <- nrow(x_t)
  n_s <- nrow(x_s)
  
  #ищем и выводим K, k
  K <- get_Kij(x_s, eps)
  print("K matrix:")
  print(K)
  k <- get_ki(x_s, x_t, eps)
  print("k matrix:")
  print(k)
  
  #просим по графику ввести границы лямбды
  
  borders <- c(-3.5, 4)
  
  newgrid <- 10^seq(borders[1], borders[2], length = m)
  
  #опорный вектор найдем так:
  
  betta_ref <- solve.QP(K, k, K, k, meq = 0, factorized = FALSE)
  
  indexarray <- data.frame()
  bettaarray <- data.frame()
  lambdanames <- vector()
  
  for (i in 1:m)
  {
    j <- newgrid[i]
    myans <- optim(betta_ref$solution, minimize, j = j, k = k, K = K, n_s = n_s)$par
    indexes <- order(myans, decreasing = T)
    bettaarray <- rbind(bettaarray, sort(myans, decreasing = T))
    indexarray <-rbind(indexarray, indexes)
    lambdanames <- append(lambdanames, paste("l", i, sep="_"))
  }
  indexarray <- t(indexarray)
  bettaarray <- t(bettaarray)
  
  colnames(indexarray) <-  colnames(bettaarray) <- lambdanames
  rownames(indexarray) <- rownames(bettaarray) <- rep("b_i", times = n_s)
  
  answer <<- bettaarray
  answer_indexes <<- indexarray
  
  #рисуем зависимость лямбда от K и k
  
  grid <- 10^seq(10, -5, length = 100)
  lasso <- glmnet(K, t(k), alpha = 0.3, lambda = grid)
  X11()
  plot(lasso, xvar = "lambda", label = T, lwd = 2, ylim = c(-3, 6))
  
  while(names(dev.cur()) !='null device') Sys.sleep(1)
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

path1 <- "C:\\source_domain.csv"
path2 <- "C:\\target_domain.csv"

source <- as.matrix(read.csv(path1, row.names = 1))
target <- as.matrix(read.csv(path2, row.names = 1))

KMM(target, source, 0.1)

write.xlsx(answer, file = "C:\\KMM_results.xlsx", sheetName="Beta")
write.xlsx(answer_indexes, file = "C:\\KMM_results.xlsx", sheetName="Indexes", append = T)
