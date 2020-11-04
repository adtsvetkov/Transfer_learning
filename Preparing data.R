#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#необходимо подготовить датасет исходных данных и целевых по алгоритму Friedman из статьи Chen et al

create_names <- function(size_x)
{
  names_vec <- vector()
  for (i in 1:size_x)
  { 
    names_vec <- append(names_vec, paste("x", i, sep="_"))
  }
  names_vec <- append(names_vec, "y")
  return(names_vec)
}

Friedman_function <- function(x, a, b, c)
{
  N <- rnorm(1, 0, 1)
  return (a[1]*10*sin(pi*(b[1]*x[1] + c[1])*(b[2]*x[2] + c[2])) + a[2]*20*(b[3]*x[3] + c[3] - 0.5)^2 + a[3]*10*(b[4]*x[4] + c[4]) + a[4]*5*(b[5]*x[5] + c[5]) + N)
}

#creates one target domain instance

target_instance <- function(size_p, size_x)
{
  vec <- vector()
  
  a <- rep(1, times = size_p)
  b <- rep(1, times = size_p)
  c <- rep(0, times = size_p)
  
  x <- runif(size_x, 0, 1)
  
  vec <- append(vec, x)
  
  vec <- append(vec, Friedman_function(x, a, b, c))
  
  return(vec)
}

#creates one source domain instance

source_instance <- function(size_p, size_x)
{
  vec <- vector()
  
  a <- rnorm(size_p, 1, 0.1)
  b <- rnorm(size_p, 1, 0.1)
  c <- rnorm(size_p, 0, 0.05)
  
  x <- runif(size_x, 0, 1)
  
  vec <- append(vec, x)
  
  vec <- append(vec, Friedman_function(x, a, b, c))
  
  return(vec)
}

#creates source data matrix

source_data <- function(size_p, size_x, size_s, names)
{
  matr <- matrix(, ncol = size_x+1, nrow = 0)
  colnames(matr) <- names
  for(i in 1:size_s)
  {
    matr <- rbind(matr, source_instance(size_p, size_x))
  }
  return(matr)
}

#creates target data matrix

target_data <- function(size_p, size_x, size_t, names)
{
  matr <- matrix(, ncol = size_x+1, nrow = 0)
  colnames(matr) <- names
  for(i in 1:size_t)
  {
    matr <- rbind(matr, target_instance(size_p, size_x))
  }
  return(matr)
}

if (length(args) == 0) 
{
  t <- 2
  size_t <- 10
  size_p <- 5
  size_x <- 10
  path1 <- "C:\\source_domain.csv"
  path2 <- "C:\\target_domain.csv"
} else if (length(args)==6) {
  t <- args[1]
  size_t <- args[2]
  size_p <- args[3]
  size_x <- args[4]
  path1 <- args[5]
  path2 <- args[6]
} else if (length(args)!=6)
{
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

size_s = t*size_t

names <- create_names(size_x)

source <- source_data(size_p, size_x, size_s, names)
target <- target_data(size_p, size_x, size_t, names)

write.csv(source, path1)
write.csv(target, path2)
