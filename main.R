require(xlsx)
require(rgp)
require(stringi)
require(xts)
require(DEoptim)

scanfiles <- function(path1, path2, path3)
{
  source <- as.matrix(read.csv(path1, row.names = 1))
  target <- as.matrix(read.csv(path2, row.names = 1))
  beta <- as.matrix(read.xlsx(path3, sheetName = "Beta")[, -1])
  indexes <-
    as.matrix(read.xlsx(path3, sheetName = "Indexes")[, -1])
  return(list(
    source = source,
    target = target,
    beta = beta,
    indexes = indexes
  ))
}

generate_names <- function(cols, rows)
{
  colsnames <- vector()
  rowsnames <- vector()
  for (i in 1:cols)
    colsnames <- append(colsnames, paste("b", i, sep = "_"))
  for (j in 1:rows)
    rowsnames <- append(rowsnames, paste("t", j, sep = "_"))
  return(list(columns = colsnames, rows = rowsnames))
}

check_function <- function(f, x_s)
{
  for (i in 1:nrow(x_s))
  {
    xlist <- as.list(x_s[i, ][-ncol(x_s)])
    res <- do.call(f, xlist)
    if (res == "NaN" || res == "Inf")
      return(F)
  }
  return(T)
}

generate_trees <- function(func, vars, p, x_s)
{
  funclist <- list()
  functionSet1 <- functionSet(list = func)
  inputVariableSet1 <- inputVariableSet(list = vars)
  constantFactorySet1 <- constantFactorySet(function()
    runif(1))
  
  repeat {
    myfunc <-
      randfunc(functionSet1,
               inputVariableSet1,
               constantFactorySet1,
               maxdepth = 10)
    as_str <- paste(deparse(body(myfunc)), collapse = "")
    if (all(stri_detect_fixed(as_str, vars)) &&
        check_function(myfunc, x_s))
      funclist <- append(funclist, myfunc)
    if (length(funclist) == p)
      break
  }
  return(funclist)
}

WMSE <- function(data, trees, weight, indexes)
{
  x_s <- data$source
  wmse <- matrix(, nrow = length(trees), ncol = ncol(weight))
  for (t in 1:length(trees))
  {
    # w - индекс вектора бет
    for (w in 1:ncol(weight))
    {
      # взяли конкретные бета для конкретной функции
      beta <- weight[, w]
      index <- indexes[, w]
      summa <- vector()
      for (i in 1:length(index))
      {
        # номер экземпляра
        instance <- index[i]
        # аргументы функции этого экземпляра
        args <- x_s[instance, ][-length(x_s[instance, ])]
        # результаты вызова функции
        funcres <- do.call(trees[[t]], as.list(args))
        # посчитанная ошибка
        err <-
          beta[i] * ((funcres - x_s[instance, length(x_s[instance])]) ^ 2)
        summa <- append(summa, err)
      }
      wmse[t, w] <- sum(summa)
    }
  }
  name <- generate_names(ncol(wmse), nrow(wmse))
  colnames(wmse) <- name$columns
  rownames(wmse) <- name$rows
  return(wmse)
}

get_best <- function(best)
{
  besttrees <- list()
  bestbetas <- vector()
  for (b in best)
  {
    besttrees <- append(besttrees, b$func)
    bestbetas <- append(bestbetas, list(b$beta))
  }
  return(list(btree = besttrees, bbeta = bestbetas))
}

findbest <- function(wmse, beta, index, trees)
{
  best <- list()
  for (t in 1:length(trees))
  {
    vec <- wmse[t, ]
    ind <- match(min(vec), vec)
    instance <-
      list(
        fvi = min(vec),
        beta = beta[, ind],
        func = trees[[t]],
        ind = index[, ind]
      )
    best <- append(best, list(instance))
  }
  best <-
    best[order(sapply(best, `[[`, "fvi"))][-seq(ncol(beta) + 1, length(trees))]
  return(best)
}

mse_func <- function(tree)
{
  res <- vector()
  for (j in 1:nrow(data$target))
  {
    args <- data$target[j, ][-length(data$target[j, ])]
    resfunc <- do.call(tree, as.list(args))
    result <-
      (resfunc - data$target[j, length(data$target[j, ])]) ^ 2
    res <- append(res, result)
  }
  return(mean(res))
}

tree_breeding <- function(getbest,
                          len,
                          func,
                          vars)
{
  j = 0
  newtrees <- list()
  functionSet1 <- functionSet(list = func)
  inputVariableSet1 <- inputVariableSet(list = vars)
  constantFactorySet1 <- constantFactorySet(function()
    runif(1))
  m <- length(getbest$btree)
  repeat {
    j = j + 1
    tryCatch({
      GP <- geneticProgramming(
        mse_func,
        population = getbest$btree,
        functionSet = functionSet1,
        populationSize = j,
        inputVariables = inputVariableSet1,
        constantSet = constantFactorySet1,
        verbose = F
      )$population
    },
    error = function(cond)
    {
      GP = NULL
    },
    warning = function(cond)
    {
      GP = NULL
    })
    if (!is.null(GP))
    {
      individuals <- unique(GP)
      if (j == m)
        j = 0
      for (i in individuals)
      {
        as_str <- paste(deparse(body(i)), collapse = "")
        if (all(stri_detect_fixed(as_str, vars)) &&
            check_function(i, data$source))
          newtrees <- append(newtrees, i)
      }
    }
    if (length(newtrees) >= len)
      break
  }
  return(first(newtrees, len))
}

wmse_breeding <- function(beta, tree, x_s)
{
  summa <- vector()
  for (i in 1:nrow(x_s))
  {
    # аргументы функции этого экземпляра
    args <- x_s[i,][-length(x_s[i,])]
    # результаты вызова функции
    funcres <- do.call(tree, as.list(args))
    # посчитанная ошибка
    err <-
      beta[i] * ((funcres - x_s[i, length(x_s[i])]) ^ 2)
    summa <- append(summa, err)
  }
  return(sum(summa))
}

beta_breeding <- function(trees, x_s, maxbeta)
{
  indexarray <- data.frame()
  bettaarray <- data.frame()
  lambdanames <- vector()
  low <- rep(0, nrow(x_s))
  up <- rep(maxbeta, nrow(x_s))
  
  for (i in 1:length(trees))
  {
    myans <-
      DEoptim(
        wmse_breeding,
        low,
        up,
        #control = DEoptim.control(itermax = 1, trace = F),
        control = DEoptim.control(trace = F),
        tree = trees[[i]],
        x_s = x_s
      )$optim$bestmem
    indexes <- order(myans, decreasing = T)
    bettaarray <- rbind(bettaarray, sort(myans, decreasing = T))
    indexarray <- rbind(indexarray, indexes)
    lambdanames <- append(lambdanames, paste("l", i, sep = "_"))
  }
  indexarray <- t(indexarray)
  bettaarray <- t(bettaarray)
  
  colnames(indexarray) <-  colnames(bettaarray) <- lambdanames
  rownames(indexarray) <-
    rownames(bettaarray) <- rep("b_i", times = nrow(x_s))
  
  return(list(beta = as.matrix(bettaarray), index = as.matrix(indexarray)))
}

ITGP <- function(data,
                 p,
                 vars,
                 func,
                 trees, 
                 steps = 10,
                 maxbeta = 10)
{
  beta <- data$beta
  indexes <- data$indexes
  envir <- list()
  
  # сгенерировали деревья
  current_steps = 0
  
  repeat {
    current_steps = current_steps + 1
    # идем оцениваем деревья на сорсах
    wmse <- WMSE(data, trees, beta, indexes)
    
    # отбираем лучшие беты на каждом дереве:
    best <- findbest(wmse, beta, indexes, trees)
    getbest <- get_best(best)
    
    # генерим новое поколение бет
    
    breeded <- beta_breeding(getbest$btree, data$source, maxbeta)
    
    beta <- breeded$beta
    indexes <- breeded$index
    
    # размножаем деревья по оценкам на таргетах
    
    trees <- tree_breeding(getbest, p, func, vars)
    
    if (current_steps %% 5 == 0)
    {
      answer <- list()
      for (t in trees) answer <- append(answer, list(list(tree = t, err = mse_func(t))))  
      answer <-
        answer[order(sapply(answer, `[[`, "err"))]
      
      envir <- append(envir, list(list(tree = answer[[1]]$tree, error = answer[[1]]$err, iter = current_steps, time = Sys.time())))
    }
    if (current_steps == steps)
      break
  }
  
  return(envir)
}

path1 <- "source_domain.csv"
path2 <- "target_domain.csv"
path3 <- "KMM_results.xlsx"

data <- scanfiles(path1, path2, path3)

func <- list("+", "*", "-", "/")
vars <- as.list(colnames(data$source)[-ncol(data$source)])
trees <- generate_trees(func, vars, 20, data$source)

Sys.time()
env <- ITGP(data, 20, vars, func, trees, steps=100)