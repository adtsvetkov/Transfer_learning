require(xlsx)
require(rgp)
require(stringi)

scanfiles <- function(path1, path2, path3)
{
  source <- as.matrix(read.csv(path1, row.names = 1))
  target <- as.matrix(read.csv(path2, row.names = 1))
  beta <- as.matrix(read.xlsx(path3, sheetName = "Beta")[, -1])
  indexes <- as.matrix(read.xlsx(path3, sheetName = "Indexes")[, -1])
  return(list(source = source, target = target, beta = beta, indexes=indexes))
}

generate_names <- function(cols, rows)
{
  colsnames <- vector()
  rowsnames <- vector()
  for (i in 1:cols) colsnames <- append(colsnames, paste("b", i, sep = "_"))
  for (j in 1:rows) rowsnames <- append(rowsnames, paste("t", j, sep = "_"))
  return(list(columns = colsnames, rows = rowsnames))
}

check_function <- function(f, x_s)
{
  for (i in 1:nrow(x_s))
  {
    xlist <- as.list(x_s[i, ][-ncol(x_s)])
    res <- do.call(f, xlist)
    if(res == "NaN" || res == "Inf") return(F)
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
    if (all(stri_detect_fixed(as_str, vars)) && check_function(myfunc, x_s))
      funclist <- append(funclist, myfunc)
    if (length(funclist) == p)
      break
  }
  return(funclist)
}

WMSE <- function(data, trees)
{
  weight <- data$beta
  indexes <- data$indexes
  x_s <- data$source
  wmse <- matrix(, nrow=length(trees), ncol=ncol(weight))
  for(t in 1:length(trees))
  {
    # w - индекс вектора бет
    for (w in 1:ncol(weight))
    {
      # взяли конкретные бета для конкретной функции
      beta <- weight[, w]
      index <- indexes[, w]
      summa <- vector()
      for(i in 1:length(index))
      {
        # номер экземпляра
        instance <- index[i]
        # аргументы функции этого экземпляра
        args <- x_s[instance, ][-length(x_s[instance, ])]
        # результаты вызова функции
        funcres <- do.call(trees[[t]], as.list(args))
        # посчитанная ошибка
        err <- beta[i]*((funcres-x_s[instance, length(x_s[instance])])^2)
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

findbest <- function(wmse, beta, index, trees)
{
  best <- list()
  for(t in 1:length(trees))
  {
    vec <- wmse[t, ]
    ind <- match(min(vec), vec)
    instance <- list(fvi = min(vec), beta = beta[, ind], func = trees[[t]], ind = index[, ind])
    best <- append(best, list(instance))
  }
  best <- best[order(sapply(best, `[[`, "fvi"))][-seq(ncol(beta) + 1, length(trees))]
  return(best)
}

MSE <- function(x_t, best)
{
  mse <- vector()
  for(i in 1:length(best))
  {
    res <- vector()
    for(j in 1:nrow(x_t))
    {
      args <- x_t[j, ][-length(x_t[j, ])]
      resfunc <- do.call(best[[i]]$func, as.list(args))
      result <- (resfunc - x_t[j, length(x_t[j, ])])^2
      res <- append(res, result)
    }
    mse <- append(mse, mean(res))
  }
  name <- generate_names(1, length(best))
  names(mse) <- name$rows
  return(mse)
}

ITGP <- function(data, p)
{
  beta <- data$beta
  
  func <- list("+", "*", "-", "/")
  vars <- as.list(colnames(data$source)[-ncol(data$source)])
  
  # сгенерировали деревья, p = 15 штук
  trees <- generate_trees(func, vars, p, data$source)
  # идем оценивать их на сорсах
  wmse <- WMSE(data, trees)
  # отбираем лучшие беты на каждом дереве:
  best <- findbest(wmse, beta, data$indexes, trees)
  # идем оценивать на таргетах
  mse <- MSE(data$target, best)
  print(mse)
  # размножаем деревья
  
  # генерим новое поколение бет на новых деревьях
}

path1 <- "source_domain.csv"
path2 <- "target_domain.csv"
path3 <- "KMM_results.xlsx"

data <- scanfiles(path1, path2, path3)
ITGP(data, 15)