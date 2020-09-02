# creates a random coefficient matrix with a given number of species
generate_coef_mat <- function(n_spp, distribution = "rnorm", disp = 1){
  if (distribution == "rnorm") {
    IM_vals <- rnorm(n = n_spp*n_spp, mean = 0, sd = disp)
  } else if (distribution == "runif") {
    IM_vals <- runif(n = n_spp*n_spp, min = -disp, max = disp)
  } else {
    ### gives an error message instead of a warning (so that it won't try to access IM) 
    stop("generate_coef_mat must have a distribution of 'rnorm' or 'runif'" )
  }
  IM <- matrix(data = IM_vals, nrow = n_spp, ncol = n_spp)
  return(IM)
}

# takes a coefficient matrix and returns a list with mutualism indices 
get_indices <- function(mat){
  
  pos_ind <- which(mat > 0)
  pos_t_ind <- which(t(mat) > 0) 
  
  # index values that are positive in both original matrix and transpose
  mut_ind <- pos_ind[which((pos_ind) %in% (pos_t_ind))]
  
  return(mut_ind)
}

# takes a coefficient matrix and creates a matrix with only mutualisms
get_mut_mat <- function(mat){
  # get mutualism indices
  mut_ind <- get_indices(mat)
  
  # Create a blank matrix and fill with only mutualisms
  mut_mat <- matrix(data = 0, nrow = nrow(mat), ncol = nrow(mat))
  mut_mat[mut_ind] <- mat[mut_ind]
  
  # remove the values on the diagonal (in the previous steps these are counted as mutualisms)
  diag(mut_mat) <- 0
  
  return (mut_mat)
}

# generates a coefficient matrix that contains at least one mutualistic interaction
generate_coef_mat2 <- function(n_spp, distribution = "rnorm", disp = 1){
  # generate a random coefficient matrix
  mat <- generate_coef_mat(n_spp, distribution, disp)
  mut_mat <- get_mut_mat(mat)
  n_mut <- nnzero(mut_mat)
  while(n_mut == 0){
    mat <- generate_coef_mat(n_spp, distribution, disp)
    mut_mat <- get_mut_mat(mat)
    n_mut <- nnzero(mut_mat)
  }
  
  diag(mat) <- 0
  return(mat)
}

# takes a coefficient matrix and returns a matrix with only non-mutualism interactions
get_not_mut_mat <- function(mat){
  # get indices for mutualism
  mut_ind <- get_indices(mat)
  
  # index values not in the mutualism matrix
  all_ind <- c(1:nrow(mat)^2)
  not_mut_ind <- all_ind[which(!(all_ind %in% mut_ind))]
  
  # Create a blank matrix and fill with only non-mutualistic interactions
  not_mut_mat <- matrix(data = 0, nrow = nrow(mat), ncol = nrow(mat))
  not_mut_mat[not_mut_ind] <- mat[not_mut_ind]
  
  # make diagonals 0 (don't want intraspecific effects since we're using K)
  diag(not_mut_mat) <- 0 
  
  return(not_mut_mat)
}

# given desolve output, find out how many species persisted
get_num_persist <- function(out){
  num_persist <- 0
  for (sp in out[nrow(out), 2:ncol(out)]){
    if(!is.na(sp)){
      if (sp > 0){
        num_persist <- num_persist + 1
      }}}
  return(num_persist)
}

# tryCatch error handler for exponential growth
errorFunc <- function(e){
  unbounded <- integer(1)
  return(unbounded)
}

## event triggered if state variable <= 0.01 (extinction) or if state variable >= 1000 (unbounded growth)
rootfun <- function (t, n, parms) {
  n1 <- rep(NA, length(n))
  for(i in seq(length(n))){
    if( n[i] > 10){
      n1[i] <- n[i]-1000
    }
    else{
      n1[i] <- n[i] - eps
    }
  }
  return(as.vector(n1))
}


## sets state variable = 0 if extinct or stops executing and produces an error if unbounded growth is happening
eventfun <- function(t, n, parms) {
  for (i in seq(length(n))){
    if (n[i] <= eps){
      n[i] <-0 
    } 
    else if (n[i] >= 1000){
      stop("error: unbounded growth")
    }
  }
  return(as.vector(n))
}

# run ode solver with tryCatch, return deSolve object if no exp growth, return integer(1) if exp growth
odeSolve <- function (y, times, func, parms, method, rootfun, events){
  result <- tryCatch({
    ode(y = y, times = times, func = func, parms = parms, method = method,
        rootfun = rootfun, events = events)
  }, error = errorFunc)
  return(result)
}

# return true if there was exponential growth
expGrowth <- function(out){
  if(is.integer(out)){
    return(TRUE)
  }
  return(FALSE)
}

# count up extinctions and persistances (NA if there was exp growth)
#returns a list of persistances and extinctions
sumExtPer <- function(out, expGrowth, n_spp){
  if (expGrowth){
    extinct <- NA
    persist <- NA
  } else{
    persist <- get_num_persist(out)
    extinct <- n_spp - persist
  }
  return(list(persist, extinct))
}

# numerator (intraspecific benefit) mutualism model
NumeratorMut <- function(t, n, parms){
  with(as.list(c(n, parms)), {
    dns <- r*n*((k - n + mut_mat%*%n + not_mut_mat%*%n)/k)
    return(list(c(dns)))
  })
}

# denominator (equilibrium-benefit) mutualism model
DenominatorMut <- function(t, n, parms){
  with(as.list(c(n, parms)), {
    dns <- r*n*((k - n + mut_mat%*%n + not_mut_mat%*%n)/(k + mut_mat%*%n))
    return(list(c(dns)))
  })
}