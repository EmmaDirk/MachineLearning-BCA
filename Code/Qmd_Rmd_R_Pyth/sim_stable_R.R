sim_stable_R <- function(timepoints = 5, #number of timepoints to keep
                       burnin = 45, #burnin period
                       N = 500, #sample size
                       ndat = 1, #number of datasets
                       phi, #matrix of autoregressive and cross-lagged effects
                       betac_cont, #beta matrix continuous confounders
                       betac_dich = NULL, #beta matrix dichotomous confounders
                       psi = NULL, #covariance matrix of the residuals of x and y
                       intercepts = c(0, 0), #intercepts of x and y
                       meansc = 0, #mean of the continuous confounders
                       sigmaC = 1, #variance of the continuous confounders (scalar or covariance matrix)
                       seed = NULL,
                       progress = FALSE) { #optional seed
  if (!is.null(seed)) {
    set.seed(seed) #set seed
  }
  isdich <- !is.null(betac_dich) #check if dichotomous confounders are present
  ncont <- ncol(betac_cont) #number of continuous confounders
  nconf <- ncont #number of confounders
  if (isdich) { #if dichotomous confounders are present
    ndich <- ncol(betac_dich) #number of dichotomous confounders
    nconf <- ncont + ndich #total number of confounders
    beta_full <- cbind(betac_cont, betac_dich) #combine betas of continuous and dichotomous confounders
  } else{ #if dichotomous confounders are not present
    beta_full <- beta_c_cont #set beta_full to beta_c_cont
    ndich <- 0 #set ndich to 0
  }
  if (length(meansc) != nconf) { #check if meansc is a scalar
    meansc1 <- meansc #save meansc
    meansc <- rep(meansc, times = nconf) #repeat meansc
  }
  len_sigmac <- (length(sigmaC) == 1) #check if sigmaC is a scalar
  if (len_sigmac) { #if sigmaC is a scalar
    sigmaC <- diag(rep(sqrt(sigmaC), times = ncont)) #set sigmaC to a diagonal matrix
  }
  if(class(psi) == "matrix"){
    psi <- list(psi)
  }
  if(length(psi) == 1){
    psi <- rep(psi, timepoints+1)
  }
  tot_timepoints <- burnin + timepoints #total number of timepoints
  ## create variable names for x and y
  varnames <- purrr::map(.x = c(-(burnin - 1):(timepoints)),
                         function(x)
                           purrr::map2(
                             .x = c("x", "y"),
                             .y = x,
                             .f = paste0
                           )) %>%
    unlist()
  covdat <- vector("list", length = ndat) #initialize list
  meandat <- vector("list", length = ndat)
  if(progress) cli::cli_progress_bar("Simulating Data", total = ndat)
  for(d in 1:ndat){
    conf <- MASS::mvrnorm(n = N, mu = meansc, Sigma = sigmaC) #generate confounders
    for(i in (ncont + 1):nconf){
      conf[,i] <- as.numeric(conf[,i] > 0)
    }
    xy <- matrix(NA, nrow = N, ncol = (timepoints+burnin)*2) #initialize xy
    for(n in 1:N){
      xy[n, 1:2] <- beta_full%*%conf[n,] + MASS::mvrnorm(n = 1, mu = c(0,0), Sigma = psi[[1]])
      for(t in 2:(timepoints+burnin)){
        if(t <= burnin){
          index_psi <- 1
        } else {
          index_psi <- t-burnin
        }
        xy[n, ((t-1)*2 + 1):(t*2)] <- phi%*%xy[n, ((t-2)*2 + 1):((t-1)*2)] + beta_full%*%conf[n,] + MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = psi[[index_psi]])
      }
    }
    colnames(xy) <- varnames
    covdat[[d]] <- cov(xy[,-(1:(burnin*2))])
    meandat[[d]] <- apply(xy[,-(1:(burnin*2))], 2, mean)
    cli::cli_progress_update()
  }
  ## create list of parameter inputs
  parameters <- list(
    timepoints = timepoints,
    burnin = burnin,
    N = N,
    ndat = ndat,
    phi = phi,
    betac_cont = betac_cont,
    betac_dich = betac_dich,
    psi = psi,
    intercepts = intercepts,
    meansc = meansc,
    sigmaC = sigmaC,
    seed = seed
  )
  ## create list with output
  datasets <- list(parameters = parameters, #parameters
                   covmats = covdat,
                   means = meandat) #datasets
  return(datasets) #return datasets
}
