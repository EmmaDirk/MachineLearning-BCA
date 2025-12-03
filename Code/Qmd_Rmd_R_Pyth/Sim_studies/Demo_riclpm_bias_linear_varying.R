############################################################
##  Libraries
############################################################

library(mvtnorm)
library(lavaan)
library(dplyr)
library(ggplot2)
library(parallel)

############################################################
##  Simulation function with SHRINK PARAMETER
##  Now returns both df and B_list
############################################################

simulate_clpm_timevarying <- function(
  N, T,
  ax = 0.2, ay = 0.2, bx = 0.10, by = 0.10, rho = 0.10,
  k = 3, gamma_x = c(0.15,0.10,0.05), gamma_y = c(0.15,0.10,0.05),
  tau2 = rep(1,3),
  mean_change_x_vec, mean_change_y_vec,
  sd_change, seed = NULL,
  shrink_param = 0.10
){
  if(!is.null(seed)) set.seed(seed)

  A   <- matrix(c(ax,by,bx,ay), 2, 2, byrow=TRUE)
  Psi <- diag(tau2)

  confounder_strength <- function(B, Psi){
    S_U <- B %*% Psi %*% t(B)
    list(
      S_U = S_U,
      var_X = S_U[1,1],
      var_Y = S_U[2,2]
    )
  }

  enforce_confounder_cap_auto <- function(B, Psi){
    cs <- confounder_strength(B, Psi)
    vX <- cs$var_X
    vY <- cs$var_Y

    if (vX <= 1 && vY <= 1) return(B)

    s <- min(
      if (vX > 0) sqrt(1/vX) else Inf,
      if (vY > 0) sqrt(1/vY) else Inf,
      1
    )

    if (!is.finite(s) || s <= 0) stop("Bad scaling factor")

    B * s
  }

  find_c <- function(A, S_U, rho){
    f <- function(c){
      S2 <- matrix(c(1,c,c,1),2,2)
      Sdyn <- S2 - S_U
      Sig  <- Sdyn - t(A)%*%Sdyn%*%A
      (Sig[1,2] / sqrt(Sig[1,1]*Sig[2,2])) - rho
    }
    tryCatch(uniroot(f, interval=c(-0.99,0.99))$root,
             error=function(e) NA_real_)
  }

  B0 <- rbind(gamma_x, gamma_y)
  B0 <- enforce_confounder_cap_auto(B0, Psi)

  shrink <- shrink_param

  B_list <- vector("list", T)
  B_list[[1]] <- B0

  for(t in 2:T){
    step <- rbind(
      rnorm(k, mean_change_x_vec[t], sd_change),
      rnorm(k, mean_change_y_vec[t], sd_change)
    )
    B_prev <- B_list[[t-1]]
    B_new  <- B_prev + step + shrink*(B0 - B_prev)
    B_list[[t]] <- enforce_confounder_cap_auto(B_new, Psi)
  }

  Sigma_e_list <- vector("list",T)
  for(t in 1:T){
    S_U <- B_list[[t]] %*% Psi %*% t(B_list[[t]])
    c_t <- find_c(A, S_U, rho)
    if (is.na(c_t)) return(NULL)

    S2 <- matrix(c(1,c_t,c_t,1),2,2)
    Sdyn <- (S2 - S_U)
    Sdyn <- (Sdyn + t(Sdyn))/2

    Sig <- Sdyn - t(A)%*%Sdyn%*%A
    Sig <- (Sig + t(Sig))/2

    if (any(diag(Sig) <= 0) || any(is.na(Sig))) return(NULL)
    Sigma_e_list[[t]] <- Sig
  }

  U <- tryCatch(rmvnorm(N, rep(0,k), Psi), error=function(e) NULL)
  if(is.null(U)) return(NULL)

  Ddyn <- tryCatch(rmvnorm(N, c(0,0), Sigma_e_list[[1]]),
                   error=function(e) NULL)
  if(is.null(Ddyn)) return(NULL)

  df <- matrix(NA, nrow=N, ncol=2*T + k)
  colnames(df) <- c(
    paste0("x",1:T),
    paste0("y",1:T),
    paste0("c",1:k)
  )
  df[, (2*T+1):(2*T+k)] <- U

  obs1 <- Ddyn + U %*% t(B_list[[1]])
  df[,1]   <- obs1[,1]
  df[,1+T] <- obs1[,2]

  for(i in 2:T){
    Ddyn_step <- tryCatch(rmvnorm(N, sigma=Sigma_e_list[[i]]),
                          error=function(e) NULL)
    if(is.null(Ddyn_step)) return(NULL)

    Ddyn <- Ddyn %*% t(A) + Ddyn_step
    obs  <- Ddyn + U %*% t(B_list[[i]])

    df[,i]     <- obs[,1]
    df[,i + T] <- obs[,2]
  }

  list(
    df = df,
    B_list = B_list
  )
}



############################################################
##  lavaan helpers
############################################################

build_riclpm_model <- function(T){

  ri_block <- paste0(
    "rix =~ ", paste(paste0("1*x",1:T), collapse=" + "), "\n",
    "riy =~ ", paste(paste0("1*y",1:T), collapse=" + "), "\n",
    "rix ~~ rix\nriy ~~ riy\nrix ~~ riy\n"
  )

  resid_block <- paste0(
    paste(paste0("x",1:T," ~~ 0*x",1:T), collapse="; "), "\n",
    paste(paste0("y",1:T," ~~ 0*y",1:T), collapse="; "), "\n"
  )

  within_lat <- paste0(
    paste(paste0("wx",1:T," =~ 1*x",1:T), collapse="; "), "\n",
    paste(paste0("wy",1:T," =~ 1*y",1:T), collapse="; "), "\n"
  )

  orth <- paste0(
    "rix ~~ ", paste(paste0("0*wx",1:T), collapse=" + "), "\n",
    "rix ~~ ", paste(paste0("0*wy",1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(paste0("0*wx",1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(paste0("0*wy",1:T), collapse=" + "), "\n"
  )

  within_var <- paste0(
    paste(paste0("wx",1:T," ~~ wx",1:T), collapse="; "), "\n",
    paste(paste0("wy",1:T," ~~ wy",1:T), collapse="; "), "\n"
  )

  within_cov <- paste0(
    paste(paste0("wy",1:T," ~~ wx",1:T), collapse="; "), "\n"
  )

  regress <- paste(
    unlist(lapply(2:T, function(t){
      c(
        paste0("wx",t," ~ wx",t-1," + wy",t-1),
        paste0("wy",t," ~ wx",t-1," + wy",t-1)
      )
    })), collapse="\n"
  )

  means <- paste0(
    paste(paste0("x",1:T), collapse=" + "), " ~ mx*1\n",
    paste(paste0("y",1:T), collapse=" + "), " ~ my*1\n"
  )

  paste(
    ri_block,
    resid_block,
    within_lat,
    orth,
    within_var,
    within_cov,
    regress,
    means,
    sep="\n"
  )
}

fit_riclpm <- function(data, model_string){
  lavaan(model_string, data=as.data.frame(data), estimator="MLR")
}

extract_crosslags_riclpm <- function(fit, T){
  pe <- parameterEstimates(fit)
  out <- c()

  for(t in 2:T){
    out <- c(
      out,
      pe$est[pe$lhs == paste0("wy",t) & pe$rhs == paste0("wx",t-1)],
      pe$est[pe$lhs == paste0("wx",t) & pe$rhs == paste0("wy",t-1)]
    )
  }
  out
}

############################################################
##  run_one: bias, est_mean, sd_B, mean_B, NA indicator
############################################################

run_one <- function(var_value, rep_id,
                    N, T,
                    mean_change_x_vec,
                    mean_change_y_vec,
                    model_string,
                    shrink_param,
                    true_crosslag = 0.1){

  sim <- simulate_clpm_timevarying(
    N=N, T=T,
    mean_change_x_vec=mean_change_x_vec,
    mean_change_y_vec=mean_change_y_vec,
    sd_change=sqrt(var_value),
    seed = 111111 + as.integer(var_value*1e6) + rep_id,
    shrink_param = shrink_param
  )

  if (is.null(sim)) {
    return(data.frame(
      var_value = var_value,
      bias      = NA,
      var_est   = NA,
      est_mean  = NA,
      sd_B      = NA,
      mean_B    = NA,
      is_na_run = 1
    ))
  }

  df     <- sim$df
  B_list <- sim$B_list

  fit <- tryCatch(fit_riclpm(df, model_string), error=function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(
      var_value = var_value,
      bias      = NA,
      var_est   = NA,
      est_mean  = NA,
      sd_B      = NA,
      mean_B    = NA,
      is_na_run = 1
    ))
  }

  ests <- tryCatch(extract_crosslags_riclpm(fit, T), error=function(e) NA)

  if (all(is.na(ests))) {
    return(data.frame(
      var_value = var_value,
      bias      = NA,
      var_est   = NA,
      est_mean  = NA,
      sd_B      = NA,
      mean_B    = NA,
      is_na_run = 1
    ))
  }

  # bias & var of cross-lags
  bias_vals <- abs(ests - true_crosslag)
  bias      <- mean(bias_vals, na.rm=TRUE)
  var_est   <- var(ests, na.rm=TRUE)
  est_mean  <- mean(ests, na.rm=TRUE)

  # SD of B coefficients across time + mean B
  # B_list: list of length T, each 2 x k
  B_array <- simplify2array(B_list)  # 2 x k x T
  sd_array <- apply(B_array, c(1,2), sd, na.rm=TRUE)  # 2 x k
  avg_sd_B <- mean(sd_array, na.rm=TRUE)
  mean_B_all <- mean(B_array, na.rm=TRUE)

  data.frame(
    var_value = var_value,
    bias      = bias,
    var_est   = var_est,
    est_mean  = est_mean,
    sd_B      = avg_sd_B,
    mean_B    = mean_B_all,
    is_na_run = 0
  )
}

############################################################
##  Main simulation wrapper
############################################################

run_full_bias_variance <- function(
  var_list          = seq(0.001, 0.15, by=0.001),
  reps              = 200,
  N                 = 10000,
  T                 = 5,
  mean_change_x_vec = rep(0,5),
  mean_change_y_vec = rep(0,5),
  shrink_param      = 0.10,
  cores             = detectCores() - 1
){

  start_time <- Sys.time()

  tasks <- expand.grid(
    var_value = var_list,
    rep_id    = 1:reps,
    KEEP.OUT.ATTRS = FALSE
  )
  total_tasks <- nrow(tasks)

  model_string <- build_riclpm_model(T)

  cl <- makeCluster(cores)
  clusterExport(cl, c(
    "simulate_clpm_timevarying",
    "fit_riclpm",
    "extract_crosslags_riclpm",
    "run_one",
    "model_string",
    "N","T",
    "mean_change_x_vec","mean_change_y_vec",
    "shrink_param"
  ), envir=environment())

  clusterEvalQ(cl, {
    library(mvtnorm)
    library(lavaan)
  })

  pb <- txtProgressBar(min=0, max=total_tasks, style=3)
  progress_counter <- 0
  next_mark <- 0.10

  indices <- 1:total_tasks
  chunk_size <- max(1, floor(total_tasks/50))
  chunks <- split(indices, ceiling(seq_along(indices) / chunk_size))

  results_list <- vector("list", length(chunks))

  for(ci in seq_along(chunks)){
    idx_vec <- chunks[[ci]]

    res_chunk <- parLapplyLB(
      cl,
      idx_vec,
      function(i){
        row <- tasks[i,]
        run_one(
          var_value=row$var_value,
          rep_id=row$rep_id,
          N=N,
          T=T,
          mean_change_x_vec=mean_change_x_vec,
          mean_change_y_vec=mean_change_y_vec,
          model_string=model_string,
          shrink_param=shrink_param
        )
      }
    )

    results_list[[ci]] <- res_chunk

    progress_counter <- progress_counter + length(idx_vec)
    setTxtProgressBar(pb, progress_counter)

    frac <- progress_counter / total_tasks
    if (frac >= next_mark) {
      elapsed <- as.numeric(Sys.time() - start_time)
      cat(sprintf("\n[%d%% completed] elapsed: %.1f sec\n",
                  round(next_mark * 100), elapsed))
      next_mark <- next_mark + 0.10
    }
  }

  close(pb)
  stopCluster(cl)

  raw_df <- bind_rows(unlist(results_list, recursive=FALSE))

  na_runs_total <- sum(raw_df$is_na_run)

  coef_summary <- raw_df %>%
    group_by(var_value) %>%
    summarise(
      mean_estimate = mean(est_mean, na.rm=TRUE),
      mean_bias     = mean(bias, na.rm=TRUE),
      sd_B          = mean(sd_B, na.rm=TRUE),
      mean_B        = mean(mean_B, na.rm=TRUE),
      na_runs       = sum(is_na_run),
      total_runs    = n(),
      na_rate       = na_runs / total_runs
    )

  # Plot 1: bias vs variance
  p_bias <- ggplot(coef_summary, aes(x=var_value, y=mean_bias)) +
    geom_line(size=1.1) +
    geom_point(size=2) +
    ylab("Mean absolute bias") +
    xlab("Variance")

  # Plot 2: estimated coefficient vs variance
  p_coef <- ggplot(coef_summary, aes(x=var_value, y=mean_estimate)) +
    geom_line(size=1.1, color="blue") +
    geom_point(size=2, color="blue") +
    ylab("Mean estimated cross-lag") +
    xlab("Variance")

  # Plot 3: NA rate vs variance
  p_na <- ggplot(coef_summary, aes(x=var_value, y=na_rate)) +
    geom_line(size=1.1, color="red") +
    geom_point(size=2, color="red") +
    ylab("NA rate") +
    xlab("Variance")

  # Plot 4: SD(B) vs variance
  p_sdB <- ggplot(coef_summary, aes(x=var_value, y=sd_B)) +
    geom_line(size=1.1, color="darkblue") +
    geom_point(size=2, color="darkblue") +
    ylab("SD of B across time") +
    xlab("Variance")

  # Plot 5: mean(B) vs variance
  p_meanB <- ggplot(coef_summary, aes(x=var_value, y=mean_B)) +
    geom_line(size=1.1, color="darkgreen") +
    geom_point(size=2, color="darkgreen") +
    ylab("Mean B across time & entries") +
    xlab("Variance")

  total_runtime <- Sys.time() - start_time

  list(
    raw        = raw_df,
    summary    = coef_summary,
    plot_bias  = p_bias,
    plot_coef  = p_coef,
    plot_na    = p_na,
    plot_sd_B  = p_sdB,
    plot_mean_B = p_meanB,
    runtime    = total_runtime,
    na_runs    = na_runs_total
  )
}

############################################################
##  Shrinkage sweep over variances:
##  x = variance, y = NA rate, color = shrinkage
############################################################

run_shrinkage_sweep <- function(
    variance_values = seq(0.005, 0.05, by = 0.005),
    shrink_vals     = seq(0, 0.4, by = 0.1),
    reps = 80,
    N = 10000,
    T = 5,
    mean_change_x_vec = rep(0,5),
    mean_change_y_vec = rep(0,5),
    cores = detectCores() - 1
){

    results <- data.frame(
        variance   = numeric(0),
        shrink     = numeric(0),
        na_runs    = numeric(0),
        na_rate    = numeric(0),
        total_runs = numeric(0)
    )

    for (s in shrink_vals) {

        cat("\n=== Shrinkage:", s, "===\n")

        for (v in variance_values) {

            out <- tryCatch(
                run_full_bias_variance(
                    var_list          = v,
                    reps              = reps,
                    N                 = N,
                    T                 = T,
                    mean_change_x_vec = mean_change_x_vec,
                    mean_change_y_vec = mean_change_y_vec,
                    shrink_param      = s,
                    cores             = cores
                ),
                error = function(e) NULL
            )

            if (is.null(out) || is.null(out$na_runs) || length(out$na_runs) != 1) {
                na_here <- NA
            } else {
                na_here <- out$na_runs
            }

            results <- rbind(
                results,
                data.frame(
                    variance   = v,
                    shrink     = s,
                    na_runs    = na_here,
                    na_rate    = na_here / reps,
                    total_runs = reps
                )
            )
        }
    }

    p <- ggplot(results, aes(x = variance, y = na_rate, color = factor(shrink))) +
        geom_line(size=1.3) +
        geom_point(size=2) +
        labs(
            x = "Variance",
            y = "NA rate",
            color = "Shrinkage",
            title = "NA rate vs variance for different shrinkage values"
        ) +
        theme_minimal(base_size = 14)

    list(results = results, plot = p)
}

############################################################
##  EXAMPLE RUNS (commented)
############################################################

### Example 1: Basic variance sweep at one shrinkage
 out <- run_full_bias_variance(
   var_list     = seq(0.001, 0.20, by=0.001),
   reps         = 1000,
   N            = 10000,
   T            = 5,
   shrink_param = 0.20,
   cores        = 7
 )
 out$plot_bias
 out$plot_coef
 out$plot_na
 out$plot_sd_B
 out$plot_mean_B
 out$runtime
 out$summary

### Example 2: Shrinkage sweep (multiple lines)
sweep_res <- run_shrinkage_sweep(
   variance_values = seq(0.001, 0.15, by = 0.001),
   shrink_vals     = seq(0, 0.5, by = 0.1),
   reps            = 50,
   N               = 10000,
   T               = 5,
   cores           = 7
 )
 sweep_res$plot



