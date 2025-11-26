run_full_bias_variance_windows <- function(
  var_list       = seq(0.001, 0.03, by=0.001),
  reps           = 200,
  N              = 10000,
  T              = 5,
  mean_change_x_vec = c(0, 0, 0, 0, 0),
  mean_change_y_vec = c(0, 0, 0, 0, 0),
  cores          = parallel::detectCores() - 1
){

  library(mvtnorm)
  library(lavaan)
  library(tidyverse)
  library(parallel)

  # ============================================================
  # Helpers inside wrapper
  # ============================================================

  simulate_clpm_timevarying <- function(
    N, T, ax=0.2, ay=0.2, bx=0.10, by=0.10, rho=0.10,
    k=3, gamma_x=c(0.15,0.10,0.05), gamma_y=c(0.15,0.10,0.05),
    tau2=rep(1,3),
    mean_change_x_vec, mean_change_y_vec,
    sd_change, seed
  ){
    if(!is.null(seed)) set.seed(seed)

    A   <- matrix(c(ax,by,bx,ay),2,2,byrow=TRUE)
    Psi <- diag(tau2)
    B   <- rbind(gamma_x, gamma_y)

    # Time-varying confounder effects (random walk, with clipped changes)
    B_list <- vector("list",T)
    B_list[[1]] <- B

    for(t in 2:T){
      # draw changes from N(mean_change, sd_change^2)
      change_x_raw <- rnorm(3, mean_change_x_vec[t], sd_change)
      change_y_raw <- rnorm(3, mean_change_y_vec[t], sd_change)

      # variance of the change distribution
      var_change <- sd_change^2
      bound <- 2 * var_change  # ± 2 * variance

      # clamp to [-bound, bound]
      change_x <- pmax(pmin(change_x_raw,  bound), -bound)
      change_y <- pmax(pmin(change_y_raw,  bound), -bound)

      B_list[[t]] <- rbind(
        B_list[[t-1]][1,] + change_x,
        B_list[[t-1]][2,] + change_y
      )
    }

    # equation solver for c_t, robust to failures
    find_c <- function(A, S_U, rho){
      f <- function(c){
        S2   <- matrix(c(1,c,c,1),2,2)
        Sdyn <- S2 - S_U
        Sig  <- Sdyn - t(A)%*%Sdyn%*%A
        (Sig[1,2] / sqrt(Sig[1,1]*Sig[2,2])) - rho
      }
      out <- tryCatch(
        uniroot(f, interval=c(-0.99,0.99))$root,
        error = function(e) NA_real_
      )
      out
    }

    # Innovation matrices; if anything goes wrong, return NULL
    Sigma_e_list <- vector("list",T)
    for(t in 1:T){
      S_U <- B_list[[t]] %*% Psi %*% t(B_list[[t]])
      c_t <- find_c(A, S_U, rho)
      if(is.na(c_t)) return(NULL)

      S2   <- matrix(c(1,c_t,c_t,1),2,2)
      Sdyn <- S2 - S_U
      Sdyn <- (Sdyn + t(Sdyn))/2

      Sig <- Sdyn - t(A)%*%Sdyn%*%A
      Sig <- (Sig + t(Sig))/2

      # basic sanity check: variances must be positive
      if(any(diag(Sig) <= 0) || any(is.na(Sig))) return(NULL)

      Sigma_e_list[[t]] <- Sig
    }

    Sigma_e1 <- Sigma_e_list[[1]]

    # rmvnorm can also fail -> catch
    U <- tryCatch(
      rmvnorm(N, rep(0,3), Psi),
      error = function(e) NULL
    )
    if(is.null(U)) return(NULL)

    Ddyn <- tryCatch(
      rmvnorm(N, c(0,0), Sigma_e1),
      error = function(e) NULL
    )
    if(is.null(Ddyn)) return(NULL)

    df <- matrix(NA, nrow=N, ncol=2*T + 3)
    colnames(df) <- c(paste0("x",1:T),
                      paste0("y",1:T),
                      paste0("c",1:3))
    df[, (2*T+1):(2*T+3)] <- U

    obs1 <- Ddyn + U %*% t(B_list[[1]])
    df[,1]   <- obs1[,1]
    df[,1+T] <- obs1[,2]

    for(i in 2:T){
      Ddyn_step <- tryCatch(
        rmvnorm(N, sigma=Sigma_e_list[[i]]),
        error = function(e) NULL
      )
      if(is.null(Ddyn_step)) return(NULL)

      Ddyn <- Ddyn %*% t(A) + Ddyn_step
      obs  <- Ddyn + U %*% t(B_list[[i]])
      df[,i]   <- obs[,1]
      df[,i+T] <- obs[,2]
    }

    df
  }

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

  fit_riclpm <- function(data, T){
    model <- build_riclpm_model(T)
    # lavaan can also fail -> catch in run_one
    lavaan(model, data=as.data.frame(data), estimator="MLR")
  }

  extract_crosslags_riclpm <- function(fit, T){
    pe <- parameterEstimates(fit, standardized=TRUE)
    out <- c()

    for(t in 2:T){
      wx_t   <- paste0("wx",t)
      wy_t   <- paste0("wy",t)
      wx_tm1 <- paste0("wx",t-1)
      wy_tm1 <- paste0("wy",t-1)

      r_xy <- pe[pe$op=="~" & pe$lhs==wy_t & pe$rhs==wx_tm1,]
      if(nrow(r_xy)==1) out <- c(out, r_xy$est)

      r_yx <- pe[pe$op=="~" & pe$lhs==wx_t & pe$rhs==wy_tm1,]
      if(nrow(r_yx)==1) out <- c(out, r_yx$est)
    }

    out
  }

  run_one <- function(var_value, rep_id){
    df <- simulate_clpm_timevarying(
      N=N, T=T,
      mean_change_x_vec=mean_change_x_vec,
      mean_change_y_vec=mean_change_y_vec,
      sd_change=sqrt(var_value),
      seed = 111111 + as.integer(var_value*1e6) + rep_id
    )

    # if simulation blew up, return NA row
    if(is.null(df)) {
      return(data.frame(
        var_value = var_value,
        bias      = NA_real_,
        var_est   = NA_real_
      ))
    }

    fit <- tryCatch(
      fit_riclpm(df, T),
      error = function(e) NULL
    )
    if(is.null(fit)) {
      return(data.frame(
        var_value = var_value,
        bias      = NA_real_,
        var_est   = NA_real_
      ))
    }

    ests <- tryCatch(
      extract_crosslags_riclpm(fit, T),
      error = function(e) NA_real_
    )

    if(all(is.na(ests))) {
      return(data.frame(
        var_value = var_value,
        bias      = NA_real_,
        var_est   = NA_real_
      ))
    }

    avg_est <- mean(ests, na.rm=TRUE)
    bias    <- abs(0.1 - avg_est)
    var_est <- var(ests, na.rm=TRUE)

    data.frame(
      var_value = var_value,
      bias      = bias,
      var_est   = var_est
    )
  }

  # ============================================================
  # Parallel execution with chunked progress bar
  # ============================================================

  tasks <- expand.grid(
    var_value = var_list,
    rep_id    = 1:reps,
    KEEP.OUT.ATTRS = FALSE
  )

  total_tasks <- nrow(tasks)

  cl <- makeCluster(cores)

  # export everything
  clusterExport(cl, varlist=ls(environment()), envir=environment())
  clusterEvalQ(cl, {
    library(mvtnorm)
    library(lavaan)
    library(tidyverse)
  })

  # progress bar
  pb <- txtProgressBar(min = 0, max = total_tasks, style = 3)
  progress_counter <- 0

  # split into chunks so we can update progress as chunks finish
  indices    <- 1:total_tasks
  chunk_size <- max(1, floor(total_tasks / 50))  # about 50 updates
  chunks     <- split(indices, ceiling(seq_along(indices) / chunk_size))

  results_list <- vector("list", length(chunks))

  for(ci in seq_along(chunks)){
    idx_vec <- chunks[[ci]]

    res_chunk <- parLapply(
      cl,
      idx_vec,
      function(i){
        row <- tasks[i,]
        run_one(row$var_value, row$rep_id)
      }
    )

    results_list[[ci]] <- res_chunk

    progress_counter <- progress_counter + length(idx_vec)
    setTxtProgressBar(pb, progress_counter)
  }

  close(pb)
  stopCluster(cl)

  results <- unlist(results_list, recursive = FALSE)
  mc_df   <- bind_rows(results)

  # ============================================================
  # Summaries + plot
  # ============================================================

  plot_df <- mc_df %>%
    group_by(var_value) %>%
    summarise(
      bias_mean = mean(bias, na.rm=TRUE),
      bias_low  = quantile(bias, 0.025, na.rm=TRUE),
      bias_high = quantile(bias, 0.975, na.rm=TRUE)
    )

  p <- ggplot(plot_df, aes(x=var_value, y=bias_mean)) +
    geom_ribbon(aes(ymin=bias_low, ymax=bias_high), alpha=0.25) +
    geom_line(size=1.2) +
    geom_point(size=2) +
    labs(
      x = "Variance of time-varying confounder effects",
      y = "Bias (0.1 - average cross-lag)"
    )

  list(
    raw     = mc_df,
    summary = plot_df,
    plot    = p
  )
}

out <- run_full_bias_variance_windows(
  reps  = 200,
  cores = 7
)

out$plot
