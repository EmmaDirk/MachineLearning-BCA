library(pbapply)
library(parallel)

# -------------------------
# Settings: keep identical on both machines
# -------------------------
n_tasks <- 200
task_size <- 3000
n_repeats <- 50   # number of benchmark runs

# -------------------------
# Detect workers
# -------------------------
n_cores <- detectCores(logical = TRUE)
n_workers <- max(1, n_cores - 1)

cat("Detected cores :", n_cores, "\n")
cat("Using workers  :", n_workers, "\n")
cat("Tasks per run  :", n_tasks, "\n")
cat("Task size      :", task_size, "\n")
cat("Repeats        :", n_repeats, "\n\n")

# -------------------------
# CPU-heavy test task
# Replace this with your own real function if needed
# -------------------------
benchmark_task <- function(i, task_size) {
  set.seed(i)
  
  x <- matrix(rnorm(task_size * task_size / 10), ncol = task_size / 10)
  y <- crossprod(x)
  z <- eigen(y, symmetric = TRUE, only.values = TRUE)$values
  
  sum(z)
}

# -------------------------
# One benchmark run
# -------------------------
run_once <- function(run_id, n_tasks, task_size, n_workers) {
  cat("========== Run", run_id, "of", n_repeats, "==========\n")
  
  cl <- makeCluster(n_workers)
  on.exit(stopCluster(cl), add = TRUE)
  
  clusterExport(cl, varlist = c("benchmark_task", "task_size"), envir = environment())
  
  timing <- system.time({
    results <- pblapply(
      X = 1:n_tasks,
      FUN = function(i) benchmark_task(i, task_size),
      cl = cl
    )
  })
  
  elapsed <- unname(timing["elapsed"])
  
  cat("User time   :", unname(timing["user.self"]), "\n")
  cat("System time :", unname(timing["sys.self"]), "\n")
  cat("Elapsed     :", elapsed, "seconds\n\n")
  
  return(elapsed)
}

# -------------------------
# Repeat benchmark
# -------------------------
elapsed_times <- numeric(n_repeats)

for (r in seq_len(n_repeats)) {
  elapsed_times[r] <- run_once(r, n_tasks, task_size, n_workers)
}

# -------------------------
# Summary
# -------------------------
avg_time <- mean(elapsed_times)
median_time <- median(elapsed_times)
min_time <- min(elapsed_times)
max_time <- max(elapsed_times)
sd_time <- sd(elapsed_times)

cat("=====================================\n")
cat("Benchmark summary\n")
cat("=====================================\n")
cat("Elapsed times :", paste(round(elapsed_times, 3), collapse = ", "), "\n")
cat("Average       :", round(avg_time, 3), "seconds\n")
cat("Median        :", round(median_time, 3), "seconds\n")
cat("Min           :", round(min_time, 3), "seconds\n")
cat("Max           :", round(max_time, 3), "seconds\n")
cat("Std. dev.     :", round(sd_time, 3), "seconds\n")