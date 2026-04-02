# 09_model_set_helpers.R
# These helpers define how a user-supplied set of model requests is translated
# into the most efficient execution plan for the simulation study.
#
# The key idea is:
# - the user specifies a list of model specifications
# - we normalize those specifications into one common format
# - we identify which models can share the same stage-1 residualised data
# - we identify which XGB models can also share the same one-time tuning object
#
# This lets the top-level runner:
# 1) simulate one data set per replication
# 2) residualise only once per unique stage-1 recipe
# 3) fit every compatible SEM on that same prepared data
# 4) do the same sharing logic again inside each bootstrap draw
# -------------------------------------------------------------------------------------------------

# convert NULL to an empty list where helpful
null_to_empty_list <- function(x) {
  if (is.null(x)) list() else x
}


# sort exclude vectors into a stable canonical format
normalize_exclude_vector <- function(exclude) {

  if (is.null(exclude) || length(exclude) == 0) {
    return(NULL)
  }

  exclude <- as.character(exclude)
  exclude <- unique(exclude)
  exclude <- sort(exclude)

  if (length(exclude) == 0) {
    return(NULL)
  }

  exclude
}


# recursively sort named lists so equality checks do not depend on argument order
canonicalize_object <- function(x) {

  if (is.null(x)) {
    return(NULL)
  }

  # keep data frames and matrices as they are
  if (is.data.frame(x) || is.matrix(x)) {
    return(x)
  }

  # recursively canonicalize ordinary lists
  if (is.list(x)) {

    # canonicalize every element first
    x <- lapply(x, canonicalize_object)

    # if the list is named, sort by name
    if (!is.null(names(x))) {
      x <- x[order(names(x))]
    }

    return(x)
  }

  x
}


# compare two objects after canonicalization
same_canonical_object <- function(x, y) {
  identical(canonicalize_object(x), canonicalize_object(y))
}


# convenient constructor for one model specification
make_model_spec <- function(
    name = NULL,
    residualizer = c("none", "linear", "xgb"),
    sem_model = c("clpm", "riclpm", "dpm"),
    confounder_order = 1,
    exclude = NULL,
    free_loadings = FALSE,
    bootstrap_B = 50,
    xgb_tuning = NULL,
    tune_xgb = TRUE,
    xgb_tune_args = list(),
    residualizer_args = list()
) {

  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  list(
    name = name,
    residualizer = residualizer,
    sem_model = sem_model,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    xgb_tuning = xgb_tuning,
    tune_xgb = tune_xgb,
    xgb_tune_args = xgb_tune_args,
    residualizer_args = residualizer_args
  )
}


# normalize one user-supplied model specification into a stable internal format
normalize_model_spec <- function(spec, idx = NULL) {

  # allow the user to pass an unnamed list with the same fields
  if (!is.list(spec)) {
    stop("Each model specification must be a list, for example created by make_model_spec().")
  }

  # defaults
  default <- list(
    name = NULL,
    residualizer = "none",
    sem_model = "clpm",
    confounder_order = 1L,
    exclude = NULL,
    free_loadings = FALSE,
    bootstrap_B = 50L,
    xgb_tuning = NULL,
    tune_xgb = TRUE,
    xgb_tune_args = list(),
    residualizer_args = list()
  )

  # fill missing entries from defaults
  for (nm in names(default)) {
    if (is.null(spec[[nm]]) && !(nm %in% names(spec))) {
      spec[[nm]] <- default[[nm]]
    }
  }

  # match the two required choices
  spec$residualizer <- match.arg(spec$residualizer, c("none", "linear", "xgb"))
  spec$sem_model <- match.arg(spec$sem_model, c("clpm", "riclpm", "dpm"))

  # fill remaining defaults explicitly
  if (is.null(spec$name)) {
    spec$name <- paste0("model_", if (is.null(idx)) 1L else as.integer(idx))
  }

  if (!is.character(spec$name) || length(spec$name) != 1 || !nzchar(spec$name)) {
    stop("Each model specification must have a non-empty character 'name'.")
  }

  if (is.null(spec$confounder_order)) {
    spec$confounder_order <- default$confounder_order
  }
  spec$confounder_order <- as.integer(spec$confounder_order[1])

  if (is.na(spec$confounder_order) || !(spec$confounder_order %in% c(0L, 1L, 2L, 3L))) {
    stop("Each model specification must have confounder_order in {0, 1, 2, 3}.")
  }

  spec$exclude <- normalize_exclude_vector(spec$exclude)
  spec$free_loadings <- isTRUE(spec$free_loadings)

  if (is.null(spec$bootstrap_B)) {
    spec$bootstrap_B <- default$bootstrap_B
  }
  spec$bootstrap_B <- as.integer(spec$bootstrap_B[1])

  if (is.na(spec$bootstrap_B) || spec$bootstrap_B < 0L) {
    stop("Each model specification must have a non-negative integer bootstrap_B.")
  }

  spec$tune_xgb <- isTRUE(spec$tune_xgb)
  spec$xgb_tune_args <- null_to_empty_list(spec$xgb_tune_args)
  spec$residualizer_args <- null_to_empty_list(spec$residualizer_args)

  if (!is.list(spec$xgb_tune_args)) {
    stop("Each model specification must have 'xgb_tune_args' as a list.")
  }

  if (!is.list(spec$residualizer_args)) {
    stop("Each model specification must have 'residualizer_args' as a list.")
  }

  # group ids are assigned later
  spec$xgb_tuning_group_id <- NA_integer_
  spec$stage1_group_id <- NA_integer_

  spec
}


# normalize a whole list of model specifications
normalize_model_spec_list <- function(model_specs) {

  if (is.null(model_specs) || length(model_specs) == 0) {
    stop("You must supply at least one model specification.")
  }

  specs <- lapply(seq_along(model_specs), function(i) {
    normalize_model_spec(model_specs[[i]], idx = i)
  })

  spec_names <- vapply(specs, function(x) x$name, character(1))

  if (anyDuplicated(spec_names) > 0) {
    dup <- unique(spec_names[duplicated(spec_names)])
    stop("Model specification names must be unique. Duplicates: ", paste(dup, collapse = ", "))
  }

  names(specs) <- spec_names
  specs
}


# decide whether two XGB model specs can share the same one-time tuning object
same_xgb_tuning_recipe <- function(spec_a, spec_b) {

  if (spec_a$residualizer != "xgb" || spec_b$residualizer != "xgb") {
    return(FALSE)
  }

  if (!identical(spec_a$confounder_order, spec_b$confounder_order)) {
    return(FALSE)
  }

  if (!identical(spec_a$exclude, spec_b$exclude)) {
    return(FALSE)
  }

  # if one spec already brings its own tuning object and the other does not,
  # keep them separate to avoid assuming they are the same
  if (is.null(spec_a$xgb_tuning) != is.null(spec_b$xgb_tuning)) {
    return(FALSE)
  }

  # if both bring a tuning object, it must be identical
  if (!is.null(spec_a$xgb_tuning) && !is.null(spec_b$xgb_tuning)) {
    if (!same_canonical_object(spec_a$xgb_tuning, spec_b$xgb_tuning)) {
      return(FALSE)
    }
  }

  # when tuning is done internally, the requested tuning setup must match
  if (!same_canonical_object(spec_a$xgb_tune_args, spec_b$xgb_tune_args)) {
    return(FALSE)
  }

  if (!identical(spec_a$tune_xgb, spec_b$tune_xgb)) {
    return(FALSE)
  }

  TRUE
}


# assign one XGB tuning-group id to every specification
assign_xgb_tuning_group_ids <- function(model_specs) {

  specs <- model_specs
  next_group <- 1L

  for (i in seq_along(specs)) {

    if (specs[[i]]$residualizer != "xgb") {
      specs[[i]]$xgb_tuning_group_id <- NA_integer_
      next
    }

    if (!is.na(specs[[i]]$xgb_tuning_group_id)) {
      next
    }

    specs[[i]]$xgb_tuning_group_id <- next_group

    if (i < length(specs)) {
      for (j in (i + 1L):length(specs)) {
        if (same_xgb_tuning_recipe(specs[[i]], specs[[j]])) {
          specs[[j]]$xgb_tuning_group_id <- next_group
        }
      }
    }

    next_group <- next_group + 1L
  }

  specs
}


# decide whether two model specs can share the same prepared stage-1 data
same_stage1_recipe <- function(spec_a, spec_b) {

  if (spec_a$residualizer != spec_b$residualizer) {
    return(FALSE)
  }

  # no residualisation means raw renamed data, so all "none" models can share it
  if (spec_a$residualizer == "none") {
    return(TRUE)
  }

  if (!identical(spec_a$confounder_order, spec_b$confounder_order)) {
    return(FALSE)
  }

  if (!identical(spec_a$exclude, spec_b$exclude)) {
    return(FALSE)
  }

  if (!same_canonical_object(spec_a$residualizer_args, spec_b$residualizer_args)) {
    return(FALSE)
  }

  if (spec_a$residualizer == "xgb") {
    return(identical(spec_a$xgb_tuning_group_id, spec_b$xgb_tuning_group_id))
  }

  TRUE
}


# assign one stage-1 group id to every specification
assign_stage1_group_ids <- function(model_specs) {

  specs <- model_specs
  next_group <- 1L

  for (i in seq_along(specs)) {

    if (!is.na(specs[[i]]$stage1_group_id)) {
      next
    }

    specs[[i]]$stage1_group_id <- next_group

    if (i < length(specs)) {
      for (j in (i + 1L):length(specs)) {
        if (same_stage1_recipe(specs[[i]], specs[[j]])) {
          specs[[j]]$stage1_group_id <- next_group
        }
      }
    }

    next_group <- next_group + 1L
  }

  specs
}


# build a readable stage-1 group object list
build_stage1_groups <- function(model_specs) {

  group_ids <- sort(unique(vapply(model_specs, function(x) x$stage1_group_id, integer(1))))

  lapply(group_ids, function(gid) {

    members <- Filter(function(x) identical(x$stage1_group_id, gid), model_specs)

    list(
      stage1_group_id = gid,
      prototype = members[[1]],
      model_names = vapply(members, function(x) x$name, character(1))
    )
  })
}


# split a combined results frame back into the originally requested model-specific frames
split_results_by_model_name <- function(results_df, model_specs) {

  out <- setNames(vector("list", length(model_specs)),
                  vapply(model_specs, function(x) x$name, character(1)))

  for (nm in names(out)) {
    out[[nm]] <- results_df[results_df$model_name == nm, , drop = FALSE]
  }

  out
}
