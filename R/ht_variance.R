#' Calculation of Horvitz-Thompson Variance for Multistage Designs
#'
#' @description
#' Computes analytical population totals and variance estimates using the
#' classical Horvitz-Thompson method across complex multi-stage sampling designs.
#'
#' @param data A data table or data frame containing the survey sampling records.
#'   For example, `data = survey_data`.
#'
#' @param ids A formula or character vector defining multi-stage clusters.
#'   Variables must be ordered from the first sampling stage to the last.
#'   For example, `ids = ~ PSUid + SSUid` or
#'   `ids = c("PSUid", "SSUid")`.
#'
#' @param y A formula or character vector specifying the target analysis variable.
#'   For example, `y = ~ TargetVar` or
#'   `y = "TargetVar"`.
#'
#' @param strata A formula or character vector indicating stratification groups.
#'   Use `NULL` if the sampling design is not stratified.
#'   For example, `strata = ~ Year + Sector + Area + Quarter`,
#'   `strata = c("Year", "Sector", "Area", "Quarter")`, or
#'   `strata = NULL`.
#'
#' @param N A formula or character vector of population sizes per stage.
#'   Variables must be provided in the same stage order as in `ids`.
#'   For example, `N = ~ N_PSU + N_SSU` or
#'   `N = c("N_PSU", "N_SSU")`.
#'
#' @param n A formula or character vector of sample sizes per stage.
#'   Variables must be provided in the same stage order as in `ids`.
#'   For example, `n = ~ n_PSU + n_SSU` or
#'   `n = c("n_PSU", "n_SSU")`.
#'
#' @returns A data table containing analytical variances grouped by strata.
#'
#' @examples
#' \dontrun{
#' Generate a reproducible synthetic two-stage fishery sample
#' fishery_data <- generate_fishery_data(
#'   seed = 123,
#'   trip_range = c(15, 20),
#'   operation_range = c(3, 5)
#' )
#'
#' ht_variance(
#'   data   = fishery_data,
#'   ids    = ~ FTid + FOid,
#'   y      = ~ SAtotalWtLive,
#'   strata = ~ Year + Fleet + FishingArea + Quarter,
#'   N      = ~ N_FT + M_FO,
#'   n      = ~ n_FT + m_FO
#' )
#' }
#'


ht_variance <- function(data, ids = NULL, y = NULL, strata = NULL, N = NULL, n = NULL) {

  # Convert input data to a data table format
  data.table::setDT(data)

  # Helper function to convert formulas or character vectors into standard variable names
  parse_argument <- function(arg) {
    if (is.null(arg)) return(character())
    if (inherits(arg, "formula")) return(all.vars(arg))
    if (is.character(arg)) return(arg)
    stop("Type Error: Argument must be either a formula or a character vector.")
  }

  vars_ids <- parse_argument(ids)

  # Create a temporary ID if no cluster IDs are provided (single-stage design)
  if (length(vars_ids) == 0) {
    vars_ids <- ".temp_id"
    data[, .temp_id := .I]
    vars_ids <- ids <- ".temp_id"
    on.exit(if (".temp_id" %in% names(data)) data[, .temp_id := NULL], add = TRUE)
  }

  var_y       <- parse_argument(y)
  vars_strata <- parse_argument(strata)
  vars_n      <- parse_argument(n)
  vars_N      <- parse_argument(N)
  length_ids <- length(vars_ids)
  length_n   <- length(vars_n)
  length_N   <- length(vars_N)

  # Check that sample sizes (n) and population sizes (N) match in number of stages
  if (length_n != length_N) {
    stop(sprintf("Consistency Error: 'n' has %d variables but 'N' has %d variables.", length_n, length_N))
  }
  if (length_ids == 0 && length_n > 1) {
    stop(sprintf("Structure Error: No 'ids' provided, but 'n' and 'N' have %d stages.", length_n))
  }
  if (length_ids > 0 && length_ids != length_n) {
    stop(sprintf("Consistency Error: Provided %d 'ids' but %d stages for 'n' and 'N'.", length_ids, length_n))
  }

  cols_necessaries <- unique(c(vars_strata, var_y, vars_ids, vars_n, vars_N))
  missing_vars     <- cols_necessaries[!cols_necessaries %in% names(data)]

  # Check that all required variables are present in the dataset
  if (length(missing_vars) > 0) {
    stop(sprintf("Variable Error: Missing variables: %s", paste(missing_vars, collapse = ", ")))
  }

  if (length(var_y) != 1) {
    stop("Variable Error: 'y' must resolve to a single variable. Check your formula or character input.")
  }

  stages <- vector("list", length(vars_ids))

  # Iterate backwards across stages to compute variance hierarchically
  for (r in rev(seq_along(vars_ids))) {
    current_level <- vars_ids[r]
    ids_groups <- c(
      vars_strata,
      head(vars_ids, r - 1),
      head(vars_n, r - 1),
      head(vars_N, r - 1)
    )

    if (r == length(vars_ids)) {
      # Final stage calculations (sample variance and stage totals)
      stages[[r]] <- data[, {
        n <- get(vars_n[r])[1]
        N <- get(vars_N[r])[1]
        val_y <- get(var_y)
        sum_y <- sum(val_y)
        s2 <- if (n > 1) sum((val_y - (sum_y / n))^2) / (n - 1) else 0
        .(
          n = n,
          sum_y = sum_y,
          s2 = s2,
          N = N
        )
      }, by = ids_groups]

      stages[[r]][, total_y := sum_y * N / n]
      stages[[r]][, f := n / N]
      stages[[r]][, variance_total := (N^2) * (1 - f) * (s2 / n)]

    } else {
      # Intermediate stages calculations (inheriting variance from lower stages)
      stages[[r]] <- stages[[r+1]][, {
        n <- get(vars_n[r])[1]
        N <- get(vars_N[r])[1]
        val_y <- get("total_y")
        sum_y <- sum(val_y)
        s2 <- if (n > 1) sum((val_y - (sum_y / n))^2) / (n - 1) else 0
        inherited_var <- sum(variance_total)
        .(
          n = n,
          sum_y = sum_y,
          s2 = s2,
          N = N,
          inherited_var = inherited_var
        )
      }, by = ids_groups]

      stages[[r]][, total_y := sum_y * N / n]
      stages[[r]][, f := n / N]
      stages[[r]][, stage_var := (N^2) * (1 - f) * (s2 / n)]
      stages[[r]][, variance_total := stage_var + (N / n) * inherited_var]
    }
  }

  cols_to_keep <- c(vars_strata, "variance_total")
  cols_to_drop <- setdiff(names(stages[[1]]), cols_to_keep)

  # Clean up and retain only the final variance results grouped by strata
  if (length(cols_to_drop) > 0) {
    stages[[1]][, (cols_to_drop) := NULL]
  }

  return(stages[[1]])
}
