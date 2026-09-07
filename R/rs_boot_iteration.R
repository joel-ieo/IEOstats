#' Single Iteration of the Multistage Rescaled Bootstrap Algorithm
#'
#' @description
#' Executes a single bootstrap iteration using the multistage rescaled 
#' bootstrap method (Preston, 2009). Sampling weights and scale factors are adjusted 
#' across all sample stages to calculate randomized totals by strata.
#'
#' @param data A data table or data frame with the sampling records.
#' @param ids A formula or character vector defining the cluster IDs for each stage.
#' @param y A formula or character vector specifying the target variable to estimate.
#' @param strata A formula or character vector indicating the stratification variables.
#' @param N A formula or character vector for population sizes at each stage.
#' @param n A formula or character vector for sample sizes at each stage.
#' @param filter A character string with a condition to subset the data.
#' @param additional_col A formula or character vector of extra columns to keep.
#'
#' @returns A data table with the randomized bootstrap estimates (y_hat) grouped by strata.
#' 
#' @examples
#' \dontrun{
#' rs_boot_iteration(
#'   data   = fishery_data,
#'   ids    = ~ FTid + FOid,
#'   y      = ~ SAtotalWtLive,
#'   strata = ~ Year + Fleet + FishingArea + Quarter,
#'   N      = ~ N_FT + M_FO,
#'   n      = ~ n_FT + m_FO
#' )
#' }
#' 
#' @export

rs_boot_iteration <- function(data, ids = NULL, y = NULL, strata = NULL, N = NULL, n = NULL, filter = NULL, additional_col = NULL) {
  
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
  vars_add    <- parse_argument(additional_col)
  
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
  
  cols_necessaries <- unique(c(vars_strata, var_y, vars_ids, vars_n, vars_N, vars_add))
  missing_vars     <- cols_necessaries[!cols_necessaries %in% names(data)]
  
  # Check that all requested variables exist in the dataset
  if (length(missing_vars) > 0) {
    stop(sprintf("Variable Error: Missing variables: %s", paste(missing_vars, collapse = ", ")))
  }
  
  if (length(var_y) != 1) {
    stop("Variable Error: 'y' must resolve to a single variable. Check your formula or character input.")
  }
  
  # Set up the hierarchical structure using the version 2 function with strata support
  dt_processed <- mark_multistage(data[, ..cols_necessaries], vars_ids, strata)
  
  vars_delta <- paste0("delta_", vars_ids)
  missing_deltas <- vars_delta[!vars_delta %in% names(dt_processed)]
  
  if (length(missing_deltas) > 0) {
    stop(sprintf("Matrix Error: Expected delta columns missing: %s", paste(missing_deltas, collapse = ", ")))
  }
  
  # Initialize baseline tracking variables for the bootstrap weights
  dt_processed[, weight_total := 1.0]
  
  dt_processed[, cum_prod_lambda_prev := 1.0]
  dt_processed[, cum_prod_root_prev   := 1.0]
  dt_processed[, cum_prod_weight_prev := 1.0]
  
  # Loop through each sampling stage to calculate adjustment factors
  for (r in seq_len(length_n)) {
    v_N     <- vars_N[r]
    v_n     <- vars_n[r]
    v_delta <- vars_delta[r]
    
    dt_processed[, `:=`(
      N_r     = as.numeric(get(v_N)),
      n_r     = as.numeric(get(v_n)),
      delta_r = as.numeric(get(v_delta))
    )]
    
    dt_processed[, `:=`(
      n_ast_r = pmax(1, n_r %/% 2),
      f_r     = n_r / N_r,
      w_r     = N_r / n_r
    )]
    
    dt_processed[, lambda_r := data.table::fcase(
      n_r == n_ast_r, 0, 
      default = sqrt((n_ast_r * cum_prod_lambda_prev * (1 - f_r)) / (n_r - n_ast_r))
    )]
    
    dt_processed[, root_term_r := sqrt(n_r / n_ast_r) * delta_r]
    
    dt_processed[, img_formula_r := - lambda_r * cum_prod_root_prev + 
                   lambda_r * cum_prod_root_prev * (n_r / n_ast_r) * delta_r]
    
    if (r == 1) {
      dt_processed[, cum_sum_img := img_formula_r]
    } else {
      dt_processed[, cum_sum_img := cum_sum_img + img_formula_r]
    }
    
    dt_processed[, weight_term_r := w_r * (1 + cum_sum_img)]
    
    dt_processed[, weight_ast_r := cum_prod_weight_prev * weight_term_r]
    
    dt_processed[, weight_total := weight_total * weight_ast_r]
    
    dt_processed[, cum_prod_lambda_prev := cum_prod_lambda_prev * f_r]
    dt_processed[, cum_prod_root_prev   := cum_prod_root_prev * root_term_r]
    dt_processed[, cum_prod_weight_prev := cum_prod_weight_prev * (w_r / weight_term_r)]
  }
  
  # Filter the dataset if a condition is provided
  if (!is.null(filter)) {
    dt_processed <- dt_processed[eval(parse(text = filter))]
  }
  
  # Sum up the final values using the adjusted weights for each strata group
  result <- dt_processed[,
                         .(y_hat = sum(weight_total * as.numeric(get(var_y)), na.rm = TRUE)), 
                         by = vars_strata
  ]
  
  return(result)
}