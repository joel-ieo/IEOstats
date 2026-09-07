#' Estimation of Fishery Discard Variance and Totals
#'
#' @description
#' Population totals, empirical variances, and confidence intervals for fishery
#' discard data are estimated using either a multistage rescaled bootstrap 
#' (Preston, 2009) or the classical Horvitz-Thompson analytical method.
#'
#' @param data A data table or data frame containing fishery sampling records.
#' @param ids A formula or character vector defining multi-stage clusters.
#' @param y A formula or character vector specifying the target analysis variable.
#' @param strata A formula or character vector indicating stratification variables.
#' @param n A formula or character vector of sample sizes per stage.
#' @param N A formula or character vector of population sizes per stage.
#' @param B Numeric value defining the number of bootstrap replicates.
#' @param additional_col A formula or character vector of auxiliary variables.
#' @param method Character string defining the estimation method ("ht" or "rs_boot").
#' @param filter Character string with a logical condition to subset data.
#' @param plot Logical flag to render a graphical representation.
#' @param seed Integer to initialize the random number generator.
#' @param alpha Numeric value establishing the significance level.
#' @param parallel Logical flag indicating whether parallel processing is utilized.
#' @param cores Integer specifying the number of computational cores.
#'
#' @returns A data table containing estimated totals, empirical variances, and confidence intervals by strata.
#'
#' @examples
#' \dontrun{
#' strata_estimates <- estimate_fishery_variance(
#'   data              = fishery_data,
#'   ids               = ~ FTid + FOid,
#'   y                 = ~ SAtotalWtLive,
#'   strata            = ~ Year + Fleet + FishingArea + Quarter,
#'   N                 = ~ N_FT + M_FO,
#'   n                 = ~ n_FT + m_FO,
#'   method            = "rs_boot",
#'   filter            = 'SAspeCodeFAO == "MEG" & SScatchFra == "Dis"',
#'   parallel          = FALSE,
#'   B                 = 1000,
#'   additional_col    = ~ SAspeCodeFAO + SScatchFra + SAid,
#'   seed              = 123,
#'   alpha             = 0.05,
#'   plot              = TRUE
#' )
#' }
#'
#' @export 

estimate_discard_variance <- function(
    data,
    ids = NULL,
    y = NULL,
    strata = NULL,
    n = NULL,
    N = NULL,
    B = 1000,
    additional_col = NULL,
    method = c("ht", "rs_boot"),
    filter = NULL,
    plot = TRUE,
    seed = 123,
    alpha = 0.05,
    parallel = FALSE,
    cores = max(1, future::availableCores() - 1)
){
  
  time <- system.time({
    # Parameter validation for bootstrap replicates and safety limits
    if (!is.numeric(B) || length(B) != 1 || B %% 1 != 0 || B < 1) {
      stop("Variable Error: 'B' (number of bootstrap replicates) must be a single positive integer greater than or equal to 1.")
    }
    
    safe_cores_limit <- max(1, future::availableCores() - 1)
    
    if (!is.numeric(cores) || length(cores) != 1 || cores %% 1 != 0 || cores < 1) {
      stop("Variable Error: 'cores' must be a single positive integer greater than or equal to 1.")
    }
    
    if (cores > safe_cores_limit) {
      stop(paste0(
        "Variable Error: Requested ", cores, " cores, but the maximum safe limit ",
        "to prevent system freezing is ", safe_cores_limit,
        " (out of ", future::availableCores, " detected cores)."
      ))
    }
    
    if (!is.null(seed)) {
      if (!is.numeric(seed) || length(seed) != 1 || seed %% 1 != 0) {
        stop("Variable Error: 'seed' must be a single integer or NULL.")
      }
    }
    
    if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
      stop("Variable Error: 'alpha' must be a single numeric value strictly between 0 and 1.")
    }
    
    if (!is.logical(parallel) || length(parallel) != 1 || is.na(parallel)) {
      stop("Variable Error: 'parallel' must be a logical value (TRUE or FALSE).")
    }
    
    if (!is.logical(plot) || length(plot) != 1 || is.na(plot)) {
      stop("Variable Error: 'plot' must be a logical value (TRUE or FALSE).")
    }
    
    if (!is.null(filter)) {
      if (!is.character(filter) || length(filter) != 1 || is.na(filter) || nchar(filter) == 0) {
        stop("Variable Error: 'filter' must be a single non-empty character string or NULL.")
      }
    }
    
    # Data initialization and formatting as data table
    data.table::setDT(data)
    
    # Helper function to parse formulas or character vectors into standard variable names
    parse_argument <- function(arg) {
      arg_name <- deparse(substitute(arg))
      
      if (is.null(arg)) return(character())
      if (inherits(arg, "formula")) return(all.vars(arg))
      if (is.character(arg)) return(arg)
      stop(sprintf("Type Error: '%s' must be either a formula or a character vector.", arg_name))
    }
    
    vars_ids <- parse_argument(ids)
    
    # Create a temporary ID if no cluster IDs are provided (single-stage design)
    if (length(vars_ids) == 0) {
      vars_ids <- ".temp_id"
      data[, .temp_id := .I]
      vars_ids <- ids <- ".temp_id"
      on.exit(if (".temp_id" %in% names(data)) data[, .temp_id := NULL], add = TRUE)
    }
    
    var_y <- parse_argument(y)
    vars_strata <- parse_argument(strata)
    vars_n <- parse_argument(n)
    vars_N <- parse_argument(N)
    vars_additional_col <- parse_argument(additional_col)
    
    # Structural consistency verification
    all_requested_vars <- c(vars_ids, var_y, vars_strata, vars_n, vars_N, vars_additional_col)
    missing_vars <- all_requested_vars[!all_requested_vars %in% names(data)]
    
    if (length(missing_vars) > 0) {
      stop(sprintf(
        "Variable Error: The following variables do not exist in the dataset: %s",
        paste(missing_vars, collapse = ", ")
      ))
    }
    
    if (length(var_y) != 1) {
      stop("Variable Error: 'y' must resolve to a single variable. Check your formula or character input.")
    }
    
    length_ids <- length(vars_ids)
    length_n <- length(vars_n)
    length_N <- length(vars_N)
    
    if (length_n != length_N) {
      stop(sprintf(
        "Consistency Error: 'n' has %d variables but 'N' has %d variables. They must match.",
        length_n, length_N
      ))
    }
    
    if (length_ids == 0 && length_n > 1) {
      stop(sprintf(
        "Structure Error: You provided no 'ids' (Single-Stage), but 'n' and 'N' have %d stages. For multiple stages, you must specify the 'ids' for each level.",
        length_n
      ))
    }
    
    if (length_ids > 0 && length_ids != length_n) {
      stop(sprintf(
        "Consistency Error: You provided %d 'ids' but %d stages for 'n' and 'N'. Every identified stage needs its own n and N.",
        length_ids, length_n
      ))
    }
    
    data.table::setDT(data)
    
    # Apply filtering condition if provided
    data_filter <- if (!is.null(filter)) data[eval(parse(text = filter))] else data.table::copy(data)
    
    # Execute classical Horvitz-Thompson variance if selected
    if("ht" %in% method){
      var_result <- ht_variance(data_filter, ids, y, strata, N, n)
    }
    
    # Execute rescaled bootstrap iterations sequentially or in parallel
    if("rs_boot" %in% method){
      
      if (isTRUE(parallel)) {
        
        old_plan <- future::plan()
        
        future::plan(future::multisession, workers = cores)
        
        on.exit(future::plan(old_plan), add = TRUE)
        
        y_hat_boot <- future.apply::future_replicate(
          n = B,
          expr = {
            rs_boot_iteration(
              data              = data,
              ids               = ids,
              y                 = y,
              strata            = strata,
              N                 = N,
              n                 = n,
              filter            = filter,
              additional_col    = additional_col
            )
          },
          simplify    = FALSE,
          future.seed = seed
        )
      } else {
        set.seed(seed)
        y_hat_boot <- replicate(
          B,
          rs_boot_iteration(data, ids, y, strata, N, n, filter, additional_col),
          simplify = FALSE
        )
      }
      
      res_bootstrapped <- data.table::rbindlist(y_hat_boot, idcol = "iter")
      
      # Summarize bootstrap results and compute empirical confidence intervals
      var_result <- res_bootstrapped[, .(
        estimated_total   = mean(y_hat, na.rm = TRUE),
        variance_total    = var(y_hat, na.rm = TRUE),
        LC                = quantile(y_hat, 1 - alpha / 2, na.rm = TRUE),
        UC                = quantile(y_hat, alpha / 2, na.rm = TRUE)
      ), by                 = vars_strata
      ]
    }
    
    # Compute Horvitz-Thompson expansion totals by strata
    expansion_terms <- paste0("(", vars_N, " / ", vars_n, ")")
    expr_text <- paste(c(var_y, expansion_terms), collapse = " * ")
    expr_parsed <- parse(text = expr_text)
    total_result <- data_filter[, .(ht_total = sum(eval(expr_parsed), na.rm = TRUE)), by = vars_strata]
    
    final_result <- if (length(vars_strata)) merge(total_result, var_result, by = vars_strata) else cbind(total_result, var_result)
    
    # Calculate parametric confidence limits if classical method is chosen
    if("ht" == method){
      final_result[, `:=`(
        LC = ht_total - qt(1 - alpha / 2, nrow(data) - 1) * sqrt(variance_total),
        UC = ht_total + qt(1 - alpha / 2, nrow(data) - 1) * sqrt(variance_total)
      )]
    }
    
    # Graphical representation module using ggplot2
    if (isTRUE(plot)) {
      
      final_result[, `:=`(
        label_strata    = do.call(paste, c(.SD, sep = " || ")),
        strata          = paste0("label_", .I)),
        .SDcols         = vars_strata
      ]
      
      final_result[, strata := factor(strata, levels = paste0("label_", sort(as.numeric(gsub("label_", "", strata)))))]
      
      p <- ggplot2::ggplot(final_result, ggplot2::aes(x = strata, y = ht_total)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#e28743", linewidth = 1) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = LC, ymax = UC), width = 0.2, color = "#4a6572", linewidth = 0.8) +
        ggplot2::geom_point(color = "#f9aa33", size = 4) +
        ggplot2::geom_text(
          ggplot2::aes(label = label_strata),
          angle       = 90,
          position    = ggplot2::position_nudge(x = 0.25),
          hjust       = 0,
          vjust       = 0.5,
          size        = 2.5,
          color       = "gray30"
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.25))) +
        ggplot2::labs(
          title       = "ESTIMATIONS",
          subtitle    = paste0((1 - alpha) * 100, "% confidence intervals"),
          x           = "Terms",
          y           = "Estimation"
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
          panel.grid.minor    = ggplot2::element_blank(),
          panel.grid.major.x  = ggplot2::element_blank(),
          panel.grid.major.y  = ggplot2::element_line(color = "#f0f0f0", linewidth = 0.5),
          axis.line.x         = ggplot2::element_line(color = "#bdc3c7", linewidth = 0.5),
          plot.title          = ggplot2::element_text(face = "bold", color = "#2c3e50", size = 15, hjust = 0),
          plot.subtitle       = ggplot2::element_text(color = "gray40", size = 11, hjust = 0),
          axis.text.x         = ggplot2::element_text(angle = 45, hjust = 1, color = "gray40"),
          axis.text.y         = ggplot2::element_text(color = "gray40"),
          plot.margin         = ggplot2::margin(15, 15, 15, 15)
        )
      print(p)
      
    }
  })
  
  print(time)
  
  return(final_result)
}