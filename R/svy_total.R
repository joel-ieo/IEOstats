#' Estimation of Population Totals and Variances
#'
#' @description
#' Estimates population totals, variances, and confidence intervals from
#' complex survey samples. Estimation can be performed using either the
#' classical Horvitz-Thompson analytical estimator or a multistage rescaled
#' bootstrap approach.
#'
#' The function supports stratified and multistage survey sampling designs. Sampling
#' units, analysis variables, strata, sample sizes, and population sizes can be
#' specified using formulas or character vectors.
#'
#' @param data A data frame or data table containing the survey observations
#'   and sampling-design variables.
#'   For example, `data = fishery_data`.
#'
#' @param ids A one-sided formula or character vector identifying the sampling
#'   units at each stage of the sampling design. Variables must be ordered
#'   from the first sampling stage to the last.
#'   For example, `ids = ~ FTid + FOid` or
#'   `ids = c("FTid", "FOid")`.
#'
#' @param y A one-sided formula or character vector identifying the analysis
#'   variable for which the population total is estimated.
#'   For example, `y = ~ SAtotalWtLive` or
#'   `y = "SAtotalWtLive"`.
#'
#' @param strata A one-sided formula or character vector identifying the
#'   stratification variables. Use `NULL` if the sampling design is not
#'   stratified.
#'   For example, `strata = ~ Year + Fleet + FishingArea + Quarter`,
#'   `strata = c("Year", "Fleet", "FishingArea", "Quarter")`, or
#'   `strata = NULL`.
#'
#' @param n A one-sided formula or character vector identifying the observed
#'   sample-size variables at each sampling stage. Variables must be provided
#'   in the same stage order as in `ids`.
#'   For example, `n = ~ n_FT + m_FO` or
#'   `n = c("n_FT", "m_FO")`.
#'
#' @param N A one-sided formula or character vector identifying the
#'   population-size variables at each sampling stage. Variables must be
#'   provided in the same stage order as in `ids`.
#'   For example, `N = ~ N_FT + M_FO` or
#'   `N = c("N_FT", "M_FO")`.
#'
#' @param B A positive integer specifying the number of bootstrap replicates.
#'   This argument is used when `method = "rs_boot"`. Larger values generally
#'   provide more stable bootstrap variance estimates but require more
#'   computation time.
#'   For example, `B = 1000`.
#'
#' @param additional_col A one-sided formula or character vector identifying
#'   additional variables that must be retained during estimation. Use `NULL`
#'   if no additional variables are required.
#'   For example, `additional_col = ~ SAspeCodeFAO + SScatchFra`,
#'   `additional_col = c("SAspeCodeFAO", "SScatchFra")`, or
#'   `additional_col = NULL`.
#'
#' @param method A character string specifying the variance-estimation method.
#'   Available options are `"ht"` for the Horvitz-Thompson analytical method
#'   and `"rs_boot"` for the multistage rescaled bootstrap.
#'   For example, `method = "ht"` or `method = "rs_boot"`.
#'
#' @param filter An optional character string containing a logical expression
#'   used to select observations before estimation. Use `NULL` to retain all
#'   observations.
#'   For example, `filter = "SScatchFra == 'Lan'"` or `filter = NULL`.
#'
#' @param plot Logical value indicating whether a graphical representation of
#'   the estimates should be produced.
#'   For example, `plot = TRUE` or `plot = FALSE`.
#'
#' @param seed A single integer used to initialize the random-number generator
#'   for bootstrap estimation. It is used when `method = "rs_boot"` to obtain
#'   reproducible bootstrap results.
#'   For example, `seed = 123`.
#'
#' @param alpha A numeric value greater than zero and less than one specifying
#'   the significance level used to construct confidence intervals. A value
#'   of `0.05` corresponds to 95 percent confidence intervals.
#'   For example, `alpha = 0.05`.
#'
#' @param parallel Logical value indicating whether bootstrap replicates
#'   should be processed using parallel computation. This argument is
#'   relevant when `method = "rs_boot"`.
#'   For example, `parallel = TRUE` or `parallel = FALSE`.
#'
#' @param cores A positive integer specifying the number of processing cores
#'   used when `parallel = TRUE`. It should not exceed the number of cores
#'   available on the system.
#'   For example, `cores = 2`.
#'
#' @return A data table containing population-total estimates, estimated
#'   variances, and confidence-interval limits for each stratum or analysis
#'   domain.
#'
#' @details
#' When `method = "ht"`, totals and variances are estimated analytically using
#' the Horvitz-Thompson approach.
#'
#' When `method = "rs_boot"`, variance is estimated from repeated rescaled
#' bootstrap samples. In this case, `B` controls the number of bootstrap
#' replicates and `seed` can be used to obtain reproducible results.
#'
#' Variables supplied through `ids`, `n`, and `N` must follow the same
#' sampling-stage order. For example, the first variable in each argument must
#' describe the first sampling stage, the second variable must describe the
#' second stage, and so forth.
#'
#' @references
#' Cochran, W. G. (1977).
#' \emph{Sampling Techniques} (3rd ed.).
#' John Wiley & Sons.
#'
#' Preston, J. (2009).
#' Rescaled bootstrap for stratified multistage sampling.
#' \emph{Survey Methodology}, 35(2), 227--234.
#' \url{https://www150.statcan.gc.ca/n1/pub/12-001-x/2009002/article/11044-eng.pdf}
#'
#'
#' @examples
#' # Generate a reproducible synthetic two-stage fishery sample
#' fishery_data <- generate_fishery_data(
#'   seed = 123,
#'   trip_range = c(15, 20),
#'   operation_range = c(3, 5)
#' )
#'
#' # Horvitz-Thompson analytical estimation
#' estimates_ht <- svy_total(
#'   data = fishery_data,
#'   ids = ~ FTid + FOid,
#'   y = ~ SAtotalWtLive,
#'   strata = ~ Year + Fleet + FishingArea + Quarter,
#'   n = ~ n_FT + m_FO,
#'   N = ~ N_FT + M_FO,
#'   method = "ht",
#'   plot = FALSE,
#'   alpha = 0.05,
#' )
#'
#' estimates_ht
#'
#' \dontrun{
#' # Multistage rescaled bootstrap estimation
#' estimates_boot <- svy_total(
#'   data = fishery_data,
#'   ids = ~ FTid + FOid,
#'   y = ~ SAtotalWtLive,
#'   strata = ~ Year + Fleet + FishingArea + Quarter,
#'   n = ~ n_FT + m_FO,
#'   N = ~ N_FT + M_FO,
#'   B = 1000,
#'   method = "rs_boot",
#'   plot = FALSE,
#'   seed = 123,
#'   alpha = 0.05,
#'   parallel = TRUE,
#'   cores = 1
#' )
#'
#' estimates_boot
#'
#' # Estimate totals only for landed catches
#' estimates_landed <- svy_total(
#'   data = fishery_data,
#'   ids = ~ FTid + FOid,
#'   y = ~ SAtotalWtLive,
#'   strata = ~ Year + Fleet + FishingArea + Quarter,
#'   n = ~ n_FT + m_FO,
#'   N = ~ N_FT + M_FO,
#'   B = 1000,
#'   additional_col = ~ SAspeCodeFAO + SScatchFra,
#'   method = "rs_boot",
#'   filter = "SScatchFra == 'Lan'",
#'   plot = FALSE,
#'   seed = 123,
#'   alpha = 0.05,
#'   parallel = FALSE
#' )
#'
#' estimates_landed
#' }
#'
#' @export

svy_total <- function(
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

    if (length(vars_strata) == 0) {

      data[, .global_strata := "stratum"]

      vars_strata <- ".global_strata"
      strata <- ~.global_strata

      on.exit(
        {
          if (".global_strata" %in% names(data)) {
            data[, .global_strata := NULL]
          }
        },
        add = TRUE
      )
    }

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

    # Estimate classical Horvitz-Thompson variance if selected
    if("ht" %in% method){
      var_result <- ht_variance(data_filter, ids, y, strata, N, n)
    }

    # Run rescaled bootstrap iterations sequentially or in parallel
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
        LC                = quantile(y_hat, alpha / 2, na.rm = TRUE),
        UC                = quantile(y_hat, 1 - alpha / 2, na.rm = TRUE)
      ), by                 = vars_strata
      ]
    }

    # Compute Horvitz-Thompson expansion estimates by stratum
    expansion_terms <- paste0("(", vars_N, " / ", vars_n, ")")
    expr_text <- paste(c(var_y, expansion_terms), collapse = " * ")
    expr_parsed <- parse(text = expr_text)
    total_result <- data_filter[, .(ht_total = sum(eval(expr_parsed), na.rm = TRUE)), by = vars_strata]

    final_result <- if (length(vars_strata)) merge(total_result, var_result, by = vars_strata) else cbind(total_result, var_result)

    # Calculate parametric confidence limits if classical method is selected
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
