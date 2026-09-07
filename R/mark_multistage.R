#' Mark Multistage Hierarchical Units and Generate Perturbation Flags
#'
#'@description
#' Prepares the dataset for the rescaled bootstrap method by validating cluster 
#' variables and generating random hierarchical perturbation indicators (delta flags) 
#' across multi-stage sampling levels.
#'
#' @param data A data table or data frame containing sampling records.
#' @param ids A formula or character vector defining the multi-stage cluster variables.
#' @param strata A formula or character vector indicating the stratification variables.
#'
#' @returns A modified data table containing original data and delta flag columns by stratum.
#' 
#' @examples
#' \dontrun{
#' mark_multistage(
#'   data   = fishery_data,
#'   ids    = ~ FTid + FOid,
#'   strata = ~ Year + Fleet + FishingArea + Quarter
#' )
#' }
#' 
#' @export

mark_multistage <- function(data, ids = NULL, strata = NULL) {
  
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
  
  vars_strata <- parse_argument(strata)
  
  cols_necessaries <- unique(c(vars_strata, vars_ids))
  
  missing_vars <- cols_necessaries[!cols_necessaries %in% names(data)]
  
  # Check that all requested variables exist in the dataset
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "Variable Error: The following variables do not exist in the dataset: %s",
      paste(missing_vars, collapse = ", ")
    ))
  }
  
  # Initialize delta flag columns with zeros for each stage
  for (id in vars_ids) data.table::set(data, j = paste0("delta_", id), value = 0L)
  
  # Iteratively select primary and nested sampling units respecting strata groups
  for (i in seq_along(vars_ids)) {
    current_level <- vars_ids[i]
    current_flag  <- paste0("delta_", current_level)
    current_grouping <- unique(c(vars_strata, vars_ids[seq_len(i - 1)]))
    by_cols <- if (length(current_grouping) > 0) current_grouping else NULL
    
    if (i == 1) {
      # First stage: randomly select about half of the units within each stratum group
      data[, (current_flag) := {
        units <- unique(get(current_level))
        n_sel <- max(1L, length(units) %/% 2L)
        selected <- units[sample.int(length(units), n_sel, replace = FALSE)]
        
        .(as.integer(get(current_level) %in% selected))
      }, by = by_cols]
      
    } else {
      # Subsequent stages: randomly select nested units within selected parent units inside strata
      parent_level <- vars_ids[i - 1]
      parent_flag  <- paste0("delta_", parent_level)
      
      data[get(parent_flag) == 1, (current_flag) := {
        units <- unique(get(current_level))
        n_sel <- max(1L, length(units) %/% 2L)
        selected <- units[sample.int(length(units), n_sel, replace = FALSE)]
        
        .(as.integer(get(current_level) %in% selected))
      }, by = by_cols]
    }
  }
  return(data)
}