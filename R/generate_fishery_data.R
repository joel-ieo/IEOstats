#' Generate Synthetic Fishery Sampling Data
#'
#' @description
#' Creates a synthetic fishery survey dataset with a two-stage sampling
#' structure (Fishing Trips and Fishing Operations), suitable for testing
#' variance estimators such as Horvitz-Thompson and rescaled bootstrap.
#'
#' @param seed Integer value used to initialize the random number generator.
#'
#' @param trip_range Integer vector of length 2.
#'   The first element is the minimum number of sampled fishing trips per
#'   stratum and the second element is the maximum.
#'   Example:
#'   \code{c(15, 30)}.
#'
#' @param operation_range Integer vector of length 2.
#'   The first element is the minimum number of fishing operations per trip
#'   and the second element is the maximum.
#'   Example:
#'   \code{c(3, 8)}.
#'
#' @return A `data.table` with one row per sampled fishing operation and
#' the following variables:
#'
#' \describe{
#' \item{Year}{Survey year.}
#' \item{Fleet}{Fishing fleet.}
#' \item{FishingArea}{Fishing area.}
#' \item{Quarter}{Calendar quarter.}
#' \item{FTid}{Fishing-trip identifier.}
#' \item{FOid}{Fishing-operation identifier.}
#' \item{SAtotalWtLive}{Synthetic total live weight.}
#' \item{N_FT}{Number of fishing trips in the population.}
#' \item{n_FT}{Number of sampled fishing trips.}
#' \item{M_FO}{Number of fishing operations in the trip population.}
#' \item{m_FO}{Number of sampled fishing operations.}
#' \item{SAspeCodeFAO}{Synthetic FAO species code.}
#' \item{SScatchFra}{Catch fraction, discarded or landed.}
#' \item{SAid}{Sampling-record identifier.}
#' }
#'
#' @examples
#' fishery_data <- generate_fishery_data()
#'
#' fishery_data <- generate_fishery_data(
#'   seed = 999,
#'   trip_range = c(20, 40),
#'   operation_range = c(4, 10)
#' )
#'
#' @export
generate_fishery_data <- function(
    seed = 123,
    trip_range = c(15, 30),
    operation_range = c(3, 8)
) {

  # -------------------------------------------------------------------
  # ARGUMENT VALIDATION
  # -------------------------------------------------------------------

  if (
    !is.numeric(seed) ||
    length(seed) != 1L ||
    is.na(seed) ||
    !is.finite(seed) ||
    seed < 0 ||
    seed > .Machine$integer.max ||
    seed != floor(seed)
  ) {
    stop(
      paste0(
        "`seed` must be a single integer between 0 and ",
        .Machine$integer.max,
        "."
      ),
      call. = FALSE
    )
  }

  validate_range <- function(x, argument) {
    if (
      !is.numeric(x) ||
      length(x) != 2L ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x != floor(x)) ||
      any(x < 1L) ||
      x[1L] > x[2L]
    ) {
      stop(
        sprintf(
          paste0(
            "`%s` must contain two positive integers ",
            "in increasing order."
          ),
          argument
        ),
        call. = FALSE
      )
    }

    as.integer(x)
  }

  seed <- as.integer(seed)

  trip_range <- validate_range(trip_range, "trip_range")

  operation_range <- validate_range(
    operation_range,
    "operation_range"
  )

  set.seed(seed)

  #--------------------------------------------------
  # STRATA DEFINITION
  #--------------------------------------------------

  strata_def <- data.table::CJ(
    Year = c(2023, 2024),
    Fleet = c("TRAWL", "LONGLINE"),
    FishingArea = c("ATLANTIC", "MEDITERRANEAN"),
    Quarter = 1:4
  )

  #--------------------------------------------------
  # FIRST STAGE: FISHING TRIPS
  #--------------------------------------------------

  trip_list <- vector("list", nrow(strata_def))

  trip_counter <- 1

  for (i in seq_len(nrow(strata_def))) {

    n_FT <- sample(
      seq.int(
        from = trip_range[1L],
        to = trip_range[2L]
      ),
      size = 1L
    )

    N_FT <- sample(
      seq.int(
        from = max(120L, n_FT),
        to = max(250L, n_FT)
      ),
      size = 1L
    )

    trip_list[[i]] <- data.table::data.table(

      FTid = sprintf(
        "FT%04d",
        trip_counter:(trip_counter + n_FT - 1)
      ),

      Year = strata_def$Year[i],
      Fleet = strata_def$Fleet[i],
      FishingArea = strata_def$FishingArea[i],
      Quarter = strata_def$Quarter[i],

      n_FT = n_FT,
      N_FT = N_FT
    )

    trip_counter <- trip_counter + n_FT
  }

  trips <- data.table::rbindlist(trip_list)

  #--------------------------------------------------
  # SECOND STAGE: FISHING OPERATIONS
  #--------------------------------------------------

  operation_list <- vector("list", nrow(trips))

  fo_counter <- 1

  for (i in seq_len(nrow(trips))) {

    m_FO <- sample(
      seq.int(
        from = operation_range[1L],
        to = operation_range[2L]
      ),
      size = 1L
    )

    M_FO <- sample(
      seq.int(
        from = max(20L, m_FO),
        to = max(60L, m_FO)
      ),
      size = 1L
    )

    operation_list[[i]] <- data.table::data.table(

      FOid = sprintf(
        "FO%05d",
        fo_counter:(fo_counter + m_FO - 1)
      ),

      FTid = trips$FTid[i],

      m_FO = m_FO,
      M_FO = M_FO
    )

    fo_counter <- fo_counter + m_FO
  }

  operations <- data.table::rbindlist(operation_list)

  #--------------------------------------------------
  # MERGE STAGES
  #--------------------------------------------------

  fishery_data <- merge(
    operations,
    trips,
    by = "FTid",
    allow.cartesian = TRUE
  )

  #--------------------------------------------------
  # AUXILIARY VARIABLES
  #--------------------------------------------------

  fishery_data[, SAspeCodeFAO := sample(
    c("MEG", "HKE", "COD", "ANK", "PLE"),
    .N,
    replace = TRUE,
    prob = c(0.30, 0.25, 0.20, 0.15, 0.10)
  )]

  fishery_data[, SScatchFra := sample(
    c("Dis", "Lan"),
    .N,
    replace = TRUE,
    prob = c(0.35, 0.65)
  )]

  fishery_data[, SAid := sprintf(
    "SA%06d",
    seq_len(.N)
  )]

  #--------------------------------------------------
  # TARGET VARIABLE
  #--------------------------------------------------

  fishery_data[, SAtotalWtLive := stats::rgamma(
    .N,
    shape = ifelse(SScatchFra == "Dis", 3, 10),
    scale = ifelse(Fleet == "TRAWL", 22, 12)
  )]

  #--------------------------------------------------
  # COLUMN ORDER
  #--------------------------------------------------

  data.table::setcolorder(
    fishery_data,
    c(
      "Year",
      "Fleet",
      "FishingArea",
      "Quarter",
      "FTid",
      "FOid",
      "SAtotalWtLive",
      "N_FT",
      "n_FT",
      "M_FO",
      "m_FO",
      "SAspeCodeFAO",
      "SScatchFra",
      "SAid"
    )
  )

  data.table::setorder(
    fishery_data,
    Year,
    Fleet,
    FishingArea,
    Quarter,
    FTid,
    FOid
  )

  return(fishery_data)
}
