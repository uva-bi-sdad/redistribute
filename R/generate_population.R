#' Generate a Population
#'
#' Simulate a population of individuals within households, with complex relationships between
#' demographic and location features.
#'
#' The population is generated in two steps:
#'
#' \strong{First}, households are generated and placed within \code{regions}. Placement within regions
#' is determined by total distance from one or more regions selected as attraction loci. If
#' coordinates are not provided, these are first randomly generated.
#'
#' After households are placed, household incomes (of the first member) is generated
#' based on cost loci, which are then used to generate building types (where types are increasingly
#' associated with income) and then household size (based on size loci, income, and building type).
#' Renting status is then generated based on income and building type:
#' 60\% chance if income is under the mean income, and 20\% otherwise, multiplied by .8
#' if the building type is of a selected renting type, or .3 otherwise.
#'
#' \strong{Second}, individuals are generated for each household. To generate an individual,
#' first, neighbors are searched for, based on \code{n_neighbors} and \code{neighbor_range}.
#' Any neighbors are summarized: average age and income, and tabulated race.
#'
#' These then affect the first member of the household: age is first drawn from a Beta
#' distribution (with shapes of \code{1} and \code{2} if renting or \code{1.5} otherwise,
#' multiplied by \code{80}) and added to \code{18}, then adjusted toward a random value centered
#' on the average neighbor age (floored Gaussian with a standard deviation of \code{1}), and race
#' is sampled (that with the highest result of a Binomial draw with \code{n_races} trials
#' proportion of \code{neighbors * base rate} chance of success for each race group).
#'
#' Neighbors also affect the income of the second member of the household if the first
#' member's income is under the neighbor mean income (or under \code{40,000} given no neighbors);
#' in this case, the second member's income is drawn from a Gaussian distribution centered on
#' the first member's income, with a standard deviation of \code{10,000}.
#'
#' The second member's age is based on that of the first member; a floored Gaussian centered on
#' the first member's age, with a standard deviation of \code{15} if the first member's age
#' if over \code{40}, or \code{5} otherwise, trimmed to be between \code{18} and \code{90}.
#'
#' The second member's race has a 70\% chance to be that of the first member,
#' and a 30\% chance to be selected like the first member's.
#'
#' Members after the second have no income, have age randomly selected from a
#' uniform distribution between \code{0} and the first member's age minus \code{15}
#' (which is then rounded up), and have race determined by either the first or second
#' member (50\% chance).
#'
#' Sex has a 50\% chance to be \code{0} or \code{1} for all but the second member; their sex
#' has a 10\% chance to be the same as the first member's, and 90\% chance to be opposite.
#'
#' @param N Number of initial individuals to generate. Final number of individuals will be larger.
#' @param regions A vector of region IDs, a matrix of coordinates, or an \code{sf} object with
#' geometries from which coordinates can be derived. If not specified
#' (and \code{capacities} is not specified), regions similar to housing units (with a mix of
#' single and multi-family locations) will be generated.
#' @param capacities A vector with the maximum number of households for each entry in \code{regions}.
#' @param attraction_loci Number of locations selected to be centers of attractiveness,
#' which influence where households are located.
#' @param random_regions A number between \code{0} and \code{1}, which determines the proportion of
#' people who are randomly relocated, in the case that there is more capacity than households.
#' @param cost_loci Number of locations selected to be centers of cost, which influences
#' the initial income associated with households.
#' @param size_loci Number of locations selected to be centers of size, which influence
#' household sizes.
#' @param similarity_metric Name of a metric to use to calculate nearness between
#' neighbors; see \code{\link[lingmatch]{lma_simets}}.
#' @param n_neighbors Number of neighbors used to influence each new household's
#' initial age and race.
#' @param neighbor_range Minimum similarity between people to be considered
#' neighbors, between \code{0} and \code{1} (where \code{0} means unrestricted,
#' and \code{1} means same region only).
#' @param n_races Number of different race groups to sample from.
#' @param n_building_types Number of different building types to sample from.
#' @param verbose Logical; if \code{TRUE}, will show status messages.
#' @examples
#' generate_population(2)
#' @returns A list with entries for \code{params} (with initial settings), and two
#' \code{data.frames}:
#' \describe{
#'   \item{\code{households}}{
#'     \tabular{ll}{
#'       \code{household} \tab Household ID. \cr
#'       \code{region} \tab Region ID. \cr
#'       \code{head_income} \tab Income of the first household member. \cr
#'       \code{size} \tab Number of individuals in the household. \cr
#'       \code{building_type} \tab Categorical indicator of building type. \cr
#'       \code{renting} \tab Binary indicator of renting status. \cr
#'     }
#'   }
#'   \item{\code{individuals}}{
#'     \tabular{ll}{
#'       \code{household} \tab Household ID. \cr
#'       \code{person} \tab Person ID. \cr
#'       \code{neighbors} \tab Number of neighbors who bore on variables. \cr
#'       \code{age} \tab Age in years. \cr
#'       \code{sex} \tab Categorical indicator of sex. \cr
#'       \code{race} \tab Categorical indicator of race. \cr
#'       \code{income} \tab Income of individual. \cr
#'     }
#'   }
#' }
#' @export

generate_population <- function(N = 1000, regions = NULL, capacities = NULL, attraction_loci = 2,
                                random_regions = .1, cost_loci = 2, size_loci = 2,
                                similarity_metric = "euclidean", n_neighbors = 50,
                                neighbor_range = .5, n_races = 6, n_building_types = 3, verbose = FALSE) {
  gen_regions <- missing(regions)
  gen_capacities <- missing(capacities)
  if (gen_regions) {
    if (gen_capacities) {
      if (verbose) cli_alert_info("regions: sequence along {.arg N}, with housing unit capacities")
      regions <- seq_len(N)
      capacities <- rep(1, N)
      n_multi <- round(N * .1)
      capacities[sample.int(N, n_multi)] <- sample(5:300, n_multi, TRUE)
    } else {
      if (verbose) {
        cli_alert_info(
          "regions: sequence along specified {.arg capacities} with specified capacities"
        )
      }
      regions <- seq_along(capacities)
    }
  }
  if (length(dim(regions)) == 2) {
    if (length(capacities) == 1 && is.character(capacities) && capacities %in% colnames(regions)) {
      if (verbose) {
        cli_alert_info(
          "regions: specified coordinates, with capacities from {.field {capacities}} column"
        )
      }
      su <- colnames(regions) != capacities
      capacities <- regions[[capacities]]
      regions <- regions[, su, drop = FALSE]
    } else if (verbose) {
      cli_alert_info(paste(
        "regions: specified coordinates, with",
        if (gen_capacities) "{.field N / nrow(regions)}" else "specified", "capacities"
      ))
    }
    rids <- if (is.null(rownames(regions))) seq_len(nrow(regions)) else rownames(regions)
  } else {
    rids <- if (length(regions) == 1 && is.numeric(regions)) {
      if (verbose) {
        cli_alert_info(paste(
          "regions: sequence of {.field N}, with",
          if (gen_capacities) "{.field N / nrow(regions)}" else "specified", "capacities"
        ))
      }
      as.character(seq_len(regions))
    } else {
      if (verbose && !gen_regions) {
        cli_alert_info(paste(
          "regions:", if (gen_capacities) {
            "specified region, with {.field N / nrow(regions)}"
          } else {
            "specified regions and"
          }, "capacities"
        ))
      }
      regions
    }
    regions <- matrix(
      ceiling(rnorm(length(rids) * 2, 1e5, 1e4)),
      ncol = 2, dimnames = list(rids, c("X", "Y"))
    )
  }
  nr <- length(rids)
  capacities <- rep_len(if (!is.numeric(capacities)) N / nr else capacities, nr)

  if (verbose) cli_alert_info("preparing {.field {N}} households")
  household <- seq_len(N)

  # select attraction loci
  if (!is.numeric(attraction_loci) || attraction_loci < 1) {
    attraction_loci <- 1
    if (verbose) cli_alert_warning("setting {.arg attraction_loci} to {.field 1}")
  } else if (attraction_loci > nr) {
    attraction_loci <- max(1, nr - 1)
    if (verbose) {
      cli_alert_warning(
        "too many attraction loci; setting {.arg attraction_loci} to {.field {attraction_loci}}"
      )
    }
  } else if (verbose) {
    cli_alert_info("attraction loci: {.field {attraction_loci}}")
  }
  selected_attraction_loci <- sample.int(nr, attraction_loci)
  if (verbose) cli_alert_info("calculating region similarities")
  if (nr == 1) {
    space <- Matrix(1, dimnames = list(rids, rids), sparse = TRUE)
  } else {
    if (inherits(regions, c("sfc", "sf"))) regions <- st_coordinates(st_centroid(regions))
    space <- lma_simets(regions, metric = similarity_metric, symmetrical = TRUE)
  }
  if (verbose) cli_alert_info("rescaling region similarities")
  min_space <- min(space@x)
  space@x <- (space@x - min_space) / (max(space@x) - min_space)
  space_cols <- integer(ncol(space))

  # select initial locations
  space_cols[selected_attraction_loci] <- -1L
  ro <- order(space %*% space_cols)
  space_cols[selected_attraction_loci] <- 0L
  capacities <- ceiling(capacities)
  if (sum(capacities) < N) {
    if (verbose) cli_alert_info("{.arg N} is larger than total capacities; expanding to accomodate")
    capacities <- ceiling(capacities / sum(capacities) * N)
  }
  rids_expanded <- rep(ro, capacities[ro])
  region <- rids_expanded[household]
  if (anyNA(region)) cli_abort("NAs in regions")

  if (N < length(rids_expanded) && is.numeric(random_regions)) {
    # redistribute a portion
    moves <- min(length(rids_expanded) - N, N * random_regions)
    if (moves > 0) {
      if (verbose) cli_alert_info("randomly placing {.field {moves}} households in remaining regions")
      region[sample(household, moves)] <- sample(rids_expanded[-household], moves)
    }
  }

  # select head-of-household income
  if (!is.numeric(cost_loci) || cost_loci < 1) {
    cost_loci <- 1
    if (verbose) cli_alert_warning("setting {.arg cost_loci} to {.field 1}")
  } else if (cost_loci > nr) {
    cost_loci <- max(1, nr - 1)
    if (verbose) cli_alert_warning("too many cost loci; setting {.arg cost_loci} to {.field {cost_loci}}")
  } else if (verbose) {
    cli_alert_info("cost loci: {.field {cost_loci}}")
  }
  selected_cost_loci <- sample.int(nr, cost_loci)
  space_cols[selected_cost_loci] <- 1L
  region_cost <- rnorm(nr, (as.numeric(space %*% space_cols) / nr) * 3e5, 3e4)
  space_cols[selected_cost_loci] <- 0L
  head_income <- abs(as.integer(
    rnorm(N, region_cost[region], 2e3) + (rbeta(N, .1, 1) - .6) * region_cost[region]
  ))

  # selecting building type
  if (verbose) cli_alert_info("drawing building type and rental types")
  building_rent_boost <- if (n_building_types < 2) 1 else sample.int(n_building_types, ceiling(n_building_types * .3))
  building_type <- if (n_building_types < 2) rep(1, nr) else sample.int(n_building_types, length(head_income), TRUE)

  # select household sizes
  if (!is.numeric(size_loci) || size_loci < 1) {
    size_loci <- 1
    if (verbose) cli_alert_warning("setting {.arg size_loci} to {.field 1}")
  } else if (size_loci > nr) {
    size_loci <- max(1, nr - 1)
    if (verbose) cli_alert_warning("too many size loci; setting {.arg size_loci} to {.field {size_loci}}")
  } else if (verbose) {
    cli_alert_info("size loci: {.field {size_loci}}")
  }
  selected_size_loci <- sample.int(nr, size_loci)
  space_cols[selected_size_loci] <- 1L
  size_prob <- ((space %*% space_cols) / nr)[region, ] +
    rbeta(N, head_income / max(head_income) * building_type / max(building_type), 1)
  size <- rpois(N, 1 + size_prob * (max(size_prob) - size_prob) / 1) + 1L

  # selecting renting status
  if (verbose) cli_alert_info("drawing renting status")
  renting <- rbinom(N, 1, ((head_income < mean(head_income)) * .2 + .4) *
    ((building_type %in% building_rent_boost) * .3 + .5))

  # select race base rates
  if (verbose) cli_alert_info("drawing race baserates")
  race_rates <- rbeta(n_races, 1, 5)
  race_rates <- race_rates / (race_rates + max(race_rates))

  # making space sparse
  if (verbose) cli_alert_info("defining region neighbors")
  space@x[space@x < neighbor_range] <- 0
  space <- drop0(space)

  # preparing individual data
  individuals <- data.frame()
  if (verbose) {
    cli_progress_step(
      "generating individuals...",
      msg_done = "generated {.field {nrow(individuals)}} individuals"
    )
  }
  individuals <- generate_individuals(
    region, head_income, size, renting, space, n_neighbors, race_rates
  )

  list(
    params = list(
      neighbors = n_neighbors, range = neighbor_range, races_rates = race_rates,
      n_building_types = n_building_types, building_type_rent = building_rent_boost
    ),
    households = data.frame(
      household,
      region = rids[region], head_income, size, building_type, renting
    ),
    individuals = as.data.frame(individuals)
  )
}
