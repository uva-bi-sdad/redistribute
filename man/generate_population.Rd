% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_population.R
\name{generate_population}
\alias{generate_population}
\title{Generate a Population}
\usage{
generate_population(N = 1000, regions = NULL, capacities = NULL,
  region_ids = NULL, attraction_loci = 3, random_regions = 0.1,
  cost_loci = 2, size_loci = 5, similarity_metric = "euclidean",
  n_neighbors = 50, neighbor_range = 0.5, n_races = 6,
  n_building_types = 3, verbose = FALSE)
}
\arguments{
\item{N}{Number of initial individuals to generate. Final number of individuals will be larger.}

\item{regions}{A vector of region IDs, a matrix of coordinates, or an \code{sf} object with
geometries from which coordinates can be derived. If not specified
(and \code{capacities} is not specified), regions similar to housing units (with a mix of
single and multi-family locations) will be generated.}

\item{capacities}{A vector with the maximum number of households for each entry in \code{regions}.}

\item{region_ids}{A vector of unique IDs for each \code{regions}, or a column name in
\code{regions} containing IDs.}

\item{attraction_loci}{Number of locations selected to be centers of attractiveness,
which influence where households are located.}

\item{random_regions}{A number between \code{0} and \code{1}, which determines the proportion of
people who are randomly relocated, in the case that there is more capacity than households.}

\item{cost_loci}{Number of locations selected to be centers of cost, which influences
the initial income associated with households.}

\item{size_loci}{Number of locations selected to be centers of size, which influence
household sizes.}

\item{similarity_metric}{Name of a metric to use to calculate nearness between
neighbors; see \code{\link[lingmatch]{lma_simets}}.}

\item{n_neighbors}{Number of neighbors used to influence each new household's
initial age and race.}

\item{neighbor_range}{Minimum similarity between people to be considered
neighbors, between \code{0} and \code{1} (where \code{0} means unrestricted,
and \code{1} means same region only).}

\item{n_races}{Number of different race groups to sample from.}

\item{n_building_types}{Number of different building types to sample from.}

\item{verbose}{Logical; if \code{TRUE}, will show status messages.}
}
\value{
A list with entries for \code{params} (with initial settings), and two
\code{data.frames}:
\describe{
  \item{\code{households}}{
    \tabular{ll}{
      \code{household} \tab Household ID. \cr
      \code{region} \tab Region ID. \cr
      \code{head_income} \tab Income of the first household member. \cr
      \code{size} \tab Number of individuals in the household. \cr
      \code{building_type} \tab Categorical indicator of building type. \cr
      \code{renting} \tab Binary indicator of renting status. \cr
    }
  }
  \item{\code{individuals}}{
    \tabular{ll}{
      \code{household} \tab Household ID. \cr
      \code{person} \tab Person ID. \cr
      \code{neighbors} \tab Number of neighbors who bore on variables. \cr
      \code{age} \tab Age in years. \cr
      \code{sex} \tab Categorical indicator of sex. \cr
      \code{race} \tab Categorical indicator of race. \cr
      \code{income} \tab Income of individual. \cr
    }
  }
}
}
\description{
Simulate a population of individuals within households, with complex relationships between
demographic and location features.
}
\details{
The population is generated in two steps:

\strong{First}, households are generated and placed within \code{regions}. Placement within regions
is determined by total distance from one or more regions selected as attraction loci. If
coordinates are not provided, these are first randomly generated.

After households are placed, household incomes (of the first member) is generated
based on cost loci, which are then used to generate building types (where types are increasingly
associated with income) and then household size (based on size loci, income, and building type).
Renting status is then generated based on income and building type:
60\% chance if income is under the mean income, and 20\% otherwise, multiplied by .8
if the building type is of a selected renting type, or .3 otherwise.

\strong{Second}, individuals are generated for each household. To generate an individual,
first, neighbors are searched for, based on \code{n_neighbors} and \code{neighbor_range}.
Any neighbors are summarized: average age and income, and tabulated race.

These then affect the first member of the household: age is first drawn from a Beta
distribution (with shapes of \code{1} and \code{2} if renting or \code{1.5} otherwise,
multiplied by \code{80}) and added to \code{18}, then adjusted toward a random value centered
on the average neighbor age (floored Gaussian with a standard deviation of \code{1}), and race
is sampled (that with the highest result of a Binomial draw with \code{n_races} trials
proportion of \code{neighbors * base rate} chance of success for each race group).

Neighbors also affect the income of the second member of the household if the first
member's income is under the neighbor mean income (or under \code{40,000} given no neighbors);
in this case, the second member's income is drawn from a Gaussian distribution centered on
the first member's income, with a standard deviation of \code{10,000}.

The second member's age is based on that of the first member; a floored Gaussian centered on
the first member's age, with a standard deviation of \code{15} if the first member's age
if over \code{40}, or \code{5} otherwise, trimmed to be between \code{18} and \code{90}.

The second member's race has a 70\% chance to be that of the first member,
and a 30\% chance to be selected like the first member's.

Members after the second have no income, have age randomly selected from a
uniform distribution between \code{0} and the first member's age minus \code{15}
(which is then rounded up), and have race determined by either the first or second
member (50\% chance).

Sex has a 50\% chance to be \code{0} or \code{1} for all but the second member; their sex
has a 10\% chance to be the same as the first member's, and 90\% chance to be opposite.
}
\examples{
generate_population(2)
}
