% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redistribute.R
\name{redistribute}
\alias{redistribute}
\title{Redistribute Data}
\usage{
redistribute(source, target = NULL, map = list(), source_id = "GEOID",
  target_id = source_id, weight = NULL, source_variable = NULL,
  source_value = NULL, aggregate = NULL, weight_agg_method = "auto",
  outFile = NULL, overwrite = FALSE, make_intersect_map = FALSE,
  overlaps = "keep", use_all = TRUE, return_geometry = TRUE,
  return_map = FALSE, verbose = FALSE)
}
\arguments{
\item{source}{A matrix-like object you want to distribute from; usually this will be
the real or more complete dataset, and is often at a lower resolution / higher level.}

\item{target}{A matrix-like object you want to distribute to: usually this will be
the dataset you want but isn't available, and is often at a higher resolution / lower level.}

\item{map}{A list with entries named with \code{source} IDs (or aligning with those IDs),
containing vectors of associated \code{target} IDs (or indices of those IDs). Entries
can also be numeric vectors with IDs as names, which will be used to weigh the relationship.
If IDs are related by substrings (the first characters of \code{target} IDs are \code{source} IDs),
then a map can be automatically generated from them. If \code{source} and \code{target}
contain \code{sf} geometries, a map will be made with \code{\link[sf]{st_intersects}}
(\code{st_intersects(source, target)}). If an intersects map is made, and \code{source}
is being aggregated to \code{target}, and map entries contain multiple target IDs,
those entries will be weighted by their proportion of overlap with the source area.}

\item{source_id, target_id}{Name of a column in \code{source} / \code{target},
or a vector containing IDs. For \code{source}, this will default to the first column. For
\code{target}, columns will be searched through for one that appears to relate to the
source IDs, falling back to the first column.}

\item{weight}{Name of a column, or a vector containing weights (or single value to apply to all cases),
which apply to \code{target} when disaggregating, and \code{source} when aggregating.
Defaults to unit weights (all weights are 1).}

\item{source_variable, source_value}{If \code{source} is tall (with variables spread across
rows rather than columns), specifies names of columns in \code{source} containing variable names
and values for conversion.}

\item{aggregate}{Logical; if specified, will determine whether to aggregate or disaggregate
from \code{source} to \code{target}. Otherwise, this will be \code{TRUE} if there are more
\code{source} observations than \code{target} observations.}

\item{weight_agg_method}{Means of aggregating \code{weight}, in the case that target IDs contain duplicates.
Options are \code{"sum"}, \code{"average"}, or \code{"auto"} (default; which will sum if \code{weight}
is integer-like, and average otherwise).}

\item{outFile}{Path to a CSV file in which to save results.}

\item{overwrite}{Logical; if \code{TRUE}, will overwrite an existing \code{outFile}.}

\item{make_intersect_map}{Logical; if \code{TRUE}, will opt to calculate an intersect-based map
rather than an ID-based map, if both seem possible. If specified as \code{FALSE}, will
never calculate an intersect-based map.}

\item{overlaps}{If specified and not \code{TRUE} or \code{"keep"}, will assign \code{target}
entities that are mapped to multiple \code{source} entities to a single source entity. The value
determines how entities with the same weight should be assigned, between \code{"first"} (default),
\code{"last"}, and \code{"random"}.}

\item{use_all}{Logical; if \code{TRUE} (default), will redistribute map weights so they sum to 1.
Otherwise, entities may be partially weighted.}

\item{return_geometry}{Logical; if \code{FALSE}, will not set the returned \code{data.frame}'s
geometry to that of \code{target}, if it exists.}

\item{return_map}{Logical; if \code{TRUE}, will only return the map, without performing the
redistribution. Useful if you want to inspect an automatically created map, or use it in a later call.}

\item{verbose}{Logical; if \code{TRUE}, will show status messages.}
}
\value{
A \code{data.frame} with a row for each \code{target_ids} (identified by the first column,
\code{id}), and a column for each variable from \code{source}.
}
\description{
Distribute data from a source frame to a target frame.
}
\examples{
# minimal example
source <- c(a = 1, b = 2)
target <- 1:5
(redistribute(source, target, verbose = TRUE))

# multi-entity example
source <- data.frame(id = c("a", "b"), cat = c("aaa", "bbb"), num = c(1, 2))
target <- data.frame(
  id = sample(paste0(c("a", "b"), rep(1:5, 2))),
  population = sample.int(1e5, 10)
)
(redistribute(source, target, verbose = TRUE))
}
