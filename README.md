# Redistribute
An R package to redistribute data.

The data being redistributed (`source`) are those observed in a given frame (across rows
with associated IDs). These data are redistributed to new frames (`target`; different rows with IDs
mapped to those of the source).

Generally, frames represent groupings, where the lowest-level frame contains single observations
(e.g., of an individual, an individual at single time-point, or a time-point from a single
source, etc.), and the highest-level frame is a single observation of the entire population.

For example, the [U.S. Census Bureau](https://www.census.gov) releases data at different geolevels,
such as county, tract, and block group, which are increasingly lower-level (higher-resolution).
These represent the different frames between which observations might be redistributed. This could
be useful if you only have observations in one frame, but you want observations in another.

Frames may also be at roughly the same level (in that they have similar observation-group sizes),
but are arranged differently (e.g., they group individuals along different dimensions). In this
case, it might be best to redistribute data to a lower-level frame, then back up to the other
higher-level frame.

## Installation
Download R from [r-project.org](https://www.r-project.org), then install the package from an R console:

```R
# install.packages("remotes")
remotes::install_github("uva-bi-sdad/redistribute")
```

And load the package:
```R
library(redistribute)
```
