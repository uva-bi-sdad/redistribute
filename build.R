library(redistribute)

# unload to rebuild
dyn.unload(system.file("libs/x64/redistribute.dll", package = "redistribute"))

# builds documentation and site
styler::style_pkg(filetype = c("R", "Rmd"))
spelling::spell_check_package()
devtools::document()
pkgdown::build_site(lazy = TRUE)
covr::report(covr::package_coverage(quiet = FALSE), "docs/coverage.html")

# runs checks
devtools::check_win_devel()
devtools::check_rhub(interactive = FALSE)
rhub::check(platforms = "macos-highsierra-release-cran")

# releases
devtools::release()
