# This code installs the packages only available on GitHub (not on CRAN)
renv::install(c(
  "EpiModel/EpiModelHIV-p@EDP",
  "EpiModel/EpiModelHPC"
))

# Updates EpiModelHIV after pushing code to Github
renv::update("EpiModelHIV")

# Force `renv` to discover the following packages
if (FALSE) {
  library("rmarkdown")
  library("pkgload")
  library("sessioninfo")
}
