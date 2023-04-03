# Installing packages

library(remotes)
remotes::install_github("EpiModel/EpiModel")
remotes::install_github("EpiModel/EpiModelHIV-p@EDP")
remotes::install_github("EpiModel/ARTnetData")
remotes::install_github("EpiModel/ARTnet")
remotes::install_github("EpiModel/EpiModelHPC")

pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModelHIV-p")
pkgload::load_all("C:\\Users\\clchand\\OneDrive - Emory University\\EpiModel-repos\\EpiModel")

renv::update()
renv::status()
renv::snapshot()
