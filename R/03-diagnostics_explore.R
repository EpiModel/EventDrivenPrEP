##
## 02. Network Model Diagnostics: Interactive Analysis
##

# scp -r sph:/projects/epimodel/BigNets/data/input input/data/

# Setup ------------------------------------------------------------------------
suppressMessages({
  library("EpiModelHIV")
})

# Interactive Dx Analysis ------------------------------------------------------

# Load the `NETSIZE` value and the formatted `netsize_string`
# NETSIZE <- 1e4 # to override (before sourcing the file)
source("R/utils-netsize.R")

# Main
dx <- readRDS(paste0("data/input/netdx-main-", netsize_string, ".rds"))
print(dx$dx_main, digits = 2)
plot(dx$dx_main)

print(dx$dx_main_static, digits = 2)
plot(dx$dx_main_static)

# Casual
dx <- readRDS(paste0("data/input/netdx-casl-", netsize_string, ".rds"))
print(dx$dx_casl, digits = 2)
plot(dx$dx_casl)

print(dx$dx_casl_static, digits = 2)
plot(dx$dx_casl_static)

# Inst
dx <- readRDS(paste0("data/input/netdx-inst-", netsize_string, ".rds"))
print(dx$dx_inst, digits = 2)
plot(dx$dx_inst)
