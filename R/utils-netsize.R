##
##  Shared variables used all over the project
##
##  The variables should be changed here only to avoid discrepencies in the
##  different scripts
##

# Size of network to use
if (! exists("NETSIZE")) {
  NETSIZE <- 10 * 1e3
}

netsize_string <- format(NETSIZE, scientific = FALSE)
