# Set directories for script and general and K Cuts and Human Density data files, and results
SOURCE_DIR <- "C:\\Users\\dafcluster\\Desktop\\Mammoth_paleopop"
DATA_DIR <- file.path(SOURCE_DIR, "data")
K_CUTS_DIR <- file.path(SOURCE_DIR, "k_cuts")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")

# Parallel cores available on machine
parallel_cores <- 60

# Install paleopop and dependencies
loadPackage <- function(pkg, min_version = NULL){
  needs_update <- ifelse(!is.null(min_version) && compareVersion(as.character(packageVersion(pkg)), min_version) < 0, TRUE, FALSE)
  if(pkg %in% rownames(installed.packages()) == FALSE || needs_update) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)
}
loadPackage("doParallel", min_version = "1.0.14")
loadPackage("foreach", min_version = "1.4.4")
loadPackage("gdistance", min_version = "1.2.2")
loadPackage("geosphere", min_version = "1.5.7")
loadPackage("lhs", min_version = "1.0")
loadPackage("metRology", min_version = "0.9.28.1")
loadPackage("R6", min_version = "2.4.0")
loadPackage("raster", min_version = "2.8.19")
loadPackage("sf", min_version = "0.7.7")
if ("paleopop" %in% rownames(installed.packages()) == FALSE) {
  install.packages(file.path(SOURCE_DIR, "paleopop_1.0.0.tar.gz"))
}
library(paleopop)

# Read the multi-simulation manager from the full (90000) simulation runs
multi_simulator <- MultiSimulationManager$new()$read_from_rds(path = DATA_DIR)

# Change the sample data frame to the (900) models selected by the ABC validation process
multi_simulator$sample_data <- read.csv(file = file.path(DATA_DIR, "mammoth_sample_data_selected.csv"))

# Human density scenarios
#   Periods: T1 = 21-15kBP (plus burn-in); T2 = 15-11kBP; T3 = 11-5kBP (plus tail)
human_density_scenarios <- c("t123", "tX23", "tXX3", "tXXX")

# Store the original human occupancy mask
original_human_occupancy_mask <- multi_simulator$human_density_model$human_occupancy_mask

t1 = Sys.time()
for (scenario in human_density_scenarios) {

  # Alter the human density model for scenario with appropriate human occupancy mask modifications
  if (scenario == "t123") { # Scenario human density present in T1, T2 & T3 (original)
    multi_simulator$model_template$harvest <- TRUE
    multi_simulator$human_density_model$human_occupancy_mask <- original_human_occupancy_mask
  } else if (scenario == "tX23") { # Scenario human density present in T2 & T3 (no T1)
    multi_simulator$model_template$harvest <- TRUE
    mask_modifier <- array(1, c(7042, 921))
    mask_modifier[, 1:320] <- 0
    multi_simulator$human_density_model$human_occupancy_mask <- original_human_occupancy_mask*mask_modifier
  } else if (scenario == "tXX3") { # Scenario human density present in  T3 only
    multi_simulator$model_template$harvest <- TRUE
    mask_modifier <- array(1, c(7042, 921))
    mask_modifier[, 1:480] <- 0
    multi_simulator$human_density_model$human_occupancy_mask <- original_human_occupancy_mask*mask_modifier
  } else if (scenario == "tXXX") { # Scenario no human density present
    multi_simulator$model_template$harvest <- FALSE
    multi_simulator$human_density_model$human_occupancy_mask <- array(0, c(7042, 921)) # redundant
  }

  # Set the results directory for the scenario
  multi_simulator$results_dir = file.path(RESULTS_DIR, "scenarios", scenario)

  # Run the simulations for the scenario
  sim_log <- multi_simulator$run()
  print(sim_log$summary)

}
Sys.time() - t1
