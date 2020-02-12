# Set directories for source script (and package), input data, and result files
# SOURCE_DIR <- "C:\\Users\\Sean\\GoogleDrive\\ARC-Grant\\MetaPop_Toolset\\Mammoth"
SOURCE_DIR <- "C:\\Users\\dafcluster\\Desktop\\Mammoth_paleopop"
DATA_DIR <- file.path(SOURCE_DIR, "data")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")
SAVED_DIR <- file.path(SOURCE_DIR, "saved")

# Load paleopop (assuming the same machine that ran the simulations with v1.1)
library(paleopop)

# Read the multi-simulation manager from the full (120,000) simulation runs
multi_simulator <- MultiSimulationManager$new()$read_from_rds(path = file.path(SAVED_DIR, "MultiSimulationManager_v1.1.RDS"))

# Create a results model for encapsulating (and dynamically generating additional) simulation results
results_model = PaleoPopResultsModel$new(coordinates = multi_simulator$model_template$coordinates,
                                         burn_in_duration = 80,
                                         trend_interval = (841 - 18550/25):(841 - 12550/25))

# Initialize the results summary manager with the existing multi-simulation manager and results model
metrics_manager <- ResultsSummaryManager$new(simulation_manager = multi_simulator,
                                             results_model = results_model)

# Remove generative models for carrying capacities and human densities (not used here - faster without)
metrics_manager$niche_k_model <- NULL
metrics_manager$human_density_model <- NULL

# Define metrics and metric functions
metrics_manager$summary_metrics <- c("abundance_trend", "extinction_time", "fossil_distance", "fossil_sites")
metrics_manager$summary_functions <- list()

# Function (via direct value) Sen's slope of abundance across the median extinction confidence interval (defined as trend interval above)
metrics_manager$summary_functions$abundance_trend = "all$abundance_trend"

# Function for extinction time (-BP)
metrics_manager$summary_functions$extinction_time = function(results_model) {
  extinction_time <- -25*(841 - results_model$all$extirpation)
  if (is.na(extinction_time)) {
    extinction_time <- 25
  }
  return(extinction_time)
}

# Function for distance to last fossil record
metrics_manager$summary_functions$fossil_distance = function(results_model) {
  young_fossil_location <- matrix(c(-179, 71.5), nrow = 1, ncol = 2)
  if (all(is.numeric(results_model$all$extinction_location))) {
    return(round(geosphere::distGeo(results_model$all$extinction_location, young_fossil_location)/1000))
  } else {
    return(NA)
  }
}

# Function for number of fossil sites where mammoths are simulated
fossil_data <- read.csv(file.path(DATA_DIR, "Woolly_QBonly_Processed.csv"))
id_coordinates <- data.table::as.data.table(cbind(pop_index = 1:7042, metrics_manager$results_model$coordinates))
indexed_fossil_data <- as.data.frame(id_coordinates[fossil_data, on = c(x = "x", y = "y")])
metrics_manager$results_model$attached$fossil_site_indices <- indexed_fossil_data$pop_index
metrics_manager$summary_functions$fossil_sites = function(results_model) {
  population_occupancy <- +(.rowSums(results_model$occupancy, m = nrow(results_model$occupancy), n = ncol(results_model$occupancy)) > 0)
  return(sum(population_occupancy[results_model$attached$fossil_site_indices]))
}

# Generate the metric summary data and save the metrics manager
t1 = Sys.time()
gen_log <- metrics_manager$generate()
Sys.time() - t1
print(gen_log$summary)
metrics_manager$save_to_rds(path = file.path(SAVED_DIR, "ResultsMetricsManager_v1.1.RDS"))
