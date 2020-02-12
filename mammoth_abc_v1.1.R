# Set directories for source script (and package), input data, and result files
# SOURCE_DIR <- "C:\\Users\\Sean\\GoogleDrive\\ARC-Grant\\MetaPop_Toolset\\Mammoth"
SOURCE_DIR <- "C:\\Users\\dafcluster\\Desktop\\Mammoth_paleopop"
DATA_DIR <- file.path(SOURCE_DIR, "data")
SAVED_DIR <- file.path(SOURCE_DIR, "saved")

# Load paleopop (assuming the same machine that ran the simulations with v1.1)
library(paleopop)

# Read the metrics manager that generated the summary metrics
metrics_manager <- ResultsSummaryManager$new()$read_from_rds(path = file.path(SAVED_DIR, "ResultsMetricsManager_v1.1.RDS"))

# Use simulation parameters from sample data plus additional niche model parameters
sample_data <- metrics_manager$sample_data
sample_parameters <- sample_data[, c("growth_rate_max", "standard_deviation", "local_threshold", 
                                     "density", "dispersal_proportion", "dispersal_max_distance", 
                                     "harvest_max", "harvest_z", "human_density_sample")]

# Additional parameters derived from original niche model attributes
OMI_lookup <- list()
for (breadth in c(40, 50, 60, 70, 80, 90)) {
  OMI_lookup[[as.character(breadth)]] <- read.csv(file = file.path(DATA_DIR, sprintf("OMI_breadth_%s.csv", breadth)))[, c("OMI", "Tol")]
}
additional_parameters <- data.frame(OMI = array(NA, nrow(sample_data)), Tol = array(NA, nrow(sample_data)))
for (i in 1:nrow(sample_data)) {
  additional_parameters[i, ] <- OMI_lookup[[as.character(sample_data$niche_breadth[i])]][sample_data$niche_cuts[i],]
}
additional_parameters$OMI_Fac <- as.integer(Hmisc::cut2(additional_parameters$OMI, g = 1000, digits = 3))
additional_parameters$Tol_Fac <- as.integer(Hmisc::cut2(additional_parameters$Tol, g = 1000, digits = 3))

# Build an ABC validation model with simulation parameters, metrics, and targets
abc_model <- AbcValidationModel$new(simulation_parameters = cbind(sample_parameters, 
                                                                  additional_parameters[, c("OMI_Fac", "Tol_Fac")]),
                                    simulation_summary_metrics = metrics_manager$summary_metric_data[, -1],
                                    observed_metric_targets = c("abundance_trend" = -1589.12,
                                                                "extinction_time" = -4662.5,
                                                                "fossil_distance" = 0,
                                                                "fossil_sites" = 94))

# Replace NAs in fossil distance metric with maximum value
abc_model$non_finite_replacements <- list(fossil_distance = max)

# Run the ABC validation
abc_model$run_abc(tol = 0.005, # best 600 models
                  method = "neuralnet",
                  numnet = 10,
                  sizenet = 3,
                  lambda = 0.001,
                  maxit = 1000, 
                  trace = FALSE)

# Generate diagnostics (PDF)
abc_model$generate_abc_diagnostics(output_dir = SAVED_DIR)

# Save the selected simulation samples and their weights
selected_samples <- abc_model$selected_simulations
sample_data_selected <- cbind(sample_data[selected_samples$index,], weight = selected_samples$weight)
row.names(sample_data_selected) <- NULL
write.csv(sample_data_selected, file = file.path(SAVED_DIR, "mammoth_sample_data_selected_v1.1.csv"), row.names = FALSE)

# Save ABC validator
abc_model$save_to_rds(path = file.path(SAVED_DIR, "AbcValidationModel_v1.1.RDS"))
