# Set directories for source script (and package), input data, niche K cut data, and result files
SOURCE_DIR <- "C:\\Users\\dafcluster\\Desktop\\Mammoth_paleopop"
DATA_DIR <- file.path(SOURCE_DIR, "data")
K_CUTS_DIR <- file.path(SOURCE_DIR, "k_cuts")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")
SAVED_DIR <- file.path(SOURCE_DIR, "saved")

# Load paleopop (assuming the same machine that ran the full simulations)
library(paleopop)

# Read the multi-simulation manager from the full (120,000) simulation runs
multi_simulator <- MultiSimulationManager$new()$read_from_rds(path = file.path(SAVED_DIR, "MultiSimulationManager_v1.1.RDS"))

# Change the sample data frame to the (600) models selected by the ABC validation process
multi_simulator$sample_data <- read.csv(file = file.path(SAVED_DIR, "mammoth_sample_data_selected_v1.1.csv"))

# First run additional scenario simulations

# Run a scenario with no human impact/harvesting
multi_simulator_no_humans <- multi_simulator$clone()
multi_simulator_no_humans$model_template$harvest <- FALSE
multi_simulator_no_humans$results_dir = file.path(RESULTS_DIR, "scenarios_v1.1", "no_humans")
t1 = Sys.time()
sim_log <- multi_simulator_no_humans$run()
Sys.time() - t1
print(sim_log$summary)

# Run a scenario with no niche climate change (maximized carrying capacity)
multi_simulator_no_climate <- multi_simulator$clone()
multi_simulator_no_climate$results_dir = file.path(RESULTS_DIR, "scenarios_v1.1", "no_climate")
# Modify the carrying capacity generation via inheriting the generative model and overwriting the file read method
NicheKModel2 <- R6::R6Class("NicheKModel2",
                            inherit = NicheCarryingCapacityModel,
                            public = list(
                              read_file = function(param) {
                                array(apply(super$read_file(param), 1, max), c(7042, 921))
                              }
                            ))
niche_k_model <- multi_simulator_no_climate$niche_k_model
multi_simulator_no_climate$niche_k_model <- NicheKModel2$new(niche_occupancy_mask = niche_k_model$niche_occupancy_mask,
                                                             file_templates = niche_k_model$file_templates)
t1 = Sys.time()
sim_log <- multi_simulator_no_climate$run()
Sys.time() - t1
print(sim_log$summary)

# Save the scenario simulation managers
multi_simulator_no_humans$save_to_rds(path = file.path(SAVED_DIR, "MultiSimulationManager_no_humans_v1.1.RDS"))
multi_simulator_no_climate$save_to_rds(path = file.path(SAVED_DIR, "MultiSimulationManager_no_climate_v1.1.RDS"))

# Now gather results data to compare the original and scenario simulations

# Create a results model for encapsulating (and dynamically generating additional) simulation results
results_model = PaleoPopResultsModel$new(coordinates = multi_simulator$model_template$coordinates,
                                         burn_in_duration = 80)

# Initialize a results summary manager with the original multi-simulation manager and results model
results_manager <- ResultsSummaryManager$new(simulation_manager = multi_simulator,
                                             results_model = results_model)

# Define summary result matrices and functions
results_manager$summary_matrices <- c("abundance", "extirpation")
results_manager$summary_functions <- list(abundance = "all$abundance",
                                          extirpation = "extirpation")

# Build the results managers for the scenarios
results_manager_no_humans <- ResultsSummaryManager$new(simulation_manager = multi_simulator_no_humans,
                                                       results_model = results_model,
                                                       summary_matrices = results_manager$summary_matrices,
                                                       summary_functions = results_manager$summary_functions)
results_manager_no_climate <- ResultsSummaryManager$new(simulation_manager = multi_simulator_no_climate,
                                                       results_model = results_model,
                                                       summary_matrices = results_manager$summary_matrices,
                                                       summary_functions = results_manager$summary_functions)

# Generate and save the results summaries (matrices)
t1 = Sys.time()
gen_log_orig <- results_manager$generate()
print(gen_log_orig$summary)
gen_log_no_humans <- results_manager_no_humans$generate()
print(gen_log_no_humans$summary)
gen_log_no_climate <- results_manager_no_climate$generate()
print(gen_log_no_climate$summary)
Sys.time() - t1

# Save the results summaries
results_manager$save_to_rds(path = file.path(SAVED_DIR, "ResultsSummaryManager_original_v1.1.RDS"))
results_manager_no_humans$save_to_rds(path = file.path(SAVED_DIR, "ResultsSummaryManager_no_humans_v1.1.RDS"))
results_manager_no_climate$save_to_rds(path = file.path(SAVED_DIR, "ResultsSummaryManager_no_climate_v1.1.RDS"))

# Calculate the weighted averages of the abundance and extirpation (matrices)
na_replacements <- list(abundance = 0, extirpation = 842)
results_manager$calculate_summary_weighted_averages(na_replacements = na_replacements)
results_manager_no_humans$calculate_summary_weighted_averages(na_replacements = na_replacements)
results_manager_no_climate$calculate_summary_weighted_averages(na_replacements = na_replacements)

# Abundance plots (all and weighted average)
managers <- list(results_manager, results_manager_no_humans, results_manager_no_climate)
titles <- c("Original", "No Humans", "Optimal Climate")
y_max <- c(3000000, 7000000, 15000000)
time_seq <- -25*(841 - 1:841)
for (i in 1:3) {
  plot(x = time_seq, y = array(0, 841), ylim = c(0, y_max[i]), col = "white", 
       xlab = "Year -BP", ylab = "Abundance", main = titles[i])
  matplot(x = matrix(time_seq, nrow = 600, ncol = 841, byrow = TRUE),
          y = managers[[i]]$summary_matrix_list$abundance, 
          type = "l", col = rgb(0.2, 0.2, 0.2, alpha = 0.05), add = TRUE, lwd = 0.25)
  lines(x = time_seq, y = managers[[i]]$summary_matrix_weighted_averages$abundance)
}

# Extirpation map 21,000 years BP to 25 (persistent to present day)
extirpation_original <- results_manager$summary_matrix_weighted_averages$extirpation
extirpation_no_humans <- results_manager_no_humans$summary_matrix_weighted_averages$extirpation
extirpation_no_climate <- results_manager_no_climate$summary_matrix_weighted_averages$extirpation
extirpations <- cbind("Original" = -25*(841 - extirpation_original),
                      "No Humans" = -25*(841 - extirpation_no_humans),
                      "Optimal Climate" = -25*(841 - extirpation_no_climate))
sp::spplot(raster::rasterFromXYZ(cbind(results_manager$results_model$coordinates, 
                                       extirpations)))

# Weighted abundance grid/map animation
normalized_weights <- results_manager$sample_data$weight/sum(results_manager$sample_data$weight)
weighted_abundances <- array(0, c(7042, 841))
for (i in 1:nrow(results_manager$sample_data)) {
  results <- readRDS(file = file.path(results_manager$results_dir,
                                      paste0(results_manager$get_results_filename(i), ".RDS")))
  results_clone <- results_manager$results_model$new_clone(results = results)
weighted_abundances <- weighted_abundances + normalized_weights[i]*results_clone$abundance
}
animation_plots = function() {
  for (i in 1:841) { # extra cell (under *) to maintain consistent colour scales
    raster::plot(raster::rasterFromXYZ(cbind(rbind(results_manager$results_model$coordinates, c(-169.5,40.5)), 
                                             c(weighted_abundances[, i], 835))),
                 axes = FALSE)#, main = "Weighted Mean Mammoth Abundance for Selected Models")
    text(x = -170.5, y = 39.6, 
         labels = paste("*", format(25*(841 - i),big.mark = ","), "years BP"), 
         adj = c(0, NA))
  }
}
animation::ani.options(interval = 0.1, ani.width = 1280, ani.height = 600,
                       ani.dev = function(...){png(res = 96, bg = "white", type = "cairo-png", ...)})
animation::saveVideo(animation_plots(),
                     clean = TRUE,
                     video.name = "weighted_abundance_animation.mkv",
                     other.opts = "-s:v 1280x600 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p")

# Save the original results manager again
results_manager$attached$weighted_abundances <- weighted_abundances
results_manager$save_to_rds(path = file.path(SAVED_DIR, "ResultsSummaryManager_original_v1.1.RDS"))

