# Set directories for source script (and package), input data, niche K cut data, and result files
SOURCE_DIR <- "C:\\Users\\dafcluster\\Desktop\\Mammoth_paleopop"
DATA_DIR <- file.path(SOURCE_DIR, "data")
K_CUTS_DIR <- file.path(SOURCE_DIR, "k_cuts")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")
SAVED_DIR <- file.path(SOURCE_DIR, "saved")

# Parallel cores available on machine
parallel_cores <- 60

# Install paleopop v1.1 and dependencies
loadPackage <- function(pkg, min_version = NULL){
  needs_update <- ifelse(!is.null(min_version) && compareVersion(as.character(packageVersion(pkg)), min_version) < 0, TRUE, FALSE)
  if(pkg %in% rownames(installed.packages()) == FALSE || needs_update) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)
}
loadPackage("abc", min_version = "2.1")
loadPackage("doParallel", min_version = "1.0.14")
loadPackage("foreach", min_version = "1.4.4")
loadPackage("gdistance", min_version = "1.2.2")
loadPackage("geosphere", min_version = "1.5.7")
loadPackage("lhs", min_version = "1.0")
loadPackage("metRology", min_version = "0.9.28.1")
loadPackage("R6", min_version = "2.4.0")
loadPackage("raster", min_version = "2.8.19")
loadPackage("sf", min_version = "0.7.7")
loadPackage("trend", min_version = "1.1.1")
if ("paleopop" %in% rownames(installed.packages()) == FALSE ||
    compareVersion(as.character(packageVersion("paleopop")), "1.1.0") < 0) {
  install.packages(file.path(SOURCE_DIR, "paleopop_1.1.0.tar.gz"))
}
library(paleopop)

# Build the template population model (fixed parameters)
model_template = PaleoPopModel$new(duration = 921,
                                   years_per_step = 25,
                                   random_seed = 1234,
                                   transition_rate = 1.0,
                                   populations = 7042,
                                   coordinates = file.path(DATA_DIR, "mammoth_coordinates.csv"),
                                   density_dependence = "contest",
                                   occupancy_threshold = 2,
                                   dispersal_target_k_threshold = 10,
                                   harvest = TRUE,
                                   harvest_g = 0.4)
model_template$results_selection <- c("abundance")

# Build a correlation model and calculate & attach Cholesky decomposition data
correlation_model <- CorrelationModel$new(coordinates = model_template$coordinates,
                                          amplitude = 0.99, breadth = 8.0)
model_template$compact_decomposition <- correlation_model$get_compact_decomposition()

# Build a niche carrying capacity model for dynamically generating sample carrying capacities
# and initial abundances utilizing niche suitability files (and sample density)
niche_k_model <- NicheCarryingCapacityModel$new(niche_occupancy_mask = file.path(DATA_DIR, "mammoth_occupancy_mask.csv"))
niche_k_model$add_file_template("carrying_capacities",
                                path_template = file.path(K_CUTS_DIR, "NB%s_K_cut_%s.RData"),
                                path_params = c("niche_breadth", "niche_cuts"),
                                file_type = "RDS")

# Build a dispersal model and pre-calculate distance data for dynamically generating sample dispersals
dispersal_model <- DispersalModel$new(coordinates = model_template$coordinates,
                                      dispersal_function_data = file.path(DATA_DIR, "mammoth_dispersal_function_data.csv"))
dispersal_model$set_distance_classes(minimum = 100, maximum = 500, interval = 20)
dispersal_model$calculate_distance_data()

# Build a human density model and pre-calculate sampling parameters for dynamically generating sample human densities
human_density_model <- HumanDensityModel$new(carrying_capacity_mean = file.path(DATA_DIR, "mammoth_human_density_mean.RDS"),
                                             carrying_capacity_sd = file.path(DATA_DIR, "mammoth_human_density_sd.RDS"),
                                             human_occupancy_mask = file.path(DATA_DIR, "mammoth_human_occupancy_mask.RDS"))
human_density_model$calculate_sampling_parameters(mean_upper = 300, max_upper = 500)

# Generate 20,000 model and generative parameter samples for each niche breadth cut group
niche_breadth_cuts = list("40" = 609, "50" = 509, "60" = 401, "70" = 301, "80" = 201, "90" = 101)
sample_data <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = 12))
names(sample_data) <- c("sample", "standard_deviation", "growth_rate_max", "local_threshold", "harvest_max",
                        "harvest_z", "density", "niche_breadth", "niche_cuts", "dispersal_proportion",
                        "dispersal_max_distance", "human_density_sample")
for (nb in names(niche_breadth_cuts)) {
  lhs_generator <- LatinHypercubeSampler$new()
  lhs_generator$set_uniform_parameter("standard_deviation", lower = 0, upper = 0.175)
  lhs_generator$set_uniform_parameter("growth_rate_max", lower = 1.28, upper = 6.84)
  lhs_generator$set_uniform_parameter("local_threshold", lower = 0, upper = 500, decimals = 0)
  lhs_generator$set_uniform_parameter("harvest_max", lower = 0, upper = 0.35)
  lhs_generator$set_uniform_parameter("harvest_z", lower = 1, upper = 2)
  lhs_generator$set_uniform_parameter("density", lower = 625, upper = 10000) # also alias for harvest_max_n
  lhs_generator$set_class_parameter("niche_breadth", as.numeric(nb))
  lhs_generator$set_class_parameter("niche_cuts", 1:niche_breadth_cuts[[nb]])
  lhs_generator$set_uniform_parameter("dispersal_proportion", lower = 0.05, upper = 0.25)
  lhs_generator$set_uniform_parameter("dispersal_max_distance", lower = 100, upper = 500)
  lhs_generator$set_uniform_parameter("human_density_sample", lower = 0, upper = 1)
  lhs_generator$generate_samples(number = 20000, random_seed = 5678)
  sample_data <- rbind(sample_data, data.frame(sample = 1:20000, lhs_generator$sample_data))
}
write.csv(sample_data, file = file.path(SAVED_DIR, "mammoth_sample_data_v1.1.csv"), row.names = FALSE)

# Build, run and save a multi-simulation manager
multi_simulator <- MultiSimulationManager$new(model_template = model_template,
                                              sample_data = sample_data,
                                              parallel_cores = parallel_cores,
                                              niche_k_model = niche_k_model,
                                              dispersal_model = dispersal_model,
                                              human_density_model = human_density_model,
                                              results_dir = file.path(RESULTS_DIR, "full_v1.1"),
                                              results_filename_attributes = c("niche_breadth", "sample"))
t1 = Sys.time()
sim_log <- multi_simulator$run()
Sys.time() - t1
print(sim_log$summary)
multi_simulator$save_to_rds(path = file.path(SAVED_DIR, "MultiSimulationManager_v1.1.RDS"))
