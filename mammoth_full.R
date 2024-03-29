# Set directories for script and general and K Cuts and Human Density data files, and results
SOURCE_DIR <- getwd()
DATA_DIR <- file.path(SOURCE_DIR, "data")
K_CUTS_DIR <- file.path(SOURCE_DIR, "k_cuts")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")

# create results dir if needed
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# Parallel cores available on machine
# 128 max
parallel_cores <- 95

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

# The example data only works with Paleopop v 1.0.0
## There have been code-breaking changes introduced in later versions
## including the complete removal of some functions and classes.
## Here we force installation of paleopop 1.0.0
if ("paleopop" %in% rownames(installed.packages())) {
  # if installed verion is > 1.0.0 then reinstall v 1.0.0 from source
  if (as.character(packageVersion("paleopop")) != "1.0.0") {
    warning(sprintf("paleopop v %s is currently installed. Forcing installation of v1.0.0", packageVersion("paleopop")),
            immediate. = TRUE, call. = FALSE)
    install.packages(file.path(SOURCE_DIR, "paleopop_1.0.0.tar.gz"), repos = NULL, type = "SOURCE",
                     lib = .libPaths()[1]) 
  } else {
    warning("paleopop v 1.0.0 is already installed.", immediate. = TRUE, call. = FALSE)
  }
  # if not installed, install v 1.0.0
  } else if (!"paleopop" %in% rownames(installed.packages())) {
  warning("paleopop  is not currently installed. Installing paleopop v1.0.0 from source",
          immediate. = TRUE, call. = FALSE)
  install.packages(file.path(SOURCE_DIR, "paleopop_1.0.0.tar.gz"), repos = NULL, type = "SOURCE",
                   lib = .libPaths()[1])
  }
library(paleopop, verbose = TRUE)

# Build the template population model (fixed parameters)
model_template = PaleoPopModel$new(duration = 921,
                                   years_per_step = 25,
                                   transition_rate = 1.0,
                                   populations = 7042,
                                   coordinates = file.path(DATA_DIR, "mammoth_coordinates.csv"),
                                   dispersal_target_k_threshold = 10,
                                   harvest = TRUE,
                                   harvest_g = 0.4)

# export mammoth abundances, total harvested, and extirpation time
model_template$results_selection <- c("abundance", "harvested", "extirpation")

# Build a correlation model and calculate & attach Cholesky decomposition data
correlation_model <- CorrelationModel$new(coordinates = model_template$coordinates,
                                          amplitude = 0.99, breadth = 8.0)
model_template$compact_decomposition <- correlation_model$get_compact_decomposition()

# Build a niche carrying capacity model for dynamically generating sample carrying capacities
# and initial abundances utilizing niche suitability files (and sample density)
niche_k_model <- NicheCarryingCapacityModel$new()
niche_k_model$add_file_template("carrying_capacities",
                                path_template = file.path(K_CUTS_DIR, "NB%s_K_cut_%s.RData"),
                                path_params = c("niche_breadth", "niche_cuts"),
                                file_type = "RDS")

# Build a dispersal model and pre-calculate distance data for dynamically generating sample dispersals
dispersal_model <- DispersalModel$new(coordinates = model_template$coordinates,
                                      dispersal_function_data = file.path(DATA_DIR, "mammoth_dispersal_function_data.csv"),
                                      decimals = 3)
dispersal_model$set_distance_classes(minimum = 100, maximum = 500, interval = 20)
dispersal_model$calculate_distance_data()

# Build a human density model and pre-calculate sampling parameters for dynamically generating sample human densities
human_density_model <- HumanDensityModel$new(carrying_capacity_mean = file.path(DATA_DIR, "mammoth_human_density_mean.RDS"),
                                             carrying_capacity_sd = file.path(DATA_DIR, "mammoth_human_density_sd.RDS"),
                                             human_occupancy_mask = file.path(DATA_DIR, "mammoth_human_occupancy_mask.RDS"),
                                             decimals = 4)
human_density_model$calculate_sampling_parameters(mean_upper = 300, max_upper = 500)

# Generate 15,000 model and generative parameter samples for each niche breadth cut group
## 2,122 possible niche samples - we provide 10
niche_breadth_samples = list("40" = 609, "50" = 509, "60" = 401, "70" = 301, "80" = 201, "90" = 101)
sample_data <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = 12))
names(sample_data) <- c("sample", "standard_deviation", "growth_rate_max", "local_threshold", "harvest_max",
                        "harvest_z", "density", "niche_breadth", "niche_cuts", "dispersal_proportion",
                        "dispersal_max_distance", "human_density_sample")
for (nb in names(niche_breadth_samples)) {
  lhs_generator <- LatinHypercubeSampler$new()
  lhs_generator$set_uniform_parameter("standard_deviation", lower = 0, upper = 0.175)
  lhs_generator$set_uniform_parameter("growth_rate_max", lower = 1.28, upper = 6.84, decimals = 3)
  lhs_generator$set_uniform_parameter("local_threshold", lower = 0, upper = 500, decimals = 0)
  lhs_generator$set_uniform_parameter("harvest_max", lower = 0, upper = 0.35, decimals = 3)
  lhs_generator$set_uniform_parameter("harvest_z", lower = 1, upper = 2, decimals = 3)
  lhs_generator$set_uniform_parameter("density", lower = 625, upper = 10000) # also alias for harvest_max_n
  lhs_generator$set_class_parameter("niche_breadth", as.numeric(nb))
  lhs_generator$set_class_parameter("niche_cuts", 1:niche_breadth_samples[[nb]])
  lhs_generator$set_uniform_parameter("dispersal_proportion", lower = 0.05, upper = 0.25)
  lhs_generator$set_uniform_parameter("dispersal_max_distance", lower = 100, upper = 500)
  lhs_generator$set_uniform_parameter("human_density_sample", lower = 0, upper = 1)
  lhs_generator$generate_samples(number = 15000)
  sample_data <- rbind(sample_data, data.frame(sample = 1:15000, lhs_generator$sample_data))
}
# Write sample parameters to csv
# write.csv(sample_data, file = file.path(DATA_DIR, "mammoth_sample_data.csv"), row.names = FALSE)

# As we only provide a subset of niche cuts we need to subset our sample_data
# to only those cuts we provide
## will run 246 of a possible 90,000 simulations
sub_sample_data <- with(sample_data, sample_data[niche_breadth == 40 & niche_cuts %in% 1:10, ])

# Build, run and save a multi-simulation manager
{multi_simulator <- MultiSimulationManager$new(model_template = model_template,
                                              # to run ALL niche cuts replace sub_sample_data with sample_data
                                              sample_data = sub_sample_data,
                                              parallel_cores = parallel_cores,
                                              niche_k_model = niche_k_model,
                                              dispersal_model = dispersal_model,
                                              human_density_model = human_density_model,
                                              results_dir = file.path(RESULTS_DIR, "mammoth_example"),
                                              results_filename_attributes = c("niche_breadth", "niche_cuts", "sample"))
t1 = Sys.time()
sim_log <- multi_simulator$run()
Sys.time() - t1} # 3 mins on 95 cores.

# How many models ran successfully?
print(sim_log$summary)

# Save multi_simulator configuration to RDS so we can read in later for counter-factual scenarios
multi_simulator$save_to_rds(path = DATA_DIR)