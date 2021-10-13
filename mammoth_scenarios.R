# Set directories for script and general and K Cuts and Human Density data files, and results
SOURCE_DIR <- getwd()
DATA_DIR <- file.path(SOURCE_DIR, "data")
K_CUTS_DIR <- file.path(SOURCE_DIR, "k_cuts")
RESULTS_DIR <- file.path(SOURCE_DIR, "results")

# Parallel cores available on machine
parallel_cores <- 95

# Install paleopop and dependencies
loadPackage <- function(pkg, min_version = NULL){
  needs_update <- ifelse(!is.null(min_version) && compareVersion(as.character(packageVersion(pkg)), min_version) < 0, TRUE, FALSE)
  if(pkg %in% rownames(installed.packages()) == FALSE || needs_update) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)
}
{loadPackage("doParallel", min_version = "1.0.14")
loadPackage("foreach", min_version = "1.4.4")
loadPackage("gdistance", min_version = "1.2.2")
loadPackage("geosphere", min_version = "1.5.7")
loadPackage("lhs", min_version = "1.0")
loadPackage("metRology", min_version = "0.9.28.1")
loadPackage("R6", min_version = "2.4.0")
loadPackage("raster", min_version = "2.8.19")
loadPackage("sf", min_version = "0.7.7")}

# Paleopop V 1.0.0 should already be installed if mammoth_full has been run
library(paleopop)
stopifnot(packageVersion("paleopop") == "1.0.0")

library(pbapply)

# Read the multi-simulation manager from the full (90000) simulation runs
multi_simulator <- MultiSimulationManager$new()$read_from_rds(path = DATA_DIR)

# If the full 90,000 simulations had been run, and the ABC analysis had been completed
# the next line would be run, which would replace the sample_data in the multi_simulator
# object to the sample_data for the 900 ABC selected runs.
# However, we will leave it as is, and run the counter-factual scenarios on the 246 simulations
# produced in mammoth_full.R
ABC.complete <- FALSE
if (ABC.complete) {
  # Change the sample data frame to the (900) models selected by the ABC validation process
  multi_simulator$sample_data <- read.csv(file = file.path(DATA_DIR, "mammoth_sample_data_selected.csv"))
}

# Human density scenarios
#   Periods: T1 = 21-15kBP (plus burn-in); T2 = 15-11kBP; T3 = 11-5kBP (plus tail)
human_density_scenarios <- c("t123", "tX23", "tXX3", "tXXX")

# Store the original human occupancy mask
original_human_occupancy_mask <- multi_simulator$human_density_model$human_occupancy_mask

# Run the counter-factual scenarios (see manuscript for details)
## About 11 minutes on 95 cores.
{t1 = Sys.time()
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
    ## harvest is required for model to run, but set harvest values to 0
    multi_simulator$model_template$harvest <- TRUE
    multi_simulator$model_template$harvest_max <- 0
    multi_simulator$model_template$harvest_max_n <- 0
    multi_simulator$human_density_model$human_occupancy_mask <- array(0, c(7042, 921))
  }
  # Set the results directory for the scenario
  multi_simulator$results_dir = file.path(RESULTS_DIR, "scenarios", scenario)
  # Run the simulations for the scenario
  sim_log <- multi_simulator$run()
  print(sim_log$summary)
}
Sys.time() - t1}

# Read in abundances for each counter factual scenario and plots trends in abundance
cf_abunds <- lapply(human_density_scenarios, function(scen, ...) {
  # extract abundances
  out_sims <- list.files(path = file.path(RESULTS_DIR, "scenarios", scen), 
                         pattern = "_results.RDS",
                         full.names = TRUE)
  res_abund <- pbsapply(out_sims, function(result,...) {
    # extract abundances from results
    abund <- readRDS(result)$abundance
    abund <- colSums(abund)
    return(abund)
  })
  res_abund
})
names(cf_abunds) <- human_density_scenarios

# which simulations go extinct in original scenario
ext_sim <- apply(cf_abunds$t123, 2, function(x) match(0, x))
sum(+(!is.na(ext_sim)))
ext_simidx <- unname(which(!is.na(ext_sim)))

# CF tXXX
ext_simCF <- apply(cf_abunds$tXXX, 2, function(x) match(0, x))
sum(+(!is.na(ext_simCF)))
ext_simCFidx <- unname(which(!is.na(ext_simCF)))

# which sims are not extinct under tXXX, but is in t123?
(chSim <- sapply(strsplit(names(ext_sim)[which(!ext_simidx %in% ext_simCFidx)], "/"), tail, 1)[1])

# Plot the abundances for t123 and tXXX scenarios

# time for x-axis
ts <- rev(seq(from = 0, by = -25, length = nrow(cf_abunds$t123)))

{
  # plot of abundances for full scenario
  matplot(cf_abunds$t123, type = "l", xlab = "Year (ka BP)", ylab = "total abundance",
        x = ts, col = "#a1a1a165")
  # extinct sims
  matlines(ts, y = cf_abunds$t123[, ext_simidx], lty = 1, col = "black")
  lines(ts, y = floor(apply(cf_abunds$t123, 1, mean)), lty = 2, col = "blue2")
  lines(ts,  y = floor(apply(cf_abunds$t123, 1, median)), lty = 1, col = "blue2")
  # extinct sims for tXXX (no humans)
  matlines(ts, y = cf_abunds$tXXX[, ext_simCFidx], lty = 1, col = "red")
  lines(ts, y = floor(apply(cf_abunds$tXX3, 1, mean)), lty = 2, col = "orangered")
  lines(ts,  y = floor(apply(cf_abunds$tXX3, 1, median)), lty = 1, col = "orangered")
}

# Are abundances statistically different (lower) in the harvest scenario (t123)?
t.test(x = floor(apply(cf_abunds$t123, 1, mean)),
       y = floor(apply(cf_abunds$tXX3, 1, mean)),
       alternative = "less", paired = TRUE)

# Extirpation differences between scenarios

## Extirpation time rasters
cf_ext <- lapply(human_density_scenarios[c(1, 4)], function(scen) {
  message(scen)
  # only look at 1 scenario for example
  out_sims <- list.files(file.path(RESULTS_DIR, "scenarios", scen), "_results.RDS",
                         full.names = TRUE)
  out_sims <- out_sims[grepl(chSim, out_sims)]
  res_ext <- pbsapply(out_sims, function(result) {
   (readRDS(result)$extirpation)*-25
  })
  res_ext
})
names(cf_ext) <- human_density_scenarios[c(1,3)]
region_raster <- rasterFromXYZ(
  cbind(multi_simulator$model_template$coordinates, cf_ext$t123, cf_ext$tXX3),
  crs = 4326, res = 1)
names(region_raster) <- c("t123", "tXXX")

# 0 values mean never occupied after burn-in
## the example files are for a very small niche breadth
## so only cover a small geographic region
## turn 0 to NA
region_raster[region_raster == 0] <- NA
region_raster

# plot of small area
ranges <- c(min(cellStats(region_raster, min, na.rm=TRUE)),
            max(cellStats(region_raster, max, na.rm=TRUE)))
ranges
plot(region_raster, 
     col = hcl.colors(100, "Lajolla"),
     nr = 2, legend = TRUE, zlim = ranges,
     ext = c(2,131,53,84),
     addfun = function() lines(rnaturalearth::ne_coastline(50, "sp")))

## difference?
delta <- region_raster[["tXXX"]] - region_raster[["t123"]]
delta
ranges <- max(abs(c(min(cellStats(delta, min, na.rm=TRUE)),
            max(cellStats(delta, max, na.rm=TRUE)))))
plot(delta, 
     legend = TRUE, zlim = c(-ranges, ranges),
     col = hcl.colors(100, "Zissou"), ext = c(2,131,53,84),
     addfun = function() lines(rnaturalearth::ne_coastline(50, "sp")))
