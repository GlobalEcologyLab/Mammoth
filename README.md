# Process-explicit models reveal pathway to extinction for woolly mammoth using pattern-oriented validation

This respository contains the `R` code and some example data (for 10 unique niche samples) to run the Mammoth model `paleopop` simulations.

In order to run the example simulations, `paleopop v 1.0.0` must be installed. `paleopop v 1.0.0` is packaged in this respository and the code checks to make sure it is installed. If a newer version is installed, it is reverted to v 1.0.0 due to code-breaking changes introduced in later versions.

The code has a number of package dependencies listed below:

- "doParallel", min version = "1.0.14"
- "foreach", min version = "1.4.4"
- "gdistance", min version = "1.2.2"
- "geosphere", min version = "1.5.7"
- "lhs", min version = "1.0"
- "metRology", min version = "0.9.28.1"
- "R6", min version = "2.4.0"
- "raster", min version = "2.8.19"
- "sf", min version = "0.7.7"

The script checks the above dependencies are advises if a newer version of a package is required. If using `sf > 1.0.0` please set `sf_use_s2` to `FALSE`.

The scripts are designed to run in parallel across *n* sessions - please set *n* accordingly. For a guide, using 95 cores on a workstation uses ~ 90GB RAM to run the example simulations. Example simulations run in approximately 3 minutes.

There are two scripts contained in the repository. They need to be run in the following order:

- `mammoth_full.R`
- `mammoth_scenarios.R`

The first script build, runs, and saves the output from a multi-simulation manager (see `paleopop` documentation for details). The second script, imports the multi-simulation manager and then runs three counter-factual scenarios each with differing rates of human harvest.

The second set of counter-factual simulations run in ~11 minutes across 95 cores.

Generally, the second script would only be run following model selection with, for example, Approximate Bayesian Computation [ABC]. This would ensure that the counter-factual scenarios are only run on the models that best fit observations. See the main manuscript for details on our ABC model selection process.

All outputs are saved in the `results/` directory.

The second script also produces plots of temporal abundance for the example simulations, and creates maps of extirpation time and the difference in extirpation time for a scenario run with and without human harvesting.

All niche samples are available on reasonable request from Damien A Fordham (damien.fordham@adelaide.edu.au).

Additionally, scripts found in `extras/` directory can be used to recreate the sampling effort map in the supporting information of the main manuscript.

Directory structure:

```
Mammoth/
  - Data/ # data necessary to build the models
  - extras/ # additional data to create maps of fossil sampling effort
  - k_cuts/ # niche samples required for the models
  - results/ # output directory for simulation results
```

This code is released under the following licence:

GNU GENERAL PUBLIC LICENSE
   Version 3, 29 June 2007

Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
Everyone is permitted to copy and distribute verbatim copies
of this license document, but changing it is not allowed.
