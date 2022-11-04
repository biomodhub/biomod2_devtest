terraVersion <- TRUE
library(tictoc)
requireNamespace("cli")
requireNamespace("tinytest")
requireNamespace("dplyr") # sorry je n'ai pas pu m'empecher ;-)
requireNamespace("magrittr")

# to run on installed package
# library(biomod2) 
# to run on development package
# (to be adjusted depending on where the dev version is)
devtools::load_all("../biomod2/") 

tic("Test_BIOMOD_FormatingData.R")
source("R/Test_BIOMOD_FormatingData.R")
toc()

tic("Test_bm_PseudoAbsences.R")
source("R/Test_bm_PseudoAbsences.R")
toc()

tic("Test_BIOMOD_Modeling.R")
source("R/Test_BIOMOD_Modeling.R")
toc()

tic("Test_BIOMOD_EnsembleModeling.R")
source("R/Test_BIOMOD_EnsembleModeling.R")
toc()

tic("Test_bm_BinaryTransformation.R")
source("R/Test_bm_binaryTransformation.R")
toc()

tic("Test_bm_plotResponseCurves.R")
source("R/Test_bm_plotResponseCurves.R")
toc()

tic("Test_BIOMOD_Projection.R")
source("R/Test_BIOMOD_Projection.R")
toc()

tic("Test_BIOMOD_EnsembleForecasting.R")
source("R/Test_BIOMOD_EnsembleForecasting.R")
toc()

tic("Test_get_predictions.R")
source("R/Test_get_predictions.R")
toc()


tic("Test_BIOMOD_RangeSize.R")
source("R/Test_BIOMOD_RangeSize.R")
toc()

tic("Test_BIOMOD_PresenceOnly.R")
source("R/Test_BIOMOD_PresenceOnly.R")
toc()

tic("Test_BIOMOD_Parallel.R")
source("R/Test_BIOMOD_Parallel.R")
toc()

cli::cli_h1("Error Summary")

dfres <- data.frame(
  Step = c(
    "BIOMOD_Formating",
    "bm_PseudoAbsences",
    "BIOMOD_Modeling",
    "BIOMOD_EnsembleModeling",
    "bm_BinaryTransformation",
    "bm_PlotResponseCurves", 
    "BIOMOD_Projection",
    "BIOMOD_EnsembleForecasting", 
    "get_predictions(bm.proj.out)",
    "BIOMOD_RangeSize",
    "BIOMOD_PresenceOnly",
    "BIOMOD_Parallel"
  ),
  Error = c(
    Error_Formating, 
    Error_PseudoAbsences,
    Error_Modeling, 
    Error_EnsembleModeling,
    Error_BinaryTransformation, 
    Error_PlotResponseCurves,
    Error_Projection, 
    Error_EnsembleForecasting,
    Error_get_pred, 
    Error_RangeSize,
    Error_PresenceOnly,
    Error_Parallel
  )
)
print(dfres)

