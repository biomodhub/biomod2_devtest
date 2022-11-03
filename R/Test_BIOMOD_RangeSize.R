cli::cli_h1("BIOMOD_RangeSize")

Error_RangeSize <- 0

# Preparation -------------------------------------------------------------

library(terra)
library(raster)
if(terraVersion){
  data("DataSpecies")
  data("bioclim_current")
  myExpl <- terra::rast(bioclim_current)
  myExplFuture <- terra::rast(bioclim_future)
} else {
  myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
  DataSpecies <- read.csv(myFile, row.names = 1)
  myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
  myExpl <- rast(raster::stack(system.file(myFiles, package = 'biomod2')))
  myFiles <- paste0('external/bioclim/future/bio', c(3, 4, 7, 11, 12), '.grd')
  myExplFuture <- rast(raster::stack(system.file(myFiles, package = 'biomod2')))
}
myRespName <- 'GuloGulo'

## myResp ------------------------------------------------------------------
myResp <- as.numeric(DataSpecies[, myRespName])
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )

## myExpl ----------------------------------

myExpl.df <- extract(myExpl, y = myResp.SpatVector, ID = FALSE)
myExplFuture.df <- extract(myExplFuture, y = myResp.SpatVector, ID = FALSE)
myExpl.matrix <- as.matrix(myExpl.df)
myExpl.raster <- stack(myExpl)


## Load Models -------------------------------------------------------------

file.out <- paste0(myRespName, "/", myRespName, ".NoCat_NoEval_Presence-Absence.models.out")
myBiomodModelOut_noCat <- 
  try(suppressWarnings(
    get(load(file.out))
  ), silent = TRUE)

file.out <- paste0(myRespName, "/", myRespName, ".NoCat_NoEval_Presence-Absence.ensemble.models.out")
myBiomodEnsembleOut_noCat <- 
  try(suppressWarnings(
    get(load(file.out))
  ), silent = TRUE)


# Get projections ----------------------------------------------------
invisible(
  capture.output(suppressWarnings(suppressMessages({
    
    
    ## Simple model - current --------------------------------------------------
    myBiomodProj <-
      BIOMOD_Projection(
        bm.mod = myBiomodModelOut_noCat,
        proj.name = 'Current',
        new.env = myExpl,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    myBiomodProj.df <-
      BIOMOD_Projection(
        bm.mod = myBiomodModelOut_noCat,
        proj.name = 'Current.df',
        new.env = myExpl.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    threshold_table <-
      get_evaluations(myBiomodModelOut_noCat, as.data.frame = TRUE) %>% 
      dplyr::filter(Eval.metric == "TSS") %>% 
      dplyr::mutate(full.name = glue::glue("GuloGulo_{Dataset}_{Run}_{Algo}")) %>% 
      dplyr::select(full.name, Cutoff)
    rownames(threshold_table) <- threshold_table$full.name
    threshold_table <- threshold_table[,-1, drop = FALSE]
    
    ProjCurrent.nonbinary.df <- 
      get_predictions(myBiomodProj.df, 
                      as.data.frame = TRUE)
    ProjCurrent.nonbinary <- 
      get_predictions(myBiomodProj)
    
    ProjCurrent.nonbinary.raster <- stack(ProjCurrent.nonbinary)
    
    ProjCurrent.df <- 
      get_predictions(myBiomodProj.df, 
                      as.data.frame = TRUE) %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_table[get_projected_models(myBiomodProj),"Cutoff"])
    
    ProjCurrent.1.df <-
      get_predictions(myBiomodProj.df, 
                      as.data.frame = TRUE,
                      model = "GLM",
                      run.eval = "RUN1") %>% 
      bm_BinaryTransformation(
        threshold = threshold_table["GuloGulo_AllData_RUN1_GLM", "Cutoff"])
    
    ProjCurrent.all.SpatRaster <-
      get_predictions(myBiomodProj)  %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_table[get_projected_models(myBiomodProj),"Cutoff"])
    
    ProjCurrent.1.SpatRaster <-
      get_predictions(myBiomodProj,
                      model = "GLM",
                      run.eval = "RUN1") %>% 
      bm_BinaryTransformation(
        threshold = threshold_table["GuloGulo_AllData_RUN1_GLM", "Cutoff"])
    
    
    ## Simple model - future --------------------------------------------------
    myBiomodFuture <-
      BIOMOD_Projection(
        bm.mod = myBiomodModelOut_noCat,
        proj.name = 'Future',
        new.env = myExplFuture,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    myBiomodFuture.df <-
      BIOMOD_Projection(
        bm.mod = myBiomodModelOut_noCat,
        proj.name = 'Future.df',
        new.env = myExplFuture.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    
    ProjFuture.nonbinary.df <- 
      get_predictions(myBiomodFuture.df, 
                      as.data.frame = TRUE)
    ProjFuture.nonbinary <- 
      get_predictions(myBiomodProj)
    
    ProjFuture.nonbinary.raster <- stack(ProjFuture.nonbinary)
    
    ProjFuture.df <-
      get_predictions(myBiomodFuture.df, as.data.frame = TRUE) %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_table[get_projected_models(myBiomodProj),"Cutoff"])
    
    ProjFuture.1.df <- 
      get_predictions(myBiomodFuture.df, as.data.frame = TRUE,
                      model = "GLM",
                      run.eval = "RUN1") %>% 
      bm_BinaryTransformation(
        threshold = threshold_table["GuloGulo_AllData_RUN1_GLM", "Cutoff"])
    
    ProjFuture.2.df <- cbind(ProjFuture.1.df,ProjFuture.1.df)
    colnames(ProjFuture.2.df) <- c("ssp1","ssp2")
    
    ProjFuture.all.SpatRaster <-
      get_predictions(myBiomodFuture) %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_table[get_projected_models(myBiomodProj),"Cutoff"])
    
    ProjFuture.1.SpatRaster <- 
      get_predictions(myBiomodFuture,
                      model = "GLM",
                      run.eval = "RUN1") %>% 
      bm_BinaryTransformation(
        threshold = threshold_table["GuloGulo_AllData_RUN1_GLM", "Cutoff"])
    
    ProjFuture.2.SpatRaster <- c(ProjFuture.1.SpatRaster,ProjFuture.1.SpatRaster)
    # names(ProjFuture.2.SpatRaster) <- c("ssp1","ssp2")
    
    ## Ensemble model - current --------------------------------------------------
    
    myBiomodEnsembleForecast <-
      BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEnsembleOut_noCat,
        proj.name = 'Current',
        new.env = myExpl,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    myBiomodEnsembleForecast.df <-
      BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEnsembleOut_noCat,
        proj.name = 'Current.df',
        new.env = myExpl.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    threshold_ensemble_table <-
      get_evaluations(myBiomodEnsembleOut_noCat, as.data.frame = TRUE) %>% 
      dplyr::filter(Eval.metric == "TSS") %>% 
      dplyr::mutate(full.name = glue::glue("GuloGulo_{Model.name}")) %>% 
      dplyr::select(full.name, Cutoff)
    rownames(threshold_ensemble_table) <- threshold_ensemble_table$full.name
    threshold_ensemble_table <- threshold_ensemble_table[,-1, drop = FALSE]
    
    EnsembleCurrent.df <- 
      get_predictions(myBiomodEnsembleForecast.df, 
                      as.data.frame = TRUE) %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table[
            get_projected_models(myBiomodEnsembleForecast), "Cutoff"])
    
    EnsembleCurrent.1.df <- 
      get_predictions(myBiomodEnsembleForecast.df, 
                      as.data.frame = TRUE,
                      model = "EMca") %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table["GuloGulo_EMcaByTSS_mergedAlgo_mergedRun_mergedData",
                                   "Cutoff"])
    EnsembleCurrent.all.SpatRaster <-
      get_predictions(myBiomodEnsembleForecast) %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table[
            get_projected_models(myBiomodEnsembleForecast), "Cutoff"]) 
    
    EnsembleCurrent.1.SpatRaster <- 
      get_predictions(
        myBiomodEnsembleForecast,
        model = "EMca") %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table["GuloGulo_EMcaByTSS_mergedAlgo_mergedRun_mergedData",
                                   "Cutoff"])
    
    ## Ensemble model - future --------------------------------------------------
    myBiomodEnsembleFuture <-
      BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEnsembleOut_noCat,
        proj.name = 'Future',
        new.env = myExplFuture,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    myBiomodEnsembleFuture.df <-
      BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEnsembleOut_noCat,
        proj.name = 'Future.df',
        new.env = myExplFuture.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    EnsembleFuture.df <- 
      get_predictions(myBiomodEnsembleFuture.df, 
                      as.data.frame = TRUE) %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table[
            get_projected_models(myBiomodEnsembleForecast), "Cutoff"]) 
    EnsembleFuture.1.df <- 
      get_predictions(myBiomodEnsembleFuture.df, 
                      as.data.frame = TRUE,
                      model = "EMca") %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table["GuloGulo_EMcaByTSS_mergedAlgo_mergedRun_mergedData",
                                   "Cutoff"])
    
    EnsembleFuture.2.df <- cbind(EnsembleFuture.1.df,EnsembleFuture.1.df)
    colnames(EnsembleFuture.2.df) <- c("ssp1","ssp2")
    
    EnsembleFuture.all.SpatRaster <- 
      get_predictions(myBiomodEnsembleFuture) %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table[
            get_projected_models(myBiomodEnsembleForecast), "Cutoff"]) 
    
    EnsembleFuture.1.SpatRaster <- 
      get_predictions(
        myBiomodEnsembleFuture,
        model = "EMca") %>% 
      bm_BinaryTransformation(
        threshold = 
          threshold_ensemble_table["GuloGulo_EMcaByTSS_mergedAlgo_mergedRun_mergedData",
                                   "Cutoff"])
    
    EnsembleFuture.2.SpatRaster <- c(EnsembleFuture.1.SpatRaster,
                                     EnsembleFuture.1.SpatRaster)
    names(EnsembleFuture.2.SpatRaster) <- c("ssp1","ssp2")
    
  })))
)


# remove NA in df ---------------------------------------------------------

ProjCurrent.df <- ProjCurrent.df 
ProjFuture.df <- ProjFuture.df


ProjCurrent.df
ProjFuture.df

ProjCurrent.1.df
ProjFuture.1.df

ProjFuture.2.df


# Check that non binary output fails ------------------------------------------

cli::cli_h2("non binary output check")


## SpatRaster -----------------------------------------------------------

cli::cli_process_start("SpatRaster")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      # both non binary
      myTest <- tinytest::expect_error(
        myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.nonbinary,
                                        proj.future = ProjFuture.nonbinary)
      )
      stopifnot(myTest)
      
      # Current non binary
      myTest <- tinytest::expect_error(
        myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.nonbinary,
                                        proj.future = ProjFuture.all.SpatRaster)
      )
      stopifnot(myTest)
      
      # Future non binary
      myTest <- tinytest::expect_error(
        myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.all.SpatRaster,
                                        proj.future = ProjFuture.nonbinary)
      )
      stopifnot(myTest)
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame --------------------------------------------------------------

cli::cli_process_start("data.frame")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      # myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.df,
      #                                 proj.future = ProjFuture.df)
      
      # both non binary
      myTest <- tinytest::expect_error(
        myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.nonbinary.df,
                                        proj.future = ProjFuture.nonbinary.df)
      )
      stopifnot(myTest)
      
      # Current non binary
      myTest <- tinytest::expect_error(
        myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.nonbinary.df,
                                        proj.future = ProjFuture.df)
      )
      stopifnot(myTest)
      
      # Future non binary
      myTest <- tinytest::expect_error(
        myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.df,
                                        proj.future = ProjFuture.nonbinary.df)
      )
      stopifnot(myTest)
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# Simple Model ------------------------------------------------------------

cli::cli_h2("Simple Model")

## SpatRaster ------------------------------------------------------------
### all vs all ------------------------------------------------------------
cli::cli_process_start("SpatRaster all vs all")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.all.SpatRaster,
                                      proj.future = ProjFuture.all.SpatRaster)
      
      plot(myRangeSize$Diff.By.Pixel)
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### 1 Current vs 2 Future ------------------------------------------------------------

cli::cli_process_start("SpatRaster 1 Current vs 2 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.1.SpatRaster,
                                      proj.future = ProjFuture.2.SpatRaster)
      
      plot(myRangeSize$Diff.By.Pixel)
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### 1 Current vs 1 Future ------------------------------------------------------------

cli::cli_process_start("SpatRaster 1 Current vs 1 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.1.SpatRaster,
                                      proj.future = ProjFuture.1.SpatRaster)
      
      plot(myRangeSize$Diff.By.Pixel)
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      # stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## RasterLayer/RasterStack ------------------------------------------------------------
### all vs all ------------------------------------------------------------
cli::cli_process_start("RasterLayer/RasterStack all vs all")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = stack(ProjCurrent.all.SpatRaster),
                                      proj.future = stack(ProjFuture.all.SpatRaster))
      
      plot(myRangeSize$Diff.By.Pixel)
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### 1 Current vs 2 Future ------------------------------------------------------------

cli::cli_process_start("RasterLayer/RasterStack 1 Current vs 2 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = raster(ProjCurrent.1.SpatRaster),
                                      proj.future = stack(ProjFuture.2.SpatRaster))
      
      plot(myRangeSize$Diff.By.Pixel)
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### 1 Current vs 1 Future ------------------------------------------------------------

cli::cli_process_start("RasterLayer/RasterStack 1 Current vs 1 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = raster(ProjCurrent.1.SpatRaster),
                                      proj.future = raster(ProjFuture.1.SpatRaster))
      
      plot(myRangeSize$Diff.By.Pixel)
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      # stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## data.frame ------------------------------------------------------------

### all vs all ------------------------------------------------------------
cli::cli_process_start("data.frame all vs all")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.df,
                                      proj.future = ProjFuture.df)
      
      # myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      # stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      # stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      # stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### 1 Current vs 2 Future ------------------------------------------------------------

cli::cli_process_start("data.frame 1 Current vs 2 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.1.df,
                                      proj.future = ProjFuture.2.df)
      
      # myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      # stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      # stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      # stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### 1 Current vs 1 Future ------------------------------------------------------------

cli::cli_process_start("data.frame 1 Current vs 1 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.1.df,
                                      proj.future = ProjFuture.1.df)
      
      # myPlotRangeSize <- bm_PlotRangeSize(myRangeSize)
      # stopifnot(inherits(myPlotRangeSize$plot.ca, "ggplot"))
      # stopifnot(inherits(myPlotRangeSize$plot.count, "ggplot"))
      # stopifnot(inherits(myPlotRangeSize$plot.perc, "ggplot"))
      
    })))
  )
}, silent = TRUE)
dev.off()

if(inherits(this_try, "try-error")){
  Error_RangeSize <- Error_RangeSize + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
