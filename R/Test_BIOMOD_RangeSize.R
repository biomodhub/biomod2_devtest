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
    
    ProjCurrent.nonbinary.df <- get_predictions(myBiomodProj.df)
    
    ProjCurrent.nonbinary <-    get_predictions(myBiomodProj)
    ProjCurrent.nonbinary.raster <- stack(ProjCurrent.nonbinary)
    
    
    ProjCurrent.df <- 
      get_predictions(myBiomodProj.df, metric.binary = "TSS")
    
    ProjCurrent.1.df <-
      get_predictions(myBiomodProj.df,
                      algo = "GLM",
                      run = "RUN1",
                      metric.binary = "TSS")
    ProjCurrent.all.SpatRaster <-
      get_predictions(myBiomodProj, metric.binary = "TSS")
    
    ProjCurrent.1.SpatRaster <-
      get_predictions(myBiomodProj,
                      algo = "GLM",
                      run = "RUN1",
                      metric.binary = "TSS") 
    
    
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
    
    ProjFuture.nonbinary.df <-  get_predictions(myBiomodFuture.df) 
    ProjFuture.nonbinary <-  get_predictions(myBiomodProj)
    
    ProjFuture.nonbinary.raster <- stack(ProjFuture.nonbinary)
    
    ProjFuture.df <-
      get_predictions(myBiomodFuture.df, metric.binary = "TSS")
    
    ProjFuture.1.df <- 
      get_predictions(myBiomodFuture.df,
                      algo = "GLM",
                      run = "RUN1",
                      metric.binary = "TSS")
    ProjFuture.2.df <- rbind(ProjFuture.1.df,ProjFuture.1.df)
    ProjFuture.2.df$full.name <- rep(c("ssp1","ssp2"), each = nrow(ProjFuture.1.df))
    
    ProjFuture.all.SpatRaster <-
      get_predictions(myBiomodFuture, metric.binary = "TSS") 
    
    ProjFuture.1.SpatRaster <- 
      get_predictions(myBiomodFuture,
                      algo = "GLM",
                      run = "RUN1",
                      metric.binary = "TSS")
    
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
    
    EnsembleCurrent.df <- 
      get_predictions(myBiomodEnsembleForecast.df, metric.binary = "TSS") 
    
    EnsembleCurrent.1.df <- 
      get_predictions(myBiomodEnsembleForecast.df,
                      algo = "EMca",
                      metric.binary = "TSS")
    
    EnsembleCurrent.all.SpatRaster <-
      get_predictions(myBiomodEnsembleForecast, 
                      metric.binary = "TSS")
    
    EnsembleCurrent.1.SpatRaster <- 
      get_predictions(
        myBiomodEnsembleForecast,
        algo = "EMca",
        metric.binary = "TSS") 
    
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
      get_predictions(myBiomodEnsembleFuture.df, metric.binary = "TSS") 
    EnsembleFuture.1.df <- 
      get_predictions(myBiomodEnsembleFuture.df,
                      algo = "EMca", metric.binary = "TSS") 
    
    EnsembleFuture.2.df <- cbind(EnsembleFuture.1.df,EnsembleFuture.1.df)
    colnames(EnsembleFuture.2.df) <- c("ssp1","ssp2")
    
    EnsembleFuture.all.SpatRaster <- 
      get_predictions(myBiomodEnsembleFuture, metric.binary = "TSS")
    
    EnsembleFuture.1.SpatRaster <- 
      get_predictions(
        myBiomodEnsembleFuture,
        algo = "EMca",
        metric.binary = "TSS")
    
    EnsembleFuture.2.SpatRaster <- c(EnsembleFuture.1.SpatRaster,
                                     EnsembleFuture.1.SpatRaster)
    names(EnsembleFuture.2.SpatRaster) <- c("ssp1","ssp2")
    
  })))
)




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
      
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize, do.maps = FALSE, do.mean = FALSE, do.plot = FALSE)
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

cli::cli_process_start("data.frame 1 Current vs 2 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.1.df,
                                      proj.future = ProjFuture.2.df)
      
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize, do.maps = FALSE, do.mean = FALSE, do.plot = FALSE)
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

cli::cli_process_start("data.frame 1 Current vs 1 Future")
pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      
      myRangeSize <- BIOMOD_RangeSize(proj.current = ProjCurrent.1.df,
                                      proj.future = ProjFuture.1.df)
      
      myPlotRangeSize <- bm_PlotRangeSize(myRangeSize, do.maps = FALSE, do.mean = FALSE, do.plot = FALSE)
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
