cli::cli_h1("BIOMOD_Parallel")

Error_Parallel <- 0

# Preparation -------------------------------------------------------------

library(terra)
library(raster)
if(terraVersion){
  data("DataSpecies")
  data("bioclim_current")
  myExpl <- terra::rast(bioclim_current)
} else {
  myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
  DataSpecies <- read.csv(myFile, row.names = 1)
  myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
  myExpl <- rast(raster::stack(system.file(myFiles, package = 'biomod2')))
}
myRespName <- 'GuloGulo'

## myExpl.cat ----------------------------------

myExpl.cat <- myExpl
myExpl.cat[[1]] <-
  terra::classify(
    myExpl.cat[[1]],
    matrix(
      c(0,20,1,
        20,40,2,
        40,60,3,
        60,Inf,4),
      ncol = 3, byrow = TRUE
    )
  )
names(myExpl.cat)[1] <- "bio3_factor"
myExpl.cat <- terra::categories(myExpl.cat, layer = 1, 
                                data.frame(ID = c(1,2,3,4),
                                           bio3 = c("low","medium","high","very high")),
                                active = 1) 


## myResp ------------------------------------------------------------------
myResp <- as.numeric(DataSpecies[, myRespName])
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

tmpindex <- which(myResp == 1)
myResp_PO <- myResp[tmpindex]
myRespXY_PO <- myRespXY[tmpindex,]

tmpindex <- which(myResp == 0)
myResp_PO_NA <- myResp
myResp_PO_NA[tmpindex] <- NA
myRespXY_PO_NA <- myRespXY

tmpindex <- sample(which(myResp == 0), size = 500, replace = FALSE)
myResp_PA2 <- myResp
myResp_PA2[tmpindex] <- NA
myRespXY_PA2 <- myRespXY


# Parallelization ------------------------------------------------


## BIOMOD_Modeling ------------
cli::cli_process_start("BIOMOD_Modeling")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName)
      
      myBiomodModelOut <- 
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Parallel',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
          seed.val = 42,
          nb.cpu = 4)
      get_predictions(myBiomodModelOut)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Parallel <- Error_Parallel + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## BIOMOD_Projection ------------
cli::cli_process_start("BIOMOD_Projection")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      
      
      myProj <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut,
        proj.name = 'Current',
        new.env = myExpl,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE,
        nb.cpu = 4)
      get_predictions(myProj)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Parallel <- Error_Parallel + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## BIOMOD_EnsembleModeling ------------
cli::cli_process_start("BIOMOD_EnsembleModeling")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({

      
      myBiomodEM <- BIOMOD_EnsembleModeling(
        bm.mod = myBiomodModelOut,
        em.by = 'all',
        metric.select = c('TSS'),
        metric.select.thresh = c(0.7),
        var.import = 3,
        metric.eval = c('TSS', 'ROC'),
        prob.mean = TRUE,
        prob.median = TRUE,
        prob.cv = TRUE,
        prob.ci = TRUE,
        prob.ci.alpha = 0.05,
        committee.averaging = TRUE,
        prob.mean.weight = TRUE,
        prob.mean.weight.decay = 'proportional',
        seed.val = 42,
        nb.cpu = 4)
      get_predictions(myBiomodEM)
      
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Parallel <- Error_Parallel + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## BIOMOD_EnsembleForecasting ------------
cli::cli_process_start("BIOMOD_EnsembleForecasting")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myProjEM <- BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEM,
        proj.name = 'Current',
        new.env = myExpl,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE,
        nb.cpu = 4)
      get_predictions(myProjEM)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Parallel <- Error_Parallel + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}