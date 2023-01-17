cli::cli_h1("build_clamping_mask")

Error_Clamping <- 0

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
  myExpl <- raster::stack(system.file(myFiles, package = 'biomod2'))
}
myRespName <- 'GuloGulo'

## myResp ------------------------------------------------------------------
myResp <- as.numeric(DataSpecies[, myRespName])
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )
myResp.spdf <- as(myResp.SpatVector, "Spatial")

## myExpl.cat ----------------------------------

if(terraVersion){
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
} else {
  myExpl.cat <- myExpl
  bio3_true <- myExpl[[1]]
  bio3_factor <- calc(bio3_true, function(x){
    factor(
      ifelse(
        x <= 20,"low",  
        ifelse( x <= 40,"medium",
                ifelse( x <= 60, "high", "very high"))),
      levels = c("low","medium","high", "very high")
    )
  })
  # plot(bio3_factor)
  names(bio3_factor) <- "bio3"
  myExpl.cat[[1]] <- bio3_factor
  
}


## myExpl.df ---------------------------------------------------------------

if(terraVersion){
  myExpl.df <- extract(myExpl, y = myResp.SpatVector, ID = FALSE)
  myExpl.cat.df <- extract(myExpl.cat, y = myResp.SpatVector, ID = FALSE)
} else {
  myExpl.df <- extract(myExpl, y = myResp.spdf)
  myExpl.cat.df <- extract(myExpl.cat, y = myResp.spdf)
}

# build model -------------------------------------------------

this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      file.out <- paste0(myRespName, "/", myRespName, ".clamping_nocat.models.out")
      if (file.exists(file.out)) {
        myBiomodModelOut <-  get(load(file.out))
      } else {
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
            modeling.id = 'clamping_nocat',
            nb.rep = 2,
            data.split.perc = 80,
            var.import = 3,
            metric.eval = c('TSS','ROC'),
            do.full.models = FALSE,
            seed.val = 42,
            nb.cpu = 4)
      }
      
      
      file.out <- paste0(myRespName, "/", myRespName, ".clamping_cat.models.out")
      if (file.exists(file.out)) {
        myBiomodModelOut.cat <-  get(load(file.out))
      } else {
        
        myBiomodData.cat <- 
          BIOMOD_FormatingData(
            resp.var = myResp,
            expl.var = myExpl.cat,
            resp.xy = myRespXY,
            resp.name = myRespName)
        
        myBiomodModelOut.cat <- 
          BIOMOD_Modeling(
            bm.format = myBiomodData.cat,
            bm.options = BIOMOD_ModelingOptions(),
            modeling.id = 'clamping_cat',
            nb.rep = 2,
            data.split.perc = 80,
            var.import = 3,
            metric.eval = c('TSS','ROC'),
            do.full.models = FALSE,
            seed.val = 42,
            nb.cpu = 4)
      }
    }))
  )
}, silent = TRUE)


# No Categorical variables ------------
cli::cli_h2("No categorical variables")

## SpatRaster ------------

cli::cli_process_start("SpatRaster")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      MinMax <- get_formal_data(myBiomodModelOut, 'MinMax')
      
      tic()
      myMask <- .build_clamping_mask(myExpl, MinMax = MinMax)
      toc()
      
      myProj <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut,
        proj.name = 'Current',
        new.env = myExpl,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = TRUE,
        do.stack = TRUE,
        nb.cpu = 4)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Clamping <- Error_Clamping + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame ------------

cli::cli_process_start("data.frame")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      MinMax <- get_formal_data(myBiomodModelOut, 'MinMax')
      myMask <- .build_clamping_mask(myExpl.df, MinMax = MinMax)
      
      myProj <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut,
        proj.name = 'Current',
        new.env = myExpl.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = TRUE,
        do.stack = TRUE,
        nb.cpu = 4)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Clamping <- Error_Clamping + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# No Categorical variables ------------
cli::cli_h2("With categorical variables")

## SpatRaster ------------

cli::cli_process_start("SpatRaster")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      MinMax <- get_formal_data(myBiomodModelOut.cat, 'MinMax')
      myMask <- .build_clamping_mask(myExpl.cat, MinMax = MinMax)
      
      myProj <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut.cat,
        proj.name = 'Current',
        new.env = myExpl.cat,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = TRUE,
        do.stack = TRUE,
        nb.cpu = 4)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Clamping <- Error_Clamping + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame ------------

cli::cli_process_start("data.frame")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      MinMax <- get_formal_data(myBiomodModelOut.cat, 'MinMax')
      myMask <- .build_clamping_mask(myExpl.cat.df, MinMax = MinMax)
      
      myProj <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut.cat,
        proj.name = 'Current',
        new.env = myExpl.cat.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = TRUE,
        do.stack = TRUE,
        nb.cpu = 4)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Clamping <- Error_Clamping + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
