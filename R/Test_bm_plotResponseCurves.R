cli::cli_h1("bm_plotResponseCurves")

# setup data --------------------------------------------------------------
cli::cli_h2("Setup data")
Error_PlotResponseCurves <- 0

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

# myExpl.cat.raster -------------------------------------------------------
myExpl.cat.raster <- stack(myExpl)
bio3_true <- myExpl.cat.raster[[1]]
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
myExpl.cat.raster[[1]] <- bio3_factor


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


# with df -----------------------------------------------------------------

myResp <- as.numeric(DataSpecies[, myRespName])
myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )
myExpl.df <- extract(myExpl, y = myResp.SpatVector)
myExpl.cat.df <- extract(myExpl.cat, y = myResp.SpatVector, ID =FALSE)


# Single Model ; No Categorical Variables --------------------------------------

invisible(
  capture.output(suppressWarnings(suppressMessages({
    myBiomodData <- 
      BIOMOD_FormatingData(
        resp.var = myResp,
        expl.var = myExpl,
        resp.xy = myRespXY,
        resp.name = myRespName)
    
    file.out <- paste0(myRespName, "/", myRespName, ".NoCat_NoEval_Presence-Absence.models.out")
    if (file.exists(file.out)) {
      myBiomodModelOut <- get(load(file.out))
    } else {
      myBiomodModelOut <- 
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_NoEval_Presence-Absence',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
          seed.val = 42
        )
    }
  }))))


## new.env = NULL  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = NULL")


### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "mean",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "min",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "max",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## new.env = SpatRaster  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = SpatRaster")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               new.env = myExpl,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               new.env = myExpl,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               new.env = myExpl,
                                               fixed.Var = "max",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## new.env = matrix  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = matrix")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = as.matrix(myExpl.df),
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## new.env = data.frame  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = data.frame")

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.df,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               new.env = myExpl.df,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               new.env = myExpl.df,
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               new.env = myExpl.df,
                                               fixed.Var = "min",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# Single Model ; With Categorical Variables --------------------------------------

invisible(
  capture.output(suppressWarnings(suppressMessages({
    myBiomodData <- 
      BIOMOD_FormatingData(
        resp.var = myResp,
        expl.var = myExpl.cat,
        resp.xy = myRespXY,
        resp.name = myRespName)
    
    file.out <- paste0(myRespName, "/", myRespName, ".Cat_NoEval_Presence-Absence.models.out")
    if (file.exists(file.out)) {
      myBiomodModelOut <- get(load(file.out))
    } else {
      myBiomodModelOut <- 
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_NoEval_Presence-Absence',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
          seed.val = 42
        )
    }
  }))))

## new.env = NULL  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = NULL")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "mean",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "min",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               fixed.Var = "max",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## new.env = RasterStack  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = RasterStack")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.raster,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## new.env = SpatRaster  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = SpatRaster")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               new.env = myExpl.cat,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               new.env = myExpl.cat,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               new.env = myExpl.cat,
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## new.env = matrix  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = matrix")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      thistest <- tinytest::expect_error(
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = as.matrix(myExpl.cat.df),
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      )
      stopifnot(thistest)
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## new.env = data.frame  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = data.frame")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1],
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:4],
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodModelOut, 
                                               models.chosen = get_built_models(myBiomodModelOut)[1:2],
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# Ensemble Model ; No Categorical Variables --------------------------------------
cli::cli_h2("Ensemble models without Categorical Variables")

invisible(
  capture.output(suppressWarnings(suppressMessages({
    myBiomodData <- 
      BIOMOD_FormatingData(
        resp.var = myResp,
        expl.var = myExpl,
        resp.xy = myRespXY,
        resp.name = myRespName)
    
    file.out <- paste0(myRespName, "/", myRespName, ".NoCat_NoEval_Presence-Absence.ensemble.models.out")
    if (file.exists(file.out)) {
      myBiomodEnsembleModelOut <- get(load(file.out))
    } 
  }))))

## new.env = NULL  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = NULL")


### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "mean",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "min",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "max",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## new.env = SpatRaster  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = SpatRaster")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               new.env = myExpl,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               new.env = myExpl,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               new.env = myExpl,
                                               fixed.Var = "max",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## new.env = matrix  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = matrix")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = as.matrix(myExpl.df),
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## new.env = data.frame  -------------------------------------------------------------
cli::cli_h3("No Categorical Variables ; new.env = data.frame")

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.df,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               new.env = myExpl.df,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               new.env = myExpl.df,
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               new.env = myExpl.df,
                                               fixed.Var = "min",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# Ensemble Model ; With Categorical Variables --------------------------------------
cli::cli_h2("Ensemble models with Categorical Variables")

invisible(
  capture.output(suppressWarnings(suppressMessages({
    myBiomodData <- 
      BIOMOD_FormatingData(
        resp.var = myResp,
        expl.var = myExpl.cat,
        resp.xy = myRespXY,
        resp.name = myRespName)
    
    file.out <- paste0(myRespName, "/", myRespName, ".Cat_NoEval_Presence-Absence.ensemble.models.out")
    if (file.exists(file.out)) {
      myBiomodEnsembleModelOut <- get(load(file.out))
    } 
  }))))

## new.env = NULL  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = NULL")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "mean",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "min",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               fixed.Var = "max",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## new.env = SpatRaster  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = SpatRaster")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               new.env = myExpl.cat,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               new.env = myExpl.cat,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               new.env = myExpl.cat,
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## new.env = matrix  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = matrix")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      thistest <- tinytest::expect_error(
        myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                                 models.chosen = "all",
                                                 new.env = as.matrix(myExpl.cat.df),
                                                 fixed.Var = "mean",
                                                 do.bivariate = FALSE,
                                                 do.plot = FALSE,
                                                 do.progress = FALSE)
      )
      stopifnot(thistest)
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## new.env = data.frame  -------------------------------------------------------------
cli::cli_h3("With Categorical Variables ; new.env = data.frame")

### fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "median",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = max ; models.chosen = all ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = "all",
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "max",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = min ; models.chosen = 1 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1],
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "min",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE --------------
cli::cli_process_start("fixed.Var = mean ; models.chosen = 4 ; do.bivariate = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:4],
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "mean",
                                               do.bivariate = FALSE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE --------------
cli::cli_process_start("fixed.Var = median ; models.chosen = 2 ; do.bivariate = TRUE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myResponseCurve <- bm_PlotResponseCurves(myBiomodEnsembleModelOut, 
                                               models.chosen = get_built_models(myBiomodEnsembleModelOut)[1:2],
                                               new.env = myExpl.cat.df,
                                               fixed.Var = "median",
                                               do.bivariate = TRUE,
                                               do.plot = FALSE,
                                               do.progress = FALSE)
      if (any(is.na(myResponseCurve$tab$expl.val))) stop("NA in variables")
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_PlotResponseCurves <- Error_PlotResponseCurves + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
