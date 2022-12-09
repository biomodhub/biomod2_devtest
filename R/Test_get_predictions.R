cli::cli_h1("get_predictions(biomod.projection.out)")

Error_get_pred <- 0

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

myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )
myResp.spdf <- as(myResp.SpatVector, "Spatial")
myResp.SpatialPoints <- sp::SpatialPoints(myResp.spdf[1:200,])
myResp.SpatVectorPO <- vect(myResp.SpatialPoints)

## myExpl ----------------------------------

myExpl.df <- extract(myExpl, y = myResp.SpatVector, ID = FALSE)
myExpl.df1 <- myExpl.df[,1, drop = FALSE]
myExpl.matrix <- as.matrix(myExpl.df)
myExpl.matrix1 <- myExpl.matrix[,1, drop = FALSE]
myExpl.raster <- stack(myExpl)
myExpl.raster1 <- stack(myExpl[[1]])
myExpl1 <- myExpl[[1]]
tmpdf <- cbind(myRespXY,myExpl.df)
rownames(tmpdf) <- NULL
myExpl.SpatVector <- vect(
  tmpdf,
  geom = c("X_WGS84", "Y_WGS84")
)
myExpl.SpatVector$ID <- NULL
myExpl.spdf <- as(myExpl.SpatVector, "Spatial")


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
                                active = 2) 

myExpl.cat.df <- extract(myExpl.cat, y = myResp.SpatVector, ID =FALSE)
myExpl.cat.raster <- myExpl.raster
bio3_true <- myExpl.raster[[1]]
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

tmpdf <- cbind(myRespXY,myExpl.cat.df)
rownames(tmpdf) <- NULL
myExpl.cat.SpatVector <- vect(
  tmpdf,
  geom = c("X_WGS84", "Y_WGS84")
)
myExpl.cat.spdf <- as(myExpl.cat.SpatVector, "Spatial")

## Load Models -------------------------------------------------------------

file.out <- paste0(myRespName, "/", myRespName, ".NoCat_Eval_Presence-Only.models.out")
myBiomodModelOut_Eval <- 
  try(suppressWarnings(
    get(load(file.out))
  ), silent = TRUE)

model_subset <-  get_built_models(myBiomodModelOut_Eval)[
  c(1:3)+rep(c(0,11,22,33), each = 3)
  ]

try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodEnsembleOut_noCat <- BIOMOD_EnsembleModeling(
        bm.mod = myBiomodModelOut_Eval,
        models.chosen = model_subset,
        em.by = 'PA_dataset+repet',
        metric.select = c('TSS'),
        metric.select.thresh = c(0.7),
        var.import = 3,
        metric.eval = c('TSS', 'ROC'),
        prob.mean = TRUE,
        prob.median = TRUE,
        prob.cv = FALSE,
        prob.ci = FALSE,
        prob.ci.alpha = 0.05,
        committee.averaging = FALSE,
        prob.mean.weight = FALSE,
        prob.mean.weight.decay = 'proportional',
        seed.val = 42)
    })))
  )
}, silent = TRUE)



# Simple Models ------------------------------------------------

cli::cli_h2("Simple Models")

## SpatRaster ------------------------------------------------

cli::cli_h3("Predictions from SpatRaster")

try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProjOut_Eval <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut_Eval,
        proj.name = 'Current',
        new.env = myExpl,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    })))
  )
}, silent = TRUE)


### standard ------------------------------------------------

cli::cli_process_start("Standard")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### full.name = subset ------------------------------------------------

cli::cli_process_start("full.name = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      my_subset <- get_projected_models(myBiomodProjOut_Eval)[1:11]
      myPred <- get_predictions(myBiomodProjOut_Eval, full.name = my_subset)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### algo = subset ------------------------------------------------

cli::cli_process_start("algo = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myModel <- c("GLM","ANN", "SRE")
      myPred <- get_predictions(myBiomodProjOut_Eval, algo = myModel)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### run.eval = subset ------------------------------------------------

cli::cli_process_start("run.eval = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = c("RUN1","RUN2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = "RUN1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      

    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### data.set = subset ------------------------------------------------

cli::cli_process_start("data.set = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = c("PA1","PA2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = "PA1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      

    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame ------------------------------------------------

cli::cli_h3("Predictions from data.frame")

try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProjOut_Eval <- BIOMOD_Projection(
        bm.mod = myBiomodModelOut_Eval,
        proj.name = 'Current',
        new.env = myExpl.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### standard ------------------------------------------------

cli::cli_process_start("standard")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### full.name = subset ------------------------------------------------

cli::cli_process_start("full.name = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      my_subset <- get_projected_models(myBiomodProjOut_Eval)[1:10]
      myPred <- get_predictions(myBiomodProjOut_Eval, full.name = my_subset)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      

    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### algo =  subset ------------------------------------------------

cli::cli_process_start("full.name = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myModel <- c("GLM","ANN", "SRE")
      myPred <- get_predictions(myBiomodProjOut_Eval, algo = myModel)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### run.eval = subset ------------------------------------------------

cli::cli_process_start("run.eval = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = c("RUN1","RUN2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = "RUN1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### data.set = subset ------------------------------------------------

cli::cli_process_start("data.set = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = c("PA1","PA2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = "PA1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      

    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# Ensemble Models ------------------------------------------------

cli::cli_h2("Ensemble Models")

## SpatRaster ------------------------------------------------

cli::cli_h3("Predictions from SpatRaster")

try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProjOut_Eval <- BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEnsembleOut_noCat,
        proj.name = 'Current',
        new.env = myExpl,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    })))
  )
}, silent = TRUE)


### standard ------------------------------------------------


cli::cli_process_start("standard")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### full.name = subset ------------------------------------------------

cli::cli_process_start("full.name = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      my_subset <- get_projected_models(myBiomodProjOut_Eval)[c(2,6,7)]
      myPred <- get_predictions(myBiomodProjOut_Eval, full.name = my_subset)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### algo = subset ------------------------------------------------

cli::cli_process_start("full.name = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myModel <- c("EMmean")
      myPred <- get_predictions(myBiomodProjOut_Eval, algo = myModel)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      myModel <- c("EMmedian","EMmean")
      myPred <- get_predictions(myBiomodProjOut_Eval, algo = myModel)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### run.eval = subset ------------------------------------------------

cli::cli_process_start("run.eval = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = c("RUN1","RUN2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = "RUN1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = "RUN2")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### data.set = subset ------------------------------------------------

cli::cli_process_start("data.set = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = c("PA1","PA2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = "PA1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = "PA2")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "SpatRaster"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame ------------------------------------------------

cli::cli_h3("Predictions from data.frame")

try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProjOut_Eval <- BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEnsembleOut_noCat,
        proj.name = 'Current',
        new.env = myExpl.df,
        metric.binary = 'all',
        metric.filter = 'all',
        build.clamping.mask = FALSE,
        do.stack = TRUE
      )
    })))
  )
}, silent = TRUE)



### Standard ------------------------------------------------

cli::cli_process_start("Standard")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### full.name = subset ------------------------------------------------

cli::cli_process_start("full.name = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      my_subset <- get_projected_models(myBiomodProjOut_Eval)[c(2,6,7)]
      
      myPred <- get_predictions(myBiomodProjOut_Eval, full.name = my_subset)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### algo = subset ------------------------------------------------

cli::cli_process_start("full.name = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myModel <- c("EMmedian")
      myPred <- get_predictions(myBiomodProjOut_Eval, algo = myModel)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
      myModel <- c("EMmedian","EMmean")
      myPred <- get_predictions(myBiomodProjOut_Eval, algo = myModel)
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### run.eval = subset ------------------------------------------------

cli::cli_process_start("run.eval = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = c("RUN1","RUN2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = "RUN1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, run.eval = "RUN2")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### data.set = subset ------------------------------------------------

cli::cli_process_start("data.set = subset")


this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = c("PA1","PA2"))
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = "PA1")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
      myPred <- get_predictions(myBiomodProjOut_Eval, data.set = "PA2")
      stopifnot(!is.null(myPred))
      stopifnot(inherits(myPred, "data.frame"))
      
    })))
  )
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_get_pred <- Error_get_pred + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}