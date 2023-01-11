cli::cli_h1("BIOMOD_Projection")

Error_Projection <- 0

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


file.out <- paste0(myRespName, "/", myRespName, ".NoCat_NoEval_Presence-Absence.models.out")
myBiomodModelOut_noCat <- 
  try(suppressWarnings(
    get(load(file.out))
  ), silent = TRUE)

file.out <- paste0(myRespName, "/", myRespName, ".Cat_NoEval_Presence-Absence.models.out")
myBiomodModelOut_Cat <- 
  try(suppressWarnings(
    get(load(file.out))
  ), silent = TRUE)

file.out <- paste0(myRespName, "/", myRespName, ".NoCat1_NoEval_Presence-Absence.models.out")
myBiomodModelOut_noCat1 <- 
  try(suppressWarnings(
    get(load(file.out))
  ), silent = TRUE)

# No Categorical Variables ------------------------------------------------

cli::cli_h2("No Categorical Variables")

## SpatRaster ------------
cli::cli_process_start("SpatRaster")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
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
      
      get_projected_models(myBiomodProj)
      get_projected_models(myBiomodProj, run = "RUN1")
      get_projected_models(myBiomodProj, algo = "GAM")
      get_projected_models(myBiomodProj, run = "RUN1", algo = "GBM")
      get_projected_models(myBiomodProj, 
                           full.name = c("GuloGulo_allData_RUN2_MARS",
                                         "GuloGulo_allData_RUN1_MAXENT"))
      
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.binary = "TSS"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.filter = "TSS", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.binary = "ROC", algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS", algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS",
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## SpatRaster & do.stack = FALSE ------------
cli::cli_process_start("SpatRaster & do.stack = FALSE")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat,
          proj.name = 'CurrentNoStack',
          new.env = myExpl,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = FALSE
        )
      
      get_projected_models(myBiomodProj)
      get_projected_models(myBiomodProj, run = "RUN1")
      get_projected_models(myBiomodProj, algo = "GAM")
      get_projected_models(myBiomodProj, run = "RUN1", algo = "GBM")
      get_projected_models(myBiomodProj, 
                           full.name = c("GuloGulo_allData_RUN2_MARS",
                                         "GuloGulo_allData_RUN1_MAXENT"))
      
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.binary = "TSS"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.filter = "TSS", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.binary = "ROC", algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS", algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS",
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## SpatRaster1 ------------
cli::cli_process_start("SpatRaster1")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat1,
          proj.name = 'Current',
          new.env = myExpl1,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      
      get_projected_models(myBiomodProj)
      get_projected_models(myBiomodProj, run = "RUN1")
      get_projected_models(myBiomodProj, algo = "GAM")
      get_projected_models(myBiomodProj, run = "RUN1", algo = "GBM")
      get_projected_models(myBiomodProj, 
                           full.name = c("GuloGulo_allData_RUN2_MARS",
                                         "GuloGulo_allData_RUN1_MAXENT"))
      
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.binary = "TSS"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.filter = "TSS", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.binary = "ROC", algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS", algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS",
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## Raster ------------
cli::cli_process_start("Raster")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat,
          proj.name = 'Current',
          new.env = myExpl.raster,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## Raster1 ------------
cli::cli_process_start("Raster1")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat1,
          proj.name = 'Current',
          new.env = myExpl.raster1,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## Matrix ------------
cli::cli_process_start("Matrix")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat,
          proj.name = 'Current',
          new.env = myExpl.matrix,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      stopifnot(tinytest::expect_error(plot(myBiomodProj)))
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## Matrix1 ------------
cli::cli_process_start("Matrix1")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat1,
          proj.name = 'Current',
          new.env = myExpl.matrix1,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      stopifnot(tinytest::expect_error(plot(myBiomodProj)))
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## Matrix with coordinates ------------
cli::cli_process_start("Matrix with coord")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat,
          proj.name = 'Current',
          new.env = myExpl.matrix,
          new.env.xy = crds(myExpl.SpatVector),
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      
      get_projected_models(myBiomodProj)
      get_projected_models(myBiomodProj, run = "RUN1")
      get_projected_models(myBiomodProj, algo = "GAM")
      get_projected_models(myBiomodProj, run = "RUN1", algo = "GBM")
      get_projected_models(myBiomodProj, 
                           full.name = c("GuloGulo_allData_RUN2_MARS",
                                         "GuloGulo_allData_RUN1_MAXENT"))
      
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.binary = "TSS"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.filter = "TSS", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj,
             metric.binary = "ROC", algo = "GLM"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS", algo = "GLM", run = "RUN1"), "ggplot"))
      stopifnot(inherits(
        plot(myBiomodProj, 
             metric.filter = "TSS",
             full.name = c("GuloGulo_allData_RUN2_MARS",
                           "GuloGulo_allData_RUN1_MAXENT")),
        "ggplot"))
      
      dev.off()
      
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## Matrix1 with coord------------
cli::cli_process_start("Matrix1 with coord")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat1,
          proj.name = 'Current',
          new.env = myExpl.matrix1,
          new.env.xy = crds(myExpl.SpatVector),
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## data.frame ------------
cli::cli_process_start("data.frame")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat,
          proj.name = 'Current',
          new.env = myExpl.df,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      stopifnot(tinytest::expect_error(plot(myBiomodProj)))
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## data.frame1 ------------
cli::cli_process_start("data.frame1")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat1,
          proj.name = 'Current',
          new.env = myExpl.df1,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      stopifnot(tinytest::expect_error(plot(myBiomodProj)))
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame with coordinates ------------
cli::cli_process_start("data.frame with coord")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat,
          proj.name = 'Current',
          new.env = myExpl.df,
          new.env.xy = crds(myExpl.SpatVector),
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## data.frame1 with coord ------------
cli::cli_process_start("data.frame1 with coord")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_noCat1,
          proj.name = 'Current',
          new.env = myExpl.df1,
          new.env.xy = crds(myExpl.SpatVector),
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# With Categorical Variables ------------------------------------------------

cli::cli_h2("With Categorical Variables")

## SpatRaster ------------
cli::cli_process_start("SpatRaster")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_Cat,
          proj.name = 'Current',
          new.env = myExpl.cat,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## Raster ------------
cli::cli_process_start("Raster")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_Cat,
          proj.name = 'Current',
          new.env = myExpl.cat.raster,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## Matrix ------------
cli::cli_process_start("Matrix")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      thistest <- tinytest::expect_error(
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_Cat,
          proj.name = 'Current',
          new.env = myExpl.matrix,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      )
      stopifnot(thistest)
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame ------------
cli::cli_process_start("data.frame")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_Cat,
          proj.name = 'Current',
          new.env = myExpl.cat.df,
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      stopifnot(tinytest::expect_error(plot(myBiomodProj)))
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## data.frame with coordinates ------------
cli::cli_process_start("data.frame with coord")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodProj <-
        BIOMOD_Projection(
          bm.mod = myBiomodModelOut_Cat,
          proj.name = 'Current',
          new.env = myExpl.cat.df,
          new.env.xy = crds(myExpl.SpatVector),
          metric.binary = 'all',
          metric.filter = 'all',
          build.clamping.mask = FALSE,
          do.stack = TRUE
        )
      pdf(file = ".tmp.pdf", width = 20/cm(1), height = 15/cm(1))
      stopifnot(inherits(
        plot(myBiomodProj), "ggplot"))
      dev.off()
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Projection <- Error_Projection + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

graphics.off()
