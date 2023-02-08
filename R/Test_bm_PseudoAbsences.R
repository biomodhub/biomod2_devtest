cli::cli_h1("bm_PseudoAbsences")


Error_PseudoAbsences <- 0

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

## myResp ----------------------------------
myResp <- as.numeric(DataSpecies[, myRespName])
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
myResp[which(myResp == 0)] <- NA
myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )
myResp.spdf <- as(myResp.SpatVector, "Spatial")
myResp.SpatialPoints <- sp::SpatialPoints(myResp.spdf[1:200,])
myResp.SpatVectorPO <- vect(myResp.SpatialPoints)

## myExpl ----------------------------------

myExpl.df <- extract(myExpl, y = myResp.SpatVector, ID = FALSE)
myExpl.matrix <- as.matrix(myExpl.df)
myExpl.raster <- stack(myExpl)
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
                                active = 1) 

myExpl.cat.df <- extract(myExpl.cat, y = myResp.SpatVector, ID =FALSE)
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
# c(rast(myExpl.cat.raster[[1]]), rast(myExpl.cat.raster[[2]]))
tmpdf <- cbind(myRespXY,myExpl.cat.df)
rownames(tmpdf) <- NULL
myExpl.cat.SpatVector <- vect(
  tmpdf,
  geom = c("X_WGS84", "Y_WGS84")
)
myExpl.cat.spdf <- as(myExpl.cat.SpatVector, "Spatial")


## myEvalExpl.cat ----------------------------------

myEvalExpl.cat <- myExpl.cat
myEvalExpl.cat.df <- myExpl.cat.df
myEvalExpl.cat.spdf <- myExpl.cat.spdf
myEvalExpl.cat.SpatVector <- myExpl.cat.SpatVector
myEvalExpl.cat.raster <- stack(myEvalExpl.cat)


# PA.strategy = random ----------------------------------------------------

cli::cli_h2("PA.strategy = random")

## resp.var = vector -------------------------------------------------------
cli::cli_h3("resp.var = vector")


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random")
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "raster", PA = c("PA1"))
      g <- plot(myBiomodData, plot.type = "raster", do.plot = FALSE)
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      g <- plot(myBiomodData, plot.type = "raster", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ; multiple Pseudo-Absences ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random")
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "raster", PA = c("PA1"))
      g <- plot(myBiomodData, plot.type = "raster", do.plot = FALSE)
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      g <- plot(myBiomodData, plot.type = "raster", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = vector ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.matrix,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = vector ; expl.var = matrix : multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.matrix,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "raster", PA = c("PA1"))
      g <- plot(myBiomodData, plot.type = "raster", do.plot = FALSE)
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      g <- plot(myBiomodData, plot.type = "raster", plot.output = "list")
      dev.off()
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "raster", PA = c("PA1"))
      g <- plot(myBiomodData, plot.type = "raster", do.plot = FALSE)
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      g <- plot(myBiomodData, plot.type = "raster", plot.output = "list")
      dev.off()
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## resp.var = SpatialPoints ----------------------------------------------------
cli::cli_h3("resp.var = SpatialPoints")

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SPDF ----------------------
cli::cli_h3("resp.var = SPDF")

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



# PA.strategy = SRE ----------------------------------------------------

cli::cli_h2("PA.strategy = SRE")

## resp.var = vector -------------------------------------------------------
cli::cli_h3("resp.var = vector")

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = NULL")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre")
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = NULL ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre")
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre",
          PA.sre.quant = 0
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0 ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre",
          PA.sre.quant = 0
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0.1")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre",
          PA.sre.quant = 0.1
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0.1 ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre",
          PA.sre.quant = 0.1
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0.25")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre",
          PA.sre.quant = 0.25
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = vector ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.matrix,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = vector ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.matrix,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## resp.var = SpatialPoints ----------------------------------------------------
cli::cli_h3("resp.var = SpatialPoints")

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SPDF ----------------------
cli::cli_h3("resp.var = SPDF")

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster : multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "sre"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# PA.strategy = disk ----------------------------------------------------

cli::cli_h2("PA.strategy = disk")

## resp.var = vector -------------------------------------------------------
cli::cli_h3("resp.var = vector")


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatRaster ; dist.max = 50")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk",
          PA.dist.max = 50
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatRaster ; multiple Pseudo-Absences ; dist.max = 50")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk",
          PA.dist.max = 50
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatRaster ; multiple Pseudo-Absences ; dist.max = 10")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk",
          PA.dist.max = 10
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = vector ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.matrix,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = vector ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.matrix,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ; expl.var = raster ; multiple Pseudo-absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## resp.var = SpatialPoints ----------------------------------------------------
cli::cli_h3("resp.var = SpatialPoints")

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SPDF ----------------------
cli::cli_h3("resp.var = SPDF")

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = matrix ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.matrix,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = data.frame ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.df,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SPDF ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.spdf,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatVector ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.nb.absences = c(1000, 500, 500, 200),
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# Categorical + PA.strategy = random  ----------------------------

cli::cli_h2("Categorical + PA.strategy = random")

## resp.var = vector -------------------------------------------------------
cli::cli_h3("resp.var = vector")


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random")
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## resp.var = SpatialPoints ----------------------------------------------------
cli::cli_h3("resp.var = SpatialPoints")

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SPDF ----------------------
cli::cli_h3("resp.var = SPDF")

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "random"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# Categorical + PA.strategy = SRE ------------------------
# SRE should always fail with categorical variables
# 
cli::cli_h2("Categorical + PA.strategy = SRE")

## resp.var = vector -------------------------------------------------------
cli::cli_h3("resp.var = vector")

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = NULL")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre")
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre",
          PA.sre.quant = 0
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0.1")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre",
          PA.sre.quant = 0.1
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("Classic ; PA.sre.quant = 0.25")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre",
          PA.sre.quant = 0.25
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = data.frame")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SPDF")
this_try <- 
  tinytest::expect_error(
      capture.output(
        myBiomodData <- 
          BIOMOD_FormatingData(
            resp.var = myResp,
            expl.var = myExpl.cat.spdf,
            resp.xy = myRespXY,
            resp.name = myRespName,
            PA.nb.rep = 3,
            PA.nb.absences = 300,
            PA.strategy = "sre"
          )
      )
    )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SpatVector")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = raster")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## resp.var = SpatialPoints ----------------------------------------------------
cli::cli_h3("resp.var = SpatialPoints")

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ;
                       expl.var = SpatRaster")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ;
                       expl.var = raster")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SPDF ----------------------
cli::cli_h3("resp.var = SPDF")

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SpatRaster")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = data.frame")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SPDF")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SpatVector")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = raster")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SpatRaster")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = data.frame")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SPDF")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SpatVector")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = raster")
this_try <- 
  tinytest::expect_error(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "sre"
        )
    )
  )

if(!this_try){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# Categorical + PA.strategy = disk ---------------------------------------

cli::cli_h2("Categorical + PA.strategy = disk")

## resp.var = vector -------------------------------------------------------
cli::cli_h3("resp.var = vector")


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName, 
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = raster ------------
cli::cli_process_start("resp.var = vector ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## resp.var = SpatialPoints ----------------------------------------------------
cli::cli_h3("resp.var = SpatialPoints")

### resp.var = SpatialPoints ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatialPoints ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatialPoints ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatialPoints ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SPDF ----------------------
cli::cli_h3("resp.var = SPDF")

### resp.var = SPDF ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### resp.var = SPDF ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SPDF ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SPDF ; expl.var = raster ------------
cli::cli_process_start("resp.var = SPDF ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.spdf,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### resp.var = SpatVector ; expl.var = data.frame ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.df,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = SpatVector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = raster ------------
cli::cli_process_start("resp.var = SpatVector ;
                       expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 3,
          PA.nb.absences = 300,
          PA.strategy = "disk"
        )
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}




# User-defined Pseudo Absences --------------------------------------------
cli::cli_h2("User-defined Pseudo Absences")


set.seed(42)
myRespTF <- as.numeric(DataSpecies[, myRespName])

myPAtable <- data.frame(PA1 = ifelse(myRespTF == 1, TRUE, FALSE),
                        PA2 = ifelse(myRespTF == 1, TRUE, FALSE))
for (i in 1:ncol(myPAtable)){
   myPAtable[sample(which(myPAtable[, i] == FALSE), 500), i] <-  TRUE
}

# expl.var = SpatRaster --------------------------------------------
cli::cli_process_start("expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData.u <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl,
                                             resp.xy = myRespXY,
                                             resp.name = myRespName,
                                             PA.strategy = 'user.defined',
                                             PA.user.table = myPAtable)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# expl.var = SpatVector --------------------------------------------
cli::cli_process_start("expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData.u <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl.SpatVector,
                                             resp.xy = myRespXY,
                                             resp.name = myRespName,
                                             PA.strategy = 'user.defined',
                                             PA.user.table = myPAtable)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# expl.var = Categorical SpatRaster ---------------------------------
cli::cli_process_start("expl.var = Categorical SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData.u <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl.cat,
                                             resp.xy = myRespXY,
                                             resp.name = myRespName,
                                             PA.strategy = 'user.defined',
                                             PA.user.table = myPAtable)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# expl.var = Categorical SpatVector ------------------------------------
cli::cli_process_start("expl.var = Categorical SpatVector")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData.u <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl.cat.SpatVector,
                                             resp.xy = myRespXY,
                                             resp.name = myRespName,
                                             PA.strategy = 'user.defined',
                                             PA.user.table = myPAtable)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# Test for wrong parameters ----------------------------------------------------

cli::cli_h2("No absence and No pseudo-absences")

cli::cli_process_start("No absence and No PA")
this_try <- try({
  invisible(
    capture.output({
      this_test <- tinytest::expect_error({
        # should fail as there are no absences nor pseudo-absences.  
        myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp[which(myResp == 1)],
          expl.var = stack(myExpl),
          resp.xy = myRespXY[which(myResp == 1),],
          resp.name = myRespName)
      })
      stopifnot(this_test)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PseudoAbsences <- Error_PseudoAbsences + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
