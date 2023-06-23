cli::cli_h1("BIOMOD_FormatingData")

Error_Formating <- 0
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
myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )
myResp.spdf <- as(myResp.SpatVector, "Spatial")
myResp.SpatialPoints <- sp::SpatialPoints(myResp.spdf[1:200,])
myResp.SpatVectorPO <- vect(myResp.SpatialPoints)

## myExpl ----------------------------------

myExpl.df <- extract(myExpl, y = myResp.SpatVector, ID = FALSE)
myExpl.df1 <- myExpl.df[,1]
myExpl.matrix <- as.matrix(myExpl.df)
myExpl.matrix1 <- as.matrix(myExpl.df1)
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
                                active = 1) 

myExpl.cat.df <- extract(myExpl.cat, y = myResp.SpatVector, ID = FALSE)
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


## myEvalResp ----------------------------------
myEvalResp <- myResp 
myEvalRespXY <- myRespXY 
myEvalResp.SpatVector <- myResp.SpatVector 
myEvalResp.spdf <- myResp.spdf 
myEvalExpl.raster <- myExpl.raster

## myEvalExpl ----------------------------------

myEvalExpl <- myExpl
myEvalExpl.df <- myExpl.df
myEvalExpl.spdf <- myExpl.spdf
myEvalExpl.SpatVector <- myExpl.SpatVector
myEvalExpl.matrix <- myExpl.matrix
myEvalExpl.raster <- stack(myEvalExpl)

## myEvalExpl.cat ----------------------------------

myEvalExpl.cat <- myExpl.cat
myEvalExpl.cat.df <- myExpl.cat.df
myEvalExpl.cat.spdf <- myExpl.cat.spdf
myEvalExpl.cat.SpatVector <- myExpl.cat.SpatVector
myEvalExpl.cat.raster <- stack(myEvalExpl.cat)

# no eval -----------------------------------------------------------------

cli::cli_h2("No Evaluations")

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
          resp.name = myRespName)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
    }
    
    )
    
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = matrix ------------
cli::cli_process_start("resp.var = vector ; expl.var = matrix")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.matrix,
          resp.xy = myRespXY,
          resp.name = myRespName)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          # resp.xy = myRespXY,
          resp.name = myRespName)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SPDF ------------
cli::cli_process_start("resp.var = vector ; expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### resp.var = vector ; expl.var = SpatVector ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatVector")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          PA.nb.rep = 2,
          PA.strategy = 'random',
          PA.nb.absences = 500)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = SpatRaster ; multiple Pseudo-Absences ------------
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
          PA.strategy = 'random',
          PA.nb.absences = c(1000, 500, 500, 200))
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          PA.nb.rep = 2,
          PA.strategy = 'random',
          PA.nb.absences = 500)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ; multiple Pseudo-Absences ------------
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
          PA.strategy = 'random',
          PA.nb.absences = c(1000, 500, 500, 200))
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector Presence Only----------------------
cli::cli_h3("resp.var = SpatVector Presence-Only")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector Presence-Only ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVectorPO,
          expl.var = myExpl,
          resp.name = myRespName, 
          PA.nb.rep = 2,
          PA.strategy = 'random',
          PA.nb.absences = 500)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector Presence-Only ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVectorPO,
          expl.var = myExpl,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.strategy = 'random',
          PA.nb.absences = c(1000, 500, 500, 200))
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector Presence-Absence")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}




# with eval -----------------------------------------------------------------

cli::cli_h2("With Evaluations")

cli::cli_h3("expl as SpatRaster")

### eval.resp.var = vector ; eval.expl.var = SpatRaster ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.expl.var = myEvalExpl,
          eval.resp.xy = myEvalRespXY)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### eval.resp.var = vector ; eval.expl.var = NULL ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = NULL")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.resp.xy = myEvalRespXY)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### eval.resp.var = SPDF ; eval.expl.var = SPDF ------------
cli::cli_process_start("eval.resp.var = SPDF ; eval.expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp.spdf,
          eval.expl.var = myEvalExpl.spdf)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = Matrix ------------
cli::cli_process_start("eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = Matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp.SpatVector,
          eval.expl.var = myEvalExpl.matrix)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = data.frame ------------
cli::cli_process_start("eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output({
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp.SpatVector,
          eval.expl.var = myEvalExpl.df)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData)
      g <- plot(myBiomodData, plot.type = "raster")
      g <- plot(myBiomodData, plot.type = "points")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "facet")
      g <- plot(myBiomodData, plot.type = "points", plot.output = "list")
      dev.off()
      summary(myBiomodData)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_h3("expl as RasterStack")

### eval.resp.var = vector ; eval.expl.var = SpatRaster ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.expl.var = myEvalExpl.raster,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### eval.resp.var = vector ; eval.expl.var = NULL ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = NULL")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# no eval - categorical variables -----------------------------------------------------------------

cli::cli_h2("No Evaluations")

## resp.var = vector -------------------------------------------------------
cli::cli_h3("resp.var = vector")


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.df,
          # resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.spdf,
          resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.SpatVector,
          resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.raster,
          resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.strategy = 'random',
          PA.nb.absences = 500)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = SpatRaster ; multiple Pseudo-Absences ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.strategy = 'random',
          PA.nb.absences = c(1000, 500, 500, 200))
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.strategy = 'random',
          PA.nb.absences = 500)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatialPoints ; expl.var = raster ; multiple Pseudo-Absences ------------
cli::cli_process_start("resp.var = SpatialPoints ; expl.var = raster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatialPoints,
          expl.var = myExpl.cat.raster,
          resp.name = myRespName,
          PA.nb.rep = 4,
          PA.strategy = 'random',
          PA.nb.absences = c(1000, 500, 500, 200))
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.df,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.raster,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector Presence Only----------------------
cli::cli_h3("resp.var = SpatVector Presence-Only")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector Presence-Only ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVectorPO,
          expl.var = myExpl.cat,
          resp.name = myRespName, 
          PA.nb.rep = 2,
          PA.strategy = 'random',
          PA.nb.absences = 500)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = SpatVector ; expl.var = SpatRaster ; multiple Pseudo-Absences ------------
cli::cli_process_start("resp.var = SpatVector Presence-Only ; expl.var = SpatRaster ; multiple Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVectorPO,
          expl.var = myExpl.cat,
          resp.name = myRespName, 
          PA.nb.rep = 4,
          PA.strategy = 'random',
          PA.nb.absences = c(1000, 500, 500, 200))
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## resp.var = SpatVector ----------------------
cli::cli_h3("resp.var = SpatVector Presence-Absence")

### resp.var = SpatVector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = SpatVector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp.SpatVector,
          expl.var = myExpl.cat,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.df,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.spdf,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.SpatVector,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.cat.raster,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}




# with eval - categorical -----------------------------------------------------------------

cli::cli_h2("With Evaluations")

cli::cli_h3("expl as SpatRaster")

### eval.resp.var = vector ; eval.expl.var = SpatRaster ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.expl.var = myEvalExpl.cat,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### eval.resp.var = vector ; eval.expl.var = NULL ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = NULL")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### eval.resp.var = SPDF ; eval.expl.var = SPDF ------------
cli::cli_process_start("eval.resp.var = SPDF ; eval.expl.var = SPDF")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp.spdf,
          eval.expl.var = myEvalExpl.cat.spdf)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = data.frame ------------
cli::cli_process_start("eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp.SpatVector,
          eval.expl.var = myEvalExpl.cat.df)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_h3("expl as RasterStack")

### eval.resp.var = vector ; eval.expl.var = SpatRaster ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.expl.var = myEvalExpl.cat.raster,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### eval.resp.var = vector ; eval.expl.var = NULL ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = NULL")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat.raster,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# Only one variable -----------------------------------------------------------------

cli::cli_h2("Only one variable ")

## No Evaluations -------------------------------------------------------
cli::cli_h3("No Evaluations")


### resp.var = vector ; expl.var = SpatRaster ------------
cli::cli_process_start("resp.var = vector ; expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl1,
          resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### resp.var = vector ; expl.var = RasterStack ------------
cli::cli_process_start("resp.var = vector ; expl.var = RasterStack")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster1,
          resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.matrix1,
          resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
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
          expl.var = myExpl.df1,
          # resp.xy = myRespXY,
          resp.name = myRespName)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with eval -----------------------------------------------------------------

cli::cli_h3("With Evaluations")


### eval.resp.var = vector ; eval.expl.var = SpatRaster ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = SpatRaster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl1,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.expl.var = myExpl1,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### eval.resp.var = vector ; eval.expl.var = NULL ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = NULL")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl1,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = Matrix ------------
cli::cli_process_start("eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = Matrix")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl1,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp.SpatVector,
          eval.expl.var = myExpl.matrix1)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = data.frame ------------
cli::cli_process_start("eval.resp.var = SpatVector Presence-Absence ; eval.expl.var = data.frame")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl1,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp.SpatVector,
          eval.expl.var = myExpl.df1)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### eval.resp.var = vector ; eval.expl.var = SpatRaster ------------
cli::cli_process_start("eval.resp.var = vector ; eval.expl.var = raster")
this_try <- try({
  invisible(
    capture.output(
      myBiomodData <- 
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.raster1,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myEvalResp,
          eval.expl.var = myExpl.raster1,
          eval.resp.xy = myEvalRespXY)
    )
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
