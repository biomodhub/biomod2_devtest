cli::cli_h1("bm_CrossValidation")


Error_CV <- 0

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
myResp01 <- myResp
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


## BIOMOD_FormatingData --------------------------------------------
x <- try({
  capture.output({
    myBiomodData <- 
      BIOMOD_FormatingData(
        resp.var = myResp01,
        expl.var = myExpl,
        resp.xy = myRespXY,
        resp.name = myRespName)
    
    myBiomodDataPA <- 
      BIOMOD_FormatingData(
        resp.var = myResp,
        expl.var = myExpl,
        resp.xy = myRespXY,
        resp.name = myRespName, 
        PA.nb.rep = 3,
        PA.nb.absences = 300,
        PA.strategy = "random")
  })
})

# Without Pseudo-Absences ----------------------------------------------------

cli::cli_h2("Without Pseudo-Absences")


## random ------------
cli::cli_h3("strategy = 'random'")
cli::cli_process_start("strategy = 'random'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "random",
                                        nb.rep = 10,
                                        perc = 0.8,
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## strategy = 'kfold' ------------------------------------------
cli::cli_h3("strategy = 'kfold'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "kfold",
                                        k = 5, 
                                        nb.rep = 2,
                                        balance = "presences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "kfold",
                                        k = 5, 
                                        nb.rep = 2,
                                        balance = "absences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## strategy = 'env' ------------------------------------------
cli::cli_h3("strategy = 'env'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "env",
                                        k = 5, 
                                        balance = "presences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "env",
                                        k = 5, 
                                        env.var = "bio3",
                                        balance = "absences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}




## strategy = 'strat', strat = 'x' -------------------------------
cli::cli_h3("strategy = 'strat', strat = 'x'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "strat",
                                        k = 5, 
                                        balance = "presences",
                                        strat = "x",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "strat",
                                        nb.rep = 2,
                                        k = 2, 
                                        balance = "absences",
                                        strat = "x",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'strat', strat = 'y' -------------------------------
cli::cli_h3("strategy = 'strat', strat = 'y'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "strat",
                                        k = 2, 
                                        balance = "presences",
                                        strat = "y",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "strat",
                                        k = 2, 
                                        balance = "absences",
                                        strat = "y",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'strat', strat = 'both' -------------------------------
cli::cli_h3("strategy = 'strat', strat = 'both'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "strat",
                                        k = 4, 
                                        balance = "presences",
                                        strat = "both",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "strat",
                                        k = 2, 
                                        balance = "absences",
                                        strat = "both",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'block' -------------------------------
cli::cli_h3("strategy = 'block'")

this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "block",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'user.defined' ------------
cli::cli_h3("strategy = 'user.defined'")
cli::cli_process_start("strategy = 'user.defined'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines.prep <- bm_CrossValidation(myBiomodData,
                                             strategy = "env",
                                             k = 2, 
                                             # env.var = "bio3",
                                             balance = "absences",
                                             do.full.models = FALSE)
      calib.lines <- bm_CrossValidation(myBiomodData,
                                        strategy = "user.defined",
                                        user.table = calib.lines.prep)
      
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines)
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodData, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      
      dev.off()
      summary(myBiomodData, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# With Pseudo-Absences ----------------------------------------------------

cli::cli_h2("With Pseudo-Absences")

## strategy = random ------------
cli::cli_h3("strategy = 'random'")
cli::cli_process_start("strategy = 'random'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "random",
                                        nb.rep = 2,
                                        perc = 0.8,
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'kfold' ------------------------------------------
cli::cli_h3("strategy = 'kfold'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "kfold",
                                        k = 5, 
                                        nb.rep = 2,
                                        balance = "presences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'kfold' ------------------------------------------
cli::cli_h3("strategy = 'kfold'")


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "kfold",
                                        k = 5, 
                                        nb.rep = 2,
                                        balance = "absences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## strategy = 'env' ------------------------------------------
cli::cli_h3("strategy = 'env'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "env",
                                        k = 2, 
                                        env.var = "bio3",
                                        balance = "presences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "env",
                                        k = 2, 
                                        balance = "absences",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}




## strategy = 'strat', strat = 'x' -------------------------------
cli::cli_h3("strategy = 'strat', strat = 'x'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "strat",
                                        k = 5, 
                                        balance = "presences",
                                        strat = "x",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "strat",
                                        nb.rep = 2,
                                        k = 2, 
                                        balance = "absences",
                                        strat = "x",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'strat', strat = 'y' -------------------------------
cli::cli_h3("strategy = 'strat', strat = 'y'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "strat",
                                        nb.rep = 2,
                                        k = 2, 
                                        balance = "presences",
                                        strat = "y",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "strat",
                                        nb.rep = 2,
                                        k = 2, 
                                        balance = "absences",
                                        strat = "x",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'strat', strat = 'both' -------------------------------
cli::cli_h3("strategy = 'strat', strat = 'both'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presence'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "strat",
                                        k = 4, 
                                        balance = "presences",
                                        strat = "both",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### balance = 'absences' ------------
cli::cli_process_start("balance = 'absences'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "strat",
                                        k = 4, 
                                        balance = "absences",
                                        strat = "both",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## strategy = 'block' -------------------------------
cli::cli_h3("strategy = 'block'")


### balance = 'presence' ------------
cli::cli_process_start("balance = 'presences'")
this_try <- try({
  invisible(
    capture.output({
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "block",
                                        do.full.models = FALSE)
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



while(!is.null(dev.list())){
  dev.off()
}


cli::cli_h3("strategy = 'user.defined'")

## strategy = 'user.defined' ------------
cli::cli_process_start("strategy = 'user.defined'")
this_try <- try({
  invisible(
    capture.output({
      
      calib.lines.prep <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "env",
                                        k = 2, 
                                        # env.var = "bio3",
                                        balance = "absences",
                                        do.full.models = FALSE)
      calib.lines <- bm_CrossValidation(myBiomodDataPA,
                                        strategy = "user.defined",
                                        user.table = calib.lines.prep)
      
      pdf(file = ".tmp.pdf")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines)
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster")
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                run = c("RUN1","RUN2"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                run = c("RUN1","RUN2"),
                PA = c("PA1"))
      g <- plot(myBiomodDataPA, 
                calib.lines = calib.lines, 
                plot.type = "raster",
                plot.output = "list",
                PA = c("PA1","PA2"))
      
      dev.off()
      summary(myBiomodDataPA, calib.lines = calib.lines)
      
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_CV <- Error_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


