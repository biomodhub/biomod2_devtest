cli::cli_h1("BIOMOD_Modeling")

Error_Modeling_CV <- 0

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

myExpl1 <- myExpl[[1]]

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


## myExpl.df ---------------------------------------------------------------
myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )
myExpl.df <- extract(myExpl, y = myResp.SpatVector, ID = FALSE)


# cross-validation ------------------------------------------------------------------
cli::cli_h2("Crossvalidation tests")

cli::cli_h3("Generate cross-validation table")

## Without stratification --------------------------------------------------
cli::cli_process_start("No stratification")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "random",
                                       nb.rep = 3,
                                       perc = 0.7)
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "kfold",
                                       nb.rep = 2,
                                       k = 3)
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "block")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling_CV <- Error_Modeling_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## With stratification - method = 'x' ------------------------------------
cli::cli_process_start("With stratification - method = 'x'")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "strat",
                                       k = 2,
                                       balance = "presences",
                                       strat = "x")
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "strat",
                                       k = 2,
                                       balance = "absences",
                                       strat = "x")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling_CV <- Error_Modeling_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## With stratification - method = 'y' ------------------------------------
cli::cli_process_start("With stratification - method = 'y'")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "strat",
                                       k = 2,
                                       balance = "presences",
                                       strat = "y")
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "strat",
                                       k = 2,
                                       balance = "absences",
                                       strat = "y")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling_CV <- Error_Modeling_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## With stratification - method = 'both' ------------------------------------
cli::cli_process_start("With stratification - method = 'both'")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "strat",
                                       k = 2,
                                       balance = "presences",
                                       strat = "both")
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "strat",
                                       k = 2,
                                       balance = "absences",
                                       strat = "both")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling_CV <- Error_Modeling_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## With stratification - method = 'block' ------------------------------------
cli::cli_process_start("With stratification - method = 'block'")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "block")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling_CV <- Error_Modeling_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("With stratification and expl.var as data.frame")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.df,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl.df,
          eval.resp.xy = myRespXY)
      
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "env",
                                       k = 2,
                                       balance = "presences")
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "env",
                                       k = 2,
                                       balance = "absences")
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling_CV <- Error_Modeling_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

cli::cli_h3("Modeling with Cross-validation")
## BIOMOD_Modeling with Cross-validation  - data.frame ---------------------
cli::cli_process_start("Modeling with Cross-validation - data.frame")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodCV <- bm_CrossValidation(myBiomodData,
                                       strategy = "kfold",
                                       nb.rep = 2,
                                       k = 3)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          models = c("GLM","SRE","RF"),
          modeling.id = 'CV_matrix',
          CV.strategy = 'user.defined',
          CV.user.table = as.data.frame(myBiomodCV),
          OPT.strategy = 'default',
          var.import = 2,
          metric.eval = c('TSS','ROC'),
          seed.val = 42
        )
      get_predictions(myBiomodModelOut)
      get_evaluations(myBiomodModelOut)
      get_built_models(myBiomodModelOut)
      get_formal_data(myBiomodModelOut)
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling_CV <- Error_Modeling_CV + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
