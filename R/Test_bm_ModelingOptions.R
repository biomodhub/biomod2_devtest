cli::cli_h1("bm_ModelingOptions")

Error_Formating <- 0

# Preparation -----------------------------------------------------------------
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

## myResp ---------------------------------------------------------------------
myResp <- as.numeric(DataSpecies[, myRespName])
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

tmpindex <- which(myResp == 1)
myResp_PO <- myResp[tmpindex]
myRespXY_PO <- myRespXY[tmpindex,]


## myBiomodData ---------------------------------------------------------------
myBiomodData <-
  BIOMOD_FormatingData(
    resp.var = myResp,
    expl.var = myExpl,
    resp.xy = myRespXY,
    resp.name = myRespName)

myBiomodDataPA <-
  BIOMOD_FormatingData(
    resp.var = myResp_PO,
    expl.var = myExpl,
    resp.xy = myRespXY_PO,
    resp.name = myRespName,
    PA.nb.rep = 2,
    PA.nb.absences = 500,
    PA.strategy = 'random')

## myBiomodCV -----------------------------------------------------------------
myBiomodCV <- bm_CrossValidation(bm.format = myBiomodData,
                                 strategy = "random",
                                 nb.rep = 3,
                                 perc = 0.6,
                                 do.full.models = TRUE)

myBiomodCVPA <- bm_CrossValidation(bm.format = myBiomodDataPA,
                                   strategy = "random",
                                   nb.rep = 3,
                                   perc = 0.6,
                                   do.full.models = TRUE)



# strategy = default ----------------------------------------------------------
cli::cli_h2("strategy = default")


cli::cli_process_start("No val.list ; No bm.format ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , val.list = NULL, bm.format = NULL, calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No val.list ; No bm.format ; calib.lines = Presence-Absence")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , val.list = NULL, bm.format = NULL
                           , calib.lines = myBiomodCV)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## NEED CORRECTION IN expected_CVnames
# cli::cli_process_start("No val.list ; No bm.format ; calib.lines = Presence-Only")
# this_try <- try({
#   invisible(
#     capture.output({
#       myBiomodOptions <-
#         bm_ModelingOptions(data.type = "binary"
#                            , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
#                                         , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
#                            , strategy = "default"
#                            , val.list = NULL, bm.format = NULL
#                            , calib.lines = myBiomodCVPA)
#       myBiomodOptions
#       names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
#     })
#   )
# }, silent = TRUE)
# 
# if(inherits(this_try, "try-error")){
#   Error_Formating <- Error_Formating + 1
#   cli::cli_process_failed()
# } else {
#   cli::cli_process_done()
# }


## with bm.format -------------------------------
cli::cli_h3("with bm.format")

cli::cli_process_start("No val.list ; bm.format = Presence-Absence ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , val.list = NULL
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No val.list ; bm.format = Presence-Only ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , val.list = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with bm.format ; with calib.lines ------------
cli::cli_h3("with bm.format")

cli::cli_process_start("No val.list ; bm.format = Presence-Absence ; calib.lines = Presence-Absence")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , val.list = NULL
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = myBiomodCV) ## does not work with myBiomodCVPA !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No val.list ; bm.format = Presence-Only ; calib.lines = Presence-Only")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , val.list = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = myBiomodCVPA) ## works with myBiomodCV !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# strategy = bigboss ----------------------------------------------------------
cli::cli_h2("strategy = bigboss")


cli::cli_process_start("No val.list ; No bm.format ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , val.list = NULL, bm.format = NULL, calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with bm.format -------------------------------
cli::cli_h3("with bm.format")


cli::cli_process_start("No val.list ; bm.format = Presence-Absence ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , val.list = NULL
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No val.list ; bm.format = Presence-Only ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , val.list = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# strategy = user.defined -------------------------------------------
cli::cli_h2("strategy = user.defined")

cli::cli_process_start("No bm.format ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "user.defined"
                           , val.list = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
                                             , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
                           , bm.format = NULL, calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with bm.format -------------------------------
cli::cli_h3("with bm.format")

cli::cli_process_start("No val.list ; bm.format = Presence-Absence ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "user.defined"
                           , val.list = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
                                             , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No val.list ; bm.format = Presence-Only ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "user.defined"
                           , val.list = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
                                             , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_allData_allRun`
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# strategy = tuned --------------------------------------------------
cli::cli_h2("strategy = tuned")

## with bm.format -------------------------------
cli::cli_h3("with bm.format")


cli::cli_process_start("No val.list ; bm.format = Presence-Absence ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('SRE')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , val.list = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No val.list ; bm.format = Presence-Only ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('SRE')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , val.list = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Formating <- Error_Formating + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

