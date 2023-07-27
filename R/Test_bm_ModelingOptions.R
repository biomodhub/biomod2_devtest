cli::cli_h1("bm_ModelingOptions")

Error_ModelingOptions <- 0

# Preparation -----------------------------------------------------------------
library(terra)
library(raster)
if (terraVersion) {
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
set.seed(42)
mysubsample <- sample(1:2488, size = 150)
myBiomodData <-
  BIOMOD_FormatingData(
    resp.var = myResp[mysubsample],
    expl.var = myExpl,
    resp.xy = myRespXY[mysubsample,],
    resp.name = myRespName)

myBiomodDataPA <-
  BIOMOD_FormatingData(
    resp.var = myResp_PO[1:100],
    expl.var = myExpl,
    resp.xy = myRespXY_PO[1:100,],
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


cli::cli_process_start("No user.val ; No bm.format ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , user.val = NULL, bm.format = NULL, calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# should be modified with expect_error
# ## with calib.lines -----------------------------
# cli::cli_h3("with calib.lines")
# 
# cli::cli_process_start("No user.val ; No bm.format ; calib.lines = Presence-Absence")
# this_try <- try({
#   invisible(
#     capture.output({
#       myBiomodOptions <-
#         bm_ModelingOptions(data.type = "binary"
#                            , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
#                                         , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
#                            , strategy = "default"
#                            , user.val = NULL, bm.format = NULL
#                            , calib.lines = myBiomodCV)
#       myBiomodOptions
#       names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
#     })
#   )
# }, silent = TRUE)
# 
# if (inherits(this_try, "try-error")) {
#   Error_ModelingOptions <- Error_ModelingOptions + 1
#   cli::cli_process_failed()
# } else {
#   cli::cli_process_done()
# }


## with bm.format -------------------------------
cli::cli_h3("with bm.format")

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , user.val = NULL
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with bm.format ; with calib.lines ------------
cli::cli_h3("with bm.format")

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; calib.lines = Presence-Absence")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , user.val = NULL
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = myBiomodCV) ## does not work with myBiomodCVPA !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No user.val ; bm.format = Presence-Only ; calib.lines = Presence-Only")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "default"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = myBiomodCVPA) ## works with myBiomodCV !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# strategy = bigboss ----------------------------------------------------------
cli::cli_h2("strategy = bigboss")


cli::cli_process_start("No user.val ; No bm.format ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , user.val = NULL, bm.format = NULL, calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# # should be modified with expect_error
# 
# ## with calib.lines -----------------------------
# cli::cli_h3("with calib.lines")
# 
# 
# cli::cli_process_start("No user.val ; No bm.format ; calib.lines = Presence-Absence")
# this_try <- try({
#   invisible(
#     capture.output({
#       myBiomodOptions <-
#         bm_ModelingOptions(data.type = "binary"
#                            , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
#                                         , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
#                            , strategy = "bigboss"
#                            , user.val = NULL, bm.format = NULL
#                            , calib.lines = myBiomodCV)
#       myBiomodOptions
#       names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
#     })
#   )
# }, silent = TRUE)
# 
# if (inherits(this_try, "try-error")) {
#   Error_ModelingOptions <- Error_ModelingOptions + 1
#   cli::cli_process_failed()
# } else {
#   cli::cli_process_done()
# }

## with bm.format -------------------------------
cli::cli_h3("with bm.format")


cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , user.val = NULL
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with bm.format ; with calib.lines ------------
cli::cli_h3("with bm.format")


cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; calib.lines = Presence-Absence")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , user.val = NULL
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = myBiomodCV) ## does not work with myBiomodCVPA !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No user.val ; bm.format = Presence-Only ; calib.lines = Presence-Only")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "bigboss"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = myBiomodCVPA) ## works with myBiomodCV !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
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
                           , user.val = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
                                             , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
                           , bm.format = NULL, calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# # should be modified with expect_error
# ## with calib.lines -----------------------------
# cli::cli_h3("with calib.lines")
# 
# 
# cli::cli_process_start("No user.val ; No bm.format ; calib.lines = Presence-Absence")
# this_try <- try({
#   invisible(
#     capture.output({
#       myBiomodOptions <-
#         bm_ModelingOptions(data.type = "binary"
#                            , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
#                                         , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
#                            , strategy = "user.defined"
#                            , user.val = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
#                                              , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
#                            , bm.format = NULL
#                            , calib.lines = myBiomodCV) ## does not work with myBiomodCVPA !
#       myBiomodOptions
#       names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
#       myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_allData_allRun`
#       myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_allData_RUN1`
#     })
#   )
# }, silent = TRUE)
# 
# if (inherits(this_try, "try-error")) {
#   Error_ModelingOptions <- Error_ModelingOptions + 1
#   cli::cli_process_failed()
# } else {
#   cli::cli_process_done()
# }


## with bm.format -------------------------------
cli::cli_h3("with bm.format")

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "user.defined"
                           , user.val = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
                                             , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = NULL)
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "user.defined"
                           , user.val = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
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

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with bm.format ; with calib.lines ------------
cli::cli_h3("with bm.format")


cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; calib.lines = Presence-Absence")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "user.defined"
                           , user.val = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
                                             , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
                           , bm.format = myBiomodData ## useless here
                           , calib.lines = myBiomodCV) ## does not work with myBiomodCVPA !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_allData_allRun`
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_allData_RUN1`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("No user.val ; bm.format = Presence-Only ; calib.lines = Presence-Only")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                        , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "user.defined"
                           , user.val = list(SRE.binary.biomod2.bm_SRE = list('_allData_allRun' = list(quant = 0.01))
                                             , XGBOOST.binary.xgboost.xgboost = list('_allData_allRun' = list(nrounds = 10)))
                           , bm.format = myBiomodDataPA
                           , calib.lines = myBiomodCVPA) ## works with myBiomodCV !
      myBiomodOptions
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_allData_allRun`
      myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# strategy = tuned ; with bm.format--------------------------------------------------
cli::cli_h2("strategy = tuned - with bm.format")

## Presence-Absence -------------------------------
cli::cli_h3("Presence-Absence")

### ANN -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = ANN")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$ANN.binary.nnet.nnet@args.values
      names(myBiomodOptions@options$ANN.binary.nnet.nnet@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### CTA -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = CTA")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('CTA')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$CTA.binary.rpart.rpart@args.values
      names(myBiomodOptions@options$CTA.binary.rpart.rpart@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### FDA -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = FDA")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('FDA')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$FDA.binary.mda.fda@args.values
      names(myBiomodOptions@options$FDA.binary.mda.fda@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### GAM-mgcv-gam -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = GAM-mgcv-gam")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GAM.mgcv.gam')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$GAM.binary.mgcv.gam@args.values
      names(myBiomodOptions@options$GAM.binary.mgcv.gam@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### GAM-mgcv-bam -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = GAM.mgcv.bam")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GAM.mgcv.bam')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$GAM.binary.mgcv.bam@args.values
      names(myBiomodOptions@options$GAM.binary.mgcv.bam@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
### GAM-gam -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = GAM")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GAM.gam.gam')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
        myBiomodOptions@options$GAM.binary.gam.gam@args.values
      names(myBiomodOptions@options$GAM.binary.gam.gam@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### GBM -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = GBM")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GBM')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$GBM.binary.gbm.gbm@args.values
      names(myBiomodOptions@options$GBM.binary.gbm.gbm@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### GLM -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = GLM")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GLM')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$GLM.binary.stats.glm@args.values
      names(myBiomodOptions@options$GLM.binary.stats.glm@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### MARS -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = MARS")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('MARS')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      myBiomodOptions@options$MARS.binary.earth.earth@args.values
      names(myBiomodOptions@options$MARS.binary.earth.earth@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### MAXENT -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = MAXENT")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('MAXENT')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      # myBiomodOptions@options$FDA.binary.mda.fda@args.values
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### MAXNET -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = MAXNET")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('MAXNET')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      # myBiomodOptions@options$FDA.binary.mda.fda@args.values
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### RF -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = RF")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('RF')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      # myBiomodOptions@options$FDA.binary.mda.fda@args.values
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### SRE -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = SRE")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('SRE')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      # myBiomodOptions@options$FDA.binary.mda.fda@args.values
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### XGBOOST -------------------------------

cli::cli_process_start("No user.val ; bm.format = Presence-Absence ; No calib.lines ; model = XGBOOST")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('XGBOOST')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodData
                           , calib.lines = NULL)
      myBiomodOptions
      # myBiomodOptions@options$FDA.binary.mda.fda@args.values
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## Presence-Only -------------------------------

cli::cli_h3("Presence-Only")


### ANN ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = ANN")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('ANN')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### CTA ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = CTA")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('CTA')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### FDA ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = FDA")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('FDA')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### GAM ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = GAM")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GAM')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### GBM ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = GBM")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GBM')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### GLM ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = GLM")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('GLM')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### MARS ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = MARS")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('MARS')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### MAXENT ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = MAXENT")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('MAXENT')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### MAXNET ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = MAXNET")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('MAXNET')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### RF ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = RF")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('RF')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### SRE ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = SRE")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('SRE')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



### XGBOOST ---------------------------------------------------------------
cli::cli_process_start("No user.val ; bm.format = Presence-Only ; No calib.lines ; model = XGBOOST")
this_try <- try({
  invisible(
    capture.output({
      myBiomodOptions <-
        bm_ModelingOptions(data.type = "binary"
                           , models = c('XGBOOST')
                           # 'ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , strategy = "tuned"
                           , user.val = NULL
                           , bm.format = myBiomodDataPA
                           , calib.lines = NULL)
      myBiomodOptions
      # names(myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values)
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA1_allRun`
      # myBiomodOptions@options$SRE.binary.biomod2.bm_SRE@args.values$`_PA2_allRun`
    })
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_ModelingOptions <- Error_ModelingOptions + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


