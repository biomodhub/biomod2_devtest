cli::cli_h1("BIOMOD_Modeling")

Error_Modeling <- 0

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
                                active = 2) 


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


# No Categorical Variables ------------------------------------------------

cli::cli_h2("No Categorical Variables")

## No Evaluation -------------------------------------------------------
cli::cli_h3("No Evaluation")

### Presence-Absence ------------
cli::cli_process_start("Presence-Absence")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
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
          modeling.id = 'NoCat_NoEval_Presence-Absence',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only ------------
cli::cli_process_start("Presence-Only")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO,
          expl.var = myExpl,
          resp.xy = myRespXY_PO,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_NoEval_Presence-Only',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only with NA ------------
cli::cli_process_start("Presence-Only with NA")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO_NA,
          expl.var = myExpl,
          resp.xy = myRespXY_PO_NA,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_NoEval_Presence-Only_with_NA',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### Presence-Absence and Pseudo-Absences ------------
cli::cli_process_start("Presence-Absence and Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PA2,
          expl.var = myExpl,
          resp.xy = myRespXY_PA2,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 1600,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_NoEval_PA2',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## With Evaluation -------------------------------------------------------
cli::cli_h3("With Evaluation")

### Presence-Absence ------------
cli::cli_process_start("Presence-Absence")
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
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_Eval_Presence-Absence',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only ------------
cli::cli_process_start("Presence-Only")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO,
          expl.var = myExpl,
          resp.xy = myRespXY_PO,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_Eval_Presence-Only',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only with NA ------------
cli::cli_process_start("Presence-Only with NA")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO_NA,
          expl.var = myExpl,
          resp.xy = myRespXY_PO_NA,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_Eval_Presence-Only_with_NA',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### Presence-Absence and Pseudo-Absences ------------
cli::cli_process_start("Presence-Absence and Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PA2,
          expl.var = myExpl,
          resp.xy = myRespXY_PA2,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 1600,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat_Eval_PA2',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# With Categorical Variables ------------------------------------------------

cli::cli_h2("With Categorical Variables")

## No Evaluation -------------------------------------------------------
cli::cli_h3("No Evaluation")

### Presence-Absence ------------
cli::cli_process_start("Presence-Absence")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName)
      
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
      get_predictions(myBiomodModelOut)
      get_evaluations(myBiomodModelOut)
      get_built_models(myBiomodModelOut)
      get_formal_data(myBiomodModelOut)
      
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only ------------
cli::cli_process_start("Presence-Only")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO,
          expl.var = myExpl.cat,
          resp.xy = myRespXY_PO,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_NoEval_Presence-Only',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only with NA ------------
cli::cli_process_start("Presence-Only with NA")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO_NA,
          expl.var = myExpl.cat,
          resp.xy = myRespXY_PO_NA,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_NoEval_Presence-Only_with_NA',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### Presence-Absence and Pseudo-Absences ------------
cli::cli_process_start("Presence-Absence and Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PA2,
          expl.var = myExpl.cat,
          resp.xy = myRespXY_PA2,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 1600,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_NoEval_PA2',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## With Evaluation -------------------------------------------------------
cli::cli_h3("With Evaluation")

### Presence-Absence ------------
cli::cli_process_start("Presence-Absence")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl.cat,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl.cat,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_Eval_Presence-Absence',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only ------------
cli::cli_process_start("Presence-Only")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO,
          expl.var = myExpl.cat,
          resp.xy = myRespXY_PO,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl.cat,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_Eval_Presence-Only',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only with NA ------------
cli::cli_process_start("Presence-Only with NA")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO_NA,
          expl.var = myExpl.cat,
          resp.xy = myRespXY_PO_NA,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl.cat,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_Eval_Presence-Only_with_NA',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### Presence-Absence and Pseudo-Absences ------------
cli::cli_process_start("Presence-Absence and Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PA2,
          expl.var = myExpl.cat,
          resp.xy = myRespXY_PA2,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 1600,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl.cat,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Cat_Eval_PA2',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# One non-categorical Variables ------------------------------------------------

cli::cli_h2("One non-categorical Variables")

## No Evaluation -------------------------------------------------------
cli::cli_h3("No Evaluation")

### Presence-Absence ------------
cli::cli_process_start("Presence-Absence")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl1,
          resp.xy = myRespXY,
          resp.name = myRespName)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_NoEval_Presence-Absence',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only ------------
cli::cli_process_start("Presence-Only")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO,
          expl.var = myExpl1,
          resp.xy = myRespXY_PO,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_NoEval_Presence-Only',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only with NA ------------
cli::cli_process_start("Presence-Only with NA")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO_NA,
          expl.var = myExpl1,
          resp.xy = myRespXY_PO_NA,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_NoEval_Presence-Only_with_NA',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### Presence-Absence and Pseudo-Absences ------------
cli::cli_process_start("Presence-Absence and Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PA2,
          expl.var = myExpl1,
          resp.xy = myRespXY_PA2,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 1600,
          PA.strategy = 'random')
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_NoEval_PA2',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## With Evaluation -------------------------------------------------------
cli::cli_h3("With Evaluation")

### Presence-Absence ------------
cli::cli_process_start("Presence-Absence")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp,
          expl.var = myExpl1,
          resp.xy = myRespXY,
          resp.name = myRespName,
          eval.resp.var = myResp,
          eval.expl.var = myExpl1,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_Eval_Presence-Absence',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
          seed.val = 42
        )
      invisible(get_predictions(myBiomodModelOut))
      invisible(get_evaluations(myBiomodModelOut))
      invisible(get_built_models(myBiomodModelOut))
      invisible(get_formal_data(myBiomodModelOut))
      
    })))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only ------------
cli::cli_process_start("Presence-Only")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO,
          expl.var = myExpl1,
          resp.xy = myRespXY_PO,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl1,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_Eval_Presence-Only',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### Presence-Only with NA ------------
cli::cli_process_start("Presence-Only with NA")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PO_NA,
          expl.var = myExpl1,
          resp.xy = myRespXY_PO_NA,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 500,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl1,
          eval.resp.xy = myRespXY)
      
      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_Eval_Presence-Only_with_NA',
          nb.rep = 2,
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

### Presence-Absence and Pseudo-Absences ------------
cli::cli_process_start("Presence-Absence and Pseudo-Absences")
this_try <- try({
  invisible(
    capture.output(suppressWarnings(suppressMessages({
      myBiomodData <-
        BIOMOD_FormatingData(
          resp.var = myResp_PA2,
          expl.var = myExpl1,
          resp.xy = myRespXY_PA2,
          resp.name = myRespName,
          PA.nb.rep = 2,
          PA.nb.absences = 1600,
          PA.strategy = 'random',
          eval.resp.var = myResp,
          eval.expl.var = myExpl1,
          eval.resp.xy = myRespXY)

      myBiomodModelOut <-
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'NoCat1_Eval_PA2',
          nb.rep = 2,
          models = c("MAXENT.Phillips"),
          data.split.perc = 80,
          var.import = 3,
          metric.eval = c('TSS','ROC'),
          do.full.models = FALSE,
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
  Error_Modeling <- Error_Modeling + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

