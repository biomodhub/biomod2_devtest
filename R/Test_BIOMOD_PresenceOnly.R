cli::cli_h1("BIOMOD_PresenceOnly")

# setup data --------------------------------------------------------------
cli::cli_h2("Setup data")
Error_PresenceOnly <- 0

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

## myExpl.cat.raster -------------------------------------------------------
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


## with df -----------------------------------------------------------------

myResp <- as.numeric(DataSpecies[, myRespName])
myResp.SpatVector <- vect(x = data.frame("presence" = myResp, 
                                         "x" = DataSpecies[,'X_WGS84'], 
                                         "y" = DataSpecies[,'Y_WGS84']), 
                          geom = c('x', 'y') )
myExpl.df <- extract(myExpl, y = myResp.SpatVector, ID = FALSE)
myExpl.cat.df <- extract(myExpl.cat, y = myResp.SpatVector, ID = FALSE)
tmpdf <- cbind(myRespXY,myExpl.df)
rownames(tmpdf) <- NULL
myExpl.SpatVector <- vect(
  tmpdf,
  geom = c("X_WGS84", "Y_WGS84")
)
myExpl.SpatVector$ID <- NULL
myExpl.spdf <- as(myExpl.SpatVector, "Spatial")


# Load Models --------------------------------------

## One variable + pseudo absences -------------------------------------------------------

invisible(
  capture.output(suppressWarnings(suppressMessages({
    file.out <- paste0(myRespName, "/", myRespName, ".NoCat1_NoEval_Presence-Only.models.out")
    if (file.exists(file.out)) {
      myBiomodModelOut1 <- get(load(file.out))
    } 
    file.out <- paste0(myRespName, "/", myRespName, ".NoCat1_NoEval_Presence-Only.ensemble.models.out")
    if (file.exists(file.out)) {
      myBiomodEnsembleOut1 <- get(load(file.out))
    } else {
      myBiomodEnsembleOut1 <-   BIOMOD_EnsembleModeling(
        bm.mod = myBiomodModelOut1,
        models.chosen = 'all',
        em.by = 'PA',
        metric.select = c('TSS'),
        metric.select.thresh = c(0.7),
        var.import = 3,
        metric.eval = c('TSS', 'ROC'),
        em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),
        EMci.alpha = 0.05,
        EMwmean.decay = 'proportional',
        seed.val = 42)
    }
  }))))


## No Cat + No Eval -------------------------------------------------------

invisible(
  capture.output(suppressWarnings(suppressMessages({
    file.out <- paste0(myRespName, "/", myRespName, ".NoCat_NoEval_Presence-Absence.models.out")
    if (file.exists(file.out)) {
      myBiomodModelOut_noEval <- get(load(file.out))
    } 
    file.out <- paste0(myRespName, "/", myRespName, ".NoCat_NoEval_Presence-Absence.ensemble.models.out")
    if (file.exists(file.out)) {
      myBiomodEnsembleOut_noEval <- get(load(file.out))
    } 
  }))))


## No Cat + Eval -----------------------------------------------------------


invisible(
  capture.output(suppressWarnings(suppressMessages({
    file.out <- paste0(myRespName, "/", myRespName, ".NoCat_Eval_Presence-Absence.models.out")
    if (file.exists(file.out)) {
      myBiomodModelOut <- get(load(file.out))
    } 
  }))))

## Cat + Eval -----------------------------------------------------------
invisible(
  capture.output(suppressWarnings(suppressMessages({
    file.out <- paste0(myRespName, "/", myRespName, ".Cat_Eval_Presence-Absence.models.out")
    if (file.exists(file.out)) {
      myBiomodModelOut.cat <- get(load(file.out))
    } 
    
  }))))

## Ensemble no Cat + Eval, em.by = PA_dataset+repet -----------------
invisible(
  capture.output(suppressWarnings(suppressMessages({
    myBiomodData <- 
      BIOMOD_FormatingData(
        resp.var = myResp_PO,
        expl.var = myExpl,
        resp.xy = myRespXY_PO,
        resp.name = myRespName,
        eval.resp.var = myResp,
        eval.resp.xy = myRespXY,
        PA.nb.rep = 2,
        PA.nb.absences = 500, 
        PA.strategy = 'random')
    file.out <- paste0(myRespName, "/", myRespName, ".Ensemble_PresenceOnly.models.out")
    if (file.exists(file.out)) {
      myBMout <- get(load(file.out))
    } else { 
      myBMout <- 
      BIOMOD_Modeling(
        bm.format = myBiomodData,
        bm.options = BIOMOD_ModelingOptions(),
        modeling.id = 'Ensemble_PresenceOnly',
        CV.strategy = 'random',
        CV.nb.rep = 2,
        CV.perc = 0.8,
        var.import = 2,
        metric.eval = c('TSS','ROC'),
        seed.val = 42,
        nb.cpu = 6
      )
    }
    file.out <- paste0(myRespName, "/", myRespName, ".Ensemble_PresenceOnly.ensemble.models.out")
    if (file.exists(file.out)) {
      myBiomodEM <- get(load(file.out))
    }  else {
      
    
    myBiomodEM <- BIOMOD_EnsembleModeling(
      bm.mod = myBMout,
      models.chosen = 'all',
      em.by = 'PA+run',
      metric.select = c('TSS'),
      metric.select.thresh = c(0.7),
      var.import = 3,
      metric.eval = c('TSS', 'ROC'),
      em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),
      EMci.alpha = 0.05,
      EMwmean.decay = 'proportional',
      seed.val = 42)
    }
  }))))



## Ensemble no Cat + Eval, em.by = all -----------------
invisible(
  capture.output(suppressWarnings(suppressMessages({
    myBiomodData <- 
      BIOMOD_FormatingData(
        resp.var = myResp_PO,
        expl.var = myExpl,
        resp.xy = myRespXY_PO,
        resp.name = myRespName,
        eval.resp.var = myResp,
        eval.resp.xy = myRespXY,
        PA.nb.rep = 2,
        PA.nb.absences = 500, 
        PA.strategy = 'random')
    file.out <- paste0(myRespName, "/", myRespName, ".Ensemble_PresenceOnly2.models.out")
    if (file.exists(file.out)) {
      myBMout2 <- get(load(file.out))
    } else { 
      myBMout2 <- 
        BIOMOD_Modeling(
          bm.format = myBiomodData,
          bm.options = BIOMOD_ModelingOptions(),
          modeling.id = 'Ensemble_PresenceOnly2',
          CV.strategy = 'random',
          CV.nb.rep = 2,
          CV.perc = 0.8,
          var.import = 2,
          metric.eval = c('TSS','ROC'),
          seed.val = 42,
          nb.cpu = 6
        )
    }
    file.out <- paste0(myRespName, "/", myRespName, ".Ensemble_PresenceOnly2.ensemble.models.out")
    if (file.exists(file.out)) {
      myBiomodEM2 <- get(load(file.out))
    }  else {
      
      
      myBiomodEM2 <- BIOMOD_EnsembleModeling(
        bm.mod = myBMout2,
        models.chosen = 'all',
        em.by = 'all',
        metric.select = c('TSS'),
        metric.select.thresh = c(0.7),
        var.import = 3,
        metric.eval = c('TSS', 'ROC'),
        em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),
        EMci.alpha = 0.05,
        EMwmean.decay = 'proportional',
        seed.val = 42)
    }
  })))
  )




# no bg.env ------------------------------------------------
cli::cli_h2("No bg.env")

## no eval ------------------------------------------------

cli::cli_process_start("No bg.env / No Cat / No Eval")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut_noEval,
                                        bm.em = myBiomodEnsembleOut_noEval)

      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut_noEval)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.em = myBiomodEnsembleOut_noEval)
    }))
  )
}, silent = TRUE)

if (inherits(this_try, "try-error")) {
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## no cat ------------------------------------------------

cli::cli_process_start("No bg.env / No Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBMout,
                                        bm.em = myBiomodEM)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.em = myBiomodEM)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBMout)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## with cat ----------------------------------------------------------------

cli::cli_process_start("No bg.env / With Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut.cat)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# with bg.env ------------------------------------------------

## SpatRaster ----------------------------------------------------------------
cli::cli_h2("SpatRaster")

## no eval ------------------------------------------------

cli::cli_process_start("With bg.env / No Cat / No Eval")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut_noEval,
                                        bm.em = myBiomodEnsembleOut_noEval,
                                        bg.env = myExpl)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut_noEval,
                                        bg.env = myExpl)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.em = myBiomodEnsembleOut_noEval,
                                        bg.env = myExpl)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}




### no cat ------------------------------------------------

cli::cli_process_start("with bg.env / No Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBMout,
                                        bm.em = myBiomodEM,
                                        bg.env = myExpl)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.em = myBiomodEM,
                                        bg.env = myExpl)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
                                        bg.env = myExpl)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### with cat ----------------------------------------------------------------

cli::cli_process_start("with bg.env / With Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut.cat,
                                        bg.env = myExpl.cat)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



## SpatVector ----------------------------------------------------------------
cli::cli_h2("SpatVector")

### no cat ------------------------------------------------

cli::cli_process_start("with bg.env / No Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBMout,
                                        bm.em = myBiomodEM,
                                        bg.env = myExpl.SpatVector)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.em = myBiomodEM,
                                        bg.env = myExpl.SpatVector)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
                                        bg.env = myExpl.SpatVector)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## data.frame ----------------------------------------------------------------
cli::cli_h2("data.frame")

### no cat ------------------------------------------------

cli::cli_process_start("with bg.env / No Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBMout,
                                        bm.em = myBiomodEM,
                                        bg.env = myExpl.df)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.em = myBiomodEM,
                                        bg.env = myExpl.df)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
                                        bg.env = myExpl.df)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


### with cat ----------------------------------------------------------------

cli::cli_process_start("with bg.env / With Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut.cat,
                                        bg.env = myExpl.cat.df)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## numeric ----------------------------------------------------------------
cli::cli_h2("numeric")

### no cat ------------------------------------------------

cli::cli_process_start("with bg.env / No Cat")
this_try <- try({
  invisible(
    capture.output(suppressWarnings({
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut1,
                                        bm.em = myBiomodEnsembleOut1,
                                        bg.env = myExpl.df$bio3)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.em = myBiomodEnsembleOut1,
                                        bg.env = myExpl.df$bio3)
      myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut1,
                                        bg.env = myExpl.df$bio3)
    }))
  )
}, silent = TRUE)

if(inherits(this_try, "try-error")){
  Error_PresenceOnly <- Error_PresenceOnly + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}