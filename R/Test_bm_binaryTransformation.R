cli::cli_h1("bm_BinaryTransformation")
# setup data --------------------------------------------------------------
cli::cli_h2("Setup data")
Error_BinaryTransformation <- 0
set.seed(42)
data_num <- rnorm(100, 500, 100)
data_numNA <- data_num
data_numNA[sample(1:100, 20)] <- NA
threshold_1 <- 500
threshold_2 <- runif(2, min = 400, max = 600)
threshold_3 <- runif(3, min = 400, max = 600)
data_df <- data.frame("model1" = rnorm(100, 500, 100),
                      "model2" = rnorm(100, 500, 100),
                      "model3" = rnorm(100, 500, 100))
data_dfNA <- data_df
data_dfNA[sample(1:100, size = 20), 1] <- NA
data_dfNA[sample(1:100, size = 20), 2] <- NA
data_dfNA[sample(1:100, size = 20), 3] <- NA
data_df1 <- data.frame("model1" = rnorm(100, 500, 100))
data_df1NA <- data_df1
data_df1NA[sample(1:100, size = 20), 1] <- NA
data_matrix <- matrix( rnorm(100, 500, 100), ncol = 2, nrow = 500)
colnames(data_matrix) <- c("model1","model2")
data_matrixNA <- data_matrix
data_matrixNA[sample(1:ncell(data_matrix), size = 40, replace = FALSE)] <- NA
data_matrix1 <- matrix( rnorm(100, 500, 100), ncol = 1)
colnames(data_matrix1) <- c("model1")
data_matrix1NA <- data_matrix1
data_matrix1NA[sample(1:ncell(data_matrix1), size = 20, replace = FALSE)] <- NA

data(bioclim_current)
rast_thr1 <- 80
rast_thr3 <- c(50,1000,400)
data_SpatRaster <- subset(rast(bioclim_current), 1:3)

## From continuous to binary / filtered vector
# vec.d_bin <- bm_BinaryTransformation(data = vec.d, threshold = 500)
# vec.d_filt <- bm_BinaryTransformation(data = vec.d, threshold = 500, do.filtering = TRUE)
# cbind(vec.d, vec.d_bin, vec.d_filt)
# 


# data.frame  -------------------------------------------------------------


cli::cli_h2("data.frame (3 cols) method ")


## No NA -------------------------------------------------------------------
cli::cli_h3("No NA")

cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_df, threshold = threshold_3, do.filtering = FALSE) 
    stopifnot(inherits(data_BT,"data.frame"))
    stopifnot(dim(data_BT) == c(100,3))
    stopifnot(is.numeric(unlist(data_BT)))
    stopifnot(all(na.omit(unlist(data_BT)) %in% c(0,1)))
    stopifnot(names(data_BT) == names(data_df))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_df, threshold = threshold_3, do.filtering = TRUE)
    stopifnot(inherits(data_BT,"data.frame"))
    stopifnot(dim(data_BT) == c(100,3))
    stopifnot(is.numeric(unlist(data_BT)))
    stopifnot(all(na.omit(unlist(data_BT)) >= 0))
    stopifnot(names(data_BT) == names(data_df))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


## With NA -----------------------------------------------------------------

cli::cli_h3("With NA")

cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_dfNA, threshold = threshold_3, do.filtering = FALSE) 
    stopifnot(inherits(data_BT,"data.frame"))
    stopifnot(dim(data_BT) == c(100,3))
    stopifnot(is.numeric(unlist(data_BT)))
    stopifnot(all(na.omit(unlist(data_BT)) %in% c(0,1)))
    stopifnot(names(data_BT) == names(data_df))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_dfNA, threshold = threshold_3, do.filtering = TRUE) 
    stopifnot(inherits(data_BT,"data.frame"))
    stopifnot(dim(data_BT) == c(100,3))
    stopifnot(is.numeric(unlist(data_BT)))
    stopifnot(all(na.omit(unlist(data_BT)) >= 0))
    stopifnot(names(data_BT) == names(data_df))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# data.frame 1 col --------------------------------------------------------


cli::cli_h2("data.frame (1 col) method ")

## No NA -------------------------------------------------------------------

cli::cli_h2("No NA ")


cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_df1, threshold = threshold_1)
    stopifnot(inherits(data_BT, "data.frame"))
    stopifnot(nrow(data_BT) == 100)
    stopifnot(is.numeric(unlist(data_BT)))
    stopifnot(all(na.omit(unlist(data_BT)) %in% c(0,1)))
    stopifnot(names(data_BT) == names(data_df1))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible(
    data_BT <- bm_BinaryTransformation(data = data_df1, threshold = threshold_1, do.filtering = TRUE)
  )
  stopifnot(inherits(data_BT, "data.frame"))
  stopifnot(nrow(data_BT) == 100)
  stopifnot(is.numeric(unlist(data_BT)))
  stopifnot(all(na.omit(unlist(data_BT)) >= 0))
  stopifnot(names(data_BT) == names(data_df1))
  
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## With NA -------------------------------------------------------------------

cli::cli_h2("With NA ")


cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_df1NA, threshold = threshold_1)
    stopifnot(inherits(data_BT, "data.frame"))
    stopifnot(nrow(data_BT) == 100)
    stopifnot(is.numeric(unlist(data_BT)))
    stopifnot(all(na.omit(unlist(data_BT)) %in% c(0,1)))
    stopifnot(names(data_BT) == names(data_df1))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible(
    data_BT <- bm_BinaryTransformation(data = data_df1NA, threshold = threshold_1, do.filtering = TRUE)
  )
  stopifnot(inherits(data_BT, "data.frame"))
  stopifnot(nrow(data_BT) == 100)
  stopifnot(is.numeric(unlist(data_BT)))
  stopifnot(all(na.omit(unlist(data_BT)) >= 0))
  stopifnot(names(data_BT) == names(data_df1))
  
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# matrix ------------------------------------------------------------------

cli::cli_h2("matrix (3 cols) method ")


## no NA -------------------------------------------------------------------

cli::cli_h3("No NA ")

cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrix, threshold = threshold_2)
    stopifnot(is.matrix(data_BT))
    stopifnot(dim(data_BT) == c(500,2))
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) %in% c(0,1)))
    stopifnot(colnames(data_BT) == colnames(data_matrix))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrix, threshold = threshold_2, do.filtering = TRUE)
    stopifnot(is.matrix(data_BT))
    stopifnot(dim(data_BT) == c(500,2))
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) >= 0))
    stopifnot(colnames(data_BT) == colnames(data_matrix))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## with NA -------------------------------------------------------------------

cli::cli_h3("With NA ")


cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrixNA, threshold = threshold_2)
    stopifnot(is.matrix(data_BT))
    stopifnot(dim(data_BT) == c(500,2))
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) %in% c(0,1)))
    stopifnot(colnames(data_BT) == colnames(data_matrix))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrixNA, threshold = threshold_2, do.filtering = TRUE)
    stopifnot(is.matrix(data_BT))
    stopifnot(dim(data_BT) == c(500,2))
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) >= 0))
    stopifnot(colnames(data_BT) == colnames(data_matrix))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# matrix (1 col) ------------------------------------------------------------------

cli::cli_h2("matrix (1 col) method ")

## no NA -------------------------------------------------------------------

cli::cli_h3("no NA ")



cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrix1, threshold = threshold_1)
    stopifnot(inherits(data_BT, "matrix"))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) %in% c(0,1)))
    stopifnot(colnames(data_BT) == colnames(data_matrix1))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrix1, threshold = threshold_1, do.filtering = TRUE)
    stopifnot(inherits(data_BT, "matrix"))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) >= 0))
    stopifnot(colnames(data_BT) == colnames(data_matrix1))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## With NA -------------------------------------------------------------------

cli::cli_h3("With NA ")



cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrix1NA, threshold = threshold_1)
    stopifnot(inherits(data_BT, "matrix"))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) %in% c(0,1)))
    stopifnot(colnames(data_BT) == colnames(data_matrix1))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}



cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_matrix1NA, threshold = threshold_1, do.filtering = TRUE)
    stopifnot(inherits(data_BT, "matrix"))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) >= 0))
    stopifnot(colnames(data_BT) == colnames(data_matrix1))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


# numeric -----------------------------------------------------------------


cli::cli_h2("numeric method ")

## No NA -----------------------------------------------------------------
cli::cli_h3("No NA ")

cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_num, threshold = threshold_1)
    stopifnot(is.vector(data_BT))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) %in% c(0,1)))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_num, threshold = threshold_1, do.filtering = TRUE)
    stopifnot(is.vector(data_BT))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) >= 0))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

## With NA -----------------------------------------------------------------
cli::cli_h3("With NA ")

cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_numNA, threshold = threshold_1)
    stopifnot(is.vector(data_BT))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) %in% c(0,1)))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_numNA, threshold = threshold_1, do.filtering = TRUE)
    stopifnot(is.vector(data_BT))
    stopifnot(length(data_BT) == 100)
    stopifnot(is.numeric(c(data_BT)))
    stopifnot(all(na.omit(c(data_BT)) >= 0))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

# SpatRaster --------------------------------------------------------------


cli::cli_h2("SpatRaster method ")
cli::cli_h3("1 threshold ")

cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_SpatRaster, threshold = rast_thr1)
    stopifnot(inherits(data_BT,"SpatRaster"))
    stopifnot(is.numeric(data_BT[]))
    stopifnot(all(na.omit(c(data_BT[])) %in% c(0,1)))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_SpatRaster, threshold = rast_thr1, do.filtering = TRUE)
    stopifnot(inherits(data_BT,"SpatRaster"))
    stopifnot(is.numeric(data_BT[]))
    stopifnot(all(na.omit(c(data_BT[])) >= 0))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}

cli::cli_h3("3 threshold ")

cli::cli_process_start("no filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_SpatRaster, threshold = rast_thr3)
    stopifnot(inherits(data_BT,"SpatRaster"))
    stopifnot(is.numeric(data_BT[]))
    stopifnot(all(na.omit(c(data_BT[])) %in% c(0,1)))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}


cli::cli_process_start("with filtering")
this_try <- try({
  invisible({
    data_BT <- bm_BinaryTransformation(data = data_SpatRaster, threshold = rast_thr3, do.filtering = TRUE)
    stopifnot(inherits(data_BT,"SpatRaster"))
    stopifnot(is.numeric(data_BT[]))
    stopifnot(all(na.omit(c(data_BT[])) >= 0))
  })
}, silent = TRUE)


if(inherits(this_try, "try-error")){
  Error_BinaryTransformation <- Error_BinaryTransformation + 1
  cli::cli_process_failed()
} else {
  cli::cli_process_done()
}
