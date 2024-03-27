
#' Prepare training features from an RGB image
#'
#' Processes RGB (predictors) and temperature (target) raster data to generate a feature
#' dataset suitablefor training machine learning models. This function computes several spatial
#' statistics (mean and standard deviation) within multiple window sizes (3x3, 5x5,
#' 7x7, 9x9) for each RGB channel. Additionally, it applies Sobel and Laplacian filters
#' for edge detection. The function ensures that the RGB and temperature rasters are
#' compatible in terms of spatial dimensions and resolution before proceeding.
#'
#' @param rgb_rst A `SpatRaster` object containing the RGB SpatRaster data.
#' @param temp_rst A `SpatRaster` object containing the temperature SpatRaster data.
#'                  Default to NULL in case only prediction is to be performed. For
#'                  training this needs to be defined.
#'
#' @return A data frame where each row represents a pixel and columns represent the
#'         temperature (if temp is not NULL), original RGB values, spatial statistics for
#'         each RGB channel,and the results from Sobel and Laplacian edge detection.
#'         This data frame is suitable for use in machine learning model training/prediction.
#'
#' @importFrom terra values focal
#'
#' @export


prepare_train_data <- function(rgb_rst, temp_rst=NULL){

  # Check RGB data
  if(!inherits(rgb_rst,"SpatRaster")){
    stop("rgb_rst must be an object of class SpatRaster")
  }

  # Check thermal raster data
  if(!is.null(temp_rst)){

    if(!inherits(temp_rst,"SpatRaster")){
      stop("temp_rst must be an object of class SpatRaster")
    }

    if(!compare_rasters(rgb_rst, temp_rst)) {
      stop("rgb_rst and temp_rst are not equivalent for nrows, ncols, resolution or CRS!")
    }
  }

  r <- rgb_rst[[1]] / 255
  g <- rgb_rst[[2]] / 255
  b <- rgb_rst[[3]] / 255

  ## 3 x 3 ops
  f3avg_r <- terra::focal(r, w=3, fun="mean", na.rm=TRUE)
  f3std_r <- terra::focal(r, w=3, fun="sd", na.rm=TRUE)
  f3avg_g <- terra::focal(g, w=3, fun="mean", na.rm=TRUE)
  f3std_g <- terra::focal(g, w=3, fun="sd", na.rm=TRUE)
  f3avg_b <- terra::focal(b, w=3, fun="mean", na.rm=TRUE)
  f3std_b <- terra::focal(b, w=3, fun="sd", na.rm=TRUE)

  ## 5 x 5 ops
  f5avg_r <- terra::focal(r, w=5, fun="mean", na.rm=TRUE)
  f5std_r <- terra::focal(r, w=5, fun="sd", na.rm=TRUE)
  f5avg_g <- terra::focal(g, w=5, fun="mean", na.rm=TRUE)
  f5std_g <- terra::focal(g, w=5, fun="sd", na.rm=TRUE)
  f5avg_b <- terra::focal(b, w=5, fun="mean", na.rm=TRUE)
  f5std_b <- terra::focal(b, w=5, fun="sd", na.rm=TRUE)

  ## 7 x 7 ops
  f7avg_r <- terra::focal(r, w=7, fun="mean", na.rm=TRUE)
  f7std_r <- terra::focal(r, w=7, fun="sd", na.rm=TRUE)
  f7avg_g <- terra::focal(g, w=7, fun="mean", na.rm=TRUE)
  f7std_g <- terra::focal(g, w=7, fun="sd", na.rm=TRUE)
  f7avg_b <- terra::focal(b, w=7, fun="mean", na.rm=TRUE)
  f7std_b <- terra::focal(b, w=7, fun="sd", na.rm=TRUE)

  ## 9 x 9 ops
  f9avg_r <- terra::focal(r, w=9, fun="mean", na.rm=TRUE)
  f9std_r <- terra::focal(r, w=9, fun="sd", na.rm=TRUE)
  f9avg_g <- terra::focal(g, w=9, fun="mean", na.rm=TRUE)
  f9std_g <- terra::focal(g, w=9, fun="sd", na.rm=TRUE)
  f9avg_b <- terra::focal(b, w=9, fun="mean", na.rm=TRUE)
  f9std_b <- terra::focal(b, w=9, fun="sd", na.rm=TRUE)

  # Sobel Filter
  sb <- sobel_filter(rgb_rst)
  # Laplacian Filter
  lp <- laplacian_filter(rgb_rst)

  ## ------------------------------------------------------------------- ##
  ## Load the data from raster objects

  # Define variable/column names

  vn <- c("r","g","b",
    #"R_G","R_B","G_B",
    ## 3x3 data
    "f3avg_r", "f3std_r",
    "f3avg_g", "f3std_g",
    "f3avg_b", "f3std_b",
    ## 5x5 data
    "f5avg_r", "f5std_r",
    "f5avg_g", "f5std_g",
    "f5avg_b", "f5std_b",
    ## 7x7 data
    "f7avg_r", "f7std_r",
    "f7avg_g", "f7std_g",
    "f7avg_b", "f7std_b",
    ## 9x9 data
    "f9avg_r", "f9std_r",
    "f9avg_g", "f9std_g",
    "f9avg_b", "f9std_b",
    "lap_filt","sob_filt")

  if(!is.null(temp_rst)){
    # Data for training including the targt variable
    rstdf <- terra::values(c(temp_rst,
                             r, g, b,
                             ## 3x3 data
                             f3avg_r, f3std_r,
                             f3avg_g, f3std_g,
                             f3avg_b, f3std_b,

                             ## 5x5 data
                             f5avg_r, f5std_r,
                             f5avg_g, f5std_g,
                             f5avg_b, f5std_b,

                             ## 7x7 data
                             f7avg_r, f7std_r,
                             f7avg_g, f7std_g,
                             f7avg_b, f7std_b,

                             ## 9x9 data
                             f9avg_r, f9std_r,
                             f9avg_g, f9std_g,
                             f9avg_b, f9std_b,
                             lp, sb)) %>%
      na.omit() %>%
      as.data.frame()
    ## Rename the dataset
    colnames(rstdf) <- c("temp",
                         vn)
  }else{
    # Only data for prediction
    rstdf <- terra::values(c(r, g, b,
                             ## 3x3 data
                             f3avg_r, f3std_r,
                             f3avg_g, f3std_g,
                             f3avg_b, f3std_b,

                             ## 5x5 data
                             f5avg_r, f5std_r,
                             f5avg_g, f5std_g,
                             f5avg_b, f5std_b,

                             ## 7x7 data
                             f7avg_r, f7std_r,
                             f7avg_g, f7std_g,
                             f7avg_b, f7std_b,

                             ## 9x9 data
                             f9avg_r, f9std_r,
                             f9avg_g, f9std_g,
                             f9avg_b, f9std_b,
                             lp, sb)) %>%
      na.omit() %>%
      as.data.frame()
    ## Rename the dataset
    colnames(rstdf) <- vn
  }

  return(rstdf)

}

#' Predict thermal data from RGB image using Random Forest
#'
#' Enhances the resolution or predicts thermal data based on high-resolution RGB images
#' and an existing low-resolution thermal raster using a Random Forest model. The model
#' is trained on features derived from the RGB and low-resolution thermal data, optionally
#' scaled and with high correlation features removed. The function can then predict
#' thermal data aligned with a provided high-resolution RGB image, effectively
#' super-resolving the thermal image or generating new thermal data where none exists.
#'
#' @param rgb_rst An RGB `SpatRaster` object containing RGB data aligned with `temp_rst`.
#' @param temp_rst A `SpatRaster` object containing corresponding lower-resolution
#'                 thermal data.
#' @param rgb_high_res (Optional) A `SpatRaster` object containing high-resolution
#'                     RGB data for which thermal predictions are to be made.
#' @param npix The number of pixels to sample for model training.
#' @param rm_cor Logical flag indicating whether to remove highly correlated features
#'               before model training.
#' @param rm_cor_thresh Correlation threshold for removing correlated features.
#' @param get_rf_model Logical flag indicating whether to return the Random Forest
#'                     model along with the predicted temperature raster(s).
#' @param verbose Logical flag for printing progress messages.
#' @param scale_data Logical flag indicating whether to scale the feature data before
#'                   model training.
#' @param ... Additional arguments passed to the Random Forest model training function.
#'
#' @return Depending on `get_rf_model`, either a list containing the predicted
#'         temperature raster(s) and the Random Forest model, or just the predicted
#'         temperature raster(s). If `rgb_high_res` is provided, the list will include
#'         predictions for both the training resolution and the high-resolution RGB.
#'         List components include: `pred_temp_train` is the thermal image predicted
#'         for the training data (usually with lower resolution); `pred_temp_rgb_hr`
#'         is the predicted thermal data for the high-resolution RGB; `rf_mod` is the
#'         Random Forest model trained for the image.
#'
#' @importFrom ranger ranger
#' @importFrom terra values `values<-`
#' @importFrom dplyr select any_of
#' @importFrom caret findCorrelation
#' @importFrom stats cor na.omit predict
#'
#' @export


rf_thermal_from_rgb <- function(rgb_rst,
                         temp_rst,
                         rgb_high_res = NULL,
                         npix = 10000,
                         rm_cor = FALSE,
                         rm_cor_thresh = 0.98,
                         get_rf_model = TRUE,
                         verbose = TRUE,
                         scale_data = FALSE, ...){


  ## ------------------------------------------------------------------- ##

  ts <- Sys.time()

  if(verbose) cat("|-> Calculating features ...\n")

  rstdf <- prepare_train_data(rgb_rst, temp_rst)

  if(verbose) cat("Done.\n\n")

  ## ------------------------------------------------------------------- ##

  if(scale_data){

    if(verbose) cat("|-> Scaling data for training ...\n")

    rstdf <- cbind(temp = rstdf[,1], scale(rstdf[,-1])) %>%
      as.data.frame()

    if(verbose) cat("Done.\n\n")
  }


  ## ------------------------------------------------------------------- ##

  if(verbose) cat("|-> Extracting sample data for training ...\n")

  s <- sample(1:nrow(rstdf), npix, replace = FALSE)
  samp_rstdf <- rstdf[s,]

  if(verbose) cat("Done.\n\n")

  ## ------------------------------------------------------------------- ##


  if(rm_cor){

    if(verbose) cat("|-> Removing highly correlated features ...\n")

    cm <- cor(samp_rstdf[,-1], method="spearman")
    vars_to_rm <- caret::findCorrelation(cm, cutoff = rm_cor_thresh, names = TRUE)
    rstdf <- rstdf %>% dplyr::select(-dplyr::any_of(vars_to_rm))
    samp_rstdf <- samp_rstdf %>% dplyr::select(-dplyr::any_of(vars_to_rm))

    if(verbose) cat("Done.\n\n")

  }

  ## ------------------------------------------------------------------- ##

  if(verbose) cat("|-> Training Random Forest model ...\n")

  rf <- ranger::ranger(temp ~ ., data = samp_rstdf)

  if(verbose) cat("Done.\n\n")

  ## ------------------------------------------------------------------- ##

  if(verbose) cat("|-> Super-resolving the thermal image ...\n")

  pred_temp <- temp_rst

  predv <- predict(rf, data = rstdf)
  values(pred_temp) <- predv$predictions

  if(!is.null(rgb_high_res)){

    pred_temp_hr <- rgb_high_res[[1]]
    names(pred_temp_hr)<-"temp_super_res_dl"

    if(verbose) cat("|-> Calculating features for the high-resolution RGB image ...\n")

    rstdf <- prepare_train_data(rgb_high_res)

    if(verbose) cat("Done.\n\n")

    ## ------------------------------------------------------------------- ##

    if(scale_data){

      if(verbose) cat("|-> Scaling high-resolution data for prediction ...\n")

      rstdf <- scale(rstdf) %>%
        as.data.frame()

      if(verbose) cat("Done.\n\n")
    }

    predv <- predict(rf, data = rstdf)
    values(pred_temp_hr) <- predv$predictions

  }

  if(verbose) cat("Done.\n")
  if(verbose) cat("[Run Time:",pretty_duration(Sys.time() - ts),"]\n\n")

  if(!is.null(rgb_high_res)){
      out <- list(
        pred_temp_train = pred_temp,
        pred_temp_rgb_hr = pred_temp_hr
      )
  }else{
    out <- list(
      pred_temp_train = pred_temp
    )
  }

  if(get_rf_model){

    out$rf_mod <- rf

    return(
      out
    )
  }else{
    return(out)
  }
}

#' Predict thermal data from RGB images using Deep Learning
#'
#' Enhances the resolution or predicts thermal data based on high-resolution RGB images
#' and an existing low-resolution thermal raster using a deep learning model. The model
#' is trained on features derived from the RGB and low-resolution thermal data, optionally
#' scaled and with high correlation features removed. The function can then predict
#' thermal data aligned with a provided high-resolution RGB image, effectively
#' super-resolving the thermal image or generating new thermal data where none exists.
#'
#' @param rgb_rst An RGB `SpatRaster` object containing RGB data aligned with `temp_rst`.
#' @param temp_rst A `SpatRaster` object containing corresponding lower-resolution
#'                 thermal data.
#' @param rgb_high_res (Optional) A `SpatRaster` object containing high-resolution
#'                     RGB data for which thermal predictions are to be made.
#' @param npix The number of pixels to sample for model training.
#' @param rm_cor Logical flag indicating whether to remove highly correlated features
#'               before model training.
#' @param rm_cor_thresh Correlation threshold for removing correlated features.
#' @param get_dl_model Logical flag indicating whether to return the deep learning
#'                     model along with the predicted temperature raster(s).
#' @param verbose Logical flag for printing progress messages.
#' @param scale_data Logical flag indicating whether to scale the feature data before
#'                   model training.
#' @param n_epochs Number of epochs for model training.
#' @param validation_split Fraction of the data to be used as validation data.
#' @param learning_rate Learning rate for the optimizer.
#'
#' @return Depending on `get_dl_model`, either a list containing the predicted
#'         temperature raster(s) and the deep learning model, or just the predicted
#'         temperature raster(s). If `rgb_high_res` is provided, the list will include
#'         predictions for both the training resolution and the high-resolution RGB.
#'         List components include: `pred_temp_train` is the thermal image predicted
#'         for the training data (usually with lower resolution); `pred_temp_rgb_hr`
#'         is the predicted thermal data for the high-resolution RGB; `dl_mod` is the
#'         Deep Learning model trained for the image.
#'
#' @details
#' The deep learning model used for predicting thermal data from high-resolution RGB
#' images is built using the Keras library in R. This model is structured as a sequential
#' model, comprising several densely connected neural network layers with ReLU activations,
#' dropout layers for regularization, and a final dense layer with a single unit for
#' regression output. The model aims to learn complex relationships between spatial features
#' derived from RGB images and corresponding thermal data, allowing it to predict temperature
#' values for new, high-resolution RGB images.
#'
#' Model Structure
#'
#' Input Layer: The input layer is designed to accept the flattened feature vectors derived
#' from the color/spatial/texture features of RGB and thermal data.
#'
#' The model includes multiple dense (fully connected) layers with ReLU (Rectified Linear Unit)
#' activation functions. These layers are responsible for capturing nonlinear relationships in
#' the data.
#'
#' Dropout layers are inserted between dense layers to reduce the risk of overfitting by randomly
#' setting a fraction of input units to 0 at each update during training time. This helps improve
#' model generalization.
#'
#' Activity regularization layers are applied to introduce a penalty on the layer's activation,
#' further aiding in preventing overfitting and promoting simpler models.
#'
#' Output Layer: The final layer is a dense layer with a single unit, suitable for regression tasks.
#' This layer outputs the predicted temperature value for each input feature vector.
#'
#' Training Configuration
#'
#' Loss Function: Mean Squared Error (MSE), suitable for regression problems, measuring the average
#' of the squares of the errors between true and predicted values.
#'
#' Optimizer: Adam, with a specified learning rate, an algorithm for first-order gradient-based
#' optimization of stochastic objective functions. Metrics: Mean Absolute Error (MAE), providing an
#' average of absolute differences between predicted and actual values, offering an interpretation
#' of prediction accuracy.
#'
#' Epochs and Validation Split: The model is trained for a predefined number of epochs, with a portion
#' of the data reserved for validation to monitor and prevent overfitting.
#'
#'
#' @importFrom keras keras_model_sequential layer_dense layer_dropout
#'            layer_activity_regularization optimizer_adam compile fit
#' @importFrom terra values `values<-`
#' @importFrom dplyr select any_of
#' @importFrom caret findCorrelation
#' @importFrom stats cor na.omit predict
#' @importFrom rlang `.data`
#' @export


dl_thermal_from_rgb <- function(rgb_rst,
                         temp_rst,
                         rgb_high_res = NULL,
                         npix = 10000,
                         rm_cor = FALSE,
                         rm_cor_thresh = 0.98,
                         get_dl_model = TRUE,
                         verbose = TRUE,
                         scale_data = TRUE,
                         n_epochs = 80,
                         validation_split = 0.2,
                         learning_rate = 0.01){


  ## ------------------------------------------------------------------- ##

  ts <- Sys.time()

  if(verbose) cat("|-> Calculating features ...\n")

  rstdf <- prepare_train_data(rgb_rst, temp_rst)

  if(verbose) cat("Done.\n\n")

  ## ------------------------------------------------------------------- ##

  if(scale_data){
    rstdf <- cbind(temp = rstdf[,1], scale(rstdf[,-1])) %>%
      as.data.frame()
  }


  ## ------------------------------------------------------------------- ##

  if(verbose) cat("|-> Extracting sample data for training ...\n")

  s <- sample(1:nrow(rstdf), npix, replace = FALSE)
  samp_rstdf <- rstdf[s,]

  if(verbose) cat("Done.\n\n")

  ## ------------------------------------------------------------------- ##


  if(rm_cor){

    if(verbose) cat("|-> Removing highly correlated features ...\n")

    cm <- cor(samp_rstdf[,-1], method="spearman")
    vars_to_rm <- caret::findCorrelation(cm, cutoff = rm_cor_thresh, names = TRUE)
    rstdf <- rstdf %>% dplyr::select(-dplyr::any_of(vars_to_rm))
    samp_rstdf <- samp_rstdf %>% dplyr::select(-dplyr::any_of(vars_to_rm))

    if(verbose) cat("Done.\n\n")

  }

  ## ------------------------------------------------------------------- ##

  if(verbose) cat("|-> Training Deep Learning model ...\n")

  #rf <- ranger(temp ~ ., data = samp_rstdf)

  samp_rstdf <- as.matrix(samp_rstdf)
  train_data <- samp_rstdf[,-1]
  train_targets <- samp_rstdf[,"temp"]

  model <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = 256, activation = 'relu',
                       input_shape = c(ncol(train_data))) %>%
    keras::layer_dropout(rate = 0.05) %>%
    keras::layer_activity_regularization(l1 = 0.001, l2 = 0.001) %>%
    keras::layer_dense(units = 128, activation = 'relu') %>%
    keras::layer_dropout(rate = 0.025) %>%
    keras::layer_activity_regularization() %>%
    keras::layer_dense(units = 64, activation = 'relu') %>%
    keras::layer_dropout(rate = 0.01) %>%
    keras::layer_activity_regularization() %>%
    keras::layer_dense(units = 1, name = "output") # Since this is a regression problem

  # Compile the model
  model %>% keras::compile(
    optimizer = keras::optimizer_adam(learning_rate = learning_rate),
    loss = 'mse',
    metrics = c('mae')
  )

  # Train the model
  history <- model %>% keras::fit(
    train_data,
    train_targets,
    epochs = n_epochs,
    validation_split = validation_split
  )

  if(verbose) cat("Done.\n\n")

  ## ------------------------------------------------------------------- ##

  if(verbose) cat("|-> Super-resolving the thermal image ...\n")

  pred_temp <- temp_rst

  predv <- model %>% predict(as.matrix(rstdf %>% dplyr::select(-.data$temp)))

  values(pred_temp) <- as.numeric(predv[,1])

  if(!is.null(rgb_high_res)){

    pred_temp_hr <- rgb_high_res[[1]]
    names(pred_temp_hr)<-"temp_super_res_dl"

    if(verbose) cat("|-> Calculating features for the high-resolution RGB image ...\n")

    rstdf <- prepare_train_data(rgb_high_res)

    if(verbose) cat("Done.\n\n")

    ## ------------------------------------------------------------------- ##

    if(scale_data){

      if(verbose) cat("|-> Scaling high-resolution data for prediction ...\n")

      rstdf <- scale(rstdf) %>%
        as.data.frame()

      if(verbose) cat("Done.\n\n")
    }

    predv <- model %>% predict(as.matrix(rstdf))
    values(pred_temp_hr) <- as.numeric(predv[,1])

  }

  if(verbose) cat("Done.\n")
  if(verbose) cat("[Run Time:",pretty_duration(Sys.time() - ts),"]\n\n")

  if(!is.null(rgb_high_res)){
      out <- list(
        pred_temp_train = pred_temp,
        pred_temp_rgb_hr = pred_temp_hr
      )
  }else{
    out <- list(
      pred_temp_train = pred_temp
    )
  }

  if(get_dl_model){

    out$dl_mod <- model

    return(
      out
    )
  }else{
    return(out)
  }

}

