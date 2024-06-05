

#' Remove outliers from a SpatRaster based on percentile thresholds
#'
#' This function removes outliers from a raster object by setting values outside
#' specified percentile thresholds to NA. It uses the quantile function to determine
#' the lower and upper bounds based on the `pmin` and `pmax` parameters, respectively.
#' Values below the lower percentile (`pmin`) or above the upper percentile (`pmax`)
#' are considered outliers and are set to NA. This can be useful for data
#' preprocessing or visualization, especially in spatial analysis and modeling,
#' to mitigate the impact of extreme values.
#'
#' @param rast A SpatRaster object from which outliers will be removed.
#' @param pmin The lower percentile used to define the lower threshold for outliers.
#'             Values below this percentile will be set to NA. Defaults to 5.
#' @param pmax The upper percentile used to define the upper threshold for outliers.
#'             Values above this percentile will be set to NA. Defaults to 95.
#' @param na.rm Remove NA's (default TRUE)
#'
#' @return A SpatRaster object similar to the input but with values outside the
#'         specified percentile thresholds set to NA.
#'
#' @importFrom terra values
#' @importFrom stats quantile
#'
#' @export
#'

remove_outliers <- function(rast,  pmin = 5, pmax = 100, na.rm=TRUE){

  # Get raster values
  tv <- terra::values(rast)
  # Calculate quantiles for pmnin and pmax
  qts <- quantile(tv, probs = c(pmin / 100, pmax / 100),na.rm=na.rm)

  # Set intervals as NA in extracted values
  tv[(tv < qts[1]) | tv > qts[2]] <- NA

  # Make new raster
  temp_c <- rast
  terra::values(temp_c) <- tv

  return(temp_c)
}

#' Compare two SpatRaster objects
#'
#' Compares two SpatRaster objects to determine if they have the same dimensions,
#' cell size, and coordinate reference system (CRS). This function performs a
#' series of checks to assess whether the two rasters are equivalent in terms
#' of their spatial properties, which is essential for many spatial analysis
#' tasks that require alignment of raster layers.
#'
#' @param raster1 The first raster object to compare.
#' @param raster2 The second raster object to compare.
#'
#' @return Returns TRUE if both rasters have the same number of rows, same number
#'         of columns, identical cell sizes, and the same CRS. Otherwise, returns FALSE.
#'
#' @importFrom terra nrow ncol res crs
#'
#' @export
#'

compare_rasters <- function(raster1, raster2) {
  # Check number of rows
  rows_equal <- terra::nrow(raster1) == terra::nrow(raster2)

  # Check number of columns
  cols_equal <- terra::ncol(raster1) == terra::ncol(raster2)

  # Check cell size
  cellsize_equal <- all(terra::res(raster1) == terra::res(raster2))

  # Check CRS
  crs_equal <- terra::crs(raster1) == crs(raster2)

  # Return TRUE if all checks are TRUE, otherwise FALSE
  return(rows_equal && cols_equal && cellsize_equal && crs_equal)
}


#' Plot a temperature SpatRaster using a specified color palette
#'
#' Plots a raster object using a specified color palette. This function allows for
#' the visualization of raster data with various color schemes to enhance the
#' representation of spatial patterns. The raster values are scaled by `scale_temp`
#' before plotting, and a range of palettes based on the viridis and RColorBrewer
#' packages are available. This function is useful for visualizing temperature
#' data or any other continuous variable represented in raster format.
#'
#' @param rst A raster object to be plotted.
#' @param palette A character string specifying the color palette to use.
#'                Options include "magma", "plasma", "inferno" (from the viridis
#'                package), "RdYlBu", "Spectral", "PuOr" (from the RColorBrewer
#'                package), and "heat" for heat.colors. Defaults to "magma".
#' @param scale_temp A numeric value to scale the raster values by before plotting.
#'                   Defaults to 100.
#' @param ... Additional arguments passed to the plot function.
#'
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette heat.colors
#' @importFrom viridis magma plasma inferno
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'

plot_temp_rast <- function(rst, palette = "magma",scale_temp = NULL, ...){

  if(palette == "magma")
    temp_palette <- colorRampPalette(viridis::magma(12))(100)
  if(palette == "plasma")
    temp_palette <- colorRampPalette(viridis::plasma(12))(100)
  if(palette == "inferno")
    temp_palette <- colorRampPalette(viridis::inferno(12))(100)
  if(palette == "RdYlBu")
    temp_palette <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100))
  if(palette == "Spectral")
    temp_palette <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(100))
  if(palette == "PuOr")
    temp_palette <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(100))
  if(palette == "heat")
    temp_palette <- colorRampPalette(heat.colors(12))(100)

  if(is.null(scale_temp)){
    terra::plot(rst, col = temp_palette, ...)
  }else{
    terra::plot(rst / scale_temp, col = temp_palette, ...)
  }
}

#' Upsample a SpatRaster to specified width and height dimensions
#'
#' Upsamples a `SpatRaster` object to specified dimensions using a chosen
#' interpolation method. This function is designed to increase the spatial resolution
#' of an image by interpolating new pixel values in a manner consistent with the
#' specified resampling method. This is useful in contexts where image resolution
#' needs to be increased for analysis or visualization purposes.
#'
#' @param rst A `SpatRaster` object representing the image to be upsampled.
#' @param to_height The target height (number of rows) for the upsampled image.
#'                  Defaults to 1440.
#' @param to_width The target width (number of columns) for the upsampled image.
#'                 Defaults to 1080.
#' @param method The method of interpolation to use for resampling. Supported
#'               methods include "near" for nearest neighbor interpolation,
#'               among others. Defaults to "near".
#'
#' @return A `SpatRaster` object representing the upsampled image at the specified
#'         target dimensions.
#'
#' @importFrom terra rast resample
#'
#' @export


rast_upsample <- function(rst, to_height = 1440, to_width = 1080,
                         method = "near"){

  if(!inherits(rst, "SpatRaster")){
    stop("Input object in rst must be of class SpatRaster!")
  }

  # Get image dims
  dims <- dim(rst)

  # Calculate the scale factor which will be used as the resolution
  # in terra SpatRaster
  scale_factor <- to_height / dims[1]

  # Loop across all bands
  for(i in 1:dims[3]){

    # Temp raster with resolution equal to scale and same size as input
    tmp <- terra::rast(nrows = dims[1], ncols = dims[2], res = scale_factor, crs = "",
                xmin = 0, ymin = 0, xmax = to_width, ymax = to_height)

    # Set raster values
    terra::values(tmp) <- as.numeric(terra::values(rst[[i]]))

    if(i == 1){
      trr <- tmp
    }else{
      trr <- c(trr, tmp)
    }
  }

  # Resampled raster with resolution of 1 (output)
  r <- terra::rast(nrows = to_height, ncols = to_width, res = 1, crs = "",
            xmin = 0, ymin = 0,xmax = to_width, ymax = to_height)

  # Resample raster data
  upsamp_rast <- terra::resample(trr, r, method = method)
  names(upsamp_rast) <- names(rst)
  return(upsamp_rast)
}


#' Convert a cimg image object (from imager) to SpatRaster
#'
#' Converts an RGB image in `cimg` class format to a `SpatRaster` object, optionally
#' converting pixel values to decimal format. The function creates three separate
#' raster layers for the red, green, and blue channels of the image. This is useful
#' for processing or analyzing images in a spatial context using the `terra` package.
#'
#' @param img A `cimg` class object representing the image to be converted.
#' @param as_decimal A logical flag indicating whether to convert pixel values
#'                   to decimal format (TRUE) or keep them as 8-bit integers (FALSE).
#'                   Defaults to FALSE.
#'
#' @return A `SpatRaster` object with three layers corresponding to the red,
#'         green, and blue channels of the original `cimg` image. Pixel values
#'         are either in decimal format or 8-bit integer format based on the
#'         `as_decimal` parameter.
#'
#' @importFrom terra rast values
#'
#' @export
#'

cimg_to_raster <- function(img, as_decimal = FALSE){


  if(!inherits(img,"cimg")){
    stop("Image object in img is not of class cimg! Check inputs")
  }

  dims <- dim(img)

  r <- g <- b <- terra::rast(nrows = dims[2], ncols=dims[1], res=1, crs="",
            xmin=0, ymin=0,xmax=dims[1],ymax=dims[2])

  if(!as_decimal){
    terra::values(r) <- decimal_to_8bit(as.numeric(img[,,1,1]))
    terra::values(g) <- decimal_to_8bit(as.numeric(img[,,1,2]))
    terra::values(b) <- decimal_to_8bit(as.numeric(img[,,1,3]))
  }else{
    terra::values(r) <- as.numeric(img[,,1,1])
    terra::values(g) <- as.numeric(img[,,1,2])
    terra::values(b) <- as.numeric(img[,,1,3])
  }

  rst_out <- c(r,g,b)
  names(rst_out) <- c("red","green","blue")

  return(rst_out)
}


#' Convert SpatRaster to cimg Object
#'
#' Transforms a `SpatRaster` object from the `terra` package into a `cimg`
#' object from the `imager` package. This conversion is useful for image processing
#' and analysis workflows that require the use of `imager` functions.
#'
#' @param x A `SpatRaster` object to be converted.
#' @param interpolation An integer indicating the interpolation method to use when
#' rotating the image.Zero (0) for nearest-neighbor interpolation, one (1) for linear.
#'
#' @return A `cimg` object corresponding to the input `SpatRaster`.
#'
#' @importFrom terra as.array
#' @importFrom imager as.cimg imrotate mirror
#'
#' @export
#'
#'

raster_to_cimg <- function(x, interpolation = 0){

  # Check if the object is a SpatRaster
  if(!inherits(x, "SpatRaster")){
    stop("x must be a SpatRaster object from terra package")
  }

  # Convert from SpatRaster to array
  arr_img <- terra::as.array(x)

  # Use the conversion from array to cimg
  img <- suppressWarnings(
    imager::as.cimg(arr_img))

  # Because the array is transposed we need to re-align the image
  # Rotate 90 degrees
  img <- imager::imrotate(img, angle = 90, interpolation = interpolation)
  # Mirror along the x-axis
  img <- imager::mirror(img, "x")

  return(img)
}

