

#' Apply Sobel filter for edge detection on RGB SpatRaster object
#'
#' This function applies the Sobel filter to an RGB SpatRaster for edge detection.
#' It first converts the RGB raster into a grayscale image by averaging the RGB
#' values and scaling to a [0, 1] range. Then, it applies Sobel operators (Gx and Gy)
#' for horizontal and vertical edge detection, respectively. The magnitude of the
#' gradient is computed, representing the edge strength at each pixel.
#'
#' @param rgb_rast An RGB SpatRaster object on which the Sobel filter will be applied.
#'                 The raster is expected to have three layers corresponding to
#'                 the RGB channels.
#'
#' @return A SpatRaster object representing the magnitude of the edges detected in the
#'         original RGB raster. This raster has values ranging from 0 to a maximum
#'         gradient magnitude, indicating the strength of edges at each pixel.
#'
#' @importFrom terra app focal
#'
#' @export
#'

sobel_filter <- function(rgb_rast) {

  # Use a single average image from RGB and convert to decimal grayscale
  image <- terra::app(rgb_rast, fun = "mean", na.rm = TRUE) / 255

  # Define Sobel kernels for horizontal and vertical detection
  Gx <- matrix(c(-1, 0, 1,
                 -2, 0, 2,
                 -1, 0, 1), nrow = 3, byrow = TRUE)
  Gy <- matrix(c(-1, -2, -1,
                 0,  0,  0,
                 1,  2,  1), nrow = 3, byrow = TRUE)

  # Apply the Sobel operator (convolution with the kernels)
  image_gx <- terra::focal(image, w = Gx, fun = sum, na.rm = TRUE)
  image_gy <- terra::focal(image, w = Gy, fun = sum, na.rm = TRUE)

  # Calculate the gradient magnitude
  edge_magnitude <- sqrt(image_gx^2 + image_gy^2)

  # Return the edge magnitude raster
  return(edge_magnitude)
}


#' Apply a Laplacian filter for edge detection on a RGB SpatRaster object
#'
#' Applies the Laplacian filter to an RGB SpatRaster for edge detection, supporting
#' both simple and diagonal in the kernel. The function first converts
#' the RGB raster to a grayscale image by averaging the RGB values and scaling to
#' a [0, 1] range. Depending on the `type` parameter, a Laplacian kernel without
#' diagonal elements (option: 'simple') or with diagonal elements ('diagonal') is applied. This
#' highlights edges within the image by enhancing regions of rapid intensity change.
#'
#' @param rgb_rast An RGB SpatRaster object on which the Laplacian filter will be applied.
#'                 The raster is expected to have three layers corresponding to the
#'                 RGB channels. This is converted to a grayscale image before
#'                 applying the Laplacian filter.
#' @param type A character string specifying the type of Laplacian kernel to use.
#'             "simple" uses a kernel that focuses on vertical and horizontal
#'             neighbors, while "diagonal" includes diagonal neighbors as well.
#'             Defaults to "simple". Any other value will result in an error.
#'
#' @return A SpatRaster object representing the filtered image, highlighting the edges
#'         detected in the original RGB raster. This image emphasizes areas of rapid
#'         intensity change, corresponding to edges.
#'
#' @importFrom terra app focal
#'
#' @export
#'

laplacian_filter <- function(rgb_rast, type="simple") {

  # Use a single average image from RGB and convert to decimal grayscale
  image <- terra::app(rgb_rast, fun = mean, na.rm = TRUE) / 255

  if(type == "simple"){
    # Define a Laplacian kernel (without diagonal consideration)
    L <- matrix(c(0, -1, 0,
                  -1, 4, -1,
                  0, -1, 0), nrow = 3, byrow = TRUE)
  }else if(type == "diagonal"){
    # Define a Laplacian kernel (with diagonal consideration)
    L <- matrix(c(-1, -1, -1,
                  -1,  8, -1,
                  -1, -1, -1
    ), nrow = 3, byrow = TRUE)
  }else{
    stop("Invalid option in type")
  }

  # Apply the Laplacian kernel using the focal function
  # This operation highlights edges in the image
  filtered_image <- terra::focal(image, w = L, fun = sum, na.rm = TRUE)

  # Return the filtered image
  return(filtered_image)
}

