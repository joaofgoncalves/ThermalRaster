% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TR_raster_aux_funs.R
\name{cimg_to_raster}
\alias{cimg_to_raster}
\title{Convert a cimg image object (from imager) to SpatRaster}
\usage{
cimg_to_raster(img, as_decimal = FALSE)
}
\arguments{
\item{img}{A `cimg` class object representing the image to be converted.}

\item{as_decimal}{A logical flag indicating whether to convert pixel values
to decimal format (TRUE) or keep them as 8-bit integers (FALSE).
Defaults to FALSE.}
}
\value{
A `SpatRaster` object with three layers corresponding to the red,
        green, and blue channels of the original `cimg` image. Pixel values
        are either in decimal format or 8-bit integer format based on the
        `as_decimal` parameter.
}
\description{
Converts an RGB image in `cimg` class format to a `SpatRaster` object, optionally
converting pixel values to decimal format. The function creates three separate
raster layers for the red, green, and blue channels of the image. This is useful
for processing or analyzing images in a spatial context using the `terra` package.
}
