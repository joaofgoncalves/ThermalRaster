% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TR_proc_Flir_thermal_data.R
\name{convert_all_to_thermal}
\alias{convert_all_to_thermal}
\title{Convert all FLIR radiometric images in folder to temperature SpatRasters}
\usage{
convert_all_to_thermal(
  input_folder,
  recursive = FALSE,
  output_folder = NULL,
  input_file_type = "jpg",
  exiftool_path,
  ...
)
}
\arguments{
\item{input_folder}{The folder path containing the input radiometric images.}

\item{recursive}{Recursively list all files in `input_folder`.}

\item{output_folder}{The folder path where the output temperature raster files
will be saved. Must be different from the `input_folder`.}

\item{input_file_type}{The file extension of the input images (default "jpg").
This parameter allows for filtering the files to be processed.}

\item{exiftool_path}{Path to the EXIFtool executable, required for metadata
extraction from the FLIR images.}

\item{...}{Other parameters passed to `flir_raw_to_thermal` function.}
}
\value{
The function does not return a value but writes the output temperature
        raster files to the `output_folder`. Files will be placed in `input_folder`
        with a suffix _temp in (Geo)TIFF format if `output_folder` is NULL. Or placed
        in the defined `output_folder`.
}
\description{
Batch processes radiometric images from a specified input folder, converting each
to a temperature raster using metadata for calibration. This function is designed
to work with FLIR radiometric JPEG images, utilizing an external EXIFtool for
metadata extraction and a specified function for the conversion process. The
output temperature rasters are saved in a designated output folder.
}
