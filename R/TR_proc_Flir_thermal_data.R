
#' Convert FLIR raw thermal image to a temperature SpatRaster object
#'
#' Converts raw thermal imagery from FLIR cameras into a temperature raster using
#' metadata extracted with EXIFtool. This function reads a FLIR JPEG image, extracts
#' necessary calibration parameters from its metadata, and computes the temperature
#' values for each pixel. The temperature data can be scaled by a factor and written
#' to a file as a `SpatRaster` object. This is useful for thermal image analysis,
#' allowing for the quantitative assessment of temperature distributions in the captured
#' scene.
#'
#' @param img_path Path to the input FLIR image file.
#' @param exiftool_path Path to the EXIFtool executable, used for extracting metadata
#'                      from the FLIR image.
#' @param out_path (Optional) Output file path (usually with .tif extension) where the
#'                 temperature raster file will be saved.
#' @param scale_factor (Optional) A numeric value to scale the temperature values.
#'                     Useful for converting temperature to a specific range or format
#'                     such as conversion to integer numbers for saving space.
#'                     If NULL, no scaling is applied.
#' @param datatype Data type for saving the raster file. Supported types include
#'                 "INT1S", "INT2S", "INT4S", "INT1U", "INT2U", "INT4U", "FLT4S",
#'                 and "FLT8S". Defaults to "FLT4S".
#' @param overwrite Logical flag indicating whether to overwrite the output raster file.
#'                  Defaults to TRUE.
#'
#' @return A `SpatRaster` object representing the temperature distribution in the
#'         input FLIR image. If `out_path` is TRUE, the raster is also saved
#'         to the specified output folder usually as a (Geo)TIFF file.
#'
#' @importFrom terra rast writeRaster
#' @importFrom Thermimage readflirJPG flirsettings raw2temp
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'

flir_raw_to_thermal <- function(img_path,
                                exiftool_path,
                                out_path = NULL,
                                scale_factor = NULL,
                                datatype = "FLT4S",
                                #write_raster = TRUE,
                                overwrite = TRUE){

  if(!dir.exists(exiftool_path)){
    stop("Could not find the EXIFtool executable!")
  }

  if(!file.exists(img_path)){
    stop("Could not find the input image file: \n", img_path,"!\n")
  }

  # if(write_raster && is.null(out_folder)){
  #   stop("out_folder cannot be NULL when writing data to file (write_raster=TRUE)")
  # }

  exiftool_path <- check_trailing_slash(exiftool_path)

  img  <- Thermimage::readflirJPG(img_path, exiftoolpath = exiftool_path)
  metadata <- Thermimage::flirsettings(img_path, exiftoolpath = exiftool_path, camvals="")

  # Extract needed parameters from Flir image metadata
  #
  ObjectEmissivity <-  metadata$Info$Emissivity              # Image Saved Emissivity - should be ~0.95 or 0.96
  #dateOriginal <-      metadata$Dates$DateTimeOriginal       # Original date/time extracted from file
  #dateModif <-   metadata$Dates$FileModificationDateTime     # Modification date/time extracted from file
  PlanckR1 <-    metadata$Info$PlanckR1                      # Planck R1 constant for camera
  PlanckB <-     metadata$Info$PlanckB                       # Planck B constant for camera
  PlanckF <-     metadata$Info$PlanckF                       # Planck F constant for camera
  PlanckO <-     metadata$Info$PlanckO                       # Planck O constant for camera
  PlanckR2 <-    metadata$Info$PlanckR2                      # Planck R2 constant for camera
  OD <-          metadata$Info$ObjectDistance                # object distance in metres
  FD <-          metadata$Info$FocusDistance                 # focus distance in metres
  ReflT <-       metadata$Info$ReflectedApparentTemperature  # Reflected apparent temperature
  AtmosT <-      metadata$Info$AtmosphericTemperature        # Atmospheric temperature
  IRWinT <-      metadata$Info$IRWindowTemperature           # IR Window Temperature
  IRWinTran <-   metadata$Info$IRWindowTransmission          # IR Window transparency
  RH <-          metadata$Info$RelativeHumidity              # Relative Humidity
  h <-           metadata$Info$RawThermalImageHeight         # sensor height (i.e. image height)
  w <-           metadata$Info$RawThermalImageWidth          # sensor width (i.e. image width)

  # Make temperature image
  temperature <- Thermimage::raw2temp(img, ObjectEmissivity, OD, ReflT, AtmosT, IRWinT, IRWinTran, RH,
                                      PlanckR1, PlanckB, PlanckF, PlanckO, PlanckR2)

  dims <- dim(temperature)

  # Transpose the output temperature matrix, convert to numeric to allow conversion
  # to a SpatRaster object
  if(!is.null(scale_factor)){
    # multiply by the scale factor and round the values to zero decimal plates
    v <- round(as.numeric(t(temperature))*scale_factor)
  }else{
    v <- as.numeric(t(temperature))
  }

  # Generate a SpatRaster object to hold the data
  out_raster <- terra::rast(nrows = dims[1], ncols = dims[2], vals = v, crs = NA,
                           resolution = 1, xmin = 0, xmax = dims[2], ymin = 0, ymax = dims[1])

  # Write SpatRaster data to file
  if(!is.null(out_path)){

    # if(!dir.exists(out_folder)){
    #   stop("Could not find the output dir! Check input path in out_folder")
    # }

    # Assemble output file path
    #out_fpath <- paste(out_folder,"/",file_path_sans_ext(basename(img_path)),"_temp.tif",sep="")
    terra::writeRaster(out_raster,
                       out_path,
                       overwrite = overwrite,
                       datatype = datatype)
  }

  return(out_raster)
}

#' Convert all FLIR radiometric images in folder to temperature SpatRasters
#'
#' Batch processes radiometric images from a specified input folder, converting each
#' to a temperature raster using metadata for calibration. This function is designed
#' to work with FLIR radiometric JPEG images, utilizing an external EXIFtool for
#' metadata extraction and a specified function for the conversion process. The
#' output temperature rasters are saved in a designated output folder.
#'
#' @param input_folder The folder path containing the input radiometric images.
#' @param recursive Recursively list all files in `input_folder`.
#' @param output_folder The folder path where the output temperature raster files
#'                      will be saved. Must be different from the `input_folder`.
#' @param input_file_type The file extension of the input images (default "jpg").
#'                        This parameter allows for filtering the files to be processed.
#' @param exiftool_path Path to the EXIFtool executable, required for metadata
#'                      extraction from the FLIR images.
#' @param ... Other parameters passed to `flir_raw_to_thermal` function.
#'
#' @return The function does not return a value but writes the output temperature
#'         raster files to the `output_folder`. Files will be placed in `input_folder`
#'         with a suffix _temp in (Geo)TIFF format if `output_folder` is NULL. Or placed
#'         in the defined `output_folder`.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
convert_all_to_thermal <- function(input_folder, recursive = FALSE, output_folder=NULL,
                                   input_file_type = "jpg", exiftool_path, ...){

  # if(input_folder == output_folder){
  #   stop("Input and output folders must be different!")
  # }
  if(!dir.exists(input_folder)){
    stop("input_folder does not exists!")
  }
  if(!dir.exists(output_folder)){
    stop("output_folder does not exists!")
  }

  fl <- list.files(input_folder,pattern = paste(".",input_file_type,"$",sep=""),
                   full.names = TRUE, recursive = recursive)
  #print(fl)

  cat("\n\nConvert radiometric image files to temperature:\n\n")
  pb <- utils::txtProgressBar(min = 1, max = length(fl), style = 3)

  for(i in 1:length(fl)){

    if(is.null(output_folder)){
      out_path <- paste(check_trailing_slash(input_folder),"/",
                        file_path_sans_ext(basename(fl[i])),"_temp.tif",sep="")
    }else{
      out_path <- paste(check_trailing_slash(output_folder),"/",
                        file_path_sans_ext(basename(fl[i])),"_temp.tif",sep="")
    }

    flir_raw_to_thermal(img_path = fl[i],
                          out_path = out_path,
                          exiftool_path = exiftool_path, ...)

    utils::setTxtProgressBar(pb,i)
  }
}

