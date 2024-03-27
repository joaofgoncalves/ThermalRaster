
#' Convert FLIR camera thumbnail image to SpatRaster
#'
#' Converts a thumbnail image extracted from a FLIR camera's radiometric JPEG
#' to a SpatRster object.
#'
#' @param img_path Path to the FLIR image.
#' @param as_decimal A logical flag indicating whether to convert the pixel values
#'                   to decimal format. Defaults to FALSE.
#'
#' @return A `SpatRaster` object representing the FLIR thumbnail image.
#'
#' @importFrom imager load.image
#'
#' @export
#'

flir_thumbnail_to_rast <- function(img_path, as_decimal = FALSE){

  if(!is.character(img_path)){
    stop("img_path must be a character/string object")
  }

  if(!file.exists(img_path)){
    stop("File in img_path was not found. Check the input.")
  }

  img <- try(imager::load.image(img_path))

  if(inherits(img,"try-error")){
    stop("Could not read the image data")
  }

  rst_out <- cimg_to_raster(img, as_decimal = as_decimal)

  return(rst_out)
}


#' Extract a full uncropped RGB image from FLIR Radiometric JPEG
#'
#' Extracts the full RGB image embedded within a FLIR radiometric JPEG file,
#' saving it to a specified output path. This function uses EXIFtool to extract
#' the embedded RGB image. Additionally, there's an option to load the extracted
#' image directly into R as either an `imager` package `cimg` object or a
#' `terra` package `SpatRaster` object, depending on the `load_as` parameter.
#'
#' @param img_path Path to the input FLIR radiometric JPEG file.
#' @param out_path Path where the extracted full RGB image will be saved.
#' @param replace Logical flag indicating whether an existing output file should
#'                be replaced. Defaults to FALSE.
#' @param exiftool_path Path to the EXIFtool executable.
#' @param load_as (Optional) Specifies how to load the extracted image into R: "cimg"
#'                for an `imager` package object, "SpatRaster" for a `terra` package
#'                object, or NULL to only save the file without loading. Defaults to NULL.
#'
#' @return Returns TRUE if the extraction and optional loading were successful, FALSE
#'         otherwise. If `load_as` is specified, returns the loaded image object.
#'
#' @importFrom imager load.image
#' @importFrom terra rast
#'
#'

flir_get_full_rgb <- function(img_path, out_path, replace = FALSE,
                              exiftool_path, load_as = NULL){


  if(!is.character(img_path)){
    stop("img_path must be a character/string object")
  }

  if(!file.exists(img_path)){
    stop("File in img_path was not found. Check the input.")
  }

  if(file.exists(out_path) && !replace){
    stop("File already exists but replace option option is set to FALSE.
         Do you want to replace this file?")
  }

  exiftool_path <- check_trailing_slash(exiftool_path)

  cmd <- paste0(exiftool_path,"exiftool -b -EmbeddedImage ",
                img_path," > ",out_path)

  out <- try(shell(cmd))

  if(inherits(out,"try-error")){
    message("An error occurred while extracting the image with EXIFtools")
    return(FALSE)
  }else{

    if(is.null(load_as)){
      return(TRUE)
    }else{
      if(load_as == "cimg"){
        return(imager::load.image(out_path))
      }else if(load_as == "SpatRaster"){
        return(suppressWarnings(terra::rast(out_path)))
      }else{
        stop("Invalid option in load_as!")
      }
    }
  }
}

#' Converts a FLIR image into an RGB SpatRaster
#'
#' Extracts the full RGB image from a FLIR radiometric JPEG, optionally crops it
#' based on FLIR metadata to match the thermal image's dimensions, and converts
#' it to a `SpatRaster` object. This function can save the extracted and cropped
#' RGB image as a TIFF file. Cropping is useful to align the RGB image with the
#' corresponding thermal image, enabling direct comparison or fusion of data from
#' both sources.
#'
#' @param img_path Path to the input FLIR radiometric JPEG file.
#' @param exiftool_path Path to the EXIFtool executable, including trailing slash.
#' @param out_path Optional. The location/path where to save the full RGB image. IF `NULL`
#' the directory of the input will be used with a file name suffix _full_RGB.jpg. This image
#' will be created by EXIFtool by extracting the RGB image from metadata. If a file
#' path is specified it must be of jpg extension.
#' @param crop Logical flag indicating whether to crop the extracted RGB image to
#'             match the dimensions of the thermal image. Defaults to TRUE.
#' @param out_path_crop Path to where the cropped imaged will be saved. This will use
#' `terra::writeRaster` to write the image data.
#' @param overwrite Overwrite image data in `terra::writeRaster`? Default: TRUE
#'
#' @return A `SpatRaster` object representing the extracted image. If `crop=TRUE` the cropped
#'         RGB image will be returned instead of the full RGB image.
#'
#' @importFrom terra writeRaster rast
#' @importFrom tools file_path_sans_ext
#' @importFrom Thermimage flirsettings
#' @importFrom imager resize as.cimg
#'
#' @export
#'

flir_rgb_to_rast <- function(img_path, exiftool_path, out_path = NULL,
                             crop = TRUE, out_path_crop = NULL, overwrite = TRUE){

  # Assemble the output paths for the full RGB image and
  # the cropped image (if crop = TRUE)
  if(crop){
    if(is.null(out_path)){
      fn <- tools::file_path_sans_ext(img_path)
      output_img_path_crop <- paste0(fn,"_crop_RGB.tif")
      output_img_path_full <- paste0(fn,"_full_RGB.jpg")
    }else{
      # fn <- tools::file_path_sans_ext(basename(img_path))
      # output_img_path_crop <- paste0(check_trailing_slash(output_folder)
      #                                ,fn,"_crop_RGB.tif")
      # output_img_path_full <- paste0(check_trailing_slash(output_folder)
      #                                ,fn,"_full_RGB.tif")

      output_img_path_full <- out_path
      output_img_path_crop <- out_path_crop
    }
  }else{
    if(is.null(out_path)){
      # fn <- tools::file_path_sans_ext(img_path)
      # output_img_path_full <- paste0(fn,"_full_RGB.jpg")
      fn <- tools::file_path_sans_ext(img_path)
      output_img_path_full <- paste0(fn,"_full_RGB.jpg")
    }else{
      # fn <- tools::file_path_sans_ext(basename(img_path))
      # output_img_path_crop <- paste0(check_trailing_slash(output_folder)
      #                                ,fn,"_crop_RGB.tif")
      # output_img_path_full <- paste0(check_trailing_slash(output_folder)
      #                                ,fn,"_full_RGB.jpg")
      output_img_path_full <- out_path
      output_img_path_crop <- out_path_crop
    }
  }

  exiftool_path <- check_trailing_slash(exiftool_path)

  img <- flir_get_full_rgb(img_path,
                           output_img_path_full,
                           replace = TRUE,
                           exiftool_path,
                           load_as = "cimg")

  if(crop){

    metadata <- try(Thermimage::flirsettings(img_path,
                                             exiftoolpath = exiftool_path,
                                             camvals=""))

    if(inherits(metadata,"try-error")){
      stop("An error occurred while getting metadata from the Flir image in img_path!")
    }

    # Extract the metadata parameters needed to perform the crop
    scale <- metadata$Info$Real2IR
    offset_x <- metadata$Info$OffsetX
    offset_y <- metadata$Info$OffsetY

    # Output size of the cropped image
    width <- metadata$Info$RawThermalImageWidth
    height <- metadata$Info$RawThermalImageHeight

    # Output size of the scaled image
    scale_x <- floor(width * scale)
    scale_y <- floor(height * scale)

    # Resize the original RGB image into the scaled version
    # By default uses near-neigbour interpolation
    img_resized <- try(imager::resize(img,
                              size_x = scale_x,
                              size_y = scale_y,
                              interpolation_type = 1))

    if(inherits(img_resized,"try-error")){
      stop("An error occurred while performing image resize with imager package!")
    }

    #plot(img_resized)

    # Calculate the offset image center
    ct = c(floor(scale_x / 2), floor(scale_y / 2)) # Raw image center/ without offsets
    cx = ct[1] + offset_x * scale # Image center with scaled offsets for x
    cy = ct[2] + offset_y * scale # Image center with scaled offsets for y

    # Based on the image center calculate the image ranges
    # which will be used to crop
    crop_ranges <- c(
      floor(cx - width / 2) + 1,
      floor(cx + width / 2),
      floor(cy - height / 2) + 1,
      floor(cy + height / 2)
    )

    img_cropped <-
      suppressWarnings(
        suppressMessages(
          imager::as.cimg(img_resized[crop_ranges[1]:crop_ranges[2],  # Along the x-axis
                              crop_ranges[3]:crop_ranges[4],,]) # Along the y-axis
        ))

    # Convert the image into a SpatRaster Object to be used by terra
    rst_out <- cimg_to_raster(img_cropped, as_decimal = FALSE)

    terra::writeRaster(rst_out, filename = output_img_path_crop,
                         NAflag = NA, datatype = 'INT1U', overwrite = overwrite)

    return(rst_out)
  }else{

    # Return the output full image
    return(cimg_to_raster(img, as_decimal = FALSE))
  }
}

