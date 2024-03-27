
#' Extract the name of the last directory in a path
#'
#' Given a file path, this function extracts the name of the last directory
#' in that path. The path is split based on the "/" separator, and the last
#' component is returned, which is expected to be the directory name.
#'
#' @param x A character string representing the file path.
#'
#' @return A character string containing the name of the last directory
#'         in the given path.
#'
#'
#' @export


get_sample_dir_name <- function(x){
  fn <- unlist(strsplit(x,"/"),use.names = FALSE)
  return(fn[length(fn)])
}

#' Convert a decimal value to an 8-Bit Number
#'
#' Takes a decimal value (assumed to be between 0 and 1) and converts it
#' to an 8-bit number, ranging from 0 to 255. The input is scaled by 255
#' and rounded to produce the 8-bit representation.
#'
#' @param decimalValue A numeric value between 0 and 1 representing the decimal
#'                     value to be converted. Handling of values outside this
#'                     range is not specified.
#'
#' @return An integer between 0 and 255 representing the 8-bit
#'         version of the input decimal value.
#'
#'
#' @export


decimal_to_8bit <- function(decimalValue) {
  # Convert the decimal value to an 8-bit number (0 to 255)
  bitValue <- round(decimalValue * 255)
  return(bitValue)
}


#' Format duration (difftime) into a readable string
#'
#' Converts a duration object into a human-readable string format, showing
#' the duration in minutes and seconds. This function takes a `difftime`
#' object, converts it into seconds, and then formats it as "X minutes Y seconds".
#' This can be particularly useful for displaying time durations in a user-friendly manner.
#'
#' @param duration A `difftime` object representing the duration to be formatted.
#'                 The duration should be a positive value.
#'
#' @return A character string describing the duration in terms of minutes and
#'         seconds, formatted as "X minutes Y seconds".
#'
#' @export
#'

pretty_duration <- function(duration) {

  # Convert the difftime into seconds
  total_seconds <- as.numeric(duration, units = "secs")

  minutes <- floor(total_seconds / 60)
  seconds <- round(total_seconds %% 60)
  return(paste(minutes, "minutes", seconds, "seconds"))
}

#' Ensure trailing slash in directory path
#'
#' Checks if a given file path ends with a slash ("/"). If not, it appends one to
#' the path. This function is useful for preparing file paths for operations that
#' require directory paths to end with a slash, ensuring consistency and preventing
#' errors in file path manipulations.
#'
#' @param input_path A character string representing the file or directory path.
#'
#' @return The input path with a trailing slash, if it was not already present.
#'
#' @export

check_trailing_slash <- function(input_path) {
  if (substring(input_path, nchar(input_path), nchar(input_path)) != "/") {
    input_path <- paste0(input_path, "/")
    return(input_path)
  }else{
    return(input_path)
  }
}

