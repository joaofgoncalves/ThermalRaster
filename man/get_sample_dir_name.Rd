% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TR_aux_funs.R
\name{get_sample_dir_name}
\alias{get_sample_dir_name}
\title{Extract the name of the last directory in a path}
\usage{
get_sample_dir_name(x)
}
\arguments{
\item{x}{A character string representing the file path.}
}
\value{
A character string containing the name of the last directory
        in the given path.
}
\description{
Given a file path, this function extracts the name of the last directory
in that path. The path is split based on the "/" separator, and the last
component is returned, which is expected to be the directory name.
}
