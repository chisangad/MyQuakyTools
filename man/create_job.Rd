% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myFunctions.R
\name{create_job}
\alias{create_job}
\title{Create batch job script}
\usage{
create_job(
  batch.cmd = "qsub",
  time = 1,
  nodes = 1,
  cores = 2,
  mem = "4gb",
  process.name = "job",
  params = NULL,
  path.script,
  email.results = F
)
}
\value{
Returns a character string of the batch submission script
}
\description{
Function to create a batch submission script
}
