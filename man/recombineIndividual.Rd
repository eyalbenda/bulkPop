% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IndividualFunctions.R
\name{recombineIndividual}
\alias{recombineIndividual}
\title{Simulate a meiotic recombination in an individual
Workhorse function that takes an individual and performs recombinations based on the genetic map. shouldn't usually be called directly in most cases.}
\usage{
recombineIndividual(individual, genome, Nrecombs = 1)
}
\arguments{
\item{individual}{an individual to recombine}

\item{genome}{a genome object corresponding to the individual}

\item{Nrecombs}{number of recombinations per chromosome}
}
\value{
a new recombined individual
}
\description{
Simulate a meiotic recombination in an individual
Workhorse function that takes an individual and performs recombinations based on the genetic map. shouldn't usually be called directly in most cases.
}
