% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Incompatibilities.R
\name{zeelPeelSelection}
\alias{zeelPeelSelection}
\title{Zeel peel selection
A sleection based on the peel-1/zeel-1 incompatibility from Seidel et al 2008.}
\usage{
zeelPeelSelection(genome, parent1, parent2, sexes, chosenFrom, zeelPeel)
}
\arguments{
\item{genome}{A genome object. corresponds to both parents}

\item{parent1}{the first parent of the individual. The function assumes the parents have already undergone recombination (meiosis)}

\item{parent2}{the second parent of the individual. The function assumes the parents have already undergone recombination (meiosis)}

\item{sexes}{A vector of length 2 with the sexes of the parents (0 is male, 1 is female)}

\item{chosenFrom}{which haplotype was chosen from each parent.}

\item{zeelPeel}{a list with two elements. marker: the marker which is being selected. kill: which genotype is killed (0 or 1)}
}
\description{
Zeel peel selection
A sleection based on the peel-1/zeel-1 incompatibility from Seidel et al 2008.
}
