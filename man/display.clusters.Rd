% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Analysis.R
\name{display.clusters}
\alias{display.clusters}
\title{Display Clusters}
\usage{
display.clusters(clust, Net, zoom = NA, type = 1, cols = c("yellow",
  "black", "red", "green"), col = cols[type + 1], ...)
}
\arguments{
\item{clust}{a matrix which is the output1 of \code{cluster.network}}

\item{Net}{an igraph object}

\item{zoom}{an integer to define the size of plot or NA, to plot all clusters}

\item{type}{0 = "yellow", 1 = "black", 2="red", 3 = "green"}

\item{cols}{a list of 4 colors if you want to change the default colors}

\item{col}{color to use with zoom}

\item{...}{Additional plotting parameters}
}
\description{
Given a cluster matrix, see \code{cluster.network} and a network
display the cluster matrix with the following colors for the position (i,j):
}
\details{
- Yellow if the protein and protein i doesn't belongs to cluster of protein j and
proteins don't interact in the network

-  Black if the protein and protein i belongs to cluster of protein j but
proteins don't interact in the network

- Red if the protein and protein i doesn't belongs to cluster of protein j but
proteins interact in the network

-  Green if the protein and protein i belongs to cluster of protein j and
proteins interact in the network
}

