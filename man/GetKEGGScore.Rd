% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimilarityScoreFunctions.R
\name{GetKEGGScore}
\alias{GetKEGGScore}
\title{KEGG data score similarity between a new compound and connecitivity map compounds}
\usage{
GetKEGGScore(newDrugInfo_comp, CMAPData)
}
\arguments{
\item{newDrugInfo_comp}{list object containing the information of a new compound}

\item{CMAPData}{Drug Info CMAP Class object containing the CMAP compounds information}
}
\value{
list object containing three data.frame with the similarity scores based on the Maps, Path and Target Kegg drug information
}
\description{
Compute the KEGG data score between the selected compound and the compounds of CMAP
}
