#' @title An S4 class to represent a Drug information set
#' 
#' @description
#' This class represents Drug information set in a S4 DrugInfoCMAP object.
#' 
#' @field APsdf APset. It should contains the atom pairs for each CMAP drug
#' @field GeneExpr data.frame. It should contain the gene expression profile for each CMAP drug
#' @field pmid list. It should contain the set of PMID in which the compound name occured.
#' @field target character. It is an array containing the targets of each CMAP drug
#' @field path character. It is an array containing the path of each CMAP drug
#' @field maps character. It is an array containing the maps of each CMAP drug
#' @field ATC data frame. It is a data.frame containing the ATC code divided by each drug classification level.
#' 
#' @export DrugInfoCMAP
#' @exportClass DrugInfoCMAP
DrugInfoCMAP = setClass('DrugInfoCMAP',
                     slots=c(APsdf='list', GeneExpr='data.frame', 
                             pmid='list', Target='character', Path='character', 
                             Maps = 'character', ATC= 'data.frame')
)
