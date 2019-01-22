#' Predictive Model data and CMAP similarity data
#'
#' Load Predictive Model data and CMAP compound pair similarity  
#' @param SimilairtyScoreList object obtained by GetAllSimilarityScoresCMAP function
#' @return list object containing the Predicitve Model and a data.frame containing all CMAP chemical similarity scores and the phenotypic predicted score
#' 
#' @export
#' 
LoadModelData <- function(){
  
  
  load('./data/CMAP_DrugSimilarityPrediction.rda')
  load('./data/PredictiveModel.rda')
  return(list('PredictiveModel'= PredictiveModel, 'CMAPSimilarityScores'= CMAP_SimilarityScores))
}


#' Predicition of new drug similarity scores
#'
#'
#' @param SimilairtyScoreList object obtained by GetAllSimilarityScoresCMAP function
#' @return InputDataAndPredicted data.frame object that contain the predicted score and the input data passed to the Predictive Model.
#' @export
#' 
PredictNewDrugSimilarity <- function(SimilairtyScoreList, ModelData){
  
  i=0
  for (feature in names(SimilairtyScoreList)){
      print(feature)
        i=i+1
        if (i==1){
          DrugInputData = SimilairtyScoreList[[feature]]
          DrugInputData$CompPair = rownames(DrugInputData)
        } else {
          DrugInputDataFeature = SimilairtyScoreList[[feature]]
          DrugInputDataFeature$CompPair = rownames(DrugInputDataFeature)
          DrugInputData = merge(DrugInputData, DrugInputDataFeature[,c('CompPair', feature)], by='CompPair', all=TRUE)
        }
  } 
  
  
  PredicitveModel = ModelData$PredictiveModel
  CMAP_InputScores = ModelData$CMAPSimilarityScores
  CMAP_InputScores$CompPair = CMAP_InputScores$Pair
  # mean_value to fill in the missing value in the input data that have to be passed to predict function  
  
  ALL_InputData = rbind(DrugInputData, CMAP_InputScores[,c('CompPair', as.character(names(SimilairtyScoreList)))])
  mean_value = colMeans(ALL_InputData[,names(SimilairtyScoreList)], na.rm=TRUE)
  
  rownames(DrugInputData) = DrugInputData$CompPair
  DrugInputData$CompPair = NULL
  
  
  DrugInputData_fillin = DrugInputData
  # Prediction Computation of the new Compound Pair 
  for (feature in names(mean_value)){
    DrugInputData_fillin[which(is.na(DrugInputData_fillin[,feature])),feature] <- mean_value[feature]
    
  }


  
  new_scores = predict(PredicitveModel$finalModel, DrugInputData_fillin)
  colnames(new_scores) = 'PredictedScore'
  new_scores = as.data.frame(new_scores)
 
  InputDataAndPredicted = merge(DrugInputData, new_scores, by='row.names')
  InputDataAndPredicted = InputDataAndPredicted[order(-InputDataAndPredicted$PredictedScore),,drop=FALSE]
  return(InputDataAndPredicted)
  
}
