
### Phenotypic Drug Similarity Prediction (PDSP)



Here we propose a model to predict phenotypic drug similairity. It was developed applying a Supervised Machine Learning approach, where the "target variable" (drug similarity) was obtained by extracting the Bioactivity data from PubChem database and by computing a phenotypic drug similarity, and the input data consisted of seven chemical features.
The seven chemical features are the data which the user has to provide as input in the configuration file (.conf).
The starting drug dataset is the set of small molecules provided by Connectvity Map (CMAP) ( Lamb et al. 2006). From the same resource the gene expression data was extracted, which is another input feature of the model. 
In the following table, all features used to develop the model are summarized together with the download date.

| Feature| Compound Number|  Source|Download Date|
|------|---------|----------|---------|
|Gene Expression|1309|Connectivity Map (Broad Institute)|Performed in 2006|
|Chemical Structure|1267|PubChem|Oct 2017|
|ATC|761|KEGG drug|Jan 2018|
|Target| 607 out of 878|KEGG drug|Jan 2018|
|Maps|469 out of 878|KEGG drug|Jan 2018|
|Path|486 out of 878|KEGG drug|Jan 2018|
|PMID|1200|PubMed|March 2018| 
|Bioassay|1231|PubChem Bioassay|Nov 2017|


The input features are the Gene Expression profile for a specific drug of interest, if it is available; the Chemical structure; the target, biological and synthesis pathway; the lieterature co-occurrence; the Anatomical therapeutic chemical classification (ATC) code. 
For each chemical feature a drug similarity score was computed yielding eight feature-based similarity scores. 
In the following table the metrics used to compute the similarity score are summarized:

|Feature|Compound Pair Number|Data Type|Similarity score method|
|-------|---------|-----------|--------|
|Gene Expression|856086|log2fc|GE profile|correlation|
|Chemical Structure|802011|Atom Pair descriptors|Tanimoto Coefficient|
|ATC|289180|character code|Sharing atc level|
|Target|348418|Target list|Jaccard Index|
|Maps|314878|Map list|Jaccard Index|
|Path|308367|Pathway list|Jaccard Index|
|PMID|719400|PMID list|Jaccard Index|
|Bioassay|59914 (filtered from 757065)|Bioactivity binary array|Manual Score|


The Bioassay data were just used to develop the Predictive model applying a neural network algorithm as machine learning method.
The Predcited model was developed using a small subset of compound pairs as training dataset, then the predicted model was applied to all other compound pairs to predict phenotyipc similarity score. The already predicted drug similarities are provided in data folder and they will be uploaded during novel prediction or they can be just indipendently visualized.
The code provided below allows to predict the phenotypic similarity score between the new added compound(s) and the CMAP compounds, and among new added compouds provided by the user, if they are more than one.


### Installation

The `PDSP` depends on some packages: ChemmineR, ChemmineOB, caret, stringr, easyPubMed. Please make sure that they are installed.

The latest version of `PDSP` can be installed from Github. To install from Gitbub, you need the devtools package.

``` r
install.packages('devtools')
devtools::install_github("Elisa89m/PDSP")
```

Then you should have the package installed.



### Configuration file and input data

The configuration file must have the following fields:

- Trade/Common name: The name of the drug
- CID: The Chemical Identifier of PubChem
- Target: a list of targets delimited by '|||'. The target name must be written in the KEGG format (e.g. BCR-ABL [HSA:25] [KO:K06619])
- Maps: The synthesis pathway identifer. The map must be written in the KEGG format (e.g. map07045)
- Path: The drug involved pathway delimited by '|||'. The pathway name must be written in the KEGG format  (e.g.hsa04010||hsa04370||hsa05200)
- ATC: The Anatomical therapeutic chemical classification provided by World Health Organization.

An example of configuration file is provided as data_example.conf
The Gene Expression data should be provided for a better prediction, but it could be omitted. Generally, GEO database provides many drug perturbation experiment, an example is the dataset GSE116436 that contains several drug perturbations in different cancer cell lines. The differential gene expression analysis could be performed using the GEO2R extention provide by GEO itself. The user must arrange the results keeping only the logFC column from table results and saving the file in GeneExpr folder. Three different drug perturbations are provided as an example in GeneExpr folder. 



### Example Prediction Workflow


#### Initialization workspace and loading data

To start the prediction at first you need to define the name of configuration file (i.e. the name of file that you should have previuosly compiled) and the name of output folder, which will be labelled with the current date.

``` r
filename_info = 'data_example.conf'
OutputFolderName = 'Test'
```

The `UploadAllData` function executes the reading of the configuration file and builds a R list object which will contain the new added drugs and the CMAP compound chemical information. The new added drug data are collected in "InfoNewDrug" item of the list and for each new drug there are the following fileds:
- CID
- Target
- Map
- Path
- PMID
- APsdf
- PMID
- Gene Expression

see the readDrugInfo function for more detail of output values.

``` r
DrugList_Info = UploadAllData(filename_info, OutputFolderName)
```


#### Chemical score similarity compuation

The `GetAllSimilarityScoreCMAP` function allows to compute all similarity scores for each features. The returned object ScoreList is a R list containing the seven chemical feature, which will be used as input to the Predictive Model. 

``` r
ScoreList = GetAllSimilarityScoresCMAP(DrugList_Info$CMAP,DrugList_Info$InfoNewDrug)
```

#### Loading predictive model and novel drug prediction


`LoadModelData` function loads the Predictive Model


``` r
ModelData = LoadModelData()
```

Predictive Model is extracted from ModelData list and passed to the `PredictNewDrugSimilarity` function, which requires also the ScoreList obtained in the previous step.


``` r
# Load Predictive Model
ModelResult = ModelData$PredictiveModel
# Run drug similarity prediction 
ResultPredicted  = PredictNewDrugSimilarity(ScoreList, ModelData)
```

``` r
# Change the name of the first column
colnames(ResultPredicted)[1] = 'Pairs'
# Get and split the compound pair column to separate the drug name in two different columns
pairs = as.character(ResultPredicted$Pairs)
pairs_split = strsplit( pairs,'_')
compA = unlist(lapply(pairs_split, '[[', 1))
compB = unlist(lapply(pairs_split, '[[', 2))

# Add the compound columns of each pair 
ResultPredicted$CompA = compA
ResultPredicted$CompB = compB
# Change the order of resulting table to make more clear.
ResultPredicted = ResultPredicted[,c(10:11,1,2:9)]
# The scores were rounded before to write in file
ResultPredicted[,4:11] = round(ResultPredicted[,4:11],4)
```

``` r
# The table will be saved in a file in in the output folder with the name `ResultPredicted.tab`.
write.table(ResultPredicted, paste0(output_folder,'ResultPredicted.tab'), sep= '\t', row.names=FALSE,quote=FALSE)

```

The `ResultPredicted` table contains the compound pair predicted scores obtained by the combination of new input drug and the CMAP drug dataset. To allow further consideration the similaritie scores based on other chemical features are reported in the same table.



### Inferred phenotypic scores (IPS) visualization of CMAP dataset 

In case the user doesn't want provide new drug to be predicted, (s)he can visualize the Predicted score of the CAMP dataset.


``` r
# Load Model Data which contains the CMAPSimilarityScores
ModelData = LoadModelData()
CMAPSimilarityScore = ModelData$CMAPSimilarityScores
colnames(CMAPSimilarityScore )[2] = 'PredictedScore'
CMAPSimilarityScore = CMAPSimilarityScore[,c(11:12,1:10)]

```



