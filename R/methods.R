#' Import packages
#' 
#' @description 
#' Import a new package and default packages
#'
#' @param package specific package to be loaded
#' @export
#' 
LoadPackage <- function(package=FALSE){
  
  if(package==FALSE){
    library(easyPubMed)
    library(stringr) 
    library(caret)
    library(ChemmineR)
    library(ChemmineOB)
    } else {
      library(package)}
  
  
  
}



#' @title Make result folder
#' 
#' @description 
#' Creation of new folder 
#'
#' @param DesiredFolderName character. Name of folder will contain the results
#' @return the output_folder path that will contain the date as suffix.
#' 
#' @export
makeResultsDirectory <- function(DesiredFolderName){
  
  dir.create('./Results/')
  day_date = Sys.Date()
  dir.create(paste0('./Results/',DesiredFolderName,'_', day_date))
  dir.create(paste0('./Results/',DesiredFolderName, '_', day_date,'/sdf/'))
  output_folder = paste0('./Results/',DesiredFolderName, '_', day_date,'/')
  return(output_folder)
  
}



#' Drug Information file reading
#' 
#' @description 
#' Read the file containing drug information provided as input by the user
#' 
#' @param filename_info name of the file containing drug information provided as input by the user (e.g. data_example.conf)
#' @return DrugInfo_list object containing Trade/Common name, CID,	Target,	Map,	Path,	ATC (when there are multiple items you have to separate by a double pipe '||')
#' 
#' @export
readDrugInfo <- function(filename_info){
  table_info <- read.table(filename_info, sep ='\t', header = TRUE)
  # Drug List initialization
  Info_list = lapply(1:5, function(x) NULL)
  names(Info_list) = colnames(table_info)[2:6]
  DrugInfo_list = lapply(1:dim(table_info)[1], function(x) Info_list)
  names(DrugInfo_list) = as.character(table_info[,1])
  for (compound in names(DrugInfo_list)){
    Info_list = DrugInfo_list[[compound]]
    for (info in names(Info_list)){
      print(info)
      row_ind = which(as.character(table_info[,1])==compound)
      if (info=='Path' | info == 'Maps' | info =='Target'){
        kegg_data = as.character(table_info[row_ind, info])
        kegg_data = unlist(strsplit(kegg_data, '[||]'))
        kegg_data <- kegg_data[which(kegg_data != "")]
        DrugInfo_list[[compound]][[info]] = kegg_data
      } else if (info == 'ATC'){
        atc_data = as.character(table_info[row_ind, info])
        atc_data = unlist(strsplit(atc_data, '-'))
        atc_data <- atc_data[which(atc_data != "")]
        DrugInfo_list[[compound]][[info]] = atc_data
      } else {
        DrugInfo_list[[compound]][[info]] = as.character(table_info[row_ind,info])
      }
    }
  }

return(DrugInfo_list)  
}


#' Get sdf file from PubChem and conversion of the sdf in Atom Pair descriptors
#'
#' Connection tu PubChem to download the sdf file
#' @param DrugInfo_list list containing the drug information contained in configuration file
#' @param result_folder character folder name that will be contained all output files (see makeResultsDirectory function)
#' @return DrugInfo_list the input list will be updated adding the Atom pair descriptor in APsdf key
#' 
#' @export
#' 
get_sdf <- function(DrugInfo_list, result_folder){
  for (compound in names(DrugInfo_list)){
    cid = DrugInfo_list[[compound]][["CID"]]
    URL_home = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'
    url = paste0(URL_home ,'cid/',cid , '/SDF')
    namefile = paste0(compound, '_', cid,'.sdf')
    # Download sdf file from PubChem
    result_folder_sdf = paste0(result_folder,'sdf/')
    page  = tryCatch({download.file(url,destfile=paste0(result_folder_sdf ,namefile))}, warning=function(e) print(' The connection or downloaded failed. Check the input CID number in data_example.txt file. If the problem persists try again later or download manually the sdf file. '),
       error = function(e) print(' The connection or downloaded failed. Check the input CID number in data_example.txt file.\n If the problem persists try again later or download manually the sdf file.'))
    
  }
  cat(' The sdf file of',compound, ' has been downloaded correctly\n')
  cat(' Atom Pair Computation\n')
  sdf_directory = paste0(result_folder, 'sdf/')
  sdf_list_file = list.files(sdf_directory )
  for (file_sdf in sdf_list_file){
    drug_sdf = read.SDFset(paste0(sdf_directory, file_sdf) )
    compound = strsplit( file_sdf,'_')[[1]][1]
    apset_drug = sdf2ap(drug_sdf)
    apset_drug@ID = compound
    DrugInfo_list[[compound]][['APsdf']] = apset_drug
  }
  
  return(DrugInfo_list)
}


#'get PMID file from PubMed
#'
# Connection tu PubMed to download the PMIDs
#'@param DrugInfo_list list containing the drug information contained in configuration file
#'@return DrugInfo_list the input list will be updated adding the PMIDs list of new compounds
#'@export
#'
get_pmid<- function(DrugInfo_list){
  
  new_comps = names(DrugInfo_list)
  fields = c('[Title]', '[Title/Abstract]','[MeSH Major Topic]','[MeSH Terms]','[MeSH Subheading]', '[Other Term]',  '[Pharmacological Action]','[Transliterated Title]','[Text Word]')
  PMID_newComps = list()
  ncomp=1
  for (compound in new_comps){
      cat(ncomp, compound, '\n')
      my_query <- compound
      query_match = str_match(my_query, pattern = "\\D")
      my_query = paste0('\"',tolower(compound),'\"')
      if (is.na(query_match)==FALSE){ 
        attempt=0
        while(attempt<5){
          my_full_query = paste(rep(my_query, 9), fields, sep= '', collapse=' OR ')
          # Connection to PubMed 
          my_entrez_id <- tryCatch({get_pubmed_ids(my_full_query)}, error=function(cond){print('Connection Error')})
          if ( is.list(my_entrez_id)){
            attempt=5
          } else { attempt=attempt+1
          Sys.sleep(1.5)}
        }
        my_abstracts_xml <- fetch_pubmed_data(my_entrez_id, retmax = 2005)
        my_PM_list <- articles_to_list(my_abstracts_xml)
        pmid_array = c()
        if (length(my_PM_list)>0){
          for (i in seq(1,length(my_PM_list))){
            curr_PM_record <- my_PM_list[[i]]
            lines = strsplit(curr_PM_record,'\n')
            pmid_row = lines[[1]][3]
            matches = str_match(pmid_row, pattern = ">(.*)<")
            pmid= matches[,2]
            #print(pmid)
            pmid_array = c(pmid_array, pmid)}
        }
        
      } else {
        
        pmid_array = c()		
      }
      ncomp=ncomp+1
      if (length(pmid_array)>2000){
        PMID_newComps[[compound]] = pmid_array[1:2000]
      } else {
        PMID_newComps[[compound]] = pmid_array}
  }
  
  for(compound in names(PMID_newComps)){
    
    DrugInfo_list[[compound]][['PMID']] = PMID_newComps[[compound]]
  }
 
  return(DrugInfo_list)
}

#'Gene Expression data upload when it is available
#'
#' Gene expression file upload in data.frame type object
#' @param DrugInfo list containing the info of new selected compounds
#' @return DrugInfo_list the input list will be updated adding the Gene expression profiles if it is available.
#' @export
#'
geneExpr_upload <- function(DrugInfo){
  
  gene_expr_files = list.files('GeneExpr')
  i=1
  for (filename in gene_expr_files){
      drugName = strsplit( gsub('.txt','',filename),'_')[[1]][[2]]
      print(drugName)
      ExprTab = read.csv(paste0('./GeneExpr/',filename), sep='\t')
      colnames(ExprTab) = drugName
      #if (i==1){
        new_drug_expr_tab = ExprTab
        DrugInfo[[drugName]]$GeneExpr = ExprTab
     # } else if (i>1){
      #  new_drug_expr_tab = merge( new_drug_expr_tab, ExprTab, by='row.names')
      #  DrugInfo[[drugName]]$GeneExpr = ExprTab
      #}
      #i=i+1 
  }
  
  return(DrugInfo)
}

#' Conncectivity map (CMAP) data and new chemical information
#'
#' Upload of Conncectivity map data in a DrugInfoCMAP class object and upload of the information of desired compounds
#' @param filename_info name of the file containing drug information provided as input by user
#' @param OutputFolderName character desired name of the output folder that will contain the results.
#' @return list object containing two other listes:
#' @return CMAPData DrugInfoCMAP class object contaning the CMAP data information
#' @return Drug_info_list list object containing the information of the new compounds provided in the configuration file (see readDrugInfo function)
#' @export
UploadAllData <- function(filename_info, OutputFolderName){
  
    LoadPackage()
    # Loading new drug data
    result_folder = makeResultsDirectory(OutputFolderName)
    Drug_info_list = readDrugInfo(filename_info)
    Drug_info_list = get_sdf(Drug_info_list, result_folder)
    Drug_info_list = get_pmid(Drug_info_list)
    Drug_info_list = geneExpr_upload(Drug_info_list)
    
    # Data CMAP 
    load('./data/apsetCMAP.rda')
    apsetCMAP_list = list(apsetCMAP)
    load('./data/ATC_Cmap.rda')
    load('./data/geneExprCMAP_1.rda')
    load('./data/geneExprCMAP_2.rda')
    load('./data/geneExprCMAP_3.rda')
    load('./data/geneExprCMAP_4.rda')
    load('./data/PMID_Cmap.rda')
    load('./data/keggDCmap.rda')
    geneExprCMAP12 = cbind(geneExprCMAP1,geneExprCMAP2)
    geneExprCMAP34 = cbind(geneExprCMAP3,geneExprCMAP4)
    geneExprCMAP = cbind(geneExprCMAP12, geneExprCMAP34)
    MAPS = unlist(split(as.character(keggDCmap$STR_MAP_NEW), keggDCmap$Compound))
    TARGET = unlist(split(as.character(keggDCmap$TARGET), keggDCmap$Compound))
    PATH = unlist(split(as.character(keggDCmap$NEWPATHWAY), keggDCmap$Compound))
    CMAPData = DrugInfoCMAP(APsdf=apsetCMAP_list, GeneExpr=geneExprCMAP, 
                                pmid=PMID_Cmap, Target= TARGET, Path=PATH, 
                                Maps = MAPS, ATC= ATC_Cmap)
    return(list('InfoNewDrug' = Drug_info_list, 'CMAP' = CMAPData))
  
}


#'Order
#'
#' Gene expression file upload in data.frame type object
#' @param CompA character. One or a list of compounds.
#' @param CompB character. List of compounds
#' @return comp_pairs 
#' @export
OrderCompound <- function(CompA,CompB){
  
  if(length(CompA)==1){
    CompA = rep(CompA, length(CompB))
  }
    order_name = CompA < CompB
    index_right = which(order_name==TRUE)
    index_other = which(order_name==FALSE)
    comp_pairs = rep('NA', length(CompA))
    comp_pairs[index_right] = paste(CompA[index_right], CompB[index_right], sep='_')
    comp_pairs[index_other] = paste(CompB[index_other], CompA[index_other], sep='_')
  
  return(comp_pairs)
}

