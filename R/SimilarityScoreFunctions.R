#' Structure Similarity compound-CMAPset
#' 
#' @description 
#' Compute the structure similarity score between the selected compound and the compounds of CMAP
#' @param apset_comp APset object containing the Atom Pair descriptors of new compound(s)
#' @param apsetCMAP APset object containing all Atom Pair descriptors of CMAP compounds
#' @return data.frame (Nx1) with similarity score between new_compound and apsetCMAP
#' @export 
AP_simCMAP <-function(apset_comp, apsetCMAP){
  
  # cmp.similarity compute the similarity between atom pairs
  AP_comp_CMAP = lapply(seq_along(cid(apsetCMAP)), function(x) cmp.similarity(apset_comp,apsetCMAP[x],mode = 1))# Tanimoto
  names(AP_comp_CMAP) = apsetCMAP@ID
  new_comp = apset_comp@ID
  comp_pairs = OrderCompound(new_comp, apsetCMAP@ID)
  StructureDf = data.frame('Structure'= unlist(AP_comp_CMAP), 'row.names' =  comp_pairs)
  
  return(StructureDf)
  
  
}


#' Structure Similarity between two new compound
#' 
#' @description 
#' Compute the structure similarity score between the selected compound and the compounds of CMAP
#' @param apset_compA  APset object containing the Atom Pair descriptors of new compound
#' @param apset_compB  APset object containing the Atom Pair descriptors of new compound 
#' @return data.frame containing the similarity score between two selected compounds
#' @export 
AP_sim2comp <-function(apset_compA, apset_compB){
   pair_score = cmp.similarity(apset_compA,apset_compB,mode = 1)
   if(apset_compA@ID < apset_compB@ID){
     comp_pair_name = paste0(apset_compA@ID,'_', apset_compB@ID)
   } else {
     comp_pair_name = paste0(apset_compB@ID, '_', apset_compA@ID)
   }
   pair_score_df = data.frame('Structure'= pair_score, 'row.names'=comp_pair_name)
  return(pair_score_df)
}

#' co-occurrence similarity between two new compound and connectivity map compounds
#' 
#' @description 
#' Compute the co-occurrence score between the selected compound and the compounds of CMAP
#' @param new_comp_name a character object containing the name associated to the new compound
#' @param PMID_newComps character array containing the PMIDs associated to the new compound
#' @param PMID_cmap list object containing the PMID associated to each compound of CMAP
#' @return data frame containing the compound pair similarity scores between the new compound(s) and all CMAP compounds
#' @export 
GetLitteratureCo_occ <-function(new_comp_name,PMID_newComps,PMID_cmap){
  
  commonPMID_newcomp_CMAP = unlist(lapply(names(PMID_cmap), function(comp) length(intersect(PMID_newComps, PMID_cmap[[comp]]))))
  names(commonPMID_newcomp_CMAP) = names(PMID_cmap)
  unionPMID_newcomp_CMAP = unlist(lapply(names(PMID_cmap), function(comp) length(unique(PMID_newComps, PMID_cmap[[comp]]))))
  names(unionPMID_newcomp_CMAP) = names(PMID_cmap)
  JI_score_newcomp_CMAP =  commonPMID_newcomp_CMAP/unionPMID_newcomp_CMAP
  comp_pairs = OrderCompound(new_comp_name,names(PMID_cmap))
  PMID_df = data.frame('PMID'= JI_score_newcomp_CMAP, 'row.names' =  comp_pairs)
  return(PMID_df)
  
}

#' Gene Expression similarity between two new compound and connecitivity map compounds
#' 
#' @description 
#' Compute the Correlation score between the selected compound and the compounds of CMAP
#' @param GeneExprNewDrug data.frame containing the gene expression of a new compound
#' @param GeneExprCMAP data.frame containing the gene expression profiles of CMAP 
#' @return data.frame object containing the gene expression correlation between the new compound and all CMAP compounds
#' @export 
GetGeneExprScore <-function(GeneExprNewDrug,GeneExprCMAP){
  
  new_comp = colnames(GeneExprNewDrug)
  common_gene = intersect(rownames(GeneExprNewDrug), rownames(GeneExprCMAP))
  CorrelationScore = unlist(lapply(colnames(GeneExprCMAP), function(comp) cor(as.numeric(GeneExprNewDrug[common_gene,new_comp]), as.numeric(GeneExprCMAP[common_gene,comp]))))
  names(CorrelationScore) = colnames(GeneExprCMAP)
  comp_pairs = OrderCompound(new_comp,colnames(GeneExprCMAP))
  CorrelationScoreScaled = (CorrelationScore+ abs(-1))/(abs(-1-1)) 
  GeneExprDF = data.frame('GeneExpr'= CorrelationScoreScaled, 'row.names' =  comp_pairs)
  return(GeneExprDF)
}


#' ATC code score similarity between a new compound and a set of ATC codes
#' 
#' Compute the ATC score between a selected compound and one or more compounds
#' @param newDrugInfo_comp list object containing the information of a new compound
#' @param ATCcode_df data.frame containing the ATC codes for one or more compounds
#' @return data.frame containing the compound pair similarity scores between the ATC code of the new compound and the ATC codes provided in ATCcode_df
#' @export 
GetATCcodeScore <-function(newDrugInfo_comp,ATCcode_df){
      result_table_freq = data.frame()
      new_comp = names(newDrugInfo_comp)
      ATCcodeNewDrug = newDrugInfo_comp[[new_comp]][['ATC']]
      new_atc_df = data.frame()
      i=0
      for (atc_code in ATCcodeNewDrug){
        i=i+1
        row_atc = data.frame('NAME'= new_comp,'L1'= substr(atc_code,1,1), 'L2'= substr(atc_code,1,3), 'L3'=substr(atc_code,1,4), 'L4'=substr(atc_code,1,5), 'L5'=atc_code, 'Number'= i)
        new_atc_df= rbind(new_atc_df, row_atc)
      }
    
      scores = list('L1'=0.2, 'L2'=0.4, 'L3'=0.6, 'L4'=0.8)
      level_type = 1
      for(level_atc in c('L1', 'L2', 'L3', 'L4')){
        atc_number= 0
        level_score_tab = data.frame()
        for(atc_codes_level in as.character(new_atc_df[,level_atc] )){
          atc_number = atc_number+1
           ATCcode_subset =  ATCcode_df[which(as.character(ATCcode_df[,level_atc] ) == atc_codes_level),]
          
           if(nrow(ATCcode_subset)>=1){
               L_tab_number = lapply(1:nrow(ATCcode_subset), function(row_ind) data.frame('CompA'=paste(ATCcode_subset$NAME[row_ind],ATCcode_subset$Number[row_ind], sep='_'),
                                                                                'CompB'=paste(new_comp, atc_number,sep='_'), 'LScore'=scores[level_atc], 'L_ext'= atc_codes_level ))
           
               L_tab = do.call("rbind", L_tab_number)
               colnames(L_tab)[3:4]=c(level_atc, paste0(level_atc, '_ext'))
               level_score_tab = rbind(level_score_tab, L_tab)
        }
        }
        if(level_type==1){
          
          final_table = level_score_tab
        } else {
          if(nrow(level_score_tab)>=1){
          final_table = merge(final_table, level_score_tab, by=c('CompA', 'CompB'), all=TRUE)
          }
        }
        level_type = level_type+1
      }
      
      if(nrow(final_table)>=1){
          compA = strsplit(as.character(final_table$CompA), '_')
          compA = unlist(lapply(compA, function(item) item[[1]][1] ))
          compB = strsplit(as.character(final_table$CompB), '_')
          compB = unlist(lapply(compB, function(item) item[[1]][1] ))
          final_table$CompPair = OrderCompound(compA, compB)
          max_value = apply(final_table[,names(scores)], 1, 'max', na.rm=TRUE)
          final_table$MAX = max_value
          result_table = aggregate(MAX~CompPair, final_table, 'sum')
          freq = as.data.frame(table(final_table$CompPair))
          colnames(freq)[1] = 'CompPair'
          result_table_freq = merge(result_table, freq, by='CompPair')
          result_table_freq$mean = result_table_freq$MAX/result_table_freq$Freq
      }
      
      
      return(result_table_freq)
}

#' KEGG data score similarity between a new compound and connecitivity map compounds
#' 
#' Compute the KEGG data score between the selected compound and the compounds of CMAP
#' @param newDrugInfo_comp list object containing the information of a new compound
#' @param CMAPData Drug Info CMAP Class object containing the CMAP compounds information
#' @return list object containing three data.frame with the similarity scores based on the Maps, Path and Target Kegg drug information 
#' @export 
GetKEGGScore <-function(newDrugInfo_comp,CMAPData){
  
  new_comp = names(newDrugInfo_comp)
  cmap_maps = CMAPData@Maps
  cmap_maps = cmap_maps[which(grepl('^map', cmap_maps))]
  cmap_maps = strsplit(cmap_maps, '\\|\\|')
  
  cmap_target  = strsplit(CMAPData@Target, '\\|\\|')
  cmap_path  = strsplit(CMAPData@Path, '\\|\\|')
  KEGG_cmap = list('Maps'=cmap_maps, 'Target'= cmap_target, 'Path'= cmap_path)
  
  results = list()
  for( Kegg_info in c('Maps','Target', 'Path')){
     
     CMAP_data = KEGG_cmap[[Kegg_info]]
     common_info_CMAP_newDrug = unlist(lapply(names(CMAP_data), function(cmap_comp) length(intersect(CMAP_data[[cmap_comp]], newDrugInfo_comp[[new_comp]][[Kegg_info]]))))
     union_info_CMAP_newDrug = unlist(lapply(names(CMAP_data), function(cmap_comp) length(unique(c(CMAP_data[[cmap_comp]], newDrugInfo_comp[[new_comp]][[Kegg_info]])))))
     JI_kegg_info_CMAP_newDrug = common_info_CMAP_newDrug/union_info_CMAP_newDrug
     comp_pairs = OrderCompound(new_comp,names(CMAP_data))
     KEGG_info_DF = data.frame( JI_kegg_info_CMAP_newDrug, 'row.names' =  comp_pairs)
     colnames(KEGG_info_DF) = Kegg_info
    results[[paste0('DF_', Kegg_info)]]= na.omit(KEGG_info_DF)
  }
  
  return(results)
  
  
}

#' All similarity score computation among new compounds
#' 
#' Compute all type similarity scores among new compounds
#' @param newDrugInfo list containing the drug information of the new selected compounds 
#' @param new_comp_A the name of one drug 
#' @param new_comp_B the name of another drug
#' @return list object containing the similairty scores based on all chemical features among two new compounds
#' @export 
GetAllSimilarityScores <- function(newDrugInfo, new_comp_A, new_comp_B){
           
  
             # Structure
            structure_score = AP_sim2comp(newDrugInfo[[new_comp_A]][['APsdf']],newDrugInfo[[new_comp_B]][['APsdf']] )
            
            # PMID
            JI_PMID = length(intersect(newDrugInfo[[new_comp_A]][['PMID']], newDrugInfo[[new_comp_B]][['PMID']]))/length(unique(c(newDrugInfo[[new_comp_A]][['PMID']],newDrugInfo[[new_comp_B]][['PMID']])))
            comp_pair = OrderCompound(new_comp_A,new_comp_B)
            JI_PMID_comps = data.frame('PMID'=JI_PMID, row.names=comp_pair)
            
            # GeneExpr
            CorrelationScore_df = data.frame()
            if ('GeneExpr' %in% names(newDrugInfo[[new_comp_A]]) & 'GeneExpr' %in% names(newDrugInfo[[new_comp_B]])) {
              gene_A = rownames(newDrugInfo[[new_comp_A]][['GeneExpr']])
              gene_B = rownames(newDrugInfo[[new_comp_B]][['GeneExpr']])
              common_genes  = intersect(gene_A, gene_B)
              CorrelationScore = cor(newDrugInfo[[new_comp_A]][['GeneExpr']][which(gene_A %in% common_genes),1], newDrugInfo[[new_comp_B]][['GeneExpr']][which(gene_A %in% common_genes),1])
              CorrelationScore_df = data.frame('GeneExpr'=CorrelationScore, row.names = comp_pair)
            }
            
            # KEGG 
            Kegg_results_new_drug = list()
            for( Kegg_info in c('Maps','Target', 'Path')){
              
              KeggData_A = newDrugInfo[[new_comp_A]][[Kegg_info]]
              KeggData_B = newDrugInfo[[new_comp_B]][[Kegg_info]]
              JI_kegg_info_newDrug = length(intersect(KeggData_A, KeggData_B))/length(unique(c(KeggData_A, KeggData_B)))
              KEGG_info_DF = data.frame( JI_kegg_info_newDrug, row.names =comp_pair)
              colnames(KEGG_info_DF) = Kegg_info
              Kegg_results_new_drug[[paste0('DF_', Kegg_info)]]= KEGG_info_DF
            }
           
            
            # ATC
            
            ATCcodeNewDrug = newDrugInfo[[new_comp_B]][['ATC']]
            new_atc_df_B = data.frame()
            i=0
            for (atc_code in ATCcodeNewDrug){
              i=i+1
              row_atc = data.frame('NAME'= new_comp_B,'L1'= substr(atc_code,1,1), 'L2'= substr(atc_code,1,3), 'L3'=substr(atc_code,1,4), 'L4'=substr(atc_code,1,5), 'L5'=atc_code, 'Number'= i)
              new_atc_df_B= rbind(new_atc_df_B, row_atc)
            }
            
            Score_twoComps = GetATCcodeScore(newDrugInfo[new_comp_A], new_atc_df_B)
          
            return(list('Structure'=structure_score, 'GeneExpr'=CorrelationScore_df, 'PMID'= JI_PMID_comps, 'Maps'=Kegg_results_new_drug$DF_Maps, 
                        'Target'=Kegg_results_new_drug$DF_Target, 'Path'=Kegg_results_new_drug$DF_Path,'ATC'= Score_twoComps))
            
            
}




#' All similarity score computation between two new compound and connecitivity map compounds
#' 
#' Compute the KEGG data score between the selected compound and the compounds of CMAP
#' @param CMAPData Drug Info CMAP Class object containing the CMAP compounds information
#' @param newDrugInfo list containing the drug information of the new selected compounds
#' @return list object containing the similairty scores based on all chemical features between new compounds and CMAP compounds
#' @export

GetAllSimilarityScoresCMAP <- function(CMAPData, newDrugInfo){
  
  StructureDF = data.frame()
  PMID_DF = data.frame()
  GeneExprDF = data.frame()
  Path_DF = data.frame()
  Maps_DF = data.frame()
  Target_DF = data.frame()
  ATC_DF = data.frame()
  apsetCMAP = CMAPData@APsdf[[1]]
  PMID_cmap = CMAPData@pmid
  for (new_comp in names(newDrugInfo)){
      print(new_comp)
      # Structure similarity Score computation
      apset_new_comp = newDrugInfo[[new_comp]][['APsdf']]
      StructureDFComp = AP_simCMAP(apset_new_comp,apsetCMAP)
      StructureDF = rbind(StructureDF,StructureDFComp)
      
      # Co-occurrence similarity Score computation
      PMID_newComps = newDrugInfo[[new_comp]]['PMID']
      PMID_DFComp =  GetLitteratureCo_occ(new_comp,PMID_newComps,PMID_cmap)
      PMID_DF = rbind(PMID_DF,PMID_DFComp)
      
      # Gene Expression Similarity Score computation
      if('GeneExpr' %in% names(newDrugInfo[[new_comp]])){
         GeneExprCMAP = CMAPData@GeneExpr
         GeneExprNewDrug = newDrugInfo[[new_comp]][['GeneExpr']]
         GeneExprDFComp = GetGeneExprScore(GeneExprNewDrug, GeneExprCMAP)
         GeneExprDF = rbind(GeneExprDF, GeneExprDFComp)
      }
      
      # KEGG data similarity Score computation
      KEGG_DF_list =  GetKEGGScore(newDrugInfo[new_comp], CMAPData)
      Path_DF = rbind(Path_DF,KEGG_DF_list$DF_Path)
      Maps_DF = rbind(Maps_DF,KEGG_DF_list$DF_Maps)
      Target_DF = rbind(Target_DF,KEGG_DF_list$DF_Target)
      
      # ATC similarity Score computation
      ATCcodeCMAP = CMAPData@ATC
      ATC_score_DF =  GetATCcodeScore(newDrugInfo[new_comp],ATCcodeCMAP)
      ATC_DF = rbind(ATC_DF, ATC_score_DF)
  }
  
        ################
        # Score Between new compounds
        new_comps= names(newDrugInfo)
        number_comps = length(new_comps)
        i=0
        for (new_comp_A in new_comps){
              i=i+1
              if(i!=number_comps){
                  for (new_comp_B in new_comps[(i+1):number_comps]){
                    
                    list_scores_AB  = GetAllSimilarityScores(newDrugInfo, new_comp_A ,new_comp_B)
                    
                    StructureDF = rbind(StructureDF,list_scores_AB$Structure )
                    PMID_DF = rbind(PMID_DF, list_scores_AB$PMID)
                    GeneExprDF = rbind(GeneExprDF, list_scores_AB$GeneExpr)
                    Path_DF = rbind(Path_DF, list_scores_AB$Path)
                    Maps_DF = rbind(Maps_DF, list_scores_AB$Maps)
                    Target_DF = rbind(Target_DF, list_scores_AB$Target)
                    ATC_DF = rbind(ATC_DF, list_scores_AB$ATC )
           
                
              }}
              
        }
  
        colnames(ATC_DF)[4] = 'ATC'
        rownames(ATC_DF) = ATC_DF$CompPair
  
  return(list('Structure'=StructureDF, 'GeneExpr'=GeneExprDF, 'PMID'= PMID_DF, 'Maps'=Maps_DF, 'Target'=Target_DF, 'Path'=Path_DF,'ATC'=ATC_DF))
}
  

  


