require(tidyverse)
require(ggthemes)
require(enrichR)
require(extrafont)
require(parallel)
require(doParallel)
require(foreach)
loadfonts()

### The function enrichment_test performs enrichment analysis starting from 3 vectors of gene names: input, annotation and background (the background is intended as the universe and, by default, is also used to filter gene names in the input and annotation) ----  


enrichment_test <- function(input_vec, # character vector, input genes
                            anno_vec, # character vector, annotation genes
                            background_vec, # character vector, background/universe genes
                            input_name="i1", # name for input list (i1 by default)
                            anno_name="a1", # name for annotation list (a1 by default)
                            anno_class="a", # name for annotation class (groups of annotations, a by default)
                            id_separator=";" # character separating overlap ids in the resulting string (; by default)
){
   
   # library(tidyverse)
   
   bg_vec <- background_vec %>% stringr::str_trim() %>% stringr::str_to_upper() %>% 
      unique() %>% stringr::str_sort() # clean background
   
   ib_vec <- input_vec %>% stringr::str_trim() %>% stringr::str_to_upper() %>% 
      unique() %>% stringr::str_sort() %>% 
      magrittr::extract(. %in% bg_vec) # clean input (overlap with bg)
   
   ab_vec <- anno_vec %>% stringr::str_trim() %>% stringr::str_to_upper() %>% 
      unique() %>% stringr::str_sort() %>% 
      magrittr::extract(. %in% bg_vec) # clean anno (overlap with bg)
   
   ov_vec <- ib_vec %>% magrittr::extract(. %in% ab_vec) # overlap between input and anno
   
   ovinput_vec <- input_vec %>% stringr::str_trim() %>% stringr::str_sort() %>% 
      unique() %>% magrittr::extract(stringr::str_to_upper(.) %in% ov_vec) # overlap, with input case
   ovinput_vec <- ovinput_vec[!duplicated(stringr::str_to_upper(ovinput_vec))] # remove of case sensitive duplicates
   
   IY<- length(ov_vec) # in input, in anno (overlap)
   IN<- length(ib_vec) - IY # in input, not in anno
   BY<- length(ab_vec) - IY # not in input, in anno
   BN<- length(bg_vec) - (IY+IN+BY) # not in input, not in anno   
   
   # fisher test, 1 tail (only enrichment)
   fisher_data <- fisher.test(matrix(c(IY, BY, IN, BN), 2, 2), alternative="greater") 
   # hypergeometric test (equivalent to Fisher 1 tail)
   #phyper(q=(IY-1),m=(BY+IY),n=(IN+BN),k=(IY+IN),lower.tail = F) 
   
   # generation of the output (a tibble with 1 row)
   out_df <- tibble("input_name"=input_name, # name of the input
                    "anno_name"=anno_name, # name of the annotation
                    "anno_class"=anno_class,
                    "overlap_size"=IY, # size of the overlap 
                    "p_value"=fisher_data$p.value, # p-value of the enrichment (Fisher test)
                    "odds_ratio"=fisher_data$estimate, # conditional MLE estimate from fisher.test function
                    "combined_score"=0, # default combined score (to match enrichr enrichment results)
                    #"odds_ratio_sample"= (IY/IN)/(BY/BN), # sample odds ratio
                    #"fold_enrichment"= (IY/(IY+IN))/((IY+BY)/(BY+IY+BY+BN)), # fold enrichment formula
                    "input_size"=IY+IN, # size of processed input 
                    "anno_size"=IY+BY, # size of processed annotation
                    "background_size"= IY+IN+BY+BN, # size of processed background
                    "overlap_input_ratio"= IY/(IY+IN), # overlap ratio with respect to input size
                    "overlap_anno_ratio"= IY/(IY+BY), # overlap ratio with respect to annotation size
                    "overlap_ids"= str_c(ovinput_vec,collapse=id_separator) # string with overlap ids
                    )
   return(out_df)
}



### The function enrichment_enrichr performs enrichment analysis with the enrichR bioconductor library, interface to the enrichR webserver ---- 

enrichment_enrichr <- function(input_vec, # vector with gene names
                               input_name="i1", # name for input list
                               dbs_vec="", # vector with database names
                               return_size_in_anno = T # use only genes with at least 1 annotation to calculate input size
){ 
   
   input_vec <- input_vec %>% stringr::str_trim() %>% stringr::str_sort() %>% 
      unique() # clean input vector
   input_vec <- input_vec[!duplicated(stringr::str_to_upper(input_vec))] # remove case duplicates
   
   #library(tidyverse)
   #library(enrichR)
   #library(plyr)
   
   # dbs_transcription <- c("ChEA_2016",
   #                        "ChEA_2022",
   #                        "ARCHS4_TFs_Coexp",
   #                        "TF_Perturbations_Followed_by_Expression",
   #                        "lncHUB_lncRNA_Co-Expression",
   #                        "TRANSFAC_and_JASPAR_PWMs",
   #                        "Epigenomics_Roadmap_HM_ChIP-seq",
   #                        "TargetScan_microRNA_2017",
   #                        "miRTarBase_2017",
   #                        "ENCODE_TF_ChIP-seq_2015",
   #                        "TF-LOF_Expression_from_GEO",
   #                        "ENCODE_Histone_Modifications_2015",
   #                        "TRRUST_Transcription_Factors_2019")
   # 
   # dbs_pathway<-c( "KEGG_2021_Human",
   #                 "KEGG_2019_Mouse",
   #                 "WikiPathways_2021_Human",
   #                 "WikiPathways_2019_Mouse",
   #                 "BioPlanet_2019",
   #                 "Reactome_2016",
   #                 "BioCarta_2016",
   #                 "Humancyc_2016",
   #                 "NCI-Nature_2016",
   #                 "Panther_2016",
   #                 "BioPlex_2017",
   #                 "huMAP",
   #                 "LINCS_L1000_Kinase_Perturbations_up",
   #                 "LINCS_L1000_Kinase_Perturbations_down",
   #                 "ARCHS4_Kinases_Coexp",
   #                 "Virus-Host_PPI_P-HIPSTer_2020",
   #                 "NURSA_Human_Endogenous_Complexome",
   #                 "Elsevier_Pathway_Collection",
   #                 "MSigDB_Hallmark_2020", # MSigDB hallmark collection
   #                 "SubCell_BarCode" # web portal to proteome subcellular localization https://www.subcellbarcode.org/
   #                 )
   # 
   # dbs_onto<-c("GO_Biological_Process_2021",
   #             "GO_Cellular_Component_2021",
   #             "GO_Molecular_Function_2021",
   #             "MGI_Mammalian_Phenotype_Level_4_2021",
   #             "Jensen_COMPARTMENTS"
   #             )
   # 
   # dbs_disease<-c( "OMIM_Disease",
   #                 "OMIM_Expanded",
   #                 "Human_Phenotype_Ontology",
   #                 "Jensen_DISEASES",
   #                 "Achilles_fitness_decrease",
   #                 "Achilles_fitness_increase",
   #                 "DSigDB",
   #                 "MSigDB_Oncogenic_Signatures",
   #                 "PheWeb_2019",
   #                 "ClinVar_2019",
   #                 "GWAS_Catalog_2019",
   #                 "DisGeNET",
   #                 "ARCHS4_IDG_Coexp",
   #                 "dbGaP",
   #                 "Rare_Diseases_GeneRIF_Gene_Lists",
   #                 "Rare_Diseases_AutoRIF_Gene_Lists",
   #                 "Rare_Diseases_AutoRIF_ARCHS4_Predictions",
   #                 "Rare_Diseases_GeneRIF_ARCHS4_Predictions",
   #                 "COVID-19_Related_Gene_Sets_2021",
   #                 "TG_GATES_2020",
   #                 "Orphanet_Augmented_2021",
   #                 "PhenGenI_Association_2021",
   #                 "HDSigDB_Human_2021",
   #                 "HDSigDB_Mouse_2021"
   #                 )
   # 
   # dbs_tissue<-c( "Jensen_TISSUES",
   #                "Human_Gene_Atlas",
   #                "Mouse_Gene_Atlas",
   #                "Allen_Brain_Atlas_up",
   #                "Allen_Brain_Atlas_down",
   #                "GTEx_Tissue_Sample_Gene_Expression_Profiles_up",
   #                "GTEx_Tissue_Sample_Gene_Expression_Profiles_down",
   #                "Cancer_Cell_Line_Encyclopedia",
   #                "NCI-60_Cancer_Cell_Lines",
   #                "Tissue_Protein_Expression_from_Human_Proteome_Map",
   #                "ARCHS4_Tissues",
   #                "ESCAPE",
   #                "CCLE_Proteomics_2020",
   #                "ProteomicsDB_2020",
   #                "SysMyo_Muscle_Gene_Sets",
   #                "Descartes_Cell_Types_and_Tissue_2021",
   #                "Allen_Brain_Atlas_10x_scRNA_2021",
   #                "Azimuth_Cell_Types_2021",
   #                "CellMarker_Augmented_2021",
   #                "PanglaoDB_Augmented_2021"
   #                )
   # 
   # dbs_misc<-c(  "Chromosome_Location_hg19",
   #               "InterPro_Domains_2019",
   #               "Pfam_Domains_2019",
   #               "Data_Acquisition_Method_Most_Popular_Genes"
   #               )
   # 
   # dbs_geo<-c("Aging_Perturbations_from_GEO_up",
   #            "Aging_Perturbations_from_GEO_down",
   #            "Disease_Perturbations_from_GEO_up",
   #            "Disease_Perturbations_from_GEO_down",
   #            "Drug_Perturbations_from_GEO_up",
   #            "Drug_Perturbations_from_GEO_down",
   #            "Single_Gene_Perturbations_from_GEO_up",
   #            "Single_Gene_Perturbations_from_GEO_down",
   #            "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
   #            "Ligand_Perturbations_from_GEO_up",
   #            "Ligand_Perturbations_from_GEO_down",
   #            "Microbe_Perturbations_from_GEO_up",
   #            "Microbe_Perturbations_from_GEO_down")
   # 
   # dbs_default<-c(dbs_transcription,dbs_pathway,dbs_onto,dbs_disease,dbs_tissue,dbs_misc,dbs_geo)
   # 
   dbs_default <- read_delim("dbs_enrichR.txt", delim = "\n", col_names = FALSE)[[1]]
   dbs <- listEnrichrDbs()
   if(!is.null(dbs_vec)){dbs_used<-intersect(dbs_vec,dbs$libraryName)
   } else {dbs_used<-intersect(dbs_default,dbs$libraryName)}
   
   enrich_list <- enrichr(input_vec, dbs_used)
   out_df <- plyr::ldply(enrich_list,.id="anno_class")
   rm(enrich_list)
   
   out_df <- out_df %>% dplyr::rename(anno_name=Term,
                                    p_value=P.value,
                                    fdr=Adjusted.P.value,
                                    odds_ratio=Odds.Ratio,
                                    combined_score=Combined.Score,
                                    enrichr_ids=Genes)  # rename columns
   
   transient <- stringr::str_split(out_df$Overlap,"/",simplify=T) # extract overlap and anno size
   out_df$overlap_size <- transient[,1] %>% as.numeric()
   out_df$anno_size <- transient[,2] %>% as.numeric()
   rm(transient)
   
   out_df <- out_df %>% dplyr::select(-c(Old.P.value,Old.Adjusted.P.value,Overlap)) # remove useless columns
   out_df$overlap_anno_ratio <- out_df$overlap_size/out_df$anno_size
   
   # sort genes and recover original gene names used for input
   transient <- out_df %>% dplyr::select(anno_class,anno_name,enrichr_ids) %>% 
      separate_rows(enrichr_ids,sep=";") %>% dplyr::arrange(anno_class,anno_name,enrichr_ids)
   
   out_df$enrichr_ids<-NULL 
   out_df$input_size <- ifelse(return_size_in_anno, 
                               length(unique(transient$enrichr_ids)), 
                               length(input_vec)) # calculate input size
   
   match_df <- data.frame("input_ids"=input_vec,
                          "enrichr_ids"=toupper(input_vec))
   transient<- left_join(transient,match_df)
   
   #any(is.na(transient$input_ids)) control, should be FALSE
   
   transient <- transient %>% dplyr::mutate(overlap_ids = ifelse(is.na(input_ids),enrichr_ids,input_ids))
   transient <- transient %>% dplyr::group_by(anno_class,anno_name) %>% 
      dplyr::summarise("overlap_ids"=str_c(overlap_ids,collapse=";")) %>% dplyr::ungroup() 
   
   out_df <- left_join(out_df,transient)
   rm(transient,match_df)
   
   out_df$overlap_input_ratio<-out_df$overlap_size/out_df$input_size
   out_df$input_name<- input_name
   
   out_df <-out_df %>% dplyr::select(input_name,anno_name,anno_class,
                                     overlap_size,p_value,fdr,odds_ratio,combined_score,
                                     input_size,anno_size,overlap_input_ratio,overlap_anno_ratio,
                                     overlap_ids) %>% as_tibble()
   # background_size is not estimated, odds ratio retrived from enrichr
   
   return(out_df)
}



### The function enrichment_lollipop creates an enrichment plot (lollipop) starting from an enrichment dataframe ---- 

enrichment_lollipop <- function(input_df, # input dataframe
                                x_col="combined_score", # column for x axis, necessary
                                y_col="anno_name", # column for y axis, necessary
                                
                                size_col="overlap_size", # column for size, optional
                                size_vec=c(1,5), # size vector c(min,max), optional
                                
                                shape_col="P<0.05", # column for shape, optional
                                shape_vec=c(16,21), # shape vector, optional

                                color_col="P<0.05", # column for color, optional
                                color_vec=c("grey25","grey50"), # color vector, optional
                                
                                fill_col="P<0.05", # column for fill color, optional
                                fill_vec=c("grey25","white"), # fill vector, optional
                                
                                text_col="overlap_size", # column for fill color, optional
                                
                                facet_col="anno_class", # column for facet, optional

                                char_max = 50, # maximum number of chars in text
                                sort_y = T # sort rows according to y?
                                ){

#library(tidyverse)
#library(ggthemes)
#library(extrafont)
#loadfonts()

# select and rename relevant columns for the plot   
   plot_df <- input_df %>% dplyr::select(x_col=all_of(x_col),
                                         y_col=all_of(y_col),
                                         size_col=any_of(size_col),
                                         shape_col=any_of(shape_col),
                                         color_col=any_of(color_col),
                                         fill_col=any_of(fill_col),
                                         text_col=any_of(text_col),
                                         facet_col=any_of(facet_col))
   
   # creation of labels, abbreviating text 
   plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
   transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
   plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
   rm(transient)

   # order rows by x, remove eventual duplicates
   if(sort_y){
      plot_df<- plot_df %>% dplyr::arrange(-x_col)
      }

   plot_df<-plot_df[!duplicated(plot_df %>% dplyr::select(any_of(c("y_col","facet_col")))),]
   plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)

# define a unique color, if absent
#if(!("col_col" %in% colnames(plot_df))){
#   plot_df$col_col<-""
   #col_vec=c("grey25")
#}

lp <- ggplot(plot_df, aes(x_col, y_col)) +
   
   geom_segment(aes(x = 0, y = y_col, xend = x_col, yend = y_col))+
   theme_tufte(base_family = "Arial Narrow")+
   scale_color_manual(values = "grey25", guide = "none")+
   theme(axis.title.y = element_blank(),axis.ticks.y= element_blank())+
   xlab(x_col)

if("size_col" %in% colnames(plot_df)){
   lp <- lp + geom_point(aes(size=size_col)) + scale_size(name=size_col,range = size_vec)
} else {
   lp <- lp + geom_point(size=2)      
}

if("shape_col" %in% colnames(plot_df)){
   lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
}

if("color_col" %in% colnames(plot_df)){
   lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col, values= color_vec) # drop = FALSE
}

if("fill_col" %in% colnames(plot_df)){
   lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
}

if("text_col" %in% colnames(plot_df)){
   lp <- lp + geom_text(aes(label=text_col, hjust=(-0.5)),fontface="italic", size=3, show_guide = F) +
      scale_x_continuous(expand = expansion(mult = c(.01, .1)))
}
  
if("facet_col" %in% colnames(plot_df)){
lp <- lp + ggforce::facet_col(vars(facet_col), scales = "free_y", space = "free") + 
   theme(strip.text = element_text(margin = margin(1,1,1,1), face="bold"))
}

return(lp)

}


### The function enrichment_b2b_lollipop creates a back to back enrichment plot (lollipop) starting from an enrichment dataframe with 2 analyses ----

enrichment_b2b_lollipop <- function(input_df, # input dataframe
                                  
                                bb_col="input_name", # column used for the b2b (2 levels)
                                x_col="combined_score", # column for x axis, necessary
                                y_col="anno_name", # column for y axis, necessary
                                
                                size_col="overlap_size", # column for size, optional
                                size_vec=c(1,5), # size vector c(min,max), optional
                                
                                shape_col="P<0.05", # column for shape, optional
                                shape_vec=c(16,21), # shape vector, optional
                                
                                color_col="input_name", # column for color, optional
                                color_vec=c("#16697A","#DB6400"), # color vector, optional
                                
                                fill_col="P<0.05", # column for fill color, optional
                                fill_vec=c("grey25","white"), # fill vector, optional
                                
                                text_col="overlap_size", # column for fill color, optional
                                
                                facet_col="anno_class", # column for facet, optional
                                
                                char_max = 50, # maximum number of chars in text
                                sort_y = T, # sort rows according to y?
                                shift_mult = 1, # space left for middle labels (by default, if left to 1, equal to the space dedicated to the lollipop on each side)
                                break_vec=seq(0,200,50) # specification of x breaks
){

plot_df <- input_df %>% dplyr::select(bb_col=all_of(bb_col),
                                      x_col=all_of(x_col),
                                      y_col=all_of(y_col),
                                            size_col=any_of(size_col),
                                            shape_col=any_of(shape_col),
                                            color_col=any_of(color_col),
                                            fill_col=any_of(fill_col),
                                            text_col=any_of(text_col),
                                            facet_col=any_of(facet_col))

# creation of labels, abbreviating text 
plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
rm(transient)

# factorize bb_col
if(!is.factor(plot_df$bb_col)){
plot_df$bb_col<-factor(plot_df$bb_col)
}

# sort y_col
if(sort_y){
   plot_df<- plot_df %>% dplyr::arrange(-x_col)
}

# remove eventual duplicates and factorize y_col 
plot_df<-plot_df[!duplicated(plot_df %>% dplyr::select(any_of(c("bb_col","y_col","facet_col")))),]
plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)

# define shift (if 1, roughly 1/3 of the plot)
shift<-(max(plot_df$x_col)/2)*shift_mult

f_left<- levels(plot_df$bb_col)[1]
f_right<- levels(plot_df$bb_col)[2]

lp<- ggplot(plot_df, aes(y = y_col))+
   geom_linerange(data = plot_df[plot_df$bb_col==f_left,],
               aes(xmin = -shift, xmax = -shift-x_col))+
   geom_linerange(data = plot_df[plot_df$bb_col==f_right,],
               aes(xmin = shift, xmax = shift+x_col))+
   geom_text(data = plot_df[plot_df$bb_col==f_left,],
          aes(y = y_col, x = 0, label = y_col),
          inherit.aes = F,family="Arial Narrow",size=3)+
   scale_x_continuous(limits =c((-shift-max(plot_df$x_col)),(shift+max(plot_df$x_col))),
                   breaks = c(rev(-break_vec)-shift, break_vec+shift),
                   labels = c(as.character(rev(break_vec)),as.character(break_vec)))+
   theme_tufte(base_family = "Arial Narrow")+
   theme(axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "top")+
   labs(x = x_col)


if("size_col" %in% colnames(plot_df)){
   lp <- lp + 
      geom_point(data = plot_df[plot_df$bb_col==f_left,],
                 aes(x = -shift-x_col,size=size_col))+
      geom_point(data = plot_df[plot_df$bb_col==f_right,],
                 aes(x = shift+x_col,size=size_col))+
      scale_size(name=size_col,range = size_vec)
} else {
   lp <- lp + 
      geom_point(data = plot_df[plot_df$bb_col==f_left,],
                 aes(x = -shift-x_col),size=2)+
      geom_point(data = plot_df[plot_df$bb_col==f_right,],
                 aes(x = shift+x_col),size=2)
}

if("text_col" %in% colnames(plot_df)){
   lp <- lp + geom_text(data = plot_df[plot_df$bb_col==f_left,],
                        aes(x = -shift-x_col, label=text_col, hjust=(1.1)),
                        fontface="italic", size=3, show_guide = F) +
      geom_text(data = plot_df[plot_df$bb_col==f_right,],
                aes(x = shift+x_col, label=text_col, hjust=(-0.1)),
                fontface="italic", size=3, show_guide = F) 
      scale_x_continuous(limits =c((-shift-max(plot_df$x_col)),(shift+max(plot_df$x_col))),
                         breaks = c(rev(-break_vec)-shift, break_vec+shift),
                         labels = c(as.character(rev(break_vec)),as.character(break_vec)),
                         expand = expansion(mult = c(.3, .3)))
}

if("shape_col" %in% colnames(plot_df)){
   lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
}

if("color_col" %in% colnames(plot_df)){
   lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col, values= color_vec) # drop = FALSE
}

if("fill_col" %in% colnames(plot_df)){
   lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
}

if("facet_col" %in% colnames(plot_df)){
   lp <- lp + ggforce::facet_col(vars(facet_col), scales = "free_y", space = "free") + 
      theme(strip.text = element_text(margin = margin(1,1,1,1), face="bold"))
}

return(lp)

}



### The function deps_b2b_lollipop creates a back to back deps plot (lollipop) starting from an deps dataframe with 2 analyses ----

deps_b2b_lollipop <- function(input_df, # input dataframe
                                    
                                    bb_col="class", # column used for the b2b (2 levels)
                                    x_col="N", # column for x axis, necessary
                                    y_col="comp", # column for y axis, necessary
                                    
                                    size_col="N", # column for size, optional
                                    size_vec=c(1,5), # size vector c(min,max), optional
                                    
                                    shape_col="class", # column for shape, optional
                                    shape_vec=c(25,24), # shape vector, optional
                                    
                                    color_col="class", # column for color, optional
                                    color_vec=c("grey60","grey40"), # color vector, optional
                                    
                                    fill_col="class", # column for fill color, optional
                                    fill_vec=c("grey60","grey40"), # fill vector, optional
                                    
                                    text_col="N", # column for fill color, optional
                                    
                                    facet_col=c("Down-regulated", "Up-regulated"), # column for facet, optional
                                    
                                    char_max = 50, # maximum number of chars in text
                                    sort_y = T, # sort rows according to y?
                                    shift_mult = 1, # space left for middle labels (by default, if left to 1, equal to the space dedicated to the lollipop on each side)
                                    break_vec=seq(0,1500,250) # specification of x breaks
){
  
  plot_df <- input_df %>% dplyr::select(bb_col=all_of(bb_col),
                                        x_col=all_of(x_col),
                                        y_col=all_of(y_col),
                                        size_col=any_of(size_col),
                                        shape_col=any_of(shape_col),
                                        color_col=any_of(color_col),
                                        fill_col=any_of(fill_col),
                                        text_col=any_of(text_col),
                                        facet_col=any_of(facet_col))
  
  # creation of labels, abbreviating text 
  plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
  transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
  plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
  rm(transient)
  
  # factorize bb_col
  if(!is.factor(plot_df$bb_col)){
    plot_df$bb_col<-factor(plot_df$bb_col)
  }
  
  # sort y_col
  if(sort_y){
    plot_df<- plot_df %>% dplyr::arrange(-x_col)
  }
  
  # remove eventual duplicates and factorize y_col 
  plot_df<-plot_df[!duplicated(plot_df %>% dplyr::select(any_of(c("bb_col","y_col","facet_col")))),]
  plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)
  
  # define shift (if 1, roughly 1/3 of the plot)
  shift<-(max(plot_df$x_col)/2)*shift_mult
  
  f_left<- levels(plot_df$bb_col)[1]
  f_right<- levels(plot_df$bb_col)[2]
  
  lp<- ggplot(plot_df, aes(y = y_col))+
    geom_linerange(data = plot_df[plot_df$bb_col==f_left,],
                   aes(xmin = -shift, xmax = -shift-x_col, size=0.75))+
    geom_linerange(data = plot_df[plot_df$bb_col==f_right,],
                   aes(xmin = shift, xmax = shift+x_col, size=0.75))+
    geom_text(data = plot_df[plot_df$bb_col==f_left,],
              aes(y = y_col, x = 0, label = y_col),
              inherit.aes = F,family=bf,size=4)+
    scale_x_continuous(limits =c((-shift-max(plot_df$x_col)),(shift+max(plot_df$x_col))),
                       breaks = c(rev(-break_vec)-shift, break_vec+shift),
                       labels = c(as.character(rev(break_vec)),as.character(break_vec)))+
    theme_tufte(base_family = bf)+
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none")+
    labs(x = "N° DEPs")
  
  
  if("size_col" %in% colnames(plot_df)){
    lp <- lp + 
      geom_point(data = plot_df[plot_df$bb_col==f_left,],
                 aes(x = -shift-x_col,size=size_col))+
      geom_point(data = plot_df[plot_df$bb_col==f_right,],
                 aes(x = shift+x_col,size=size_col))+
      scale_size(name=size_col,range = size_vec)
  } else {
    lp <- lp + 
      geom_point(data = plot_df[plot_df$bb_col==f_left,],
                 aes(x = -shift-x_col),size=2)+
      geom_point(data = plot_df[plot_df$bb_col==f_right,],
                 aes(x = shift+x_col),size=2)
  }
  
  if("text_col" %in% colnames(plot_df)){
    lp <- lp + geom_text(data = plot_df[plot_df$bb_col==f_left,],
                         aes(x = -shift, label=text_col, hjust=(-0.05)),
                         fontface="italic", size=4, show_guide = F) +
      geom_text(data = plot_df[plot_df$bb_col==f_right,],
                aes(x = shift, label=text_col, hjust=(1.1)),
                fontface="italic", size=4, show_guide = F) 
    # scale_x_continuous(limits =c((-shift-max(plot_df$x_col)),(shift+max(plot_df$x_col))),
    #                    breaks = c(rev(-break_vec)-shift, break_vec+shift),
    #                    labels = c(as.character(rev(break_vec)),as.character(break_vec)),
    #                    expand = expansion(mult = c(.3, .3)))
  }
  
  if("shape_col" %in% colnames(plot_df)){
    lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
  }
  
  if("color_col" %in% colnames(plot_df)){
    lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col, values= color_vec) # drop = FALSE
  }
  
  if("fill_col" %in% colnames(plot_df)){
    lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
  }
  
  # lp <- lp + ggforce::facet_col(vars(facet_col), scales = "free_y", space = "free") +
      # theme(strip.text = element_text(margin = margin(1,1,1,1), face="bold"))

  lp <- lp + geom_text(aes(x=shift, label=c("Up-regulated", rep("",nrow(plot_df)-1)), hjust=(-0.05), vjust=(-2), fontface=2),family=bf,size=4)
  lp <- lp + geom_text(aes(x=-shift, label=c("Down-regulated", rep("",nrow(plot_df)-1)), hjust=(1.1), vjust=(-2), fontface=2),family=bf,size=4)
  return(lp)
  
}




# The function enrichment_dotmatrix creates am enrichment matrix dotplot, starting from an enrichment dataframe with multiple analyses ----

enrichment_dotmatrix <- function(input_df, # input dataframe
                                x_col="input_name", # column for x axis, necessary
                                y_col="anno_name", # column for y axis, necessary
                                
                                size_col="odds_ratio", # column for size, optional
                                size_vec=c(1,5), # size vector c(min,max), optional
                                
                                shape_col="P<0.05", # column for shape, optional
                                shape_vec=c(16,21), # shape vector, optional
                                
                                color_col="P<0.05", # column for color, optional
                                color_vec=c("grey25","grey50"), # color vector, optional
                                
                                fill_col="P<0.05", # column for fill color, optional
                                fill_vec=c("grey25","white"), # fill vector, optional
                                
                                facet_col="anno_class", # column for facet, optional

                                char_max = 50, # maximum number of chars in text (y_axis)
                                sort_y = F, # sort x rows? (in development)
                                sort_x = F # sort y rows? according to size_col (in development)
){

#library(tidyverse)
#library(ggthemes)
#library(extrafont)
#loadfonts()
   
   # select and rename relevant columns for the plot   
   plot_df <- input_df %>% dplyr::select(x_col=all_of(x_col),
                                         y_col=all_of(y_col),
                                         size_col=any_of(size_col),
                                         shape_col=any_of(shape_col),
                                         color_col=any_of(color_col),
                                         fill_col=any_of(fill_col),
                                         facet_col=any_of(facet_col))
   
   # creation of labels, abbreviating text 
   plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
   transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
   plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
   rm(transient)
   
   # order rows by x, remove eventual duplicates
   if(sort_y | sort_x){
      plot_df<- plot_df %>% dplyr::arrange(-size_col)
   }
   
   plot_df<-plot_df[!duplicated(plot_df %>% dplyr::select(any_of(c("x_col","y_col","facet_col")))),]
   
   if(!(is.factor(plot_df$x_col))){
   plot_df$x_col<-factor(plot_df$x_col,levels=(unique(plot_df$x_col)))# factorization (x_col)
   }
   if(!(is.factor(plot_df$y_col))){
   plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)
   }


   lp <- ggplot(plot_df, aes(x_col, y_col)) +
      theme_tufte(base_family = "Arial Narrow")+
      theme(axis.title = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))  
   
   if("size_col" %in% colnames(plot_df)){
      lp <- lp + geom_point(aes(size=size_col)) + scale_size(name=size_col,range = size_vec)
   } else {
   lp <- lp + geom_point(size=2)      
   }
   
   if("shape_col" %in% colnames(plot_df)){
      lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
   }
   
   if("color_col" %in% colnames(plot_df)){
      lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col, values= color_vec) # drop = FALSE
   }
   
   if("fill_col" %in% colnames(plot_df)){
      lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
   }
   
   if("facet_col" %in% colnames(plot_df)){
      lp <- lp + ggforce::facet_col(vars(facet_col), scales = "free_y",  space = "free") +
         theme(strip.text = element_text( margin = margin(1,1,1,1)))
   }
   
   return(lp)
   
}



# The function enrichment_geneterm_network creates a gene-term network starting from an enrichment dataframe ----

enrichment_geneterm_network <- function(input_df){

#input_df <- read_tsv("CF_GOBP_networks.txt") # %>% dplyr::filter(input_name=="trAM_inter")
gene_col <- "overlap_ids"
term_col <- "anno_name"

# iMon "#84A298"
# IM "#EAB364" 
# moAM "#889C4C" 
# trAM" #B2473E"

col_vec<- c("#B2473E","black")

plot_df <- input_df %>% dplyr::select(term_col=all_of(term_col),
                                     gene_col=all_of(gene_col))



edges_df <- plot_df  %>% separate_rows(gene_col,sep=";")
edges_df$gene_col <- str_to_title(edges_df$gene_col) %>% str_trim()

plot_df$term_col <- str_split(plot_df$term_col,"\\(",simplify=T)[,1] %>% str_trim()


terms_df<- tibble(name=unique(edges_df$term_col),
                  class="term")
genes_df<- tibble(name=unique(edges_df$gene_col),
                  class="gene")
nodes_df<-bind_rows(terms_df,genes_df)

library(geomnet)
library(data.table)
library(igraph)
library(ggraph)

net_net <- fortify(as.edgedf(as.data.frame(edges_df)), nodes_df)
net_net <- subset(net_net, from_id!=to_id)
g <- graph_from_data_frame(net_net, directed=FALSE, vertices = nodes_df)
l<-ggraph(g,layout="kk") +
  geom_edge_link2(aes(edge_colour = "grey50",alpha=0.2)) +
  ggtitle("Network",subtitle = paste0("Layout: kk"))+
  scale_edge_width(range = c(0.1,1)) +
  scale_edge_colour_manual(values = "grey50") +
  geom_node_point(aes(fill=class), shape = 21,size=2.5,pch='.') +
  geom_node_label(aes(filter = class=="term", label = name, size=2),family = bf, repel=TRUE) +
  geom_node_text(aes(filter = class=="gene", label = name, size=2),family = bf, repel=TRUE,check_overlap = TRUE) +
  scale_fill_manual(values = col_vec) +
  scale_color_brewer(palette ="Dark2") +
  theme_graph() +
  theme(legend.position = "bottom") +
  guides(text=F)
l
# ggsave("tmp_graph.pdf",l, device=cairo_pdf, width = 10, height = 6, units = c("in"))
#   
# 
# nt<-ggplot(data = net_net,aes(from_id = from_id, 
#                               to_id = to_id,
#                               label= from_id,
#                               colour=class,
#                               fill=class,
#                               size = class,
#                               shape= class,
#                               alpha= class
# )) +
#    scale_colour_manual(values=col_vec,guide="none")+
#    scale_fill_manual(values=col_vec,guide="none")+
#    scale_size_manual(values = c(3, 5),guide="none")+
#    scale_shape_manual(values=c(2,6),guide="none")+
#    scale_alpha_manual(values=c(0.9,0.9))+
#    theme_void(base_family = "Arial Narrow")+
#   geom_net(
#     layout.alg = "kamadakawai", # "kamadakawai" "fruchtermanreingold"
#     fontsize=3, repel=T, #labelcolour="black",
#     curvature = 0.05, ealpha = 0.3,
#     ecolour = "grey50", linewidth=0.7,
#     directed = F, arrowsize = 0.3, arrowgap = 0.01,
#     show.legend=F)
# nt

#ggsave("trAM_Net_trial.pdf", nt, device=cairo_pdf, width = 8, height = 5, units = c("in"))

}

# enrichment_clusterprofiler ----

enrichment_clusterprofiler <- function()

   {

library(clusterProfiler)

gene_list<-c("a","i","g","z")
universe_list<-c("a","b","c","d","e","f","g","h","i","z")

anno_df<-data.frame(
   "term"=c("on","on","on","on","tw","tw","tw","tw","tw"),
   "gene"=c("a","b","c","d","e","f","g","h","i")
)

enricher(
   gene=gene_list,
   pvalueCutoff = 1,
   pAdjustMethod = "none",
   universe=universe_list,
   minGSSize = 1,
   maxGSSize = 10000,
   qvalueCutoff = 1,
   TERM2GENE=anno_df,
   TERM2NAME = NA
)@result

}



# Function to extract coverage data from bedtools coverage with option d (coverage for each nucleotide) ----


resume_bedcoverage_d<-function(dir,
                               tailcut=4) {
   
   list_input<-list.files(path=dir, pattern="bed", full.names=F)
   
   if(exists("o_tab")){
      rm(o_tab)
   }
   
   for(i in (list_input)){
      i_name<-substr(i,1,(nchar(i)-tailcut))
      i_tab<-read.table(paste0(dir,i),header=FALSE,sep="\t")
      if(!exists("o_tab")){
         o_tab<-data.frame(name=i_tab[,8])
         rownames(o_tab)<-paste(i_tab[,4],i_tab[,7],sep="_") # correct end (-1)
         colnames(o_tab)<-c(i_name)
      } else {
         o_tab[,i_name]<-i_tab[,8]
      }
   }
   return(o_tab)
}

# Function to create a data frame of gene counts from bedtools coverage files (bed format) ----


resume_bedcoverage<-function(dir, # string: directory with the bed files
                             name_cols=c(1,2,3,6), # numeric or vector: columns with feature names
                             cov_col=7  # numeric: column with coverage values
) {
   
   list_input<-fs::dir_ls(dir,glob="*.bed",type="file") # only files with .bed extension
   
   if(exists("o_tab")){
      rm(o_tab)
   }
   
   for(i in (list_input)){
      
      i_tab<-read_tsv(i,col_names=F, progress = F,show_col_types = F)
      
      s_name<-str_split(i,dir,simplify=T)[,2] 
      s_name<-str_split(s_name,"\\.",simplify=T)[,1]
      
      x_tab <- i_tab %>% unite(id_name,all_of(name_cols)) %>% dplyr::select(id_name)
      x_tab[,s_name] <- i_tab[,cov_col]
      
      rm(i_tab)
      
      if(!exists("o_tab")){
         o_tab<- x_tab
      } else {
         o_tab<- full_join(o_tab,x_tab,by="id_name")
      }
   }
   return(o_tab)
}



### Function to use edgeR function glmQLF generalized linear models (quasi likelihood function) ----


edgeRglmQLF<-function(mat=edge_f, # object of class DGEGLM
                      contro, # comparison, create with makeContrasts
                      cpm_mat=edge_n, #used to calculate average signal
                      label="", # include in label column
                      sig_thr=0.5, # signal threshold
                      sig_col="CPM", # signal column
                      fc_thr=0.5, # fold change threshold
                      pval_thr=0.05, # p-value threshold
                      pval_col="p_val", # p-value column
                      names=FALSE)
{ # include label in column names
   degs<-glmQLFTest(edge_f,contrast=contro)$table[,-3]
   # degs<-glmTreat(edge_f,contrast=contro,lfc=log2(1.2))$table[,-3] # alternative with FC threshold
   colnames(degs)<-c("log2_FC","log2_CPM","p_val")
   a_levels<-rownames(contro)[which(contro!=0)]
   a_samples<-which(cpm_mat$samples$group%in%a_levels)
   cpm_sele<-cpm(cpm_mat,log=T)[,a_samples]
   degs$log2_CPM<-apply(cpm_sele,1,function(x) mean(x))
   #degs<-exactTest(edge_c, pair=cond, dispersion=bcv^2)$table
   degs$p_adj<-p.adjust(degs$p_val, method ="BH")
   degs$CPM<-round(2^degs$log2_CPM,2)
   degs$class<-"="
   degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC>=fc_thr & degs[,pval_col]<=pval_thr),"class"]<-"+"
   degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC<=(-fc_thr) & degs[,pval_col]<=pval_thr),"class"]<-"-"
   degs$class<-as.factor(degs$class)
   degs$comp<-label
   degs$id<-rownames(degs)
   degs<-degs[,c("id","comp","log2_FC","CPM","p_val","p_adj","class")]
   if(names=="TRUE"){
      newnames<-paste(label,colnames(degs),sep="_")
      colnames(degs)<-newnames
   }
   return(degs)
}





### Function to load all the sheets of an excel file in a list ----


read_excel_allsheets <- function(filename, tibble = FALSE) 
{
   sheets <- readxl::excel_sheets(filename)
   x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, guess_max = 1048576))
   if(!tibble) x <- lapply(x, as.data.frame)
   names(x) <- sheets
   x
}



### Function to extract and save only the different part of two string ----


mf <- function(x,sep){
   xsplit = strsplit(x,split = sep)
   xdfm <- as.data.frame(do.call(rbind,xsplit))
   res <- list()
   for (i in 1:ncol(xdfm)){
      if (!all(xdfm[,i] == xdfm[1,i])){
         res[[length(res)+1]] <- as.character(xdfm[,i])
      }
   }
   res <- as.data.frame(do.call(rbind,res))
   res <- apply(res,2,function(x) paste(x,collapse="_"))
   return(res)
}


### Function for the limma analysis and DEqMS spectraCountBayes  ----


limmafnc<-function(type = "PROT",c_anno,dat_gene,psm_count_table,contro_list,expr_avgse_df,signal_thr,fc_thr, pval_thr, pval_fdr){
   # make design table
   design = model.matrix(~0+c_anno$condition) 
   colnames(design) = levels(as.factor(c_anno$condition))
   rownames(design)<-c_anno$sample
   gene_matrix = as.matrix(dat_gene)
   gene_matrix <- gene_matrix[,rownames(design)]

   #DEqMS part analysis
   colnames(psm_count_table)<-c("gene","count")
   rownames(psm_count_table) = psm_count_table$gene
   
   filt_contro_list <- list()
   for (i in 1:length(contro_list)) {
     if(all(str_extract_all(contro_list[i], "\\[a-zA-Z]+")[[1]] %in% colnames(design))){
       filt_contro_list<-c(filt_contro_list,i)
     }
   }
   contro_list<-contro_list[unlist(filt_contro_list)]
   if(length(contro_list) == 0){ stop("Error: No valid contrast design given. Check the match between the spell of Condition and contrast design.")}
   contrast =  limma::makeContrasts(contrasts=contro_list,levels=design)
   colnames(contrast)<-names(contro_list)
   
   #limma part analysis
   fit1 <- limma::lmFit(gene_matrix,design)
   fit2 <- limma::eBayes(contrasts.fit(fit1,contrasts = contrast))
   
   if(type == "PROT"){
     fit2$count = psm_count_table[rownames(fit2$coefficients),"count"]
     fit3 = DEqMS::spectraCounteBayes(fit2) # Perform analysis based on the number of PSMs
   }else if(type == "PEP"){
     fit3<-fit2
   }
   
   # extract results
   
   degs_w_df<-NULL
   degs_l_df<-NULL
   
   signal_col="log2_expr"
   fc_col="log2_FC"
   if(pval_fdr){pval_col="p_adj"}else{pval_col="p_val"}
   
   for(comp in colnames(fit3$coefficients)){
       if(type == "PROT"){
          degs_u<-tryCatch({
             DEqMS::outputResult(fit3,coef_col = comp) %>% 
                dplyr::rename("log2_FC"="logFC",
                              "log2_expr"="AveExpr",
                              "p_val"="sca.P.Value",
                              "p_adj"="sca.adj.pval",
                              "id"="gene") %>% dplyr::arrange(id)
          }, 
          error=function(cond){
             colnames(fit3$sca.t)<-colnames(fit3$coefficients)
             colnames(fit3$sca.p)<-colnames(fit3$coefficients)
             DEqMS::outputResult(fit3,coef_col = comp) %>% 
                dplyr::rename("log2_FC"="logFC",
                              "log2_expr"="AveExpr",
                              "p_val"="sca.P.Value",
                              "p_adj"="sca.adj.pval",
                              "id"="gene") %>% dplyr::arrange(id)
          })
       }else if(type == "PEP"){
         degs_u<-topTable(fit3, coef = comp, number = Inf) %>% 
           dplyr::rename("log2_FC"="logFC",
                         "log2_expr"="AveExpr",
                         "p_val"="P.Value",
                         "p_adj"="adj.P.Val") %>% dplyr::mutate("id"=rownames(topTable(fit3, coef = comp, number = Inf))) %>% dplyr::select(-c(t,B))
       }
      rownames(degs_u)<-degs_u$id
      #Recalculate the log2_expr
      name_log2_expr<-paste0(names(contrast[which(contrast[,comp] != 0),comp]),"_avg")
      mean_exprAvg_column<-rowMeans(x=expr_avgse_df[name_log2_expr], na.rm = TRUE)
      degs_u$log2_expr<-mean_exprAvg_column
      
      degs_u$comp<-comp
      degs_u$class<-"="
      degs_u[which(degs_u[,signal_col]>=signal_thr & 
                      degs_u[,fc_col]>=fc_thr & 
                      degs_u[,pval_col]<=pval_thr),"class"]<-"+"
      
      degs_u[which(degs_u[,signal_col]>=signal_thr & 
                      degs_u[,fc_col]<=(-fc_thr) & 
                      degs_u[,pval_col]<=pval_thr),"class"]<-"-"
      
      
      degs_l_df<-bind_rows(degs_l_df,degs_u %>% dplyr::select("id","comp","class","log2_FC","p_val","p_adj","log2_expr"))
      
      degs_add<- degs_u %>% dplyr::select("class","log2_FC","p_val","p_adj","log2_expr")
      colnames(degs_add)<-paste(comp,colnames(degs_add),sep="_")
      degs_w_df<-bind_cols(degs_w_df,degs_add)
      
   }
   
   rownames(degs_l_df)<-NULL
   degs_w_df <- degs_w_df %>% add_column(id=rownames(degs_w_df),.before=1)
   
   return(list("degs_l_df"=degs_l_df,"degs_w_df"=degs_w_df))
}




### Function for the enrichment of the protein in the dataset. It use the function enrichment_enrichr ----


enrichRfnc<-function(in_df, pval_fdr_enrich, pval_enrich_thr, overlap_size_enrich_thr, dbs=NULL){
   DEGs_lists<-NULL
   
   for(comp_c in unique(in_df$comp)){
      
      DEGs_lists[[paste0(comp_c,"_all")]]<-in_df %>% dplyr::filter(comp==comp_c,class!="=") %>% pull(id) %>% unique() %>% sort()
      DEGs_lists[[paste0(comp_c,"_up")]]<-in_df %>% dplyr::filter(comp==comp_c,class=="+") %>% pull(id) %>% unique() %>% sort()
      DEGs_lists[[paste0(comp_c,"_down")]]<-in_df %>% dplyr::filter(comp==comp_c,class=="-") %>% pull(id) %>% unique() %>% sort()
      
   }
   
   enrfcn <- function(a) {
     source("functions_2021.R")
        frg<-DEGs_lists[[a]]
        if(length(frg)>0){
           enrichment_enrichr(frg, input_name=a, dbs_vec = dbs)
        }
        # clu_n <- enrichment_enrichr(frg, input_name=a, dbs_vec = dbs)
        # enr_df<-rbind(enr_df,clu_n)
   }
   ncores <- min(detectCores()-2, length(unique(in_df$comp))*3)
   enr_df<-NULL
   enr_df<-do.call(rbind, mclapply(names(DEGs_lists), enrfcn, mc.cores = ncores))

   enr_df$"Significative"<-ifelse(if(pval_fdr_enrich){(enr_df$fdr<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)}else{(enr_df$p_value<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)},"TRUE","FALSE") %>% factor(levels=c("TRUE","FALSE"))
   enr_df$log2_OR<-log2(enr_df$odds_ratio)
   return(enr_df)
}

### Function for the enrichment of the protein of an universe. It use the function enrichment_enrichr ----


enrichRfnc_universe<-function(in_df, pval_fdr_enrich, pval_enrich_thr, overlap_size_enrich_thr, dbs=NULL){
  enr_df <- tryCatch({
    DEGs_lists<-NULL
    DEGs_lists[["Universe_all"]]<-in_df
    
    # enrfcn <- function(a) {
    #   source("functions_2021.R")
    #   frg<-DEGs_lists[[a]]
    #   if(length(frg)>0){
    #     enrichment_enrichr(frg, input_name=a, dbs_vec = dbs)
    #   }
    #   # clu_n <- enrichment_enrichr(frg, input_name=a, dbs_vec = dbs)
    #   # enr_df<-rbind(enr_df,clu_n)
    # }
    # ncores <- 3
    enr_df<-enrichment_enrichr(DEGs_lists[["Universe_all"]], input_name="Universe_all", dbs_vec = dbs)

    # enr_df<-do.call(rbind, mclapply(names(DEGs_lists), enrfcn, mc.cores = ncores))

    enr_df$"Significative"<-ifelse(if(pval_fdr_enrich){(enr_df$fdr<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)}else{(enr_df$p_value<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)},"TRUE","FALSE") %>% factor(levels=c("TRUE","FALSE"))
    enr_df$log2_OR<-log2(enr_df$odds_ratio)
    enr_df <- enr_df[-which(enr_df$Significative == "FALSE"),]
    enr_df
  },
  error=function(cond){
    print("\n ERROR: An error occur when connect to EnrichR. \n ")
    NULL
  })
  return(enr_df)
}


### Function to discover the communities in the network of STRING based on our genes. ----


find_communities <- function(genes, thr_score,string_gene_df){
   
   # links_selection ----
   
   dt_links <- as.data.table(string_gene_df) 
   
   remove(string_gene_df)
   dt_links[gene1 %in% genes & gene2 %in% genes]
   
   colnames(dt_links)<- c("from","to","weight")
   rownames(dt_links)<- NULL
   
   # community detection ----
   
   i_links <- dt_links %>% dplyr::filter(from %in% genes, to %in% genes, weight>thr_score)
   i_nodes <- data.frame("id"=unique(c(i_links$from,i_links$to)))
   i_net <- igraph::graph_from_data_frame(d=i_links, vertices=i_nodes, directed=F)
   i_net_2<-igraph::simplify(i_net)
   i_clp <- cluster_fast_greedy(i_net_2)
   
   i_comms_df<- data.frame("gene" = i_clp$names,
                           "comm_o" = as.factor(i_clp$membership),stringsAsFactors = F)
   
   comms_size <- i_comms_df %>% dplyr::group_by(comm_o) %>% dplyr::summarise("size"=n()) %>%
      dplyr::ungroup() %>% dplyr::arrange(-size) %>% dplyr::mutate("comm_n" = factor(order(-size))) %>% dplyr::select(comm_o,comm_n)
   i_comms_df <- left_join(i_comms_df,comms_size)
   rownames(i_comms_df)<-i_comms_df$gene
   
   i_comms_list<-list()
   for(comm_x in sort(unique(i_comms_df$comm_n))){
      i_comms_list[[comm_x]] <- i_comms_df %>% dplyr::filter(comm_n==comm_x) %>% pull(gene)
   }
   return(list("i_comms_df"=i_comms_df, "i_comms_list"=i_comms_list, "dt_links"=dt_links))
}


### Function to plot the PPI networks in multithreading ----

plot_networks<-function(g, scr_thr, bf, comp, colour_vector, bs, dirOutput_net, layouts){
  library(parallel)
  library(ggraph)
  library(ggplot2)
  library(doParallel)
  library(foreach)
  library(qpdf)
  #Multicore parallelizatio of plot
  ncores <- length(layouts)
  cluster_ext <- makeCluster(ncores, type = "SOCK")
  registerDoParallel(cl = cluster_ext)
  name_list=vector()
  name_list<-foreach (l = layouts, .packages = c("ggraph"), .combine = 'c') %dopar% {
    tryCatch({
      l_list<-ggraph(g,layout=l) +
        # geom_edge_link2(aes(edge_width=(0.05+abs(weight-scr_thr)/400), alpha = (weight/1000)-0.5),edge_colour = "grey70") +
        geom_edge_link2(aes(edge_width=weight, alpha = (weight/1000)-0.5),edge_colour = "grey70") +
        
        # ggtitle(paste0("Network of ",comp),subtitle = paste0("Layout: ",l))+
        scale_edge_width(range = c(0.1,1)) +
        geom_node_point(aes(fill = Community, size = Degree), shape = 21,pch='.') +
        geom_node_text(aes(label = name, color=Community),family = bf, repel = TRUE, size=bs*0.3, show.legend = FALSE) +
        scale_fill_manual(values = colour_vector[[comp]]) +
        # scale_edge_width_continuous(name=c("Weigth","Size","Community")) +
        scale_color_manual(values = colour_vector[[comp]]) +
        theme_bw(base_size = bs, base_family = bf) +
        theme(legend.position = "bottom", panel.grid = element_blank(),
              panel.border = element_blank(), panel.background = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
        guides(text=F, colour = guide_legend(override.aes = list(size=5)), edge_alpha="none")
      l_list$labels$edge_width<-"STRINGdb score"
      # l_list$labels$edge_alpha<-"Weight"
      l_list
      ggsave(paste0(dirOutput_net,gsub(comp, pattern = "\\/", replacement="vs"),"_",l,"_network.pdf"), l_list, device=cairo_pdf, width = 10, height = 6, units = c("in"))
      #
      print(paste0(dirOutput_net,gsub(comp, pattern = "\\/", replacement="vs"),"_",l,"_network.pdf"))
      # name_list<-append(name_list,paste0(dirOutput_net,comp,"_",l,"_network.pdf"))
    }, error=function(cond){print("Error: Error occur during the plot of network.\n")})
  }
  tryCatch({
    pdf_combine(input = name_list, output = paste0(dirOutput_net,comp,"_network.pdf"))
    unlink(name_list)
  }, error=function(cond){print("Error: Not possible combine PDFs network in single file.\n")})
  stopCluster(cluster_ext)
}


### Function to kinase study in multithreading. It identify the kinases, their activities and it draw the CORAL tree ----

kinase_act_phosr <- function(dirOutput_kinase, formule_CORAL, comp, dat_pep, deps_pep_l_df, psm_peptide_table, c_anno_phos, df){
  
  data("KinaseMotifs")
  data("PhosphoSitePlus")
  
  tmp_dat_pep <- dat_pep[rownames(psm_peptide_table[which(psm_peptide_table[, "GeneName"] 
                                                          %in% 
                                                            deps_pep_l_df[which(deps_pep_l_df$comp == comp 
                                                                                & 
                                                                                  (deps_pep_l_df$class == "+" 
                                                                                   | 
                                                                                     deps_pep_l_df$class == "-")),
                                                                          "id"]),
                                                    ]),
                         ]
  rownames(tmp_dat_pep) <- make.names(psm_peptide_table[rownames(tmp_dat_pep), "GeneName"], unique = T)
  
  tmp_dat_pep<- tmp_dat_pep[,(c_anno_phos[str_remove(c_anno_phos$condition, "_p\\b") %in% df$p[unlist(lapply(df$p, function(x) grepl(x, formule_CORAL[comp])))], "sample"])]
  ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(tmp_dat_pep)))
  colnames(ppe@assays@data$Quantification) <- c_anno_phos[str_remove(c_anno_phos$condition, "_p\\b") %in% df$p[unlist(lapply(df$p, function(x) grepl(x, formule_CORAL[comp])))], "condition"]
  GeneSymbol(ppe)<-psm_peptide_table[rownames(psm_peptide_table[which(psm_peptide_table[, "GeneName"] 
                                                                      %in% 
                                                                        deps_pep_l_df[which(deps_pep_l_df$comp == comp 
                                                                                            & 
                                                                                              (deps_pep_l_df$class == "+" 
                                                                                               | 
                                                                                                 deps_pep_l_df$class == "-")),
                                                                                      "id"]),
  ]),"GeneName"]
  Sequence(ppe)<-psm_peptide_table[rownames(psm_peptide_table[which(psm_peptide_table[, "GeneName"] 
                                                                    %in% 
                                                                      deps_pep_l_df[which(deps_pep_l_df$comp == comp 
                                                                                          & 
                                                                                            (deps_pep_l_df$class == "+" 
                                                                                             | 
                                                                                               deps_pep_l_df$class == "-")),
                                                                                    "id"]),
  ]), "Annotated Sequence"]
  
  mat <- SummarizedExperiment::assay(ppe, "Quantification")
  grps =colnames(ppe)
  
  substrate.list = PhosphoSite.human
  substrate.list=sapply(substrate.list, function(y){unlist(lapply(y, function(x){ gsub(";.+","",x)}))})
  
  mat.std <- standardise(mat)
  mat.std.geneSymbol <- mat.std
  rownames(mat.std.geneSymbol)<-str_to_upper(unlist(str_split_fixed(rownames(mat.std.geneSymbol), "\\.", n=2)[,1]))
  mat.std<-data.frame(mat.std)
  mat.std$geneSymbol <- unlist(str_split_fixed(rownames((mat.std)), "\\.", n=2)[,1])
  mat.std$id<-rownames(mat.std)
  seqs <- na.omit(ppe@Sequence)
  kssMat <- kinaseSubstrateScore(substrate.list = substrate.list,
                                 mat = mat.std.geneSymbol,
                                 seqs = seqs,
                                 numMotif = 5,
                                 numSub = 1,
                                 species = "human",
                                 verbose = T)
  set.seed(42)
  predMat <- kinaseSubstratePred(kssMat, inclusion = 5)
  
  colnames(kssMat$ksActivityMatrix) <- str_remove(c_anno_phos[str_remove(colnames(kssMat$ksActivityMatrix), "_p\\b"), "condition"], "_p\\b")
  #Fare magia per le varie formule
  tttt<- str_remove(c_anno_phos$condition, "_p\\b")
  design = model.matrix(~0+tttt)
  colnames(design) = levels(as.factor(tttt))
  rownames(design)<-str_remove(c_anno_phos$sample, "_p\\b")
  contrast =  limma::makeContrasts(contrasts=formule_CORAL[comp],levels=design)
  
  mean_kinase_activity <- lapply(rownames(contrast)[which((contrast != 0)[,1])], function(x){(rowMeans(kssMat$ksActivityMatrix[, grepl(x, colnames(kssMat$ksActivityMatrix))])[colnames(predMat)])*contrast[x,]})
  kinase_Act <- Reduce("+", mean_kinase_activity)
  
  write.table(data.frame(kinase_Act), file = paste0(dirOutput_kinase,comp,"_kinase_activity.txt"), col.names = F, quote = F)
  
  renderSvgPanZoom(comp, kinase_Act, dirOutput_kinase)
}
