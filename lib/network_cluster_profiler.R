library(clusterProfiler)
library(DOSE)
library(ReactomePA)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
de <- data_comms_df$TDP43$gene
de_2 <- mapIds(org.Hs.eg.db, de, 'ENTREZID', 'SYMBOL')
edo <- enrichGO(gene          = de_2,
                universe      = de_2,
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
View(edo)
View(edo@result)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, node_label="all",cex_label_category = 0.75, cex_label_gene = 0.5)
p1
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p2


#reactome
x <- enrichPathway(gene=de_2, pvalueCutoff = 0.05, readable=TRUE)
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, node_label="all",cex_label_category = 0.75, cex_label_gene = 0.5)
p1

library(pathfindR)
output_df <- run_pathfindR(RA_input)

#Cluster in network KEGG
#####################################################################################################################################
library(pathfindR)

P2_df <- degs_l_df[which(degs_l_df$comp=="S34F_ars"),]
input_pathfindR <- data.frame("Gene.symbol"=P2_df$id,"logFC"=P2_df$log2_FC,"adj.P.Val"=P2_df$p_val)

P2_enr_up <- enr_df[which((enr_df$input_name == "S34F_ars_up") & enr_df$anno_class == "KEGG_2021_Human" & enr_df$p_value < 0.05),
                 c("anno_name","log2_OR","p_value","overlap_ids")]
P2_enr_up$overlap_ids <- str_replace_all(P2_enr_up$overlap_ids, ";", ", ")
colnames(P2_enr_up)<-c("Term_Description","Fold_Enrichment","lowest_p","Up_regulated")
P2_enr_down <- enr_df[which((enr_df$input_name == "S34F_ars_down") & enr_df$anno_class == "KEGG_2021_Human" & enr_df$p_value < 0.05),
                 c("anno_name","log2_OR","p_value","overlap_ids")]
P2_enr_down$overlap_ids <- str_replace_all(P2_enr_down$overlap_ids, ";", ", ")
colnames(P2_enr_down)<-c("Term_Description","Fold_Enrichment","lowest_p","Down_regulated")

tmp<- merge(P2_enr_up,P2_enr_down, by=("Term_Description"), all.x=T,all.y=T)
tmp$Fold_Enrichment <- as.numeric(apply(tmp[,c("Fold_Enrichment.x","Fold_Enrichment.y")], 1, mean,na.rm=TRUE))
tmp$lowest_p <- as.numeric(apply(tmp[,c("lowest_p.x","lowest_p.y")], 1, min,na.rm=TRUE))

df_KEGG<- tmp[,c("Term_Description","Fold_Enrichment","lowest_p","Up_regulated","Down_regulated")]
df_KEGG$Up_regulated<-replace_na(df_KEGG$Up_regulated,"")
df_KEGG$Down_regulated<-replace_na(df_KEGG$Down_regulated,"")
term_gene_graph(df_KEGG, use_description = TRUE, node_size = "p_val")
#####################################################################################################################################

results_pathfindr <- run_pathfindR(input_pathfindR, 
                                   gene_sets = "KEGG", 
                                   min_gset_size = 5,
                                   pin_name_path = "STRING",
                                   output_dir = "results_analysis/Result_pathfindR")
head(results_pathfindr)
results_pathfindr_cluster<-cluster_enriched_terms(results_pathfindr)
term_gene_graph(results_pathfindr, use_description = TRUE)


