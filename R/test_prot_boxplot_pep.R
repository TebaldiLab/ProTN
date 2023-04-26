# prot_intensity_long <- proBatch::matrix_to_long(dat_pep[rownames(psm_peptide_table[psm_peptide_table$GeneName %in% prot_find,]),]) %>%
#   data.table::setnames(new = c("Protein","Condition","Intensity"))
# prot_intensity_long$Protein<-psm_peptide_table[prot_intensity_long$Protein,"GeneName"]
# prot_intensity_long <- prot_intensity_long %>% dplyr::group_by(Condition,Protein) %>% dplyr::summarise(Intensity=mean(Intensity)) %>% ungroup()
# prot_intensity_long <- unique((left_join(prot_intensity_long,
#                                          c_anno[prot_intensity_long$Condition, c("sample","condition")],
#                                          by=c("Condition"="sample")))[,-c(1)])  %>%
#   data.table::setnames(new = c("Protein","Intensity","Condition"))
# 
# expr_avgse_pep_df_2<-expr_avgse_pep_df
# expr_avgse_pep_df_2$id<-psm_peptide_table[expr_avgse_pep_df_2$id,"GeneName"]
# expr_avgse_pep_df_2<-aggregate(x = expr_avgse_pep_df_2[,2:ncol(expr_avgse_pep_df_2)], by = list(expr_avgse_pep_df_2$id), FUN = mean)
# colnames(expr_avgse_pep_df_2)[1]<-"id"
# 
# avg_df = as.data.frame(expr_avgse_pep_df_2)[,c("id",grep("avg",colnames(expr_avgse_pep_df_2), value = T))] %>% 
#   dplyr::filter(id %in% prot_find) %>% dplyr::select(c("id",grep("avg",colnames(expr_avgse_pep_df_2), value = T))[-c(1)])
# se_df = as.data.frame(expr_avgse_pep_df_2)[,c("id",grep("se",colnames(expr_avgse_pep_df_2), value = T))] %>% 
#   dplyr::filter(id %in% prot_find) %>% dplyr::select(c("id",grep("se",colnames(expr_avgse_pep_df_2), value = T))[-c(1)])
# prot_avg_long <- proBatch::matrix_to_long(avg_df) %>% data.table::setnames(new = c("Protein","Condition","avg"))
# prot_avg_long$Protein<-(as.data.frame(expr_avgse_pep_df_2)[,c("id",grep("avg",colnames(expr_avgse_pep_df_2), value = T))] %>% dplyr::filter(id %in% prot_find))$id[as.integer(prot_avg_long$Protein)]
# prot_se_long <- proBatch::matrix_to_long(se_df) %>% data.table::setnames(new = c("Protein","Condition","se"))
# prot_se_long$Protein<-(as.data.frame(expr_avgse_pep_df_2)[,c("id",grep("se",colnames(expr_avgse_pep_df_2), value = T))] %>% dplyr::filter(id %in% prot_find))$id[as.integer(prot_se_long$Protein)]
# prot_avg_se_long <- cbind(prot_avg_long,prot_se_long[3])
# prot_avg_se_long$Condition <- str_replace(prot_avg_se_long$Condition, "_avg", "")
# 
# # save(prot_intensity_long, cc, bs, bf,cc, prot_avg_se_long, dat_gene, prot_find, c_anno, file = "tmp_2.RData")
# g<-ggplot(prot_avg_se_long,aes(Condition,avg,fill=Condition,colour=Condition))+
#   geom_crossbar(aes(ymin=avg,ymax=avg),position = "dodge",width=.8,alpha=.9,fatten=1.5)+
#   geom_errorbar(aes(ymin=(avg-se), ymax=(avg+se)), width=.4,position=position_dodge(),show.legend=F,alpha=.8)+
#   geom_quasirandom(data=prot_intensity_long, aes(Condition,Intensity), alpha=.7,width=.1,shape=16,size=0.11*bs)+
#   scale_fill_manual(name="Condition",values=cc[sort(unique(prot_intensity_long$Condition))]) +
#   scale_colour_manual(name="Condition",values=cc[sort(unique(prot_intensity_long$Condition))])+
#   # scale_y_continuous(limits=c(0,NA),expand = expand_scale(mult = c(.1, .25)))+
#   theme_bw(base_size = bs, base_family = bf) +
#   theme(axis.title.x=element_blank()) +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, colour=cc[sort(unique(prot_intensity_long$Condition))]))+
#   theme(panel.grid.major.x=element_blank(),
#         panel.grid.minor.y=element_blank())+
#   theme(legend.text = element_text(size = 0.7*bs),
#         legend.key.size = unit((0.015*bs),"in"),
#         legend.position="none",
#         legend.title=element_blank(),
#         legend.background = element_rect(fill = NA))+
#   
#   theme(strip.text=element_text(colour="white",face="bold"))+
#   theme(panel.border=element_rect(colour=c("grey40"),size=0.03*bs))+
#   theme(strip.background=element_rect(fill="grey40",colour="grey40",size=0.03*bs))+
#   theme(plot.title = element_text(hjust = 0.5))+
#   facet_wrap(~Protein, scales = "free",ncol = if(length(prot_find)>4){round(length(prot_find)/1.9)}else{4})+
#   labs(y="Abundance")
# print(g)
