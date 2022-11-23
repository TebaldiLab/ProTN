#################### network -----
# load necessary packages, functions and data ----
setwd("C:/Users/emmab/Desktop/TDP43 stuff/Human_2018")

library(tidyverse)
library(ggrepel) # graphics
library(ggbeeswarm) # graphics
library(readxl)
library(RColorBrewer) # Color palettes 
library(scales) # heatmap colors scales

library(geomnet) # network plots
library(igraph) # Network analysis
library(data.table)

library(extrafont) # plots with arial font (and cairo_pdf library)
loadfonts()
bs<-11      # point sizes in ggplot
bf<-"Arial" # font family in ggplot

# load("rdata/string_11_gene_mm.RData")
# load("rdata/string_11_gene_hs.RData") # string annotation

info_string<-read_tsv("T:/UNITN/LM_QCB/Proteomics_internship/general_pipeline_PD_files/lib/network/9606.protein.info.v11.5.txt")[,c(1,2)]
links_string<-read.delim("T:/UNITN/LM_QCB/Proteomics_internship/general_pipeline_PD_files/lib/network/9606.protein.links.v11.5.txt", sep=" ")
links_string$protein1 <- info_string$preferred_name[match(links_string$protein1,info_string$`#string_protein_id`)]
links_string$protein2 <- info_string$preferred_name[match(links_string$protein2,info_string$`#string_protein_id`)]
string_gene_df<-data_frame("gene1"=links_string$protein1,"gene2"=links_string$protein2,"weigth"=links_string$combined_score)
  
#if use human data change from geneA, geneB to gene1,gene2
#if use human data change from gene1, gene2 to geneA,geneB

#source(paste0(utils_folder, "diff_analysis_functions.R"))
find_communities <- function(genes, thr_score){
  # links_selection ----
  
  dt_links <- as.data.table(string_gene_df) 
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

# select list of genes ---- 
#degs_enrich <- fread("degs_enrich/degs_synergy.txt", header = TRUE)
#degs_enrich<-fread("common_degs_404.txt")
# load("rdata/axonal/dt_long_1220.RData")
# ax_down_free<-dt_long_1220[analysis=="degs" & fract=="free" & comp=="ax"]
# ax_down_poly<-dt_long_1220[analysis=="degs" & fract=="poly" & comp=="ax"]
# degs_enrich<-ax_down_poly

#Load MY DATA
degs_hum <- read_tsv("T:/UNITN/LM_QCB/Proteomics_internship/astro/outputResults/sent_Basso/20220331131431_output_ev/ev_DEGs_v2.txt")
degs_hum<-degs_hum[,c(1,3,4)]


#human
# degs_hum <- read_excel("data/Human_Degs.xlsx")
# degs_hum<-degs_hum[,c(2,4,7)]
names(degs_hum)<-c("external_gene_name","log2FC","pvalue")
all_degs<-as.data.table(degs_hum)
all_degs[,class:=ifelse(log2FC>0.75 & pvalue<0.05, "+", 
                        ifelse(log2FC< -0.75 & pvalue<0.05, "-","=" ))]
degs_enrich<-all_degs
# degs_ip <- fread("degs_ip/degs_ip.txt", header = TRUE)
# degs_input <- fread("degs_input/degs_input.txt", header = TRUE)
sel_class <- "+-"
if(sel_class == "+"){g_sel <- degs_enrich[class == "+", external_gene_name]; labplot <- "up"}
if(sel_class == "-"){g_sel <- degs_enrich[class == "-", gene_symbol]; labplot <- "down"}
if(sel_class == "+-"){g_sel <- degs_enrich[class != "=", external_gene_name]; labplot <- "all"}

#labplot<-"Common_404_Degs"
#g_sel<-degs_enrich$x
# params ----
thr_score <- 150 # define max strength of interaction
scr_thr <- 850 # visual
links <- 5

# find communities ----
g_sel<-toupper(g_sel)
comm <- find_communities(g_sel, thr_score)
i_comms_df <- comm[["i_comms_df"]]
i_comms_list <- comm[["i_comms_list"]]
dt_links <- comm[["dt_links"]]

length(i_comms_list)

colour_vector <- c(c(brewer.pal(n = 8, name = "Dark2")), "#680000", "#ae0001", "#eb8c00", "#680000", "#001080", "#999999", "#434343")
colour_vector <- colour_vector[1:length(i_comms_list)]

nc <- ggplot(i_comms_df, aes(x=comm_n, fill=comm_n)) +
  geom_bar(alpha=0.9, colour="white", width=0.8)+
  scale_x_discrete(drop=FALSE)+
  scale_fill_manual(values=colour_vector,drop=FALSE)+
  theme_bw(base_size = 25) +
  labs(x = "Community", y = "# genes") +
  theme(legend.position = "none", panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

nc
# ggsave(paste0("figures_mouse/enrich_degs_", labplot, "_", thr_score, ".png"), nc, device="png", width = 5, height = 5, units = c("in"))

# Filter genes and links ----
net_genes <- g_sel
net_edges <- subset(dt_links, from%in%net_genes & to%in%net_genes & weight>scr_thr) #select genes with weight>scr_thr

gene_links<-table(c(net_edges$from,net_edges$to))
filt_genes<-names(gene_links)[which(gene_links > links)] #select genes with link > links
net_edges<-subset(net_edges,from %in% filt_genes & to %in% filt_genes) 

gene_links<-table(c(net_edges$from,net_edges$to))
filt_genes<-names(gene_links)[which(gene_links > 2)] #remove isolated genes (leaves)
net_edges<-subset(net_edges,from %in% filt_genes & to %in% filt_genes)

net_edges$group_to<-"PP"

# Vertices df ----

gene_vertices<-data.frame(label=filt_genes)
rownames(gene_vertices)<-gene_vertices$label

gene_vertices$class<-i_comms_df[as.character(gene_vertices$label),"comm_n"] #class=communities


gene_vertices$n_size<-"A" #column that define the size of the dots
#gene_vertices[which(rownames(gene_vertices) %in% specific),"n_size"]<-"B"

gene_vertices$type<-"A" # can be used for shape
#gene_vertices[which(rownames(gene_vertices)%in%sma_genes),"type"]<-"B"

ggplot(gene_vertices, aes(x=class, fill=class)) +
  geom_bar(alpha=0.9,colour="white", width=0.8)+
  scale_x_discrete(drop=FALSE)+
  scale_fill_manual(values=colour_vector, drop=FALSE)+
  theme_classic()

dim(gene_vertices)

# plot with DiagrammeR----
net_net <- fortify(as.edgedf(as.data.frame(net_edges)), gene_vertices)
uniqueN(net_net$from_id)

#fruchtermanreingold
#kamadakawai --> algortim for dots organization 

# plot with ggraph and graphlayouts----

library(networkdata)
library(igraph)
library(ggraph)
library(graphlayouts)

node_list_from<-as.data.table(unique(net_net[,c("from_id")]))
node_list_to<-as.data.table(unique(net_net[,c("to_id")]))
colnames(node_list_from)<-"node"
colnames(node_list_to)<-"node"
node_list<-as.data.table(rbind((node_list_from),(node_list_to)))
node_list<-node_list[!duplicated(node_list$node)]
node_list[,id:=seq(1, nrow(node_list))]
node_list<-merge(node_list,net_net[,c("from_id","class")],by.x="node",by.y="from_id")

node_list<-unique(node_list)

colnames(node_list)<-c("name","id","class")

edge_list<-net_net[,c("from_id","to_id","weight")]
colnames(edge_list)<-c("from","to","weight")
edge_list<-as.data.table(edge_list)
edge_list<-na.omit(edge_list)
g<- graph_from_data_frame(edge_list, directed=TRUE, vertices=node_list)


normalise <- function(x, from = range(x), to = c(0, 1)) {
  x <- (x - from[1]) / (from[2] - from[1])
  if (!identical(to, c(0, 1))) {
    x <- x * (to[2] - to[1]) + to[1]
  }
  x
}

V(g)$degree <- normalise(V(g)$degree, to = c(3, 11))
pdf(file = "T:/UNITN/LM_QCB/Proteomics_internship/general_pipeline_PD_files/lib/network/graph.pdf")
for (l in c("bipartite","star","circle","nicely","dh","gem","graphopt","grid","mds","sphere","randomly","fr","kk","drl","lgl")){
tryCatch({
  print(l)
plot<-ggraph(g,layout=l) +
  # geom_edge_link(
  #   aes(end_cap = circle(node2.degree, "pt"),
  #       edge_width=(0.05+abs(weight-scr_thr)/800),alpha=0.1),
  #   edge_colour = "grey50",
  #   arrow = arrow(
  #     angle = 10,
  #     length = unit(0.15, "inches"),
  #     ends = "last",
  #     type = "closed"
  #   )
  # ) +
geom_edge_link2(aes(edge_width=(0.05+abs(weight-scr_thr)/400)),edge_colour = "grey50",alpha=0.2) +
  ggtitle(l)+
  scale_edge_width(range = c(0.1,1)) +
  
  geom_node_point(aes(fill = class), shape = 21, size = 3,pch='.') +
  geom_node_text(aes(label = name,size=9,color=class),
                 family = "Futura", repel = TRUE
  ) +
  scale_fill_manual(values = colour_vector) +
  scale_color_brewer(palette ="Dark2") +
  #scale_size(range = c(2, 5), guide = "none") +
  theme_graph() +
  theme(legend.position = "bottom")
print(plot)
}, error=function(cond){print(l)})
}
dev.off()

ggsave(paste0("T:/UNITN/LM_QCB/Proteomics_internship/general_pipeline_PD_files/lib/network/graph", labplot, "_", thr_score, ".pdf"), plot, device=cairo_pdf, width = 20, height = 20, units = c("in"))
