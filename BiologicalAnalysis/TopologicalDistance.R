options(stringsAsFactors = FALSE)
library(igraph)
library(foreach)
library(doParallel)
library(ggplot2)

# Setting work directory first
work.dir <- '~/Cancer_Resilience/Program/Code/Data/'
setwd(work.dir)

# Register parallel computing
cl <- makeCluster(15)
registerDoParallel(cl)

# Specify the type of cancer to be analyzed (BRCA, LUAD or LUSC)
cancer.type <- 'BRCA'

# Inputting data
load('./Others/fullNetwork.rda') # Loading directed gene regulated network
expT <- read.csv(paste('./expT_', cancer.type, '.csv', sep = ''), row.names = 1)
resT_beta_topgene <- read.csv(paste('./Intermediate/Resilience_Centrality_Results/resT_beta_', 
                                    cancer.type, '.csv', sep = ''), header = F)
CGC <- read.csv('./Others/cancer_gene_census.csv')
if (cancer.type == 'BRCA') {
  CGC.gene <- CGC[grep('breast', CGC$Tumour.Types.Somatic.), 1]
} else if (cancer.type == 'LUSC') {
  CGC.gene <- CGC[grep('Lung SCC|lung|NSCLC|lung cancer|lung carcinoma|Lung SSC', CGC$Tumour.Types.Somatic.), 1]
} else if (cancer.type == 'LUAD') {
  CGC.gene <- CGC[grep('lung adenocarcinoma|lung|NSCLC|lung cancer|lung carcinoma', CGC$Tumour.Types.Somatic.), 1]
}

# Getting iKGs
th <- 200
get_top_gene <- function(res, gene_name, th){
  res <- as.vector(rank(res))
  return(gene_name[res <= th])
}
top_gene_beta <- lapply(resT_beta_topgene, get_top_gene, gene_name, th)
gene.franq <- as.data.frame(table(unlist(top_gene_beta)))
topgene <- gene.franq$Var1[order(gene.franq$Freq, decreasing = T)[1:200]]

# Construction directed gene regulatory network and getting max commected component
gene.regulation.net <- graph_from_adjacency_matrix(t(originalUnweighted), mode = 'directed')
connected.components <- components(gene.regulation.net, 'weak')
max.component <- induced_subgraph(gene.regulation.net, which(connected.components$membership == 
                                                               which(connected.components$csize ==
                                                                       max(connected.components$csize))))

# Calculation topological avarage distance from CGC genes to iKGs
num <- length(which(igraph::vertex_attr(max.component,'name') %in% topgene))
sp <- distances(max.component, which(vertex_attr(max.component,'name') %in% CGC.gene),
                which(vertex_attr(max.component,'name') %in% topgene), 'out')
sp[is.infinite(sp)] <- 50
ave.distance <- mean(sp)

# Constructing null distribution
custom.function = function(x) {
  sp <- igraph::distances(max.component, which(igraph::vertex_attr(max.component,'name') %in% CGC.gene),
                          sample(1:igraph::gorder(max.component), size = num), 'out')
  sp[is.infinite(sp)] <- 50
  return(mean(sp))
}
ave.dis.random <- foreach(x = rep(num, 10000), .combine=c) %dopar% custom.function(x)
sum(ave.dis.random < ave.distance) / 10000 # P-value
