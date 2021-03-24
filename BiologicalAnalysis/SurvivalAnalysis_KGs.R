options(stringsAsFactors = FALSE)
library(ggplot2)
library(survival)
library(survminer)

# Setting work directory first
work.dir <- '~/Cancer_Resilience/Program/Code/Data/'
setwd(work.dir)

# Specify the type of cancer to be analyzed (BRCA, LUAD or LUSC) and 
# corresponding number of cluster (4 for BRCA, 6 for LUAD, 3 for LUSC)
cancer.type <- 'BRCA'
num_cluster <- 4

# Functions about Spectral Clustering
discretisationEigenVectorData <- function(eigenVector) {
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  return(Y)
}

discretisation <- function(eigenVectors) {
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = discretisationEigenVectorData(eigenVectors %*% R)
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
  }
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

spectralClustering <- function(affinity, K, type = 3) {
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps 
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  } else if (type == 2) {
    Di = diag(1 / d)
    NL = Di %*% L
  } else if(type == 3) {
    Di = diag(1 / sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)
  return(labels)
}


# Input data
expT <- read.csv(paste('./expT_', cancer.type, '.csv', sep = ''), row.names = 1)
gene_name <- row.names(expT)
resT_beta_topgene <- read.csv(paste('./Intermediate/Resilience_Centrality_Results/resT_beta_', 
                                    cancer.type, '.csv', sep = ''), header = F)
if (cancer.type == 'BRCA') {
  clinical.info <- read.delim('./Clinical_info/nationwidechildrens.org_clinical_patient_brca.txt')
  clinical.info <- clinical.info[c(-1, -2),c(2, 14, 15, 16)] #BRCA
} else if (cancer.type == 'LUSC') {
  clinical.info <- read.delim('./Clinical_info/nationwidechildrens.org_clinical_patient_lusc.txt')
  clinical.info <- clinical.info[c(-1, -2),c(2, 50, 28, 33)] #LUSC
} else if (cancer.type == 'LUAD') {
  clinical.info <- read.delim('./Clinical_info/nationwidechildrens.org_clinical_patient_luad.txt')
  clinical.info <- clinical.info[c(-1, -2),c(2, 26, 30, 33)] #LUAD
}

# Selecting top 200 resilience centrality genes for each tumor sample
th <- 200
get_top_gene <- function(res, gene_name, th){
  res <- as.vector(rank(res))
  return(gene_name[res <= th])
}
top_gene_beta <- lapply(resT_beta_topgene, get_top_gene, gene_name, th)

# Calculation similarity matrix based on jaccard index
jaccard_matrix <- matrix(0, length(top_gene_beta), length(top_gene_beta))
for (i in 1:length(top_gene_beta)) {
  for (j in 1:length(top_gene_beta)) {
    jaccard_matrix[i, j] <- length(intersect(top_gene_beta[[i]], top_gene_beta[[j]])) / 
      length(union(top_gene_beta[[i]], top_gene_beta[[j]]))
  }
}
diag(jaccard_matrix) <- 0

# Clustering the samples based on spectral clustering
samples.name.in.resT <- colnames(expT)
samples.name.in.resT <- gsub('\\.', '-', samples.name.in.resT)
samples.name.in.clinical <- paste(clinical.info$bcr_patient_barcode, '-01', sep = '')
samples.template <- samples.name.in.clinical[samples.name.in.clinical %in% 
                                               samples.name.in.resT]
cl <- spectralClustering(jaccard_matrix, num_cluster)
cat_num <- data.frame(sample = samples.name.in.resT, cluster = as.vector(cl))
cat_num <- cat_num[match(samples.template, samples.name.in.resT),]
clinical.info <- clinical.info[match(samples.template, samples.name.in.clinical),]
clinical.info[clinical.info$last_contact_days_to == '[Not Available]',3] <- 
  clinical.info[clinical.info$last_contact_days_to == '[Not Available]',4]
clinical.info[clinical.info$vital_status == 'Alive', 2] <- 0
clinical.info[clinical.info$vital_status == 'Dead', 2] <- 1
clinical.info$cluster <- cat_num$cluster
clinical.info$vital_status <- as.numeric(clinical.info$vital_status)
clinical.info$last_contact_days_to <- as.numeric(clinical.info$last_contact_days_to)

# Survival analysis
fit <- survfit(Surv(last_contact_days_to, vital_status) ~ cluster, data = clinical.info)
ggsurvplot(fit, pval = TRUE, conf.int = F, linetype = "strata", surv.median.line = "hv", 
           ggtheme = theme_bw())
