options(stringsAsFactors = FALSE)
library(ggplot2)

# Setting work directory
work.dir <- '~/Cancer_Resilience/Program/Code/Data/'
setwd(work.dir)

outliers <- function(data) {
  quantile.value <- quantile(data, c(0.25,0.75))
  bounds <- 1.5 * (quantile.value[2] - quantile.value[1])
  return(!((data < quantile.value[1] - bounds) | (data > quantile.value[2] + bounds)))
}

species <- c('BRCA', 'LUAD', 'LUSC')
for (cancer.type in species) {
  resN_x <- read.csv(paste('./Intermediate/Parameters_Results/x_effN_', cancer.type, 
                           '.csv', sep = ''), header = F)
  resT_x <- read.csv(paste('./Intermediate/Parameters_Results/x_effT_', cancer.type, 
                           '.csv', sep = ''), header = F)
  resN_beta <- read.csv(paste('./Intermediate/Parameters_Results/beta_effN_', cancer.type, 
                              '.csv', sep = ''), header = F)
  resT_beta <- read.csv(paste('./Intermediate/Parameters_Results/beta_effT_', cancer.type, 
                              '.csv', sep = ''), header = F)
  expT <- read.csv(paste('~/Cancer_Resilience/Program/matlab/all_exp/expT_', cancer.type, 
                         '.csv', sep = ''),row.names = 1)
  
  clinical.info <- read.delim(paste('~/Cancer_Resilience/Data/TCGA/nationwidechildrens.org_clinical_patient_',
                                    tolower(cancer.type), '.txt', sep = ''))
  if (cancer.type == 'BRCA') {
    clinical.info <- clinical.info[clinical.info$vital_status == 'Dead', 
                                   c(2, 6, 14, 16, 21, 38, 39, 40, 41)]  # BRCA
  } else if (cancer.type == 'LUSC') {
    clinical.info <- clinical.info[clinical.info$vital_status == 'Dead', 
                                   c(2, 23, 25, 27, 29, 33, 50, 59)]  # LUSC
  } else if (cancer.type == 'LUAD') {
    clinical.info <- clinical.info[clinical.info$vital_status == 'Dead', 
                                   c(2, 24, 26, 33, 85, 86, 87, 88, 89, 90, 92)]   # LUAD
  }
  
  samples.name.in.resT <- colnames(expT)
  samples.name.in.resT <- gsub('\\.', '-', samples.name.in.resT)
  samples.name.in.clinical <- paste(clinical.info$bcr_patient_barcode, '-01', sep = '')
  samples.template <- samples.name.in.clinical[samples.name.in.clinical %in% samples.name.in.resT]
  sur.time <- as.numeric(clinical.info$death_days_to[match(samples.template, samples.name.in.clinical)])
  beta_effs <- resT_x[match(samples.template, samples.name.in.resT), 1]
  if (cancer.type == 'BRCA') {
    dead.samples <- data.frame(sur.time, beta_effs)
    dead.samples$age <- as.numeric(clinical.info$age_at_diagnosis
                                   [match(samples.template, samples.name.in.clinical)])
    dead.samples$clinicalps <- clinical.info$ajcc_pathologic_tumor_stage[match(samples.template, 
                                                                               samples.name.in.clinical)]
    dead.samples$clinicalpm <- clinical.info$ajcc_metastasis_pathologic_pm[match(samples.template, 
                                                                                 samples.name.in.clinical)]
    dead.samples$clinicalpt <- clinical.info$ajcc_tumor_pathologic_pt[match(samples.template, 
                                                                            samples.name.in.clinical)]
  } else {
    dead.samples <- data.frame(sur.time, beta_effs)
    dead.samples$age <- as.numeric(clinical.info$age_at_initial_pathologic_diagnosis
                                   [match(samples.template, samples.name.in.clinical)])
    dead.samples$clinicalpt <- clinical.info$ajcc_tumor_pathologic_pt[match(samples.template, 
                                                                          samples.name.in.clinical)]
    dead.samples$clinicalps <- clinical.info$ajcc_pathologic_tumor_stage[match(samples.template, 
                                                                              samples.name.in.clinical)]
    dead.samples$clinicalpm <- clinical.info$ajcc_metastasis_pathologic_pm[match(samples.template, 
                                                                                 samples.name.in.clinical)]
  }
  
  dead.samples <- dead.samples[!(is.na(dead.samples$sur.time)), ]
  dead.samples <- dead.samples[outliers(log(dead.samples$sur.time)), ]
  dead.samples <- dead.samples[!(dead.samples$clinicalps %in% c('[Discrepancy]','[Not Available]','Stage X')),]
  dead.samples$Cancer.Type <- rep(cancer.type, dim(dead.samples)[1])
  if (cancer.type == 'BRCA') res <- dead.samples
  else res <- rbind(res, dead.samples)
}

# Plotting Fig. 4A-4C
ggplot(res, aes(y = log(beta_effs), x = as.factor(Cancer.Type), fill = Cancer.Type)) + 
  geom_boxplot(outlier.colour = 'gray60') + theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 11), legend.position = "none", 
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values=c("#FF8C00", "#9932CC", "#87CEFA")) +
  theme(legend.position = "none") + xlab('') + ylab('x_eff') +
  scale_x_discrete(labels = c('BRCA', 'LUAD', 'LUSC'))
ggplot(res, aes(y = log(sur.time), x = as.factor(Cancer.Type), fill = Cancer.Type)) + 
  geom_boxplot(outlier.colour = 'gray60') + theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 11), legend.position = "none", 
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values=c("#FF8C00", "#9932CC", "#87CEFA")) +
  theme(legend.position = "none") + xlab('') + ylab('survival time') +
  scale_x_discrete(labels = c('BRCA', 'LUAD', 'LUSC'))
ggplot(res, aes(y = log(sur.time), x = log(beta_effs))) + 
  geom_point(aes(colour = Cancer.Type)) + 
  geom_smooth(method = 'lm', se = F) + theme_bw() + 
  theme(panel.grid =element_blank()) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 11), legend.position = "none", 
        legend.text = element_text(size = 10)) +
  xlab('x_eff') + ylab('survival time') + ylim(2, 10) +
  scale_colour_manual(name = "Samples info", values=c("#FF8C00", "#9932CC", "#87CEFA"))

# Plotting Fig. 4D
BRCA.dead.samples <- res[res$Cancer.Type == 'BRCA', ]
lm.BRCA <- lm(log(sur.time) ~ age + clinicalps + clinicalpm, data = BRCA.dead.samples)
BRCA.dead.samples$adj.sur.time <- lm.BRCA$residuals
print(cor.test((BRCA.dead.samples$sur.time),log(BRCA.dead.samples$beta_effs)))
ggplot(BRCA.dead.samples, aes(y = (adj.sur.time), x = log(beta_effs))) + 
  geom_point(aes(colour = Cancer.Type)) + 
  geom_smooth(method = 'lm', se = F) + theme_bw() + 
  theme(panel.grid =element_blank()) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 11), legend.position = "none", 
        legend.text = element_text(size = 10)) +
  xlab('x_eff') + ylab('survival time') + 
  scale_colour_manual(name = "Samples info", values=c("#FF8C00", "#9932CC", "#87CEFA"))
