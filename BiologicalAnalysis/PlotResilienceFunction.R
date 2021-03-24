options(stringsAsFactors = FALSE)
library(ggplot2)

# Setting work directory first
work.dir <- '~/Cancer_Resilience/Program/Code/Data/'
setwd(work.dir)

cancer.type <- 'LUSC'
decision.sur <- 4.3  # Setting decision surface (see Table S3)

resN_x <- read.csv(paste('./Intermediate/Parameters_Results/x_effN_', cancer.type, 
                         '.csv', sep = ''), header = F)
resT_x <- read.csv(paste('./Intermediate/Parameters_Results/x_effT_', cancer.type, 
                         '.csv', sep = ''), header = F)
resN_beta <- read.csv(paste('./Intermediate/Parameters_Results/beta_effN_', cancer.type, 
                            '.csv', sep = ''), header = F)
resT_beta <- read.csv(paste('./Intermediate/Parameters_Results/beta_effT_', cancer.type, 
                            '.csv', sep = ''), header = F)
res_s <- data.frame(x_eff = c(resN_x$V1, resT_x$V1), beta_eff = c(resN_beta$V1, resT_beta$V1), 
                    clinical = c(rep('Normal', dim(resN_x)[1]), rep('Tumor', dim(resT_x)[1])))
maxlocs <- read.csv(paste('./Intermediate/Landscape_Results/maxlocs_', cancer.type, 
                          '.csv', sep = ''), header = F)
maxlocs <- maxlocs[maxlocs$V1 >= min(log(resN_beta)) & maxlocs$V1 <= max(log(resT_beta)),]
maxlocs <- maxlocs[maxlocs$V2 >= decision.sur - 1 & maxlocs$V2 < decision.sur + 1,]
minlocs <- read.csv(paste('./Intermediate/Landscape_Results/minlocs_', cancer.type, 
                          '.csv', sep = ''), header = F)
minlocs <- minlocs[minlocs$V1 >= min(log(resN_beta)) & minlocs$V1 <= max(log(resT_beta)),]

# Plotting resilience function
ggplot() + 
  geom_point(data = res_s, aes(x = log(beta_eff), y = log(x_eff), colour = clinical), 
             size = 0.7, alpha = 1) +
  scale_y_continuous(breaks=seq(2, 10, 2)) +
  scale_colour_manual(values=c("#228B22","#FF4040"), name="Sample Group") + theme_bw() + 
  theme(panel.grid =element_blank()) +
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 11), legend.position = "none", 
        legend.text = element_text(size = 10)) +
  geom_smooth(data = minlocs[minlocs$V2 < decision.sur,], aes(x = V1, y = V2), color = '#228B22', size = 1.5, se = F) + 
  geom_smooth(data = minlocs[minlocs$V2 > decision.sur,], aes(x = V1, y = V2), color = '#FF4040', size = 1.5, se = F) + 
  geom_smooth(data = maxlocs, aes(x = V1, y = V2), color = 'black', 
              linetype = "dashed", size = 1.5, se = F) +
  xlab('beta_eff') + ylab('x_eff') +
  labs(title = paste('Resilience function of ', cancer.type, sep = ''))

# Plotting frequency distribution histogram of status parameter at beta = th +/- c. 
# We set th = 7.8/8.2/8.7/9.2 respectively 
th <- 7.8
c <- 0.2
fre <- res_s[(log(res_s$beta_eff) > th - c) & (log(res_s$beta_eff) < th + c), ]
ggplot() + geom_histogram(data = fre, aes(log(x_eff), ..count.. / sum(..count..), fill = clinical), 
                                color = 'white', bins = 15) +
  scale_x_continuous(breaks=seq(0, 8, 1)) + theme_bw() + theme(panel.grid =element_blank()) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12), 
        legend.position = "none") + 
  scale_fill_manual(values=c("#228B22","#FF4040"), name="Sample Group") + xlim(0, 9) +
  xlab('') + ylab('')
