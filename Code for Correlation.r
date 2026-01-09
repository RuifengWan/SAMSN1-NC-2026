## ICGC HCC US cohort分析 ##
## 2023.3.13 LHY ##
library(survival)
library(survminer)
library(openxlsx)
library(ggthemes)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(org.Hs.eg.db)
library(stringi)
library(clusterProfiler)
library(ggpubr)
library(ggbreak)
library(corrplot)
library(Hmisc)
library(cowplot)
library(patchwork)


#### US cohort ####
specimen_US <- read_delim(file = "./ICGC/US/specimen.LIHC-US.tsv",delim = "\t",col_names = T)
specimen_simplify_US <- specimen_US %>% mutate(specimen_type=ifelse(specimen_type=="Primary tumour - solid tissue","T","N")) %>%
  select(1,5,7)

clinical_US <- read_delim(file = "./ICGC/US/donor.LIHC-US.tsv",delim = "\t",col_names = T)
exp_US <- read_delim(file ="./ICGC/US/exp_seq.LIHC-US.tsv",delim = "\t",col_names = T ) %>% .[,c(1,3,8,9)]
head(exp_US)
# A tibble: 6 × 4
#   icgc_donor_id icgc_specimen_id gene_id           normalized_read_count
#   <chr>         <chr>            <chr>             <dbl>
# 1 DO23018       SP49531          FAM129B           0.00000117 
# 2 DO23018       SP49531          FAM128A           0.000140   
# 3 DO23018       SP49531          FAM128B           0.000426   
# 4 DO23018       SP49531          FAM129A           0.000000231
# 5 DO23018       SP49531          FAM129C           0.000000143
# 6 DO23018       SP49531          FAM131A           0.00000432
expwide_US <- dcast(exp_US, icgc_specimen_id + icgc_donor_id ~ gene_id, fun.aggregate = mean)

# CK5cor1 CK5 vs NK inh marker #
{
  CK5cor1 <- expwide_US %>% select(SAMSN1,KLRC1,HAVCR2,TIGIT,KLRB1,LILRB1,LAG3)
  
  {
    p_list <- list()
    
    for(i in 2:ncol(CK5cor1)){
      res <- rcorr(CK5cor1$SAMSN1, CK5cor1[,i], type = "spearman")
      p_value <- signif(res$P[1,2],2)
      cor_value <- round(res$r[1,2],2)
      data_new <- CK5cor1[,c(1,i)]
      colnames(data_new) <- c("SAMSN1","y")
      
      p <- ggplot(data_new, aes(x = log2(SAMSN1*10000000+1), y = log2(y*10000000+1)))+
        geom_point(color = "#988d7b")+
        geom_smooth(method = "lm", formula = y ~ x,
                    fill = "#b2e7fa", color = "#00aeef", alpha = 0.8)+
        theme_bw()+
        ylab(colnames(CK5cor1)[i])+
        xlab("SAMSN1")+
        theme(panel.grid = element_blank(),
              axis.title = element_text( face = "bold.italic"),
              plot.title = element_text(hjust = 0.5, size = 10)
        )+
        labs(title = paste0("r =", cor_value, "  ",
                            ifelse(p_value <= 0.0001,"q < 0.0001",paste0("q = ",p_value))))
      p_list[[i-1]] <- p
    }
    p <- plot_grid(plotlist = p_list, align = "h", 
                   ncol = 4)
    ggsave("./CK5_Cor1.pdf", plot = p, height = 6,width = 12, limitsize = F)
  }
}
