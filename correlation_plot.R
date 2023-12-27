rm(list=ls())

library(dplyr)
library(readxl)
library(Hmisc)
library(xlsx)
library(tidyr)

data_dir = "..."

pancancer = paste(data_dir, "input_file", sep = "") #read pancancer matrix file


# Tex abundance compared among pan-cancer
g=ggplot(pancancer, aes(x=reorder(cancertype, Tex, FUN=median), y=Tex))+geom_boxplot(colour='#707070')+ geom_jitter(position=position_jitter(0.2), colour = '#707070')+ylab("Tex sig.score")+xlab("Cancer Type")

c<-g+theme_classic()+theme(axis.title.x=element_blank(),
                           axis.ticks.x = element_blank(),axis.text=element_text(size=14, color = "#000000"),
                           axis.title=element_text(size=14))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(c, file = paste(data_dir, "tex_sigscore_across_pancancer_wpancreas.png", sep = ""), width = 9, height = 6, dpi= 300, units = "in", device = "png")

pancancer$group_tex <- ifelse(pancancer$Tex > quantile(pancancer$Tex, 2/3), 'High',  ifelse(pancancer$Tex < quantile(pancancer$Tex, 1/3), 'Low', 'Medium'))


# correlation analysis

scat_cor<-ggscatter(pancancer, x = "Tex", y = "TMB", # genes correlate with sig.score
                    color = "darkgray", fill="darkgray",shape = 21, size = 2,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient.
                    cor.coeff.args = list(method = "spearman", label.sep = "\n")
)+scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

g <- scat_cor + theme_classic() + theme(axis.text=element_text(size=14, color = "#000000"),
                                        axis.title=element_text(size=14)) + xlab("Tex sig.score") +ylab("Tumor mutation burden")
ggsave(g, file = paste(data_dir, "tex_tmb_corr_pancancer.png", sep=""), width = 4, height = 4, dpi= 300, units = "in", device = "png")



#R value dotplot across cancer types
#to-do: R value calculation script

ggplot(tex_corr,aes(x=cancer_type, y = variable, color = r)) + 
    geom_point(size=8) +theme_classic()+ scale_color_viridis(option = 'B',name = 'R value', limits = c(-1,1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text=element_text(size=14, color = "#000000"),
          axis.title=element_text(size=14)) +
    ylab('') + xlab('') 



# correlation analysis between Tex and other major parameters

data_df <- read.csv("C:/Users/wguo/Documents/temp_work_athome/tex_group_mut_across_all_cancer.csv", row.names = 1)
mean_df <- data_df %>% 
  group_by(cancer_type) %>% 
  summarize(sig_score=mean(sig_score),
            mut_total=mean(mut_total, na.rm = T),
            CD8A=mean(CD8A, na.rm = T),
            gpvalue=mean(gpvalue, na.rm = T))
mean_df$X<-"TCGA_MEAN"
mean_df$group="Mean"

mean_df=mean_df[order(mean_df$sig_score),]
facet_order=mean_df$cancer_type


plot_df<-rbind(mean_df, data_df)

Plot <- ggplot(gath_df, aes(x = reorder(X, gpvalue), y = value, color = group, shape=group, size=group)) + 
  geom_point(alpha=0.6) +
  facet_grid(variable~cancertype, scales = "free_y", switch="y",
             labeller = labeller(variable=c("Tex"="Tex signature", "TMB"="TMB (log10)", "CD8A"="CD8A"))) + 
  scale_color_manual(values=c("black", "firebrick", "grey", "dodgerblue")) +
  scale_shape_manual(values = c(18,16,16,16)) +
  scale_size_manual(values = c(4,0.3,0.3,0.3)) +
  scale_x_discrete(expand=c(0.2,0.2)) +
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Tex group\n(Pan-cancer)", shape = "Tex group\n(Pan-cancer)", size = "Tex group\n(Pan-cancer)",
       y = "Tex signature score", x = "TCGA cohort") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text=element_text(angle=60),
        strip.placement = "outside",
        strip.background = element_blank())