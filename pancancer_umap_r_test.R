library(tidyverse)
library(palmerpenguins)
library(umap)

theme_set(theme_bw(18))

result=read.csv('Y:/Joyce/TCGA_pancancer/archetype_pca/pancancer_wtRSR_hypoxia_clst6/XCELL_kmean_pca_umap_results.csv')
pancancer=read.csv('Y:/Joyce/TCGA_pancancer/all_extend_gene_signature_summ_pancancer_apoptosis.csv')

result$pid=result$X
View(result)
sub_pan=subset(pancancer, pid %in% result$pid)
result=subset(result, pid %in% sub_pan$pid)

merge=merge(result, sub_pan, by = 'pid')

result <- result %>% 
  drop_na() %>%
  
  mutate(ID=row_number()) 

result_meta=result %>% select(kmean_pca,pid, ID)

set.seed(142)
umap_fit <- result %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  
  mutate(ID=row_number())%>%
  inner_join(result_meta, by="ID")

colnames(umap_df)[1:2]=c('UMAP1','UMAP2')



umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, color = kmean_pca))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")
#ggsave("UMAP_plot_example1.png")
umap_df %>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2, color = cluster))+
    geom_point()+
    labs(x = "UMAP1",
         y = "UMAP2")+scale_color_manual(values=c("royalblue4","lightsalmon1","firebrick4","orangered2","#9FCAE6","cornflowerblue"))

umap_df %>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2, color = Tex))+
    geom_point()+
    labs(x = "UMAP1",
         y = "UMAP2")+scale_colour_gradientn(colours = c("blue","gray","red"))

    umap_df %>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2, color = STING))+
    geom_point()+
    labs(x = "UMAP1",
         y = "UMAP2")+scale_colour_viridis_c(option = 'turbo')

color=rev(heat.colors(10, alpha = 0.7))
umap_df %>%
ggplot(aes(x = UMAP1, 
           y = UMAP2, color = STING))+
geom_point()+
labs(x = "UMAP1",
     y = "UMAP2")+scale_colour_gradientn(colours= color)+facet_wrap(~cluster) +theme_bw()

ggsave('Y:/Joyce/TCGA_pancancer/organize_cluster/umap_batf3_facet.png', dpi=300, device = 'png', height = 4, width =7 )

+scale_colour_discrete(name='Cancer type')

ggsave('umap_r.png', dpi=300, device = 'png', height = 6, width = 9)