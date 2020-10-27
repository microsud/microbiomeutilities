


calculate_beta_stats <- function(){
  
}
library(tibble)
library(microbiomeutilities)
data("zackular2014")
ps <- zackular2014
ps <- microbiome::transform(ps, "compositional")
dist.method= "bray"
ps.dist <- distance(ps, "bray")

group <- "DiseaseState"
grp <- sym(group)
get_metaf<- function(x){
  df <- meta(x) %>% 
     rownames_to_column("Var2")
}
ps.meta <- get_metaf(ps)

ps.meta$group_var <- ps.meta[,group]

dist_melt <- reshape2::melt(
  as.matrix(ps.dist),
  varnames = c("S1", "S2")) %>% 
  mutate_if(is.factor, as.character) %>%
  left_join(ps.meta, by = c("S1" = "Var2")) %>% 
  filter(S1 != S2)
   


dist_within <- dist_melt %>% 
  group_by(S2, group_var) %>% 
  summarise(within.mean.dist=mean(value)) 

ggplot(dist_within, 
       aes(group_var, within.mean.dist)) + 
  geom_boxplot()


dist_between <- dist_melt %>% 
  left_join(ps.meta,
            by = c("S2" = "Var2"), 
            suffix = c("_1", "_2")) %>% 
  filter(group_var_1 != group_var_2) %>% 
  mutate(comparison= paste0(group_var_1, " vs ", group_var_2)) %>% 
  group_by(S2, comparison) %>% 
  summarise(between.mean.dist=mean(value))

ggplot(dist_between, 
       aes(comparison, between.mean.dist)) + 
  geom_boxplot() + geom_point() + 
  theme(axis.text.x = element_text(angle=90))

head(dist_between)

table(dist_between$S2)
head(dist_between)




head(dist_within)
