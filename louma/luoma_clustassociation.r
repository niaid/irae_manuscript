suppressMessages(library(tidyverse)) 
suppressMessages(library(here)) 

figpath = here('louma/binomial_freq_model_irae_baseline/figures/')
datapath = here('louma/binomial_freq_model_irae_baseline/data/')
dir.create(figpath, recursive = TRUE)
dir.create(datapath, recursive = TRUE )

#### Init.
s = readRDS(file = here('louma/generated_data/cd3_louma_Seurat_processed.rds'))
md = s@meta.data; rm(s); gc()

# count subclusters by sample 
d = md %>%
  filter(IRAE %in% c('colitis','control')) %>% 
  group_by(IRAE, subjectid, cluster = seurat_clusters) %>% 
  tally() %>% 
  ungroup()

d$IRAE = factor(d$IRAE, levels = c('control', 'colitis'))
tot = d %>% group_by(subjectid) %>% summarize(total = sum(n))
d = d %>%
  mutate(total = as.numeric(plyr::mapvalues(x = subjectid, from = tot$subjectid, to = tot$total))) %>%
  ungroup() %>%  
  mutate(ot = total - n) %>% 
  mutate(subjectid = factor(subjectid))
saveRDS(d, file = paste0(datapath, 'd.rds'))


#### fit aggregated binomial mixed model 
f4 = n/total ~ IRAE + cluster + IRAE:cluster + (1|subjectid) 
m4 = lme4::glmer(formula = f4, family = 'binomial',  data = d, weights = total)

# post fit checks 
emmeans::emmip(emm4, ~ IRAE | cluster)
emm4 = emmeans::emmeans(m4, specs = ~ IRAE | cluster)
plot(emm4)

cont = emmeans::contrast(object = emm4, method = "revpairwise")
link_result = 
  cont %>% 
  summary(infer = c(TRUE, TRUE), type = 'unlink') %>% 
  rbind() %>% 
  as.data.frame() %>% 
  mutate(sc = cluster) %>% 
  separate(sc, into = c("cluster","cnumber"), sep = "_")


source('util/MattPMutils.r')
link_result$col= ifelse(log(link_result$odds.ratio) > 0, yes = 'c1', no =  'c2')
cd8.clusters = c('1','4','5','7','10','13','14','18')
cs = link_result %>%
  mutate(cd8 = ifelse(cluster %in% cd8.clusters, 'CD8','N')) %>% 
  filter(cd8 == 'CD8') %>% 
  mutate(cluster = paste(cd8, cluster))

# plot 
p = ggplot(aes(x = log(odds.ratio), y = reorder(cluster, odds.ratio), fill = col  ), data = cs) + 
  theme_bw() + 
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = c( col.alpha('#f3270a', 0.5), col.alpha('#2a579b', 0.5)  )) + 
  ylab('scRNAseq cluster') + 
  xlim(c(-4.5,4.5)) +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  labs(x = 'Log Odds Colitis vs No Colitis', title = '')  + 
  theme(axis.text.y=element_text( hjust=0.1, color = 'black'))
p
ggsave(p, filename = paste0(figpath,'CD8colitis.v.nocolitis.freq.pdf'),width = 3, height = 2.2)







  