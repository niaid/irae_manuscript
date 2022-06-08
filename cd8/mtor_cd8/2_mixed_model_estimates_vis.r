# mixed model estimates of mtor sigs 
suppressMessages(library(tidyverse))
suppressMessages(library(here))

figpath = here("cd8/mtor_cd8/figuresV4//")
datapath = here("cd8/mtor_cd8/generated_dataV4//")
source('util/MattPMutils.r')

mm_res.m1 = readRDS(file = here('cd8/mtor_cd8/generated_dataV4/mm_res.m1.rds'))
data.table::fwrite(mm_res.m1,file = paste0(datapath, 'mm_res.m1.txt'), sep = '\t')
d = mm_res.m1
d$module[d$module == 'TEMRA MTOR'] = 'TM sig'
d$celltype = str_replace_all(string = d$celltype, pattern = '_' , replacement = ' ')
d$cm = paste(d$celltype, d$module,sep = ' ')


p = 
  ggplot(d, aes(x = estimatetime0_group2vs1, y = cm)) + 
  theme_bw() +
  theme(axis.text.y = element_text(color = 'black')) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_point(shape = 18, size = 2, color = col.alpha('red', 0.9)) +
  geom_segment(aes(x = (estimatetime0_group2vs1 + -1*std.errortime0_group2vs1),
                   xend = estimatetime0_group2vs1 + 1*std.errortime0_group2vs1,
                   yend = cm),color = col.alpha('red', 0.5), size = 1.3) +
  geom_point(data = d, aes(x = estimatetime1vs0, y = cm), size = 2, shape = 18, color = '#e2a359') +
  geom_segment(aes(x = (estimatetime1vs0 + -1*std.errortime1vs0),
                   xend = estimatetime1vs0 + 1*std.errortime1vs0,
                   yend = cm), color = col.alpha('#e2a359', 0.5),
               size = 1.3) + 
  xlab('single cell mixed effects \n model contrast effect size ') + ylab('') + 
  theme(axis.title.x = element_text(size = 9.4)) +
  ylab('')
ggsave(p, filename = paste0(figpath,'estimates.pdf'), width = 6, height = 10)  


# subset. 
ds = d %>% filter(celltype %in% c("CD8Tcell TEMRA", "CD8Tcell TEM","CD8Tcell naive.CD27neg"))

p = 
  ggplot(data = ds,aes(x = estimatetime0_group2vs1, y = cm)) + 
  theme_bw() +
  theme(axis.text.y = element_text(color = 'black')) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_point(shape = 18, size = 2, color = col.alpha('red', 0.9)) +
  geom_segment(aes(x = (estimatetime0_group2vs1 + -1*std.errortime0_group2vs1),
                   xend = estimatetime0_group2vs1 + 1*std.errortime0_group2vs1,
                   yend = cm),color = col.alpha('red', 0.5), size = 1.3) +
  geom_point(data = ds, aes(x = estimatetime1vs0, y = cm), size = 2, shape = 18, color = '#e2a359') +
  geom_segment(aes(x = (estimatetime1vs0 + -1*std.errortime1vs0),
                   xend = estimatetime1vs0 + 1*std.errortime1vs0,
                   yend = cm), color = col.alpha('#e2a359', 0.5),
               size = 1.3) + 
  xlab('single cell mixed effects \n model contrast effect size ') + ylab('') + 
  theme(axis.title.x = element_text(size = 9.4))
ggsave(p, filename = paste0(figpath,'estimates.sub.pdf'), width = 4.2, height = 3)  


sessionInfo()

# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] here_0.1        forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2    
# [8] tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.4       cellranger_1.1.0 pillar_1.5.0     compiler_3.6.1   dbplyr_1.4.2     tools_3.6.1     
# [7] digest_0.6.27    jsonlite_1.6.1   lubridate_1.7.4  lifecycle_1.0.0  gtable_0.3.0     pkgconfig_2.0.3 
# [13] rlang_0.4.10     reprex_0.3.0     cli_2.3.1        rstudioapi_0.11  DBI_1.1.0        haven_2.2.0     
# [19] withr_2.4.0      xml2_1.3.2       httr_1.4.1       fs_1.3.2         generics_0.0.2   vctrs_0.3.6     
# [25] hms_0.5.3        rprojroot_1.3-2  grid_3.6.1       tidyselect_1.1.0 glue_1.4.2       R6_2.4.1        
# [31] fansi_0.4.2      readxl_1.3.1     farver_2.0.3     modelr_0.1.6     magrittr_2.0.2   backports_1.2.1 
# [37] scales_1.1.0     ellipsis_0.3.1   rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 labeling_0.3    
# [43] utf8_1.1.4       stringi_1.4.6    munsell_0.5.0    broom_0.7.4      crayon_1.4.1 