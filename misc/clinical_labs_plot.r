# clinical data 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
figpath = here('misc/figures/'); dir.create(figpath)

# read clin data 
d = data.table::fread(file = here('data/starting_data/clinical/CITEseq_CPK_irAE_sorted_forMatt.csv'))

# add irae info 
d$irae = ifelse(d$ID %in% c(1:5), 'IRAE', "NO IRAE")

# format 
d2 = d %>%
  filter(Lab_ID == "Creatine Kinase" ) %>% 
  separate(Normal.Range, into = c('min_normal', 'max_normal'), sep = '-') %>%
  mutate(max_normal = as.numeric(max_normal)) %>% 
  group_by(ID) %>% 
  summarize(max_normal = max_normal, 
            Timing_StartToIRAE = Timing_StartToIRAE,
            Timing_ToStartDate = Timing_ToStartDate) %>% 
  mutate(Timing_StartToIRAE = ifelse(Timing_StartToIRAE < 0 , NA, Timing_StartToIRAE))

# plot CK 
d_ = d %>% filter(Lab_ID == "Creatine Kinase" )
p = 
  ggplot(d_, aes(y = Lab_value, x = Timing_ToStartDate, color = irae)) +
  labs(x = 'time (days) from Avelumab start date ', y = 'Creatine Kinase') + 
  scale_color_manual(values = c('firebrick', 'dodgerblue')) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        legend.position = c(0.92, 0.2), 
        axis.text.x = element_text(size = 5)) + 
  geom_line(size= 1.2) +  
  facet_wrap(~ID, nrow = 2, scales = 'free_x') + 
  geom_hline(data = d2 , mapping = aes(yintercept = max_normal),
             color = 'grey48', linetype = 'dashed') + 
  geom_vline(xintercept = 0, size = 0.5, color = 'black')
ggsave(p, filename = paste0(figpath, 'CK.pdf'), width = 5.3, height = 2.6)

# alt 
p = p + facet_wrap(~ID, nrow = 1, scales = 'free_x') + 
  theme(strip.background = element_blank(), 
        legend.position = 'right', 
        axis.text.x = element_text(size = 5)) + 
  xlim(c(-20, 200))
ggsave(p, filename = paste0(figpath, 'wide_CK.pdf'), width = 9.3, height = 1.6)

sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] here_0.1        forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3    
# [6] readr_1.3.1     tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.4        cellranger_1.1.0  pillar_1.5.0      compiler_3.6.1    dbplyr_1.4.2     
# [6] tools_3.6.1       digest_0.6.27     jsonlite_1.6.1    lubridate_1.7.4   lifecycle_1.0.0  
# [11] gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10      reprex_0.3.0      cli_2.3.1        
# [16] rstudioapi_0.11   DBI_1.1.0         haven_2.2.0       withr_2.4.0       xml2_1.3.2       
# [21] httr_1.4.1        fs_1.3.2          generics_0.0.2    vctrs_0.3.6       hms_0.5.3        
# [26] rprojroot_1.3-2   grid_3.6.1        tidyselect_1.1.0  data.table_1.12.8 glue_1.4.2       
# [31] R6_2.4.1          fansi_0.4.2       readxl_1.3.1      farver_2.0.3      modelr_0.1.6     
# [36] magrittr_2.0.2    backports_1.2.1   scales_1.1.0      ellipsis_0.3.1    rvest_0.3.5      
# [41] assertthat_0.2.1  colorspace_1.4-1  labeling_0.3      utf8_1.1.4        stringi_1.4.6    
# [46] munsell_0.5.0     broom_0.7.4       crayon_1.4.1     