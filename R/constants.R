# Setup environment
# source('~/ukb_exomes/R/constants.R')
packages = c('GGally', 'reshape2', 'data.table', 'scales', 'Matrix', 'ggplot2', 'extrafont', 'gridExtra', 'grDevices', 'grid',
             'RCurl',  'tidyverse', 'devtools', 'broom', 'plotly', 'slackr', 'magrittr', 'gapminder', 'readr', 
             'purrr', 'skimr', 'gganimate', 'gghighlight', 'plotROC', 'naniar', 'BiocManager', 'cowplot', 'corrplot', 'corrr', 'ggridges', 'RColorBrewer', 
             'ggpubr',  'pbapply', 'RMySQL', 'egg', 'ggwordcloud', 'patchwork', 'ggthemes',  'ggrepel', 'LowMACA', 'dplyr')

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
}

source('~/ukb_exomes/R/constants.R')
pop_colors['mid']  <- pop_colors['mde'] 
pop_names['mid']  <- pop_names['mde']
pop_names['ami']  <- 'Amish'
pop_colors['ami']  <- "#FFC0CB" 
pop_names['afr']  <- "African" 

themes = theme(plot.title = element_text(hjust = 0.5, color = 'Black', size = 10, face = 'bold'),
               axis.text = element_text(color = 'Black', size = 10), 
               axis.title = element_text(color = 'Black', size = 15, face = 'bold'), 
               legend.title = element_text(color = 'Black', size = 12, face = 'bold'), 
               legend.text = element_text(color = 'Black', size = 12), 
               legend.position = 'top', legend.box = 'vertical', 
               strip.text = element_text(color = 'Black', size = 10), 
               strip.background = element_rect( color = "black", size=0.5, linetype="solid") )
library(R.utils)
figure_path <- '~/Dropbox (Partners HealthCare)/contamination/rev_figures/'
data_path <- '~/Dropbox (Partners HealthCare)/contamination/data/'

# functions
format_ref_af <- function(data){
  data <- data %>%
    mutate(ref_AF = if_else(ref_AF=='None', ref_AF, paste0('ref-AF: (', ref_AF, '%, ',  100-as.numeric(ref_AF),  '%)'))) %>%
    mutate(ref_AF = factor(ref_AF, levels = c( 'None', 'ref-AF: (0%, 100%)', 'ref-AF: (1%, 99%)', 'ref-AF: (5%, 95%)', 'ref-AF: (10%, 90%)', 'ref-AF: (20%, 80%)')))
  return(data)
}
convert_to_long_format <- function(data, version){
  data_long <- data %>%
    select(-(8:13)) %>%
    pivot_longer(., 2:7) %>%
    mutate(ref_AF = str_split(name, '_') %>% map_chr(., 7)) %>%
    format_ref_af(.)
  if(version == 'v3'){
    data_long <- data_long %>%
      mutate(label = factor(label, labels = c('Broad (gnomAD)', 'Sanger (HGDP)', 'Other'), levels = c('broad', 'hgdp', 'external')),
             contamination = if_else(!is.na(recomputed_contamination), recomputed_contamination, contamination),) 
  }else(
    data_long <- data_long %>%
      mutate(label = if_else(label, 'Broad', 'Other'))
  )
  return(data_long)
}

get_n_hom_var_table <- function(data){
  n_var <- data %>%
    select(8:13)%>%
    pivot_longer(., 1:6) %>%
    mutate(ref_AF = str_split(name, '_') %>% map_chr(., 5)) %>%
    format_ref_af(.) %>%
    group_by(ref_AF) %>%
    dplyr::summarize(mean_n_var = mean(value))
  return(n_var)
}

get_age_bins <- function(age){
  bins = case_when(
    age <= 10  ~ '0-10', 
    age <= 20  ~ '10-20', 
    age <= 30  ~ '20-30', 
    age <= 40  ~ '30-40', 
    age <= 50  ~ '40-50', 
    age <= 60  ~ '50-60', 
    age <= 70  ~ '60-70', 
    age <= 80  ~ '70-80', 
    age > 80  ~ '> 80', 
  )
  return(bins)
}
age_bins <- c('0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '> 80')
dp_names <- c('no' = 'No DP filters', '10' = '10 < DP < 100', '20' = '20 < DP < 100', '_all'='')

# Full figure function
get_figure_charr_freemix_by_ref_AF <- function(data, version, var_type, DP, save=TRUE, filter=TRUE, sub_ref_AF=NULL, height=4, width=6){
  long_data <- convert_to_long_format(data, version)
  if(version == 'v3'){
    var_cnt_data <- data %>% 
      get_n_hom_var_table()
    data <- rbind(
      long_data %>% filter(label == 'Broad (gnomAD)'),
      long_data %>% filter(label ==  'Sanger (HGDP)')) %>%
      filter(!release)
  }else{
    var_cnt_data <- data  %>%
      get_n_hom_var_table()
    data <- rbind(
      long_data %>% filter(label == 'Broad'),
      long_data %>% filter(label ==  'Other')) 
  }
  if(filter){
    data <- data %>% filter(ref_AF %in% sub_ref_AF)
    var_cnt_data <- var_cnt_data %>% filter(ref_AF %in% sub_ref_AF)
  }
  figure <- data %>%
    ggplot + aes(x = contamination, y = value, color = label) +
    # labs(x = 'Freemix Score', y = paste0('CHARR (',dp_names[DP], ')'), color = 'Sequencing site') +
    labs(x = 'Freemix Score', y = paste0('CHARR Score'), color = 'Sequencing site') +
    geom_point(alpha=0.8, size = 1) + 
    geom_abline(slope = 1, intercept = 0, lty=2) +
    scale_color_manual(values = c('#045494', 'gray')) +
    themes +  theme_classic() + 
    theme(axis.title = element_text(family = 'Arial', size = 15, face = 'bold'),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          axis.text = element_text(family = 'Arial', size = 12),
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_wrap(.~ref_AF, scale = 'free') 
  print(var_cnt_data)
  # +
  #   geom_label(data = var_cnt_data, aes(label = paste('Mean(N_hom_var) = ', as.character(round(mean_n_var)))),
  #              x = if_else(version == 'v3', 0.05, 0.025), y = Inf, vjust = 1,  color = 'black', size = 3, family = 'Arial') +
  #   annotate(geom='text', x=if_else(version == 'v3', 0.095, 0.045), y=if_else(version == 'v3', 0.08, 0.035),label='y=x', family = 'Arial', size=6)
  
  if(save){
    png(paste0(figure_path, 'S2_gnomad_', version,'_contam_dp', DP, '_', var_type,'.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

# Figure 3 function 
fig3_n_way_comparison <- function(name, save=TRUE, height = 3.5, width = 5.5){
  callrate_filtered_charr <- mixed_data_n <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples_callrate_filtered_charr.csv'), sep='\t')%>%
    mutate(callrate_filtered = if_else(str_detect(AF_type, '\\('), 'Callrate filtered - Yes', 'Callrate filtered - No'),
           main_type = if_else(str_detect(AF_type, '\\('), str_split(AF_type, ' \\(') %>% map_chr(., 1), AF_type)) 
  callrate_filtered_full <- callrate_filtered_charr %>%
    select(1, 4:5, 7, charr=freemix_score) %>%
    mutate(AF_type = 'Freemix Score') %>%
    distinct() %>%
    rbind(callrate_filtered_charr %>% select(1, 4:5,7, 2,3)) %>%
    mutate(callrate_filtered = if_else(str_detect(AF_type, '\\('), T, F),
           main_type = if_else(str_detect(AF_type, '\\('), str_split(AF_type, ' \\(') %>% map_chr(., 1), AF_type)) %>% 
    mutate(contam_rate_percent = paste0(contam_rate*100, '%'))
  callrate_filtered_sum <- callrate_filtered_full %>% 
    filter(AF_type %in% c('Freemix Score', 'local AF (Callrate filtered)')) %>%
    mutate(AF_type = if_else(AF_type == 'local AF (Callrate filtered)', 'CHARR Score', AF_type)) %>%
    group_by(AF_type, contam_rate) %>%
    dplyr::summarize(
      std = sd(charr),
      mean = mean(charr)
    ) %>%
    mutate(y_max = mean + std,
           y_min = mean - std)
  figure <- callrate_filtered_sum %>%
    ggplot + 
    geom_pointrange(aes( x=contam_rate, 
                         y=mean,
                         color=AF_type,
                         fill=AF_type,
                         ymax = y_max,
                         ymin=y_min), size = 0.25)+
    geom_abline(slope = 1, intercept = 0, lty =2) + 
    scale_color_manual(breaks=c('CHARR Score', 'Freemix Score'), labels=c('CHARR Score', 'Freemix Score'), values = c( '#D95F02', '#1B9E77')) +
    scale_y_continuous(label = scales::percent_format(accuracy=0.1), breaks = c(0.005, 0.01, 0.02, 0.05, 0.1), limit = c(0, 0.11)) +
    scale_x_continuous(label = scales::percent_format(accuracy=0.1), breaks = c(0.005, 0.01, 0.02, 0.05, 0.1)) +
    labs(x='True Contamination Rate', y='Contamination estimate', color=NULL, fill = NULL, alpha = 'Callrate filtered') +
    themes + 
    theme(legend.position = c(0.2, 0.8))
  if(save){
    png(paste0(figure_path, name, '_hgdp_n_way_mixed_callrate_filtered_charr.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS1_n_hom_var_RR_above_0 <- function(save,  height = 4, width = 6){
  hom_var <- read_tsv(paste0(data_path, 'gnomad_v3_n_hom_var_final.tsv'))%>%
    mutate(label = if_else(hgdp, 'hgdp', label),
           contamination = if_else(label == 'broad', contamination/100, contamination),
           pass = if_else(c < 0.005, 'CHARR score < 0.5%', 'CHARR score > 0.5%')) %>%
    mutate(label = factor(label, labels = c('Broad (gnomAD)', 'Sanger (HGDP)', 'Other'), levels = c('broad', 'hgdp', 'external')),
           contamination = if_else(!is.na(recomputed_contamination), recomputed_contamination, contamination),) 
  figure <- hom_var %>%
    ggplot + aes(x = n_hom_var_RR_above_0, color=pass, fill=pass) + 
    labs(x = 'Number of homozygous variants with RR > 0',  y=NULL, color = NULL, fill = NULL) +
    geom_histogram(bins = 50, alpha=0.5) + 
    scale_x_log10(label = comma)+ 
    scale_color_brewer(palette = 'Dark2') + 
    scale_fill_brewer(palette = 'Dark2') +
    themes + theme_classic() +   
    theme(axis.title = element_text(family = 'Arial', size = 12),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = c(0.8,0.8),
          axis.text = element_text(family = 'Arial', size = 12),
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) 
  if(save){
    png(paste0(figure_path, 'S1_gnomad_v3_contam_n_hom_var_RR_above_0.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off() 
  }
  return(figure) 
}

# correlation figure function
get_corr_df <- function(long_data, version){
    data <- long_data %>%
      filter(!is.na(label)) %>%
      group_by(label, ref_AF) %>% 
      dplyr::summarize(tidy(cor.test(contamination, value)))
  return(data)
}

get_figure_corr_by_ref_AF <- function(data, version, var_type, DP, save=TRUE){
  if(version == 'v3'){
    data <- data %>%
      filter(label != 'Other')
  }else{
    # data <- rbind(
    #   long_data %>% filter(label == 'Broad'),
    #   long_data %>% filter(label ==  'Other')) %>%
    #   filter(releasable & ((ref_AF == 'ref-AF: (0%, 100%)') | value < 0.06))
    # colors <- c('#045494', 'gray')
    # var_lab <- 0.025
    # annt_x <- 0.045
    # annt_y <- 0.035
  }
  
  figure <- data %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    labs(title = paste0('CHARR', dp_names[DP],' ~ Freemix score (', toupper(var_type), ')'), x = 'Reference allele frequency', y = "Pearson's correlation", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = conf.low, ymax = conf.high, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_wrap(~dp_filter)+
    annotate(geom='text', x=Inf, y=-Inf,label=expression(p.values<10^-100), family = 'Arial', size=4, hjust = 1, vjust=-0.8)
  
  if(save){
    png(paste0(figure_path, 'gnomad_', version, '_contam_dp', DP, '_', var_type, '_correlation_for_supplement.png'), height = 4, width = 10, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

# Linear models
get_lm_df1 <- function(long_data, version='v3', intercept=T){
    if(intercept){
      lm_df1 <- long_data %>%
        filter(!is.na(label)) %>%
        group_by(label, ref_AF) %>%
        do(tidy(lm(value ~ contamination, data = .)))
    }else{
      lm_df1 <- long_data %>%
        filter(!is.na(label)) %>%
        group_by(label, ref_AF) %>%
        do(tidy(lm(value ~ contamination + 0, data = .)))
    }
  return(lm_df1)
}

get_lm_df2 <- function(long_data, version='v3', intercept=T){
    if(intercept){
      lm_df2 <- long_data %>%
        filter(!is.na(label)) %>%
        group_by(label, ref_AF) %>%
        do(glance(lm(value ~ contamination, data = .)))
    }else{
      lm_df2 <- long_data %>%
        filter(!is.na(label)) %>%
        group_by(label, ref_AF) %>%
        do(glance(lm(value ~ contamination + 0, data = .)))
    }
  return(lm_df2)
}

get_lm_df <- function(lm_df1, lm_df2){
  lm_df <- lm_df1 %>%
    select(1:5, 8) %>%
    rbind(., lm_df2 %>% select(1:3, 15) %>% mutate(term = 'R-squared', std.error = NA) %>% select(label, ref_AF, term, estimate = r.squared, std.error, dp_filter)) %>%
    mutate(term = factor(if_else(term=='(Intercept)', 'Intercept', if_else(term=='contamination', 'Slope', term)), levels = c('Intercept', 'Slope', 'R-squared'))) %>%
    mutate(ref_line = if_else(term == 'Intercept', 0, 1))
  return(lm_df)
}

figS3A_get_figure_lm_by_ref_AF <- function(data, version, var_type, DP, name, intercept=TRUE, save=TRUE){
  figure <-data %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    labs(x=NULL, y = "Estimates", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    geom_hline(aes(yintercept = ref_line), lty=2) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(term~dp_filter, scale = 'free')
  if(save){
    suffix = if_else(intercept, '', '_no_intercept')
    png(paste0(figure_path, name, '_gnomad_', version,'_contam_dp', DP, '_', var_type,'_linear_model', suffix,'_for_supplement.png'), 
        height = if_else(intercept, 6, 4), 
        width = 10, 
        units = 'in', 
        res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS3B_get_figure_lm_by_ref_AF <- function(data, version, var_type, DP, name, intercept=TRUE, save=TRUE){
  figure <-data %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    labs(x = 'Reference allele frequency', y = "Estimates", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    geom_hline(aes(yintercept = ref_line), lty=2) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(term~dp_filter, scale = 'free')
  if(save){
    suffix = if_else(intercept, '', '_no_intercept')
    png(paste0(figure_path, name, '_gnomad_', version,'_contam_dp', DP, '_', var_type,'_linear_model', suffix,'_for_supplement.png'), 
        height = if_else(intercept, 6, 4), 
        width = 10, 
        units = 'in', 
        res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS3_get_figure_lm_by_dp_ref_AF <- function(long_data_no_dp, long_data_dp10, long_data_dp20, version, var_type, name, intercept=TRUE, save=TRUE){
  lm_df1 <- as.data.frame(rbind(get_lm_df1(long_data_no_dp, version, intercept = intercept), 
                                get_lm_df1(long_data_dp10, version, intercept = intercept), 
                                get_lm_df1(long_data_dp20, version, intercept = intercept))) %>% 
    mutate(dp_filter = factor(rep(c('No DP filter', '10 < DP < 100', '20 < DP < 100'), each=if_else(version=='v3', if_else(intercept, 24, 12), if_else(intercept, 24, 12))), 
                              levels =c('No DP filter', '10 < DP < 100', '20 < DP < 100') ))%>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None')))
  lm_df2 <- as.data.frame(rbind(get_lm_df2(long_data_no_dp, version, intercept = intercept), 
                                get_lm_df2(long_data_dp10, version, intercept = intercept), 
                                get_lm_df2(long_data_dp20, version, intercept = intercept))) %>% 
    mutate(dp_filter = factor(rep(c('No DP filter', '10 < DP < 100', '20 < DP < 100'), each=if_else(version=='v3', 12, 6)), 
                              levels =c('No DP filter', '10 < DP < 100', '20 < DP < 100') ))%>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None')))
  lm_df <- get_lm_df(lm_df1, lm_df2)
  if(name == 'S3A'){
    figure <- figS3A_get_figure_lm_by_ref_AF(lm_df, version, var_type, '_all', name, intercept, save)
  }else if(name == 'S3B'){
    figure <- figS3B_get_figure_lm_by_ref_AF(lm_df, version, var_type, '_all', name, intercept, save)
  }
  return(figure)
}

figS4_v3_snp_indel_all_lm <- function(long_data, version, name, save, height = 6, width = 10){
  lm_df1 <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, var_type) %>%
    do(tidy(lm(value ~ contamination, data = .)))
  lm_df2 <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, var_type) %>%
    do(glance(lm(value ~ contamination, data = .)))
  lm_df <- lm_df1 %>%
    select(1:6) %>%
    rbind(., lm_df2 %>% select(1:4) %>% mutate(term = 'R-squared', std.error = NA) %>% select(label, ref_AF,var_type, term, estimate = r.squared, std.error)) %>%
    mutate(term = factor(if_else(term=='(Intercept)', 'Intercept', if_else(term=='contamination', 'Slope', term)), levels = c('Intercept', 'Slope', 'R-squared'))) %>%
    mutate(ref_line = if_else(term == 'Intercept', 0, 1))
  figure <-lm_df %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    # labs(title = paste0('CHARR (20 < DP < 100) ~ Freemix score'), x = 'Reference allele frequency', y = "LR estimates", color = 'Sequencing sites') +
    labs(title = paste0('CHARR score ~ Freemix score'), x = 'Reference allele frequency', y = "Estimates", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    geom_hline(aes(yintercept = ref_line), lty=2) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(term~var_type, scale = 'free')
  if(save){
    png(paste0(figure_path, name, '_gnomad_', version,'_contam_dp20_linear_model.png'), 
        height = height, 
        width = width, 
        units = 'in', 
        res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS5_v2_ref_af_source <- function(long_data, var_cnt_data, name, save=TRUE, height = 6, width = 9){
  data <- rbind(
    long_data %>% filter(label == 'Broad'),
    long_data %>% filter(label ==  'Other'))
  figure <- data %>%
    filter(ref_AF != 'ref-AF: (5%, 95%)') %>%
    ggplot + aes(x = contamination, y = value, color = label) +
    # labs(x = 'Freemix Score', y = paste0('CHARR (20 < DP < 100)'), color = 'Sequencing site') +
    labs(x = 'Freemix Score', y = paste0('CHARR Score'), color = 'Sequencing site') +
    geom_point(alpha=0.8, size = 1) + 
    geom_abline(slope = 1, intercept = 0, lty=2) +
    scale_color_manual(values = c('#045494', 'gray')) +
    themes +  theme_classic() + 
    theme(axis.title = element_text(family = 'Arial', size = 15, face = 'bold'),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position='top',
          axis.text = element_text(family = 'Arial', size = 12),
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(AF_source~ref_AF, scale = 'free') 
  # +
  # geom_label(data = var_cnt_data %>% filter(ref_AF != 'ref-AF: (5%, 95%)'), aes(label = paste('Mean(N_hom_var) = ', as.character(round(mean_n_var)))),x = 0.025, y = Inf, vjust = 1,  color = 'black', size = 3, family = 'Arial') +
  # annotate(geom='text', x=0.045, y=0.035,label='y=x', family = 'Arial', size=6)
  if(save){
    png(paste0(figure_path, name, '_gnomad_v2_contam_dp20_ref_af_source.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS6_v2_ref_af_source_lm <- function(long_data, name, save=TRUE, height = 6, width = 10){
  lm_df1 <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, AF_source) %>%
    do(tidy(lm(value ~ contamination, data = .)))
  lm_df2 <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, AF_source) %>%
    do(glance(lm(value ~ contamination, data = .)))
  lm_df <- lm_df1 %>%
    select(1:6) %>%
    rbind(., lm_df2 %>% select(1:4) %>% mutate(term = 'R-squared', std.error = NA) %>% select(label, ref_AF, AF_source, term, estimate = r.squared, std.error)) %>%
    mutate(term = factor(if_else(term=='(Intercept)', 'Intercept', if_else(term=='contamination', 'Slope', term)), levels = c('Intercept', 'Slope', 'R-squared')))%>%
    mutate(ref_line = if_else(term == 'Intercept', 0, 1))
  figure <-lm_df %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    # labs(title = paste0('CHARR (20 < DP < 100) ~ Freemix score'), x = 'Reference allele frequency', y = "LR estimates", color = 'Sequencing sites') +
    labs(title = paste0('CHARR score ~ Freemix score'), x = 'Reference allele frequency', y = "Estimates", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    geom_hline(aes(yintercept = ref_line), lty=2) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(term~AF_source, scale = 'free')
  if(save){
    png(paste0(figure_path, name, '_gnomad_v2_contam_dp20_ref_af_source_linear_model.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS7_v3_remove_age_gene <- function(long_data, var_cnt_data, name, save=TRUE, height = 6, width = 9){
  data <- rbind(
    long_data %>% filter(label == 'Broad (gnomAD)'),
    long_data %>% filter(label ==  'Sanger (HGDP)'))
  figure <- data %>%
    filter(ref_AF != 'ref-AF: (5%, 95%)') %>%
    ggplot + aes(x = contamination, y = value, color = label) +
    # labs(x = 'Freemix Score', y = paste0('CHARR (20 < DP < 100)'), color = 'Sequencing site') +
    labs(x = 'Freemix Score', y = paste0('CHARR Score'), color = 'Sequencing site') +
    geom_point(alpha=0.8, size = 1) + 
    geom_abline(slope = 1, intercept = 0, lty=2) +
    scale_color_manual(values = c('#045494', 'gray')) +
    themes +  theme_classic() + 
    theme(axis.title = element_text(family = 'Arial', size = 15, face = 'bold'),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position='top',
          axis.text = element_text(family = 'Arial', size = 12),
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(source~ref_AF, scale = 'free') 
  print(var_cnt_data)
  # +
  #   geom_label(data = var_cnt_data %>% filter(ref_AF != 'ref-AF: (5%, 95%)'), aes(label = paste('Mean(N_hom_var) = ', as.character(round(mean_n_var)))),x = 0.05, y = Inf, vjust = 1,  color = 'black', size = 3, family = 'Arial') +
  #   annotate(geom='text', x=0.095, y=0.08,label='y=x', family = 'Arial', size=6)
  
  if(save){
    png(paste0(figure_path, name, '_gnomad_v3_contam_dp20_age_gene.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS8_v3_remove_age_gene_lm <- function(long_data, name, save=TRUE, height = 6, width = 10){
  lm_df1 <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, source) %>%
    do(tidy(lm(value ~ contamination, data = .)))
  lm_df2 <- long_data %>%
    filter( !is.na(label)) %>%
    group_by(label, ref_AF, source) %>%
    do(glance(lm(value ~ contamination, data = .)))
  lm_df <- lm_df1 %>%
    select(1:6) %>%
    rbind(., lm_df2 %>% select(1:4) %>% mutate(term = 'R-squared', std.error = NA) %>% select(label, ref_AF, source, term, estimate = r.squared, std.error)) %>%
    mutate(term = factor(if_else(term=='(Intercept)', 'Intercept', if_else(term=='contamination', 'Slope', term)), levels = c('Intercept', 'Slope', 'R-squared'))) %>%
    mutate(ref_line = if_else(term == 'Intercept', 0, 1))
  figure <-lm_df %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    # labs(title = paste0('CHARR (20 < DP < 100) ~ Freemix score'), x = 'Reference allele frequency', y = "LR estimates", color = 'Sequencing sites') +
    labs(title = paste0('CHARR score ~ Freemix score'), x = 'Reference allele frequency', y = "Estimates", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    geom_hline(aes(yintercept = ref_line), lty=2) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(term~source, scale = 'free')
  if(save){
    png(paste0(figure_path, name, '_gnomad_v3_contam_dp20_age_gene_linear_model.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS9_v2_v3_age_vs_charr <- function(v2_long_data, v3_long_data, name, save=TRUE, height = 4, width = 7.5){
  v2_long_data <- v2_long_data %>%
    filter(ref_AF == 'ref-AF: (10%, 90%)') %>%
    mutate(age_bin = factor(get_age_bins(age), levels = age_bins),
           version = 'gnomAD v2 (WES)') %>%
    select(value, ref_AF, age_bin, version)
  v3_long_data <- v3_long_data %>%
    filter(ref_AF == 'ref-AF: (10%, 90%)') %>%
    mutate(age_bin = factor(get_age_bins(age), levels = age_bins),
           version = 'gnomAD v3 (WGS)') %>%
    select(value, ref_AF, age_bin, version)
  print('v2:')
  print(table(v2_long_data[v2_long_data$ref_AF == 'ref-AF: (10%, 90%)','age_bin']))
  print('v3:')
  print(table(v3_long_data[v3_long_data$ref_AF == 'ref-AF: (10%, 90%)','age_bin']))
  dp20 <- rbind(v3_long_data, v2_long_data) 
  figure <- dp20 %>%
    filter(!is.na(age_bin)) %>%
    filter(ref_AF == 'ref-AF: (10%, 90%)') %>%
    ggplot + aes(x = age_bin, y = value) +
    # labs(x = 'Age', y = 'CHARR (20 < DP < 100)', color = 'Age') +
    labs(x = 'Age', y = 'CHARR Score', color = 'Age') +
    geom_boxplot() + 
    coord_flip() +
    scale_y_log10() +
    themes +  theme_classic() +
    theme(axis.title = element_text(family = 'Arial', size = 15, face = 'bold'),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(.~version, scale = 'free')
  if(save){
    png(paste0(figure_path, name, '_gnomad_contam_dp20_age_vs_charr.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS10_hgdp_freemix_before_after <-function(hgdp, name, save=TRUE, height = 4, width = 6){
  figure <- hgdp %>%
    filter(!is.na(recomputed_contamination)) %>%
    ggplot + aes(x = contamination, y = recomputed_contamination, color = label) +
    labs(x = 'Old', y = 'Recomputed', color = 'Sequencing site') +
    geom_point(alpha=0.8, size = 1) + 
    geom_abline(slope = 1, intercept = 0, lty=2, size = 1) +
    scale_color_manual(values = c('#045494', 'orange', 'gray') ) +
    themes +  theme_classic() + 
    theme(axis.title = element_text(family = 'Arial', size = 15, face = 'bold'),
          legend.position='None',
          axis.text = element_text(family = 'Arial', size = 12),
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) 
  # +
  #   annotate(geom='text', x=0.05, y=0.045, label='y=x', family = 'Arial', size=6)
  if(save){
    png(paste0(figure_path, name, '_gnomad_v3_hgdp_verifybamID_comparison.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS11_hgdp_charr_freemix_before_after <-function(hgdp, name, save=TRUE, height = 4, width = 7.5){
  pop_colors['Broad (gnomAD)'] <- '#D6CFC7'
  pop_colors['mid']  <- pop_colors['mde'] 
  pop_names['mid']  <- pop_names['mde']
  pop_colors <- pop_colors[c('afr', 'amr', 'eas', 'mid', 'nfe', 'oth', 'sas', 'Broad (gnomAD)')]
  hgdp <- rbind(
    hgdp %>% filter(label ==  'Broad (gnomAD)'),
    hgdp %>% filter(label ==  'Sanger (HGDP)')
  ) %>% mutate(Version = 'Old', contamination=contamination) %>%
    rbind(
      rbind(
        hgdp %>% filter(label ==  'Broad (gnomAD)'),
        hgdp %>% filter(label ==  'Sanger (HGDP)')
      )%>%
        mutate(Version = 'Recomputed', contamination = if_else(!is.na(recomputed_contamination), recomputed_contamination, contamination))
    )
  figure <- hgdp %>%
    mutate(pop = if_else(label == 'Broad (gnomAD)', 'Broad (gnomAD)', pop)) %>%
    ggplot + aes(x = contamination, y = mean_AB_snp_biallelic_af_adjust_10, pch = label, color = pop) +
    # labs(x = 'Freemix Score', y = 'CHARR (20 < DP < 100)', pch = 'Sequencing site', color = 'Population') +
    labs(x = 'Freemix Score', y = 'CHARR Score', pch = 'Sequencing site', color = 'Population') +
    geom_point(alpha=0.8, size = 2) + 
    geom_abline(slope = 1, intercept = 0, lty=2, size = 1) +
    scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'oth', 'sas'), labels = c('AFR', 'AMR', 'EAS', 'MID', 'NFE', 'OTH', 'SAS')) +
    themes +  theme_classic() + 
    theme(axis.title = element_text(family = 'Arial', size = 15, face = 'bold'),
          legend.text = element_text(family = 'Arial', size = 10),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top', legend.box="vertical",
          axis.text = element_text(family = 'Arial', size = 12),
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) +
    # annotate(geom='text', x=0.045, y=0.05, label='y=x', family = 'Arial', size=6) + 
    facet_grid(~Version) 
  if(save){
    png(paste0(figure_path, name, '_gnomad_v3_hgdp_verifybamID_comparison.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

fig12_hgdp_sample_selection <- function(name, save=TRUE, height = 4, width = 7.5){
  id <- c('0281', '0021', '0105', '0118', '0341', '0768', '1096', '1190', '1243', '0714', '1371', '0529', '0669', '1396', '1385', '0689', '0688', '0610', '0726', '0690', '0931', '0944', '1092', '1406', '0905', '0875', '0703', '1009', '1050', '0845')
  hgdp_id <- paste0('HGDP0', id)
  sub <- read_tsv(paste0(data_path,'gnomad_v3_contamination_estimate_full_snp_bi_minimum_vds_af_adjust_dp20_100.tsv')) %>% 
    select(s, pop, recomputed_contamination, mean_AB_snp_biallelic_af_adjust_10) %>%
    filter(s %in% hgdp_id)
  figure1 <- sub %>%
    ggplot + 
    aes(x = pop, 
        # y = recomputed_contamination,
        y = mean_AB_snp_biallelic_af_adjust_10,
        color = pop) + 
    labs(x = 'Population', y='CHARR Score') +
    geom_point() + 
    scale_y_log10() +
    scale_x_discrete(labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'))) + 
    scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'), labels = pop_names[c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth')]) + 
    themes + 
    theme(legend.position = 'None')
  figure2 <- sub %>%
    ggplot + 
    aes(x = pop, 
        y = recomputed_contamination,
        color = pop) + 
    labs(x = 'Population', y='Freemix Score') +
    geom_point() + 
    scale_y_log10() +
    scale_x_discrete(labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'))) + 
    scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth')) + 
    themes + 
    theme(legend.position = 'None')
  figure = ggpubr::ggarrange(figure1, figure2, ncol = 2,
                             # labels = c('CHARR', 'Freemix Score'),
                             # common.legend = TRUE,
                             # vjust = 0, hjust = 0, 
                             font.label = list(size = 13, color = "black", face = "bold", family = 'Arial')
  )
  if(save){
    png(paste0(figure_path, name, '_selected_samples_for_simulation_analysis_freemix.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
  
}

figS13_sample_decontamination <- function(name, save=TRUE, height = 4, width = 7.5){
  decontam_data <- read.csv(paste0(data_path, 'charr_simulation_hgdp_decontaminated_full_30_samples.csv'), sep='\t')
  decontam_long <- decontam_data %>% 
    select(s, pop, vbid=recomputed_contamination, charr=mean_AB_snp_biallelic_af_adjust_10) %>% 
    rbind(decontam_data %>% select(s, pop, vbid=contam_free_vbid, charr=contam_free_charr)) %>% 
    mutate(type=factor(rep(c('Original', 'Decontaminated'), each=30), levels=c('Original', 'Decontaminated')))
  figure <- decontam_long %>% 
    ggplot + aes(x=vbid, y=charr, color=pop, pch=type) +
    geom_point(size=2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x='Freemix Score', y='CHARR Score', color='Population', pch='Type') + 
    scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'), labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'))) +
    themes + 
    theme(legend.position = 'right')
  
  if(save){
    png(paste0(figure_path, name, '_hgdp_decontaminated_30_samples_freemix_score_vs_charr_decontaminated_vs_original.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()  
  }
  return(figure)
}

figS14_2_way_color_by_truth <- function(name, save=TRUE, height = 5, width = 7.5, type='gnomAD'){
  mixed_data_2 <- read.csv(paste0(data_path, 'charr_simulation_hgdp_2_way_mixed_samples_full_150_samples.csv'), sep='\t')
  if(type != 'gnomAD'){
    mixed_data_2 <- mixed_data_2 %>%
      mutate(charr = charr_local)
  }
  p <- mixed_data_2 %>% 
    ggplot + aes(x=freemix_score, 
                 y=charr, 
                 color=factor(contam_rate, levels=c(0.005, 0.01, 0.02, 0.05, 0.1), labels = c('0.5%', '1.0%', '2.0%', '5.0%', '10.0%'))) +
    geom_point(size=2 )+
    geom_abline(slope = 1, intercept = 0, lty =2) + 
    scale_color_manual(values = c('#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695')) +
    scale_y_continuous(label = scales::percent_format(accuracy=0.1)) +
    scale_x_continuous(label = scales::percent_format(accuracy=0.1)) +
    labs(x='Freemix Score', y='CHARR Score', color='True Contamination Rate') + themes
  if(save){
    png(paste0(figure_path, name, '_hgdp_2_way_mixed_150_samples_', type ,'_AF.png'), height = height, width = width, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  return(p)
}

figS15_two_way_mixing_het_hom_ratio <- function(name, save=TRUE, height = 4, width = 7.5){
  contam_rates <- c('0.5%', '1%', '2%', '5%', '10%')
  names(contam_rates) <- c(0.005, 0.01, 0.02, 0.05, 0.1)
  mixed_data <- read.csv(paste0(data_path, 'charr_simulation_hgdp_2_way_mixed_samples_full_150_samples.csv'), sep='\t')
  hom_var <- read.csv(paste0(data_path, 'charr_simulation_hgdp_2_way_mixed_samples_full_150_samples_hom_var_summary.csv'), sep='\t')
  hom_var <- hom_var %>% 
    merge(mixed_data, by= c("original", "contaminant", "contam_rate"))
  figure <- hom_var %>%
    ggplot + aes(y=n_het_var/n_hom_var, x=pop_original, color=pop_original, group=interaction(pop_original,contam_rate)) +
    geom_boxplot() +
    labs(x = NULL, y='Het/Hom Ratio', color = 'Population (original)')+
    themes + 
    # scale_color_manual(values = c('#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695')) +
    theme(legend.position = 'top',
          legend.background = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_x_discrete(labels= NULL) + 
    scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'), labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth')))+ 
    facet_grid(~contam_rate, scales = 'free', labeller = labeller(contam_rate = contam_rates)) 
  if(save){
    png(paste0(figure_path, name, '_hgdp_2_way_mixed_150_samples_het_hom_ratio_contaminant_pop.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS16_n_way_mixing <- function(name, save=TRUE, height = 4, width = 7.5){
  callrate_filtered_charr <- mixed_data_n <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples_callrate_filtered_charr.csv'), sep='\t')%>%
    mutate(callrate_filtered = if_else(str_detect(AF_type, '\\('), 'Callrate filtered - Yes', 'Callrate filtered - No'),
           main_type = if_else(str_detect(AF_type, '\\('), str_split(AF_type, ' \\(') %>% map_chr(., 1), AF_type)) 
  figure <- callrate_filtered_charr %>% 
    filter(main_type %in% c('gnomAD AF', 'local AF')) %>%
    ggplot + aes(x=freemix_score, 
                 y=charr, 
                 color=factor(contam_rate, levels=c(0.005, 0.01, 0.02, 0.05, 0.1), labels = c('0.5%', '1.0%', '2.0%', '5.0%', '10.0%'))) +
    geom_point(size=2 )+
    geom_abline(slope = 1, intercept = 0, lty =2) + 
    theme_classic()+ 
    scale_color_manual(values = c('#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695')) +
    scale_y_continuous(label = scales::percent_format(accuracy=0.1), breaks = c(0.005, 0.01, 0.02, 0.05, 0.1), minor_breaks = NULL, limits = c(0, 0.13) ) +
    scale_x_continuous(label = scales::percent_format(accuracy=0.1), breaks = c(0.005, 0.01, 0.02, 0.05, 0.1), minor_breaks = NULL, limits = c(0, 0.13) ) +
    labs(x='Freemix Score', y='CHARR Score', color='True Contamination Rate') + themes +
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=90, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size=12),
          panel.grid.major = element_line(linetype="solid",size=0.3),
          legend.position = 'top',
          legend.text = element_text(size = 12)) +
    facet_grid(callrate_filtered ~ main_type) 
  if(save){
    png(paste0(figure_path, name, '_hgdp_n_way_mixed_callrate_filtered_charr_vs_vbid.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS17_n_way_box_plot <- function(name, save=TRUE, height = 3.5, width = 7){
  # callrate_filtered_charr <- mixed_data_n <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples_callrate_filtered_charr.csv'), sep='\t')%>%
  callrate_filtered_charr <- mixed_data_n <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples_callrate_filtered_charr_gnomad_0.csv'), sep='\t')%>%
    mutate(callrate_filtered = if_else(str_detect(AF_type, '\\('), 'Callrate filtered - Yes', 'Callrate filtered - No'),
           main_type = if_else(str_detect(AF_type, '\\('), str_split(AF_type, ' \\(') %>% map_chr(., 1), AF_type)) 
  callrate_filtered_full <- callrate_filtered_charr %>%
    select(1, 4:5, 7, charr=freemix_score) %>%
    mutate(AF_type = 'Freemix Score') %>%
    distinct() %>%
    rbind(callrate_filtered_charr %>% select(1, 4:5,7, 2,3)) %>%
    mutate(callrate_filtered = if_else(str_detect(AF_type, '\\('), T, F),
           main_type = if_else(str_detect(AF_type, '\\('), str_split(AF_type, ' \\(') %>% map_chr(., 1), AF_type)) %>% 
    mutate(contam_rate_percent = paste0(contam_rate*100, '%'))
  contam_rates <- c('0.5%', '1%', '2%', '5%', '10%')
  names(contam_rates) <- c(0.005, 0.01, 0.02, 0.05, 0.1)
  figure <- callrate_filtered_full  %>%
    filter(main_type %in% c('gnomAD AF', 'local AF', 'Freemix Score')) %>%
    filter(AF_type != 'gnomAD AF (Callrate filtered)') %>%
    mutate(contam_rate_percent = factor(contam_rate, levels = c(0.005, 0.01, 0.02, 0.05, 0.1),
                                        labels = c('0.5%', '1%', '2%', '5%', '10%'))) %>%
    mutate(AF_type = if_else(AF_type == 'local AF (Callrate filtered)', 'local AF: 100% Callrate', AF_type)) %>%
    mutate(AF_type = if_else(AF_type == 'Freemix Score', AF_type, paste0('CHARR Score \n (', AF_type, ')'))) %>%
    mutate(AF_type = factor(AF_type, levels = c('Freemix Score', 'CHARR Score \n (gnomAD AF)', 'CHARR Score \n (local AF)', 'CHARR Score \n (local AF: 100% Callrate)'))) %>%
    ggplot + aes(y=charr, 
                 x=AF_type, 
                 color=AF_type,
                 fill = AF_type,
                 group=AF_type) +
    geom_boxplot(alpha = 0.5) +
    labs(x=NULL, y ='Contamination estimate', alpha = 'Callrate filtered', color = NULL, fill= NULL) + 
    scale_y_continuous(label = scales::percent_format(accuracy=0.1)) +
    scale_color_brewer(palette='Dark2') +
    scale_fill_brewer(palette='Dark2') +
    facet_wrap(~contam_rate, ncol=3, scales = 'free', labeller = labeller(contam_rate = contam_rates)) +
    geom_hline(data = callrate_filtered_full, aes(yintercept = contam_rate), lty=2) +
    themes + 
    theme(legend.position = c(.83,.27), 
          legend.text = element_text(size=10), 
          legend.margin = margin(0,0,0,0, unit="cm"),
          legend.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.key.size = unit(1.1, 'cm')) 
  if(save){
    png(paste0(figure_path, name, '_hgdp_n_way_mixed_callrate_charr_boplot.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS18_2_way_vs_n_way <- function(name, save=TRUE, height = 5, width = 7.5, type = 'gnomAD'){
  mixed_data_2 <- read.csv(paste0(data_path, 'charr_simulation_hgdp_2_way_mixed_samples_full_150_samples.csv'), sep='\t')
  mixed_data_n <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples.csv'), sep='\t')
  mixed_full <- mixed_data_2 %>%
    select(-c('contaminant', 'pop_contaminant')) %>%
    rbind(mixed_data_n) %>%
    mutate(result_type = c(rep('2_way', 150), rep('n_way', 150)))
  if(type != 'gnomAD'){
    mixed_full <- mixed_full %>%
      mutate(charr = charr_local)
  }
  p <- mixed_full %>% 
    ggplot + aes(x=freemix_score, 
                 y=charr, 
                 color=result_type,
                 pch=result_type) +
    geom_point(size=2) +
    scale_color_brewer(palette = 'Dark2') +
    scale_y_continuous(label = scales::percent_format(accuracy=0.1)) +
    scale_x_continuous(label = scales::percent_format(accuracy=0.1)) +
    labs(x='Freemix Score', y='CHARR Score', pch='Mixing type', color = 'Mixing type') + 
    # scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'), labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'))) +
    
    facet_wrap(~contam_rate, ncol=3, scales = 'free') +
    geom_hline(data = mixed_full, aes(yintercept = contam_rate), lty=2) +
    geom_vline(data = mixed_full, aes(xintercept = contam_rate), lty=2) + 
    themes + 
    theme(legend.position = c(.8,.4), legend.margin = margin(0,0,0,0, unit="cm"), legend.background = element_blank()) 
  if(save){
    png(paste0(figure_path, name, '_hgdp_2_way_vs_n_way_mixed_',type,'_AF.png'), height = height, width = width, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  return(p)
}


#################################################### ARCHIVED #############################################################
fig3_n_way_comparison_archived <- function(name, save=TRUE, height = 5, width = 7.5, type='gnomAD'){
  mixed_data_n <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples.csv'), sep='\t')
  library(dplyr)
  print(mixed_data_n %>%
          dplyr::summarize(tidy(cor.test(freemix_score, charr))))
  if(type != 'gnomAD'){
    mixed_data_n <- mixed_data_n %>%
      mutate(charr= charr_local)
  }
  p <- mixed_data_n %>% 
    ggplot + aes(x=freemix_score, 
                 y=charr, 
                 color=factor(contam_rate, levels=c(0.005, 0.01, 0.02, 0.05, 0.1), labels = c('0.5%', '1.0%', '2.0%', '5.0%', '10.0%'))) +
    geom_point(size=2 )+
    geom_abline(slope = 1, intercept = 0, lty =2) + 
    scale_color_manual(values = c('#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695')) +
    scale_y_continuous(label = scales::percent_format(accuracy=0.1)) +
    scale_x_continuous(label = scales::percent_format(accuracy=0.1)) +
    labs(x='Freemix Score', y='CHARR', color='True Contamination Rate') 
  if(save){
    png(paste0(figure_path, name, '_hgdp_n_way_mixed_150_samples_', type,'_AF.png'), height = height, width = width, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  return(p)
}

figS8_v3_snp_indel_all_archived <- function(data, version, name, DP, save=TRUE, filter=TRUE, sub_ref_AF=NULL, height=4, width=6){
  get_n_hom_var_table <- function(data){
    n_var <- data %>%
      select(8:13, 22)%>%
      pivot_longer(., 1:6) %>%
      mutate(ref_AF = str_split(name, '_') %>% map_chr(., 5)) %>%
      format_ref_af(.) %>%
      group_by(ref_AF, var_type) %>%
      dplyr::summarize(mean_n_var = mean(value))
    return(n_var)
  }
  get_figure_charr_freemix_by_ref_AF <- function(data, version, name, save, filter, sub_ref_AF, height, width){
    long_data <- convert_to_long_format(data, version)
    if(version == 'v3'){
      var_cnt_data <- data %>% 
        get_n_hom_var_table()
      data <- rbind(
        long_data %>% filter(label == 'Broad (gnomAD)'),
        long_data %>% filter(label ==  'Sanger (HGDP)')) %>%
        filter(!release)
    }else{
      var_cnt_data <- data %>%
        get_n_hom_var_table()
      data <- rbind(
        long_data %>% filter(label == 'Broad'),
        long_data %>% filter(label ==  'Other'))
    }
    if(filter){
      data <- data %>% filter(ref_AF %in% sub_ref_AF)
      var_cnt_data <- var_cnt_data %>% filter(ref_AF %in% sub_ref_AF)
    }
    figure <- data %>%
      ggplot + aes(x = contamination, y = value, color = label) +
      # labs(x = 'Freemix Score', y = paste0('CHARR  (20 < DP < 100) '), color = 'Sequencing site') +
      labs(x = 'Freemix Score', y = paste0('CHARR'), color = 'Sequencing site') +
      geom_point(alpha=0.8, size = 1) + 
      geom_abline(slope = 1, intercept = 0, lty=2) +
      scale_color_manual(values = c('#045494', 'gray')) +
      themes +  theme_classic() + 
      theme(axis.title = element_text(family = 'Arial', size = 15, face = 'bold'),
            legend.text = element_text(family = 'Arial', size = 12),
            legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
            legend.position = 'top',
            axis.text = element_text(family = 'Arial', size = 12),
            strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
      facet_grid(~var_type, scale = 'free') 
    # +
    #   geom_label(data = var_cnt_data, aes(label = paste('Mean(N_hom_var) = ', as.character(round(mean_n_var)))),
    #              x = if_else(version == 'v3', 0.05, 0.025), y = Inf, vjust = 1,  color = 'black', size = 3, family = 'Arial') +
    #   annotate(geom='text', x=if_else(version == 'v3', 0.095, 0.045), y=if_else(version == 'v3', 0.08, 0.035),label='y=x', family = 'Arial', size=6)
    print(var_cnt_data)
    
    if(save){
      png(paste0(figure_path, name, '_gnomad_', version,'_contam_dp20_vs_freemix.png'), height = height, width = width, units = 'in', res = 300)
      print(figure)
      dev.off()
    }
    return(figure)
  }
  get_figure_charr_freemix_by_ref_AF(data, version, name, save, filter, sub_ref_AF, height, width)
}

figS10_v3_snp_indel_all_corr_archived <- function(long_data, version, name, save, height=4, width=10){
  data <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, var_type) %>% 
    dplyr::summarize(tidy(cor.test(contamination, value)))
  
  figure <- data %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    # labs(title = paste0('CHARR (20 < DP < 100)  ~ Freemix score'), x = 'Reference allele frequency', y = "Pearson's correlation", color = 'Sequencing sites') +
    labs(title = paste0('CHARR ~ Freemix score'), x = 'Reference allele frequency', y = "Pearson's correlation", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = conf.low, ymax = conf.high, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_wrap(~var_type)+
    annotate(geom='text', x=Inf, y=-Inf,label=expression(p.values<10^-100), family = 'Arial', size=4, hjust = 1, vjust=-0.8)
  
  if(save){
    png(paste0(figure_path, name,'_gnomad_', version, '_contam_dp20_correlation.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS13_v2_ref_af_source_corr_archived <- function(long_data, name, save=TRUE, height = 4, width = 7.5){
  data <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, AF_source) %>% 
    dplyr::summarize(tidy(cor.test(contamination, value)))
  figure <- data %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    # labs(title = paste0('CHARR (20 < DP < 100) ~ Freemix score'), x = 'Reference allele frequency', y = "Pearson's correlation", color = 'Sequencing sites') +
    labs(title = paste0('CHARR ~ Freemix score'), x = 'Reference allele frequency', y = "Pearson's correlation", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = conf.low, ymax = conf.high, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(~AF_source, scale = 'free') +
    annotate(geom='text', x=Inf, y=-Inf,label=expression(p.values<10^-100), family = 'Arial', size=4, hjust = 1, vjust=-0.8)
  if(save){
    png(paste0(figure_path, name, '_gnomad_v2_contam_dp20_ref_af_source_correlation.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}

figS16_v3_remove_age_gene_corr_archived <- function(long_data, name, save=TRUE, height = 4, width = 7.5){
  data <- long_data %>%
    filter(!is.na(label)) %>%
    group_by(label, ref_AF, source) %>% 
    dplyr::summarize(tidy(cor.test(contamination, value)))
  figure <- data %>%
    filter(!(ref_AF %in% c('ref-AF: (0%, 100%)', 'None'))) %>%
    mutate(ref_AF = factor(str_replace(ref_AF, 'ref-AF: ', ''), levels = c('(1%, 99%)', '(5%, 95%)', '(10%, 90%)', '(20%, 80%)'))) %>%
    ggplot + 
    # labs(title = paste0('CHARR (20 < DP < 100) ~ Freemix score'), x = 'Reference allele frequency', y = "Pearson's correlation", color = 'Sequencing sites') +
    labs(title = paste0('CHARR ~ Freemix score'), x = 'Reference allele frequency', y = "Pearson's correlation", color = 'Sequencing sites') +
    geom_pointrange(aes(x = ref_AF, y = estimate, ymin = conf.low, ymax = conf.high, group = label, color = label)) + 
    geom_line(aes(x = ref_AF, y = estimate, group = label, color = label)) +
    scale_color_manual(values = c('#045494', 'gray') ) +
    themes +  theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, family = 'Arial', size = 15, face = 'bold'),  
          axis.title = element_text(family = 'Arial', size = 13, face = 'bold'),
          axis.text = element_text(family = 'Arial', size = 12),
          axis.text.x = element_text(family = 'Arial', size = 10),
          legend.text = element_text(family = 'Arial', size = 12),
          legend.title = element_text(family = 'Arial', size = 12, face = 'bold'),
          legend.position = 'top',
          strip.text.x = element_text(family = 'Arial', size = 10, face = 'bold')) + 
    facet_grid(~source, scale = 'free') +
    annotate(geom='text', x=Inf, y=-Inf,label=expression(p.values<10^-100), family = 'Arial', size=4, hjust = 1, vjust=-0.8)
  if(save){
    png(paste0(figure_path, name, '_gnomad_v3_contam_dp20_age_gene_correlation.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}


figS14_two_way_mixing_archived <- function(name, save=TRUE, height = 4, width = 7.5){
  mixed_data <- read.csv(paste0(data_path, 'charr_simulation_hgdp_2_way_mixed_samples_full_150_samples.csv'), sep='\t')
  mixed_long <- mixed_data %>% 
    select(1:3, estimator=freemix_score, 7:8) %>% 
    rbind(mixed_data %>% select(1:3, estimator=charr, 7:8)) %>% 
    mutate(type=factor(rep(c('VBID', 'CHARR'), each=150), level=c('VBID', 'CHARR'))) 
  figure <- ggplot(mixed_long, aes(x = contam_rate, y = estimator, color=type, group=interaction(type, contam_rate))) + 
    geom_boxplot() +
    labs(x='True Contamination Rate', y='Contamination Estimator', color='Type') + 
    # scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'), labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'))) + 
    scale_color_brewer(palette = 'Dark2') +
    geom_abline(slope=1, intercept=0, lty=2)
  if(save){
    png(paste0(figure_path, name, '_hgdp_2_way_mixed_150_samples_freemix_score_vs_charr.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off() 
  }
  return(figure)
}

figS14_two_way_mixing_archived <- function(name, save=TRUE, height = 4, width = 7.5){
  mixed_data <- read.csv(paste0(data_path, 'charr_simulation_hgdp_2_way_mixed_samples_full_150_samples.csv'), sep='\t')
  mixed_long <- mixed_data %>% 
    select(1:3, estimator=freemix_score, 7:8) %>% 
    rbind(mixed_data %>% select(1:3, estimator=charr, 7:8)) %>% 
    mutate(type=factor(rep(c('VBID', 'CHARR'), each=150), levels = c('VBID', 'CHARR'))) 
  figure <- ggplot(mixed_long, aes(x = type, y = estimator, color=type,group= type)) + geom_boxplot() +
    labs(x=NULL, y='Contamination level', color='Type') + 
    # scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'), labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'))) + 
    scale_color_brewer(palette = 'Dark2') +
    facet_grid(~contam_rate) + 
    geom_hline(data = mixed_long, aes(yintercept = contam_rate), lty=2)
  if(save){
    png(paste0(figure_path, name, '_hgdp_2_way_mixed_150_samples_freemix_score_vs_charr_color_by_contaminant.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off() 
  }
  return(figure)
}


figS17_n_way_mixing_archived <- function(name, save=TRUE, height = 4, width = 7.5, type ='gnomAD'){
  mixed_data_n <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples.csv'), sep='\t')
  if(type != 'gnomAD'){
    mixed_data_n <- mixed_data_n %>%
      mutate(charr = charr_local)
  }
  mixed_long <- mixed_data_n %>% 
    select(1:2, estimator=freemix_score, 6) %>% 
    rbind(mixed_data_n %>% select(1:2, estimator=charr, 6)) %>% 
    mutate(type=factor(rep(c('VBID', 'CHARR'), each=150))) 
  figure <- ggplot(mixed_long, aes(x = type, y = estimator, color=pop_original,group=interaction(pop_original, type))) + geom_boxplot() +
    labs(x=NULL, y='Contamination level', color='Population (original)') + 
    scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'), labels= toupper(c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas','oth'))) + 
    facet_grid(~contam_rate) + 
    geom_hline(data = mixed_long, aes(yintercept = contam_rate), lty=2)
  if(save){
    png(paste0(figure_path, name, '_hgdp_n_way_mixed_150_samples_freemix_score_vs_charr_',type,'_AF.png'), height = height, width = width, units = 'in', res = 300)
    print(figure)
    dev.off() 
  }
  return(figure)
}
# figS17_n_way_mixing('S17', save=save, height = 4, width = 7.5, type = 'gnomAD')
# figS17_n_way_mixing('S17', save=save, height = 4, width = 7.5, type = 'sample')

figS18_2_way_color_by_truth_archived <- function(name, save=TRUE, height = 5, width = 7.5, type='gnomAD'){
  mixed_data_2 <- read.csv(paste0(data_path, 'charr_simulation_hgdp_n_way_mixed_samples_full_150_samples.csv'), sep='\t')
  if(type != 'gnomAD'){
    mixed_data_2 <- mixed_data_2 %>%
      mutate(charr = charr_local)
  }
  p <- mixed_data_2 %>% 
    ggplot + aes(x=freemix_score, 
                 y=charr, 
                 color=factor(contam_rate, levels=c(0.005, 0.01, 0.02, 0.05, 0.1), labels = c('0.5%', '1.0%', '2.0%', '5.0%', '10.0%'))) +
    geom_point(size=2 )+
    geom_abline(slope = 1, intercept = 0, lty =2) + 
    scale_color_manual(values = c('#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695')) +
    scale_y_continuous(label = scales::percent_format(accuracy=0.1)) +
    scale_x_continuous(label = scales::percent_format(accuracy=0.1)) +
    labs(x='Freemix Score', y='CHARR', color='True Contamination Rate') 
  if(save){
    png(paste0(figure_path, name, '_hgdp_2_way_mixed_150_samples_', type ,'_AF.png'), height = height, width = width, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  return(p)
}

# figure1 = figS18_2_way_color_by_truth('S18', save=save, height = 4, width = 6, type = 'gnomAD')
# figure2 = figS18_2_way_color_by_truth('S18', save=save, height = 4, width = 6, type = 'sample')
# figure = ggpubr::ggarrange(figure1, figure2, ncol = 2, labels = c('(A): gnomAD AF', '(B): local AF'), common.legend = TRUE, vjust = 0.5, hjust=0,
#                            #vjust = 0, hjust = 0, 
#                            font.label = list(size = 13, color = "black", face = "bold", family = 'Arial'))
# figure
# png(paste0(figure_path, 'S18_hgdp_n_way_mixed_150_samples.png'), height = 4, width = 8, units = 'in', res = 300)
# print(figure)
# dev.off()