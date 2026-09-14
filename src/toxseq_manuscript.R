# STARTUP #### 

overall_start_time = Sys.time()
tell_user = function(...) { cat(file=stderr(), paste0(...)); flush.console() }

# DEPENDENCIES ####

tell_user('Loading required packages...')

options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/toxseq')
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(janitor)
  library(openxlsx)
  library(magick)
  library(smoother)
  library(fgsea)
  library(msigdbr)
  select = dplyr::select
  summarize = dplyr::summarize
})


# OUTPUT STREAMS #### 

tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T
write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/data-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}




# FUNCTIONS ####

tell_user('done.\nDefining functions...')


plot_sf = function(sf, anno_dpi=NULL, xlims=c(0,200), ylims=c(0,1.05),
                   xbreaks=seq(0,500,100), xminor=seq(0,500,10),
                   xlabel='days post-inoculation', legend_pos='bottomleft',
                   anno_cex=0.15) {
  plot(sf, col=sf$color, lwd=1.5, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbreaks, labels=NA, tck=-0.05)
  axis(side=1, at=xbreaks, lwd=0, labels=xbreaks, line=-0.5)
  if (!is.null(xminor)) axis(side=1, at=xminor, labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text=xlabel, cex=0.8)
  axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2)
  mtext(side=2, line=2.75, text='survival', cex=0.8)
  par(xpd=T)
  if (!is.null(anno_dpi)) {
    points(x=c(anno_dpi), y=rep(max(ylims),length(anno_dpi)), pch=25, bg='#000000', cex=anno_cex)
  }
  # if (legend_pos == 'right') {
  #   par(xpd=NA)
  #   legend(x=xlims[2] + diff(xlims)*0.03, y=ylims[2], bty='n', cex=.8, lwd=2,
  #          legend=surv_meta$disp2, col=surv_meta$color, text.col=surv_meta$color)
  # } else {
  #   legend(x=legend_pos, y=1, bty='n', cex=.8, lwd=2, legend=surv_meta$disp2,
  #          col=surv_meta$color, text.col=surv_meta$color)
  # }
  par(xpd=F)
}



clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}

rank_uniq = function(x) match(x, sort(unique(x)))

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}

lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot



format_p = function(p) {
  formatted_p = paste0(' = ',formatC(signif(p,2), format='g', digits=2))
  formatted_p[signif(p,2) < 0.1] = paste0(' = ',formatC(signif(p[signif(p,2) < 0.1],2), format='f', digits=3))
  formatted_p[signif(p,2) < 0.01] = paste0(' = ',formatC(signif(p[signif(p,2) < 0.1],2), format='e', digits=1))
  formatted_p[signif(p,2) == 1] = ' = 1.0'
  formatted_p[signif(p,2) < 1e-15] = ' < 1.0e-15'
  return (formatted_p)
}


meansd = function(x, digits=1) {
  paste0(formatC(mean(x, na.rm=T),format='f',digits=digits),'±',formatC(sd(x, na.rm=T),format='f',digits=digits))
}

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x > 0 & signed, '+', ''),formatC(100*x,format='f',digits=digits),'%'))


clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

# one fewer parameter - do clipdist symmetrically around 0
clipsym = function(x, absx) {
  return (pmin(absx, pmax(-absx, x)))
}


# Hypergeometric over-representation of a gene list in the toxic
# transcriptional signature (g_tox_genes). Accepts a character vector,
# a tibble/data.frame whose first column is gene symbols, or the path
# to a one-column text file. direction='up'/'down' restricts the
# signature to up- or down-regulated tox genes; 'any' uses both.
test_set_enrichment = function(genes, name = '',
                               direction = c('any','up','down'),
                               background = NULL) {
  direction = match.arg(direction)
  if (is.character(genes) && length(genes) == 1 && file.exists(genes)) {
    genes = readLines(genes) %>% trimws()
    genes = genes[nzchar(genes)]
  } else if (is.data.frame(genes)) {
    genes = genes %>% pull(1) %>% as.character()
  }
  genes = unique(genes[!is.na(genes) & nzchar(genes)])

  signature = switch(direction,
                     up   = g_tox_genes$gene[g_tox_genes$sign > 0],
                     down = g_tox_genes$gene[g_tox_genes$sign < 0],
                     any  = g_tox_genes$gene)
  if (is.null(background)) background = volcano_data$gene
  background = unique(background)

  genes_in_universe     = intersect(genes, background)
  signature_in_universe = intersect(signature, background)
  overlap               = intersect(genes_in_universe, signature_in_universe)

  N = length(background)
  K = length(signature_in_universe)
  n = length(genes_in_universe)
  k = length(overlap)
  expected = K * n / N
  log2_fold_enrich = log2((k + 0.5) / (expected + 0.5))
  pval = phyper(k - 1, K, N - K, n, lower.tail = FALSE)

  tibble(set_name = name, direction = direction,
         n_input = length(genes), n_input_in_universe = n,
         n_signature = K, n_overlap = k,
         expected = expected, log2_fold_enrich = log2_fold_enrich, pval = pval)
}


volcano = function(tbl,
                   clipwidth=5, 
                   maxy=10, 
                   title='', 
                   colorcolname='color',
                   pvalcolname='pval',
                   log2fccolname='log2fc') {

  if (colorcolname %in% colnames(tbl)) {
    colors = tbl[colorcolname] %>% pull()
  } else {
    colors = '#000000'
  }
  xlims = c(-clipwidth, clipwidth)
  ylims = c(0, maxy*1.05)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  points(x=clipdist(tbl[[log2fccolname]], -clipwidth, clipwidth), y=pmin(maxy,-log10(tbl[pvalcolname] %>% pull())), pch=20,
       col=colors)
  axis(side=1, at=-clipwidth:clipwidth, tck=-0.025, labels=NA)
  axis(side=1, at=-clipwidth:clipwidth, lwd=0, line=-0.75, cex.axis=0.8)
  abline(v=0, lwd=.125)
  abline(h=-log10(0.05), lwd=.125)
  axis(side=2, tck=-0.025, labels=NA)
  axis(side=2, line=-0.5, lwd=0, cex.axis=0.8, las=2)
  mtext(side=1, line=1.6, text='L2FC', cex=0.7)
  mtext(side=2, line=2.0, text='-log10(P)', cex=0.7)
  mtext(side=3, line=0, text=title, cex=0.7)
  
}



# DATA ####

tell_user('done.\nReading in data...')

# metadata
treatments        = read_tsv('analytic/treatments.tsv', col_types=cols())
groups            = read_tsv('analytic/groups.tsv', col_types=cols())

# DESeq2 differential expression
contrasts         = read_tsv('analytic/deseq2_contrasts.tsv.gz', col_types=cols())
dge_score_genes   = read_tsv('analytic/legacy_dge_score_genes.tsv', col_types=cols())

# cell-type gene lists
astrocyte_genes      = read_csv('analytic/astrocyte_genes.csv', col_types=cols()) %>% clean_names()
microglia_genes      = read_csv('analytic/microglia_genes.csv', col_types=cols()) %>% clean_names()
oligodendrocyte_genes = read_csv('analytic/oligodendrocyte_genes.csv', col_types=cols()) %>% clean_names()

# survival cohorts
l_cohort          = read_tsv('analytic/L_cohort.tsv', col_types=cols())
dp_cohort         = read_tsv('analytic/DP_cohort.tsv', col_types=cols())
prp190327_survival = read_tsv('analytic/PRP190327_survival.tsv', col_types=cols())
aso2_doseresp_survival = read_tsv('analytic/aso2_doseresp_survival.tsv', col_types=cols())
toxmit_survival   = read_tsv('analytic/toxmit_survival.tsv', col_types=cols())

# NfL / biomarkers
prp190327_nfl     = read_tsv('analytic/PRP190327_nfl.tsv', col_types=cols())
nfl_timepoints    = read_tsv('analytic/nfl_timepoints.tsv', col_types=cols())

# MRI
mri_timeseries_tox = read_tsv('analytic/mri_timeseries_tox.tsv', col_types=cols())
mri_timeseries_histo = read_tsv('analytic/time_series_histo_scores.tsv', col_types=cols())

# qPCR / pharmacology
mouse_500ug_1wk_qpcr_pk    = read_tsv('analytic/mouse_500ug_1wk_qpcr_pk.tsv', col_types=cols())
mouse_500ug_1wk_fob_qpcr   = read_tsv('analytic/mouse_500ug_1wk_fob_qpcr.tsv', col_types=cols())
mouse_500ug_3xicv_10wk     = read_tsv('analytic/mouse_500ug_3xicv_10wk.tsv', col_types=cols())
mouse_axd_4wk            = read_tsv('analytic/mouse_axd_4wk.tsv', col_types=cols())
rat_3mg_1wk_qpcr           = read_tsv('analytic/rat_3mg_1wk_qpcr.tsv', col_types=cols())
rat_1mg_3xicv_4wk          = read_tsv('analytic/rat_1mg_3xicv_4wk.tsv', col_types=cols())
rat_3mg_8wk                = read_tsv('analytic/rat_3mg_8wk.tsv', col_types=cols())
dp_qpcr                    = read_tsv('analytic/DP_qpcr.tsv', col_types=cols())
prp230201_qpcr             = read_tsv('analytic/PRP230201_dose_response_qpcr.tsv', col_types=cols())
prp250210_qpcr             = read_tsv('analytic/PRP250210_qpcr.tsv', col_types=cols())

# potency / dose-response
prp240918_potency          = read_tsv('analytic/PRP240918_modified_aso2_potency.tsv', col_types=cols())
bjab_ccl2                  = read_tsv('analytic/bjab_ccl2_log10_fold_changes.tsv', col_types=cols())
itc_all_asos               = read_tsv('analytic/itc_all_asos.tsv', col_types=cols())
s2_readout_labels          = read_tsv('analytic/s2_readout_labels.tsv', col_types=cols())

# DGE scoring (from dan_contrib)
umap_data                  = read_tsv('analytic/umap_data.tsv',       col_types=cols())
volcano_data               = read_tsv('analytic/volcano_data.tsv',    col_types=cols()) %>% clean_names()
dge_scores_24              = read_tsv('analytic/dge_scores_24.tsv',   col_types=cols())
contrasts_0432             = read_tsv('analytic/all_contrasts_0432.tsv.gz', col_types=cols())

# The public data files identify ASOs by the anonymized 'aso' column (e.g. ASO1, ASOG);
# restore the internal column names this script expects.
for (v in c('treatments','groups','prp190327_survival','toxmit_survival','prp190327_nfl',
            'mri_timeseries_tox','mouse_500ug_3xicv_10wk','mouse_axd_4wk','rat_1mg_3xicv_4wk',
            'prp230201_qpcr','prp250210_qpcr','prp240918_potency','umap_data','dge_scores_24')) {
  assign(v, get(v) %>% rename(ionis_id = aso))
}
for (v in c('mouse_500ug_1wk_fob_qpcr','rat_3mg_8wk')) {
  assign(v, get(v) %>% rename(isis_no = aso))
}

volcano_data %>%
  filter(padj < 0.05 & base_mean > 10 & abs(log2fold_change) > 1) %>%
  mutate(sign = log2fold_change / abs(log2fold_change)) %>%
  select(gene, sign, log2fold_change, padj, base_mean) -> g_tox_genes


# reference data
lens                       = read_tsv('analytic/gene_lengths.tsv', col_types=cols()) # gene lengths from mm10 GTF file

# Aguzzi paper supplemental data
sn20_supp2                 = read_tsv('analytic/sorce-nuvolone-2020-table-s2.tsv',   col_types=cols())

# Gene lists
homology                   = read_tsv('analytic/jax_mouse_human_orthologs.tsv', col_types=cols()) %>% clean_names()
fazal2019                  = read_tsv('analytic/fazal-2019-table-s3.tsv', col_types=cols()) %>% clean_names()
engel2022                  = read_tsv('analytic/engel-2022-table-s3-halo-fibrillarin.tsv', col_types=cols()) %>% clean_names()
nakaya2018                 = read_tsv('analytic/nakaya-maragkakis-2018-table-s6.tsv', col_types=cols()) %>% clean_names()

## Processing cell type enrichemnt datasets ####
astro_color = '#F0A804' # astrocyte
micro_color = '#CD00CD' # microglia
odc_color = '#37FDFC' # oligodendrocytes

astrocyte_genes %>%
  mutate(celltype = '3_astrocyte', color=astro_color) %>%
  select(-feature_id) %>%
  rename(gene=feature_name, meanexp = astrocyte_average, l2fc = astrocyte_log2_fold_change, pval = astrocyte_p_value) %>%
  filter(pval < 0.01, meanexp >= 1, l2fc >= 2) -> astros

microglia_genes %>%
  mutate(celltype = '2_microglia', color=micro_color) %>%
  select(-feature_id) %>%
  rename(gene=feature_name, meanexp = microglia_average, l2fc = microglia_log2_fold_change, pval = microglia_p_value) %>%
  filter(pval < 0.01, meanexp >= 1, l2fc >= 2) -> micros

oligodendrocyte_genes %>%
  mutate(celltype = '1_oligodendrocyte', color=odc_color) %>%
  select(-feature_id) %>%
  rename(gene=feature_name, meanexp = oligodendrocyte_average, l2fc = oligodendrocyte_log2_fold_change, pval = oligodendrocyte_p_value) %>%
  filter(pval < 0.01, meanexp >= 1, l2fc >= 2) -> odcs

rbind(astros, micros, odcs) %>%
  arrange(desc(l2fc)) %>% group_by(gene) %>% slice(1) %>% # when a gene is in >1 list, pick the one with highest fold enrichment
  select(gene, meanexp, l2fc, pval, celltype, color) -> cell_type_enrichments

write_supp_table(cell_type_enrichments, 'Genes designated as being enriched in astrocytes, microglia, or oligodendrocytes.')

cell_type_enrichments %>%
  group_by(celltype) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() -> cell_type_enrichments_smry

write_supp_table(cell_type_enrichments_smry, 'Count of genes designated as being enriched in astrocytes, microglia, or oligodendrocytes.')

## Processing gene lists ####
homology %>%
  select(db_class_key, common_organism_name, symbol) %>%
  mutate(common_organism_name = gsub('mouse, laboratory','mouse', common_organism_name)) %>%
  pivot_wider(id_cols=db_class_key, names_from=common_organism_name, values_from=symbol, values_fn=list) %>%
  unnest(mouse) %>%
  unnest(human) -> orthologs

fazal2019 %>%
  filter(common_gene != 'NA') %>%
  filter(source=='Atlas_Analysis') %>%
  filter(nucleolus_log2fc > 0.75) %>% # this is the cutoff in the Fazal 2019 paper. they say it's 324 "orphan" (i.e. previously unassigned to a compartment). if you add this line you can replicate the 324 count: # filter(nucleolus_orphan==1) %>%
  distinct(common_gene) -> fazal_nucleolus_genes_human

fazal_nucleolus_genes_human %>%
  inner_join(orthologs, by=c('common_gene'='human')) %>%
  select(mouse) %>%
  rename(gene=mouse) -> fazal_nucleolus_genes_mouse

engel2022 %>%
  filter(padj < 0.05, log2fc >= 0.5) %>%
  distinct(gene) -> engel_nucleolus_genes_human 

engel_nucleolus_genes_human %>%
  inner_join(orthologs, by=c('gene'='human')) %>%
  select(mouse) %>%
  rename(gene=mouse) -> engel_nucleolus_genes_mouse

nakaya2018 %>%
  select(gene) -> nakaya_fus_genes


# DISPLAY ITEMS ####


## Table 1. ASOs used in this study #### 
treatments %>%
  filter(include=='yes' & !is.na(sequence)) %>%
  select(display_name, target, sequence, citation) -> table_1


## Table S1. Which data were collected for which ASO ####
expt_ids = list(
  'L cohort survival'             = l_cohort              %>% filter(include) %>% select(ionis_id = treatment) %>% distinct(),
  'PRP190327 survival'            = prp190327_survival    %>% distinct(ionis_id),
  'ASO 2 dose-response survival'  = tibble(ionis_id = 'ASO2'),
  'ToxMit survival'               = toxmit_survival       %>% distinct(ionis_id),
  'PRP190327 NfL'                 = prp190327_nfl         %>% distinct(ionis_id),
  'MRI'                           = mri_timeseries_tox    %>% distinct(ionis_id),
  'Mouse 500ug 1wk PK qPCR'      = mouse_500ug_1wk_qpcr_pk  %>% select(ionis_id = treatment) %>% distinct(),
  'Mouse 500ug 1wk FOB qPCR'     = mouse_500ug_1wk_fob_qpcr %>% mutate(ionis_id = as.character(isis_no)) %>% distinct(ionis_id),
  'Mouse 500ug 3xICV 10wk'       = mouse_500ug_3xicv_10wk   %>% distinct(ionis_id),
  'Mouse 4xICV 4wk'              = mouse_axd_4wk          %>% distinct(ionis_id),
  'Rat 3mg 1wk qPCR'             = rat_3mg_1wk_qpcr         %>% select(ionis_id = treatment) %>% distinct(),
  'Rat 1mg 3xICV 4wk'            = rat_1mg_3xicv_4wk        %>% distinct(ionis_id),
  'Rat 3mg 8wk'                  = rat_3mg_8wk              %>% mutate(ionis_id = as.character(isis_no)) %>% distinct(ionis_id),
  'PRP230201 dose-response qPCR' = prp230201_qpcr           %>% distinct(ionis_id),
  'PRP250210 qPCR'               = prp250210_qpcr           %>% distinct(ionis_id),
  'ASO 2 modified potency'       = prp240918_potency        %>% distinct(ionis_id),
  'BJAB CCL22'                   = bjab_ccl2               %>% select(ionis_id = treatment) %>% distinct(),
  'RNA-seq'                      = groups                  %>% distinct(ionis_id)
)

map_dfr(names(expt_ids), ~expt_ids[[.x]] %>% mutate(ionis_id = as.character(ionis_id), experiment = .x)) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = experiment, values_from = present, values_fill = FALSE) %>%
  right_join(treatments %>% filter(include == 'yes') %>% select(ionis_id, display_name),
             by = 'ionis_id') %>%
  select(display_name, all_of(names(expt_ids))) %>%
  arrange(display_name) -> table_s1

write_supp_table(table_s1, 'Figure S2: Data types collected per ASO')

clipcopy(table_s1)
# clipcopy(table_1)

## Table 2 ####

prp190327_survival %>%
  mutate(event=1) %>%
  inner_join(treatments %>% select(ionis_id, display_name), by='ionis_id') %>%
  select(display_name, dpi, event) %>%
  mutate(display_name = factor(display_name, levels=rev(sort(unique(display_name))))) -> prp190327_survival_proc
cox_model = coxph(Surv(dpi, event) ~ display_name, data=prp190327_survival_proc)
summary(cox_model)$coefficients %>%
  as_tibble(rownames='var') %>%
  mutate(display_name = gsub('display_name','',var)) %>%
  select(-var) %>%
  relocate(display_name) %>%
  rename(exp_coef = `exp(coef)`, se_coef = `se(coef)`, pval = `Pr(>|z|)`) -> cox_ps

pbs_mean = mean(prp190327_survival$dpi[prp190327_survival$ionis_id=='PBS'])
prp190327_survival_proc %>%
  group_by(display_name) %>%
  summarize(.groups='keep',
            n=n(),
            survival_dpi_median = median(dpi),
            survival_dpi_value = mean(dpi),
            survival_dpi_meansd = meansd(dpi),
            survival_dpi_delta_value = mean(dpi/pbs_mean)-1,
            survival_dpi_delta_disp = percent(mean(dpi/pbs_mean)-1,signed=T),
            survival_dpt_value = mean(dpi-120),
            survival_dpt_meansd = meansd(dpi-120)) %>%
  ungroup() %>%
  left_join(cox_ps, by='display_name') %>%
  mutate(pval = case_when(display_name=='PBS' ~ 1, TRUE ~ pval)) %>%
  arrange(desc(survival_dpt_value)) %>%
  mutate(label = case_when(survival_dpi_delta_value < -0.03 ~ 'toxic',
                           display_name=='PBS' ~ 'PBS',
                           TRUE ~ 'tolerated')) -> table2_data

table2_data %>%
  mutate(survival_dpt_meansd = case_when(survival_dpi_delta_value > -0.03 ~ '',
                                         TRUE ~ survival_dpt_meansd)) %>%
  select(display_name, n, survival_dpi_meansd, survival_dpi_delta_disp, survival_dpt_meansd, pval) -> table_2

write_supp_table(table2_data, 'Figure 1C: Summary of survival statistics for ASOs A-H.')

toxtol_colors = tribble(
  ~label, ~color,
  'toxic', '#FF6347',
  'tolerated', '#6495ED',
  'PBS', '#BBBBBB'
)
tol_col = toxtol_colors$color[toxtol_colors$label=='tolerated']
tox_col = toxtol_colors$color[toxtol_colors$label=='toxic']

mouse_500ug_1wk_qpcr_pk %>%
  inner_join(treatments %>% select(ionis_id,display_name), by=c('treatment'='ionis_id')) %>%
  group_by(display_name) %>%
  filter(analyte != 'PK') %>%
  mutate(value = value / 100) %>%
  ungroup() %>%
  group_by(display_name, analyte) %>%
  summarize(.groups='keep',
            mean_value = mean(value),
            mean_percent = percent(mean(value)),
            sd_percent = percent(sd(value)),
            meansd_percent = paste0(mean_percent,'±', sd_percent), 
            bold_value = mean_value > 1.5) %>%
  ungroup() -> mouse_500ug_1wk_qpcr_long

mouse_500ug_1wk_qpcr_long %>%
  pivot_wider(id_cols = display_name, names_from = analyte, values_from=c(meansd_percent, bold_value)) -> mouse_500ug_1wk_qpcr_wide
  
bjab_ccl2 %>%
  inner_join(treatments %>% select(ionis_id,display_name), by=c('treatment'='ionis_id')) %>%
  rename(bjab_value = value) %>%
  mutate(bold_value_bjab = bjab_value > 2) -> bjab_ccl2_wide
  



## Figure 1.  #### 
tell_user('done.\nCreating Figure 1...')
if ('figure-1'=='figure-1') {
  
  resx=600
  png('display_items/figure-1.png',width=6.5*resx,height=6.5*resx,res=resx)

  # Row 1: A (dose-resp survival), B (MRI timeseries)
  # Row 2: C (prp190327 survival), D (CCL2), E (ITC Kd)
  # Row 3: F (NfL), G (PK)
  layout_matrix = matrix(c(1,1,2,2,
                           3,3,4,5,
                           6,6,7,7), nrow=3, byrow=T)
  layout(layout_matrix)
  panel = 1

  condition_colors = c(toxic='tomato', tolerated='cornflowerblue', PBS='gray60')

  ### A. dose-responsiveness ####
  aso2_doseresp_survival %>%
    mutate(event=1) -> aso2_doseresp_survival_proc

  aso2_dr_meta = tibble(
    dose = sort(unique(aso2_doseresp_survival$dose)),
    color = colorRampPalette(colors=c('#D7D7D7','#DDDD01','#990199'))(6)
  )

  par(mar=c(3,4,3,1))
  xlims = c(0, 250)
  xbigs = seq(0,500,100)
  xats = seq(0,500,10)
  sf = survfit(Surv(dpi, event) ~ dose, data=aso2_doseresp_survival_proc)
  sf$dose = gsub('dose=','',names(sf$strata))
  sf$color = aso2_dr_meta$color[match(sf$dose, aso2_dr_meta$dose)]
  plot(sf, xlim = xlims, col=sf$color, lwd=2, ylim=c(0,1.05), xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, labels=xbigs, line=-0.5)
  if (!is.null(xats)) axis(side=1, at=xats, labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation', cex=0.8)
  axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2)
  mtext(side=2, line=2.75, text='survival', cex=0.8)
  par(xpd=T)
  text(x=0, y=1.05, labels='prion\ninoculation', pos=3, cex=0.7)
  points(0, 1.03, pch=25, bg='black', col='black', cex=0.8)
  points(x=90, y=1.05, pch=25, col='#000000', bg='#000000')
  text(x=90, y=1.05, labels='ASO', pos=3, cex=0.7)
  legend(x=1, y=0.99, aso2_dr_meta$dose, col=aso2_dr_meta$color, lwd=2, bty='n', cex=0.6,
         title = 'dose (µg)', title.col='#000000')
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1, adj = -0.15, line = 0.5); panel = panel + 1

  #### Fig 1A summary table ####
  aso2_doseresp_survival_proc %>%
    group_by(dose) %>%
    summarise(
      n        = n(),
      median   = median(dpi),
      mean     = mean(dpi),
      sd       = sd(dpi),
      .groups  = 'drop'
    ) %>%
    mutate(prop_increase_vs_0 = median / median[dose == 0] - 1) %>%
    mutate(pval_vs_0 = sapply(dose, function(d) {
      if (d == 0) return(NA_real_)
      survdiff(Surv(dpi, event) ~ dose,
               data = aso2_doseresp_survival_proc %>% filter(dose %in% c(0, d)))$chisq %>%
        { pchisq(., df = 1, lower.tail = FALSE) }
    })) %>%
    rename(dose_ug                  = dose,
           median_survival_dpi      = median,
           mean_survival_dpi        = mean,
           sd_survival_dpi          = sd,
           median_prop_increase_vs_0 = prop_increase_vs_0,
           pval_vs_0_logrank        = pval_vs_0) -> fig1a_surv_table

  write_supp_table(fig1a_surv_table, 'Figure 1A: ASO dose-response survival summary')
  write_supp_table(aso2_doseresp_survival %>% select(dose_ug = dose, survival_dpi = dpi),
                   'Figure 1A: ASO 2 dose-response per-animal survival')



  ### B. MRI timeseries ####
  mri_timeseries_tox %>%
    pivot_longer(cols = !ionis_id, names_to = 'timepoint_label', values_to = 'value') %>%
    mutate(dpi = as.integer(gsub(' dpi', '', timepoint_label))) %>%
    filter(value != '\u2014') %>%
    group_by(ionis_id) %>%
    summarise(first_toxic_dpi = min(dpi), .groups = 'drop') -> mri_first_toxic_dpi

  mri_timeseries_tox %>%
    mutate(label = case_when(
      ionis_id == 'PBS'        ~ 'PBS',
      `123 dpi` != '\u2014'   ~ 'toxic',
      TRUE                     ~ 'tolerated'
    )) %>%
    left_join(treatments %>% select(ionis_id, display_name), by='ionis_id') %>%
    left_join(mri_first_toxic_dpi, by = 'ionis_id') %>%
    arrange(factor(label, levels=c('PBS','tolerated','toxic')),
            first_toxic_dpi,
            display_name) -> mri_meta_raw

  n_pbs_mri = sum(mri_meta_raw$label == 'PBS')
  n_tol_mri = sum(mri_meta_raw$label == 'tolerated')
  n_tox_mri = sum(mri_meta_raw$label == 'toxic')

  mri_meta_raw %>%
    mutate(y = c(
      seq_len(n_pbs_mri),
      seq_len(n_tol_mri) + n_pbs_mri + 1L,
      seq_len(n_tox_mri) + n_pbs_mri + n_tol_mri + 2L
    ),
    y = max(y) + 1L - y) -> mri_meta

  mri_timeseries_tox %>%
    pivot_longer(cols=!ionis_id, names_to='timepoint_label', values_to='value') %>%
    mutate(
      dpi        = as.integer(gsub(' dpi', '', timepoint_label)),
      x          = match(dpi, c(-14L, 30L, 60L, 90L, 123L)),
      has_number = (value != '\u2014'),
      cell_color = if_else(has_number,
                           toxtol_colors$color[toxtol_colors$label=='toxic'],
                           toxtol_colors$color[toxtol_colors$label=='tolerated'])
    ) %>%
    left_join(mri_meta %>% select(ionis_id, y), by='ionis_id') -> mri_long

  ci_alpha  = 0.5
  xdpis     = c(-14, 30, 60, 90, 123)
  half_w    = 0.42
  half_h    = 0.42
  xlims_mri = c(0.5, 5.5)
  ylims_mri = c(0.5, max(mri_meta$y) + 0.5)

  par(mar=c(3, 5, 3, 1))
  plot(NA, NA, xlim=xlims_mri, ylim=ylims_mri, axes=F, ann=F, xaxs='i', yaxs='i')
  mtext(side=1, at=1:5, line=0.4, text=xdpis, cex=0.7)
  mtext(side=1, line=1.6, text='treatment dpi', cex=0.7)
  mtext(side=2, at=mri_meta$y,
       text=mri_meta$display_name, las=2, cex=0.5, line=1.8)
  
  rect(xleft   = mri_long$x - half_w,
       xright  = mri_long$x + half_w,
       ybottom = mri_long$y - half_h,
       ytop    = mri_long$y + half_h,
       col     = alpha(mri_long$cell_color, ci_alpha),
       border  = NA)
  mri_num = mri_long %>% filter(has_number)
  text(x=mri_num$x, y=mri_num$y, labels=mri_num$value, cex=0.8, col='black')
  mri_meta %>%
    left_join(treatments %>% select(ionis_id, targets_prnp), by='ionis_id') %>%
    mutate(targets_prnp = coalesce(targets_prnp, 'no')) -> mri_meta_prnp
  mtext(side=2, at=mri_meta_prnp$y, line=0.3, las=2,
        text=mri_meta_prnp$targets_prnp, cex=0.55, col='#666666')
  mtext(side=2, at=max(mri_meta_prnp$y+1), line=0.3, las=2,
        text='targets\nPrnp', cex=0.45, col='#666666')
  y_tol_range = range(mri_meta$y[mri_meta$label %in% c('PBS', 'tolerated')])
  y_tox_range = range(mri_meta$y[mri_meta$label == 'toxic'])
  overhang = 0.4
  axis(side=2, at=y_tol_range + overhang*c(-1,1), labels=NA, lwd=1, lwd.ticks=1, tck=0.03, line=4)
  mtext(side=2, at=mean(y_tol_range), text='tolerated', col=tol_col, line=4.2, cex=0.65)
  axis(side=2, at=y_tox_range + overhang*c(-1,1), labels=NA, lwd=1, lwd.ticks=1, tck=0.03, line=4)
  mtext(side=2, at=mean(y_tox_range), text='toxic', col=tox_col, line=4.2,  cex=0.65)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(mri_meta_raw %>%
                     select(display_name, label, `-14 dpi`, `30 dpi`, `60 dpi`, `90 dpi`, `123 dpi`, first_toxic_dpi),
                   'Figure 1B: MRI hyperintensity scores per ASO and timepoint')

  ### C. prp190327 survival ####
  treatments %>%
    select(ionis_id, display_name) %>%
    inner_join(table2_data %>% select(display_name, label, survival_dpi_delta_value),
               by='display_name') %>%
    arrange(desc(survival_dpi_delta_value)) %>%
    group_by(label) %>%
    mutate(group_rank = row_number()) %>%
    ungroup() -> prp_surv_sorted

  n_tox_surv = sum(prp_surv_sorted$label == 'toxic')
  n_tol_surv = sum(prp_surv_sorted$label == 'tolerated')
  tox_pal = colorRampPalette(c('#FFB1A3', '#FF6347'))(max(n_tox_surv, 1))
  tol_pal = colorRampPalette(c('#B2CAF6', '#6495ED'))(max(n_tol_surv, 1))

  prp_surv_sorted %>%
    mutate(color = case_when(
      label == 'PBS'       ~ '#000000',
      label == 'toxic'     ~ tox_pal[group_rank],
      label == 'tolerated' ~ tol_pal[group_rank]
    ),
    color = if_else(ionis_id == 'ASO6', '#41B6C4', color)) -> prp190327_surv_meta

  prp190327_survival %>%
    mutate(event = 1) %>%
    inner_join(prp190327_surv_meta %>% select(ionis_id, color), by='ionis_id') -> prp190327_surv_data

  surv_color_meta = prp190327_surv_data %>%
    distinct(ionis_id, color) %>%
    arrange(ionis_id)

  sf = survfit(Surv(dpi, event) ~ ionis_id, data=prp190327_surv_data)
  sf$color = surv_color_meta$color[match(gsub('ionis_id=', '', names(sf$strata)), surv_color_meta$ionis_id)]

  par(mar=c(3,4,3,1))
  plot_sf(sf, xlims=c(0,260), xbreaks=seq(0,300,50), xminor=seq(0,300,10))
  par(xpd=T)
  text(x=0, y=1.05, labels='prion\ninoculation', pos=3, cex=0.7)
  points(0, 1.03, pch=25, bg='black', col='black', cex=0.8)
  text(x=120, y=1.05, labels='ASO', pos=3, cex=0.7)
  points(120, 1.03, pch=25, bg='black', col='black', cex=0.8)
  par(xpd=F)
  prp190327_surv_meta %>%
    inner_join(table2_data %>% select(display_name, survival_dpi_delta_disp),
               by='display_name') %>%
    arrange(desc(survival_dpi_delta_value)) %>%
    mutate(legend_label = paste0(display_name, ' (', survival_dpi_delta_disp, ')')) -> surv_legend
  surv_legend %>%
    left_join(treatments %>% select(ionis_id, targets_prnp), by='ionis_id') %>%
    mutate(targets_prnp = coalesce(targets_prnp, 'no')) -> surv_legend_prnp
  pu       = par('usr')  # c(x1, x2, y1, y2) in plot units
  px_range = pu[2] - pu[1]
  par(xpd=T)
  leg = legend('bottomleft', inset=c(0.06, 0.02),
               legend=surv_legend_prnp$legend_label, col=surv_legend_prnp$color,
               lwd=1.5, bty='n', cex=0.65)
  dy        = abs(leg$text$y[1] - leg$text$y[2]) * 0.5
  tranche_x = leg$rect$left - px_range * 0.012
  label_x   = tranche_x     - px_range * 0.014
  tol_ys    = leg$text$y[surv_legend_prnp$label == 'tolerated']
  tox_ys    = leg$text$y[surv_legend_prnp$label == 'toxic']
  segments(x0=tranche_x, x1=tranche_x,
           y0=min(tol_ys)-dy, y1=max(tol_ys)+dy, col=tol_col, lwd=2)
  text(label_x, mean(tol_ys), 'tolerated', srt=90, col=tol_col, cex=0.45, adj=c(0.5, 0.5))
  segments(x0=tranche_x, x1=tranche_x,
           y0=min(tox_ys)-dy, y1=max(tox_ys)+dy, col=tox_col, lwd=2)
  text(label_x, mean(tox_ys), 'toxic', srt=90, col=tox_col, cex=0.45, adj=c(0.5, 0.5))
  prnp_x = leg$rect$left + leg$rect$w + px_range * 0.01
  text(prnp_x, leg$text$y, surv_legend_prnp$targets_prnp, adj=0, cex=0.5, col='#666666')
  text(prnp_x, max(leg$text$y)+0.1, 'targets\nPrnp', adj=0, cex=0.5, col='#666666')
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1, adj = -0.15, line = 0.5); panel = panel + 1



  write_supp_table(prp190327_surv_data %>%
                     left_join(treatments %>% select(ionis_id, display_name), by='ionis_id') %>%
                     select(display_name, survival_dpi = dpi, event),
                   'Figure 1C: PRP190327 per-animal survival')

  ### D. CCL2 ####
  # combine survival-based labels (table2_data) with MRI-based labels for Prnp ASOs
  mri_timeseries_tox %>%
    mutate(mri_label = case_when(
      ionis_id == 'PBS'      ~ 'PBS',
      `123 dpi` != '—'  ~ 'toxic',
      TRUE                   ~ 'tolerated'
    )) %>%
    select(ionis_id, mri_label) -> pk_mri_labels

  treatments %>%
    select(ionis_id, display_name) %>%
    left_join(table2_data %>% select(display_name, label), by='display_name') %>%
    left_join(pk_mri_labels, by='ionis_id') %>%
    mutate(label = coalesce(label, mri_label)) %>%
    filter(!is.na(label), label != 'PBS') %>%
    left_join(toxtol_colors, by='label') %>%
    semi_join(bjab_ccl2 %>% filter(!is.na(value)), by=c('ionis_id'='treatment')) %>%
    arrange(label, display_name) %>%
    mutate(y = c(seq_len(sum(label=='tolerated')),
                 seq_len(sum(label=='toxic')) + sum(label=='tolerated') + 1),
           y = max(y) + 1 - y) -> ccl2_meta

  bjab_ccl2 %>%
    inner_join(ccl2_meta, by=c('treatment'='ionis_id')) -> ccl2_plot

  xlims    = c(0, ceiling(max(ccl2_plot$value)))
  ylims    = c(0.5, max(ccl2_meta$y) + 0.5)
  barwidth = 0.6
  xats     = pretty(xlims)
  tol_span = range(ccl2_meta$y[ccl2_meta$label=='tolerated']) + c(-0.5, 0.5)
  tox_span = range(ccl2_meta$y[ccl2_meta$label=='toxic'])     + c(-0.5, 0.5)

  par(mar=c(3, 2, 3, 4))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.05, labels=NA)
  axis(side=1, at=xats, labels=xats, lwd=0, line=-0.5, cex.axis=0.8)
  mtext(side=1, line=1.8, text='CCL2 (log10 fold change)', cex=0.7)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  axis(side=2, at=ccl2_meta$y, labels=ccl2_meta$display_name,
       las=2, tck=0, cex.axis=0.65)
  abline(v=0, lty=3)
  rect(xleft=0, xright=ccl2_plot$value,
       ybottom=ccl2_plot$y - barwidth/2, ytop=ccl2_plot$y + barwidth/2,
       col=ccl2_plot$color, border=NA)
  axis(side=4, at=tol_span, labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='tolerated'])
  axis(side=4, at=tox_span, labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='toxic'])
  mtext(side=4, at=mean(tol_span), line=0.1, las=2, text='tolerated', col=toxtol_colors$color[toxtol_colors$label=='tolerated'], cex=0.6)
  mtext(side=4, at=mean(tox_span), line=0.1, las=2, text='toxic',     col=toxtol_colors$color[toxtol_colors$label=='toxic'],     cex=0.6)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(ccl2_plot %>% select(display_name, label, ccl2_log10_fold_change = value),
                   'Figure 1D: BJAB CCL2 induction per ASO')

  ### E. ITC Kd ####
  # map ITC ASO names to the toxic/tolerated designations from mri_meta (MRI-based labels);
  # drop rows flagged as bad curve fits and controls with no tox/tol designation (e.g. control ASO 0)
  itc_all_asos %>%
    filter(is.na(notes) | notes != 'bad curve fit') %>%
    mutate(display_name = gsub('active |control ', '', aso)) %>%
    inner_join(mri_meta %>% select(display_name, label), by='display_name') %>%
    mutate(kd_nm = 1e9 / ka) -> itc_kd     # KD = 1/ka, expressed in nM

  itc_kd %>%
    group_by(display_name, label) %>%
    summarize(.groups='drop',
              n    = n(),
              mean = mean(kd_nm),
              l95  = lower(kd_nm),
              u95  = upper(kd_nm)) %>%
    left_join(toxtol_colors, by='label') %>%
    arrange(label, mean) %>%
    mutate(y = c(seq_len(sum(label=='tolerated')),
                 seq_len(sum(label=='toxic')) + sum(label=='tolerated') + 1),
           y = max(y) + 1 - y) -> itc_smry

  itc_kd %>%
    inner_join(itc_smry %>% select(display_name, y, color), by='display_name') -> itc_kd_plot

  xlims    = c(0, ceiling(max(itc_kd_plot$kd_nm) / 5) * 5)
  ylims    = c(0.5, max(itc_smry$y) + 0.5)
  xats     = pretty(xlims)
  barwidth = 0.3
  tol_span = range(itc_smry$y[itc_smry$label=='tolerated']) + c(-0.5, 0.5)
  tox_span = range(itc_smry$y[itc_smry$label=='toxic'])     + c(-0.5, 0.5)

  par(mar=c(3, 4, 3, 4))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats, tck=-0.05, labels=NA)
  axis(side=1, at=xats, labels=xats, lwd=0, line=-0.5, cex.axis=0.8)
  mtext(side=1, line=1.8, text=expression('K'[D]*' (nM)'), cex=0.8)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  axis(side=2, at=itc_smry$y, labels=itc_smry$display_name,
       las=2, tck=0, cex.axis=0.75)
  points(itc_kd_plot$kd_nm, itc_kd_plot$y, col=itc_kd_plot$color, bg='white', pch=21, cex=0.7)
  segments(x0=itc_smry$mean, y0=itc_smry$y - barwidth/2, y1=itc_smry$y + barwidth/2,
           lwd=2, col=itc_smry$color)
  suppressWarnings(arrows(x0=itc_smry$l95, x1=itc_smry$u95, y0=itc_smry$y,
         code=3, angle=90, length=0.04, lwd=1.5, col=itc_smry$color))
  axis(side=4, at=tol_span, labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='tolerated'])
  axis(side=4, at=tox_span, labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='toxic'])
  mtext(side=4, at=mean(tol_span), line=0.1, las=2, text='tolerated', col=toxtol_colors$color[toxtol_colors$label=='tolerated'], cex=0.6)
  mtext(side=4, at=mean(tox_span), line=0.1, las=2, text='toxic',     col=toxtol_colors$color[toxtol_colors$label=='toxic'],     cex=0.6)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(itc_kd %>% select(display_name, label, ka, kd_nm),
                   'Figure 1E: ITC binding affinity (KD) per measurement')

  ### F. NfL ####
  # x positions: two timepoints per group, gap of 1 unit between groups
  nfl_x_meta = tibble(
    ionis_id = rep(c('PBS','ASOA','ASOE','ASOD','ASOF'), each=2),
    plot_dpi = rep(c(119, 127), times=5),
    x        = c(1,2, 4,5, 7,8, 10,11, 13,14)
  )
  nfl_group_meta = nfl_x_meta %>%
    group_by(ionis_id) %>%
    summarize(x_center = mean(x), .groups='drop') %>%
    left_join(treatments %>% select(ionis_id, display_name), by='ionis_id')

  treatments %>%
    select(ionis_id, display_name) %>%
    inner_join(table2_data %>% select(display_name, label), by='display_name') %>%
    left_join(toxtol_colors, by='label') -> nfl_toxtol

  prp190327_nfl %>%
    inner_join(nfl_timepoints, by = c('timepoint'='description')) %>%
    select(-x) %>%
    inner_join(nfl_toxtol %>% select(ionis_id, display_name, color), by='ionis_id') %>%
    filter(plot_dpi %in% c(119, 127)) %>%
    inner_join(nfl_x_meta, by = c('ionis_id', 'plot_dpi')) -> nfl_plot

  nfl_plot %>%
    group_by(ionis_id, plot_dpi, x, color) %>%
    summarize(.groups='keep',
              n    = n(),
              mean = mean(nfl),
              l95  = lower(nfl),
              u95  = upper(nfl)) %>%
    ungroup() -> nfl_smry

  nfl_smry %>%
    left_join(
      nfl_smry %>% filter(plot_dpi == 119) %>% select(ionis_id, mean_119 = mean),
      by = 'ionis_id'
    ) %>%
    mutate(prop_change_119_to_127 = if_else(plot_dpi == 127, mean / mean_119 - 1, NA_real_)) %>%
    select(-mean_119) -> nfl_smry

  write_supp_table(nfl_smry %>%
                     left_join(treatments %>% select(ionis_id, display_name), by='ionis_id') %>%
                     select(display_name, plot_dpi, n, mean, l95, u95, prop_change_119_to_127),
                   'Figure 1F: NfL summary by ASO and timepoint')

  par(mar=c(3, 4, 3, 1))
  ylims = c(0, ceiling(max(nfl_plot$nfl) / 500) * 500)
  yats  = seq(0, ylims[2], 500)
  plot(NA, NA, xlim=c(0.5, 14.5), ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=c(0.5, 14.5), labels=NA, lwd.ticks=0)
  mtext(side=1, at=nfl_x_meta$x,          line=0.02, text=nfl_x_meta$plot_dpi,       cex=0.5)
  mtext(side=1, at=nfl_x_meta$x,          line=0.50, text=rep('dpi',nrow(nfl_x_meta)),       cex=0.5)
  mtext(side=1, at=nfl_group_meta$x_center, line=1.6, text=nfl_group_meta$display_name, cex=0.7)
  axis(side=2, at=yats, tck=-0.03, labels=NA)
  axis(side=2, at=yats, labels=formatC(yats, big.mark=','), lwd=0, line=-0.5, las=2, cex.axis=0.8)
  mtext(side=2, line=2.25, text='NfL (pg/mL)', cex=0.8)
  points(nfl_plot$x, nfl_plot$nfl, col=nfl_plot$color, bg='white', pch=21)
  barwidth = 0.3
  segments(x0=nfl_smry$x - barwidth/2, x1=nfl_smry$x + barwidth/2,
           y0=nfl_smry$mean, lwd=2, col=nfl_smry$color)
  arrows(x0=nfl_smry$x, y0=nfl_smry$l95, y1=nfl_smry$u95,
         code=3, angle=90, length=0.05, lwd=1.5, col=nfl_smry$color)
  axis(side=3, at=c(3.5,  8.5), labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='tolerated'])
  axis(side=3, at=c(9.5, 14.5), labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='toxic'])
  mtext(side=3, at=6,  line=0.2, text='tolerated', col=toxtol_colors$color[toxtol_colors$label=='tolerated'], cex=0.7)
  mtext(side=3, at=12, line=0.2, text='toxic',     col=toxtol_colors$color[toxtol_colors$label=='toxic'],     cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  ### G. PK ####
  # combine survival-based labels (table2_data) with MRI-based labels for Prnp ASOs
  # (pk_mri_labels computed above in panel D)
  treatments %>%
    select(ionis_id, display_name) %>%
    left_join(table2_data %>% select(display_name, label), by='display_name') %>%
    left_join(pk_mri_labels, by='ionis_id') %>%
    mutate(label = coalesce(label, mri_label)) %>%
    filter(!is.na(label), label != 'PBS') %>%
    left_join(toxtol_colors, by='label') %>%
    semi_join(mouse_500ug_1wk_qpcr_pk %>% filter(analyte == 'PK'), by=c('ionis_id'='treatment')) %>%
    arrange(label, display_name) %>%
    mutate(x = c(seq_len(sum(label=='tolerated')),
                 seq_len(sum(label=='toxic')) + sum(label=='tolerated') + 1)) -> pk_meta

  mouse_500ug_1wk_qpcr_pk %>%
    filter(analyte == 'PK') %>%
    inner_join(pk_meta, by=c('treatment'='ionis_id')) -> pk_plot

  pk_plot %>%
    group_by(treatment, display_name, x, color) %>%
    summarize(.groups='keep',
              n    = n(),
              mean = mean(value),
              l95  = lower(value),
              u95  = upper(value)) %>%
    ungroup() -> pk_smry

  xlims  = c(0.5, max(pk_meta$x) + 0.5)
  ylims  = c(0, ceiling(max(pk_plot$value) / 5) * 5)
  yats   = pretty(ylims)
  par(mar=c(3, 4, 3, 1))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=1, line=-0.5, at=pk_meta$x, labels=pk_meta$display_name,
       las=2, tck=0, lwd=0, cex.axis=0.75)
  axis(side=2, at=yats, tck=-0.03, labels=NA)
  axis(side=2, at=yats, labels=yats, lwd=0, line=-0.5, las=2, cex.axis=0.8)
  mtext(side=2, line=2.25, text='brain ASO (\u00b5g/g)', cex=0.8)
  points(pk_plot$x, pk_plot$value, col=pk_plot$color, bg='white', pch=21, cex=0.7)
  barwidth = 0.3
  segments(x0=pk_smry$x - barwidth/2, x1=pk_smry$x + barwidth/2,
           y0=pk_smry$mean, lwd=2, col=pk_smry$color)
  suppressWarnings(arrows(x0=pk_smry$x, y0=pk_smry$l95, y1=pk_smry$u95,
         code=3, angle=90, length=0.04, lwd=1.5, col=pk_smry$color))
  tol_span = range(pk_meta$x[pk_meta$label=='tolerated']) + c(-0.5, 0.5)
  tox_span = range(pk_meta$x[pk_meta$label=='toxic'])     + c(-0.5, 0.5)
  axis(side=3, at=tol_span, labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='tolerated'])
  axis(side=3, at=tox_span, labels=NA, tck=0, lwd=1.5, col=toxtol_colors$color[toxtol_colors$label=='toxic'])
  mtext(side=3, at=mean(tol_span), line=0.2, text='tolerated', col=toxtol_colors$color[toxtol_colors$label=='tolerated'], cex=0.7)
  mtext(side=3, at=mean(tox_span), line=0.2, text='toxic',     col=toxtol_colors$color[toxtol_colors$label=='toxic'],     cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(pk_plot %>% select(display_name, label, brain_aso_ug_g = value),
                   'Figure 1G: brain ASO concentration (PK) per animal')

  silence_is_golden = dev.off()
}
### end Figure 1 ####


## Figure S1. ####
tell_user('done.\nCreating Figure S1...')
mri_timeseries_tox %>%
  pivot_longer(cols = !ionis_id, names_to = 'timepoint_label', values_to = 'value') %>%
  mutate(dosed_dpi = as.integer(gsub(' dpi', '', timepoint_label))) %>%
  mutate(label = case_when(
    value != '\u2014'   ~ 'toxic',
    TRUE                ~ 'tolerated'
  )) %>%
  mutate(harvest_dpi = case_when(label == 'toxic' ~ dosed_dpi + suppressWarnings(as.numeric(value)),
                                 label == 'tolerated' ~ 150)) %>%
  inner_join(treatments, by=c('ionis_id')) %>%
  select(display_name, dosed_dpi, harvest_dpi, label) -> timeseries_tox_designations

mri_timeseries_histo %>%
  inner_join(timeseries_tox_designations, by=c('treatment'='display_name','dosed_dpi')) -> histo


histo %>%
  inner_join(toxtol_colors, by='label') %>%
  group_by(treatment, stain, dosed_dpi, harvest_dpi, animal, label, color) %>%
  summarize(.groups='keep',
            mean_score = mean(score, na.rm=T),
            n_necrotic = sum(necrosis, na.rm=T),
            n_autolytic = sum(autolytic, na.rm=T)) %>%
  ungroup() -> histo_averages

### Plot ####
if ('figure-s1'=='figure-s1') {
  
  png('display_items/figure-s1.png', width=6.5*resx, height=6*resx, res=resx)
  layout(matrix(1:4, nrow=4))
  panel = 1
  
  resx = 600
  stains_ordered = c('ASO', 'PrP', 'GFAP', 'Iba1')
  treatment_order = c('PBS', paste('ASO', 1:8))
  
  # x-axis position for each (treatment, dosed_dpi) group
  histo_averages %>%
    mutate(treatment = factor(treatment, levels=treatment_order)) %>%
    group_by(treatment, dosed_dpi) %>%
    summarise(harvest_dpi = first(harvest_dpi), .groups='drop') %>%
    arrange(treatment, dosed_dpi) %>%
    mutate(x = row_number()) -> group_meta
  
  # Per-group summary: mean and 95% CI
  histo_averages %>%
    left_join(group_meta, by=c('treatment','dosed_dpi')) %>%
    group_by(treatment, stain, dosed_dpi, x, color) %>%
    summarise(
      .groups = 'drop',
      mean_s = mean(mean_score, na.rm=TRUE),
      l95    = pmax(0, lower(mean_score)),
      u95    = upper(mean_score)
    ) -> group_smry
  
  # Individual animal points with x positions
  histo_averages %>%
    left_join(group_meta, by=c('treatment','dosed_dpi')) -> histo_plot
  
  # Treatment-level tranche metadata for x-axis brackets
  group_meta %>%
    mutate(treatment = as.character(treatment)) %>%
    group_by(treatment) %>%
    summarise(x_lo=min(x)-0.4, x_hi=max(x)+0.4, x_mid=mean(x), .groups='drop') -> tranche_meta
  
  n_groups = nrow(group_meta)
  xlims    = c(0.5, n_groups + 0.5)
  ylims    = c(0, 3.3)
  yats     = seq(0, 3, 1)
  
  for (st in stains_ordered) {
    sub_pts  = histo_plot  %>% filter(stain == st)
    sub_smry = group_smry %>% filter(stain == st)
    
    par(mar=c(4, 4, 2, 1))
    plot(NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
    
    # Light vertical separators between treatment blocks
    abline(v=tranche_meta$x_lo[-1], lty=3, col='#DDDDDD', lwd=0.5)
    
    # Individual animal points (jittered by animal within group)
    jx = sub_pts$x + (as.numeric(factor(sub_pts$animal)) - 2.5) * 0.12
    points(jx, sub_pts$mean_score,
           pch=21, bg=sub_pts$color, col=sub_pts$color, cex=0.7)
    
    # Group mean tick
    segments(x0=sub_smry$x - 0.2, y0=sub_smry$mean_s,
             x1=sub_smry$x + 0.2, y1=sub_smry$mean_s,
             lwd=2, col=sub_smry$color)
    
    # 95% CI arrows
    suppressWarnings(arrows(x0=sub_smry$x, y0=sub_smry$l95,
           x1=sub_smry$x, y1=sub_smry$u95,
           code=3, angle=90, length=0.03, lwd=1.2, col=sub_smry$color))
    
    # y-axis
    axis(side=2, at=yats, tck=-0.04, labels=NA)
    axis(side=2, at=yats, labels=yats, lwd=0, line=-0.5, las=2, cex.axis=0.6)
    mtext(side=2, line=1, text='mean IHC score', cex=0.6)
    mtext(side=3, line=0.25, text=st, cex=0.8, font=2)
    
    # dpi labels at each x position
    axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
    mtext(side=1, line=0.15, at=group_meta$x, text=group_meta$dosed_dpi, cex=0.35)
    mtext(side=1, line=0.15, at=0, adj=1, text='dose dpi:', cex=0.35)
    mtext(side=1, line=0.85, at=group_meta$x, text=round(group_meta$harvest_dpi), cex=0.35)
    mtext(side=1, line=0.85, at=0, adj=1, text='harvest dpi:', cex=0.35)
    
    # Treatment tranche brackets
    for (i in seq_len(nrow(tranche_meta))) {
      tr = tranche_meta[i, ]
      axis(side=1, at=c(tr$x_lo, tr$x_hi), labels=NA,
           tck=0.015, lwd=1, lwd.ticks=1, line=2.0)
      mtext(side=1, at=tr$x_mid, text=tr$treatment,
            line=2.0, cex=0.5, las=1, padj=1)
    }
    
    mtext(LETTERS[panel], side=3, cex=1, adj=0.0, line=0.5); panel = panel + 1
  }
  
  silence_is_golden = dev.off()
}

write_supp_table(histo_averages %>%
                   select(treatment, stain, dosed_dpi, harvest_dpi, animal, mean_ihc_score = mean_score),
                 'Figure S1: per-animal IHC scores (ASO, PrP, GFAP, Iba1 stains) over time')
### end Figure S1 ####






## Figure S2.  ####
tell_user('done.\nCreating Figure S2...')
if ('figure-s2'=='figure-s2') {

  # Combined toxic/tolerated labels: table2_data for survival ASOs, MRI for Prnp ASOs
  mri_timeseries_tox %>%
    mutate(mri_label = case_when(
      ionis_id == 'PBS'     ~ 'PBS',
      `123 dpi` != '\u2014' ~ 'toxic',
      TRUE                  ~ 'tolerated'
    )) %>%
    select(ionis_id, mri_label) -> s2_mri_labels

  treatments %>%
    select(ionis_id, display_name) %>%
    left_join(table2_data %>% select(display_name, label), by='display_name') %>%
    left_join(s2_mri_labels, by='ionis_id') %>%
    mutate(label = coalesce(label, mri_label)) %>%
    filter(!is.na(label)) %>%
    left_join(toxtol_colors, by='label') -> s2_labels

  # Build x-position metadata for a set of ionis_ids
  build_meta = function(ionis_ids) {
    s2_labels %>%
      filter(ionis_id %in% as.character(ionis_ids), label != 'PBS') %>%
      arrange(label, display_name) %>%
      mutate(x = c(seq_len(sum(label == 'tolerated')),
                   seq_len(sum(label == 'toxic')) + sum(label == 'tolerated') + 1L))
  }

  # Plot one strip-chart panel
  plot_strip = function(df, meta, title = '', ylab = '') {
    # df has columns: ionis_id (character), value
    df2   = df %>%
      mutate(value = suppressWarnings(as.numeric(value))) %>%
      inner_join(meta %>% select(ionis_id, x, color), by = 'ionis_id') %>%
      filter(!is.na(value))
    smry  = df2 %>%
      group_by(ionis_id, x, color) %>%
      summarize(mean = mean(value, na.rm = T),
                l95  = lower(value),
                u95  = upper(value), .groups = 'drop')
    xlims = c(0.5, max(meta$x) + 0.5)
    yhi   = suppressWarnings(max(df2$value, na.rm = T))
    if (!is.finite(yhi) | yhi==0) yhi = 1
    ylims = c(0, max(pretty(yhi))* 1.05)
    yats  = pretty(c(0, yhi))
    if (grepl('fob',ylab, ignore.case=T)) {
      yats = 0:8
      ylims = c(0,8)
    } else if (grepl('calbindin',ylab, ignore.case=T)) {
      yats = 0:1
      ylims = c(0,1.2)
    } else if (grepl('iba', ylab, ignore.case=T) & all(df2$value %in% 0:4)) {
      yats = 0:4
      ylims = c(0,4)
    } else if (grepl('qPCR', title, ignore.case=T)) {
      yats = 0:8/4
      ylims = c(0,2)
    }
    
    par(mar = c(3, 2, 3, 2))
    plot(NA, NA, xlim = xlims, ylim = ylims, axes = F, ann = F, xaxs = 'i', yaxs = 'i')
    axis(side = 1, at = xlims, labels = NA, lwd.ticks = 0)
    axis(side = 1, at = meta$x, labels = meta$display_name,
         las = 2, tck = 0, lwd = 0, cex.axis = 0.5, line=-0.75)
    axis(side = 2, at = yats, tck = -0.04, labels = NA)
    yatlabels = yats
    if ( !(all(df2$value %% 1 == 0)) ) { # if non-integer, it's usually qPCR and should be percents
      yatlabels = percent(yats) 
      abline(h=1, lty=3, lwd=0.5)
    }
    axis(side = 2, at = yats, labels = yatlabels, lwd = 0, line = -0.5, las = 2, cex.axis = 0.6)
    #mtext(side = 2, line = 2.0, text = ylab, cex = 0.4)
    mtext(side = 3, line = 1, text = title, cex = 0.3, at = min(xlims) + 0.05*(max(xlims)-min(xlims)), adj=0, font = 1)
    points(df2$x, df2$value, col = df2$color, bg = 'white', pch = 21, cex = 0.3)
    segments(x0 = smry$x - 0.2, x1 = smry$x + 0.2, y0 = smry$mean, lwd = 2, col = smry$color)
    smry_ci = smry %>% filter(is.finite(l95), is.finite(u95))
    if (nrow(smry_ci) > 0)
      suppressWarnings(arrows(x0 = smry_ci$x, y0 = smry_ci$l95, y1 = smry_ci$u95,
             code = 3, angle = 90, length = 0.015, lwd = 1.2, col = smry_ci$color))
    if (any(meta$label == 'tolerated')) {
      axis(side = 3, at = range(meta$x[meta$label == 'tolerated']) + c(-0.5, 0.5),
           labels = NA, tck = 0, lwd = 1.5, col = toxtol_colors$color[toxtol_colors$label == 'tolerated'])
      mtext(side=3, cex=0.3, line=.25, at = mean(range(meta$x[meta$label == 'tolerated'])), text='tolerated', col = toxtol_colors$color[toxtol_colors$label == 'tolerated'])
    }
    if (any(meta$label == 'toxic'))
      axis(side = 3, at = range(meta$x[meta$label == 'toxic']) + c(-0.5, 0.5),
           labels = NA, tck = 0, lwd = 1.5, col = toxtol_colors$color[toxtol_colors$label == 'toxic'])
    mtext(side=3, cex=0.3, line=.25, at = mean(range(meta$x[meta$label == 'toxic'])), text='toxic', col = toxtol_colors$color[toxtol_colors$label == 'toxic'])
    
    return (xlims)
  }

  # Normalised datasets: list of (label, dataset_key, data) where data has ionis_id/readout/value
  # dataset_key matches the 'dataset' column in analytic/s2_readout_labels.tsv
  s2_datasets = list(
    list(label = 'Mouse 500 µg 1x 1wk', dataset_key = 'Mouse 500ug 1x 1wk qPCR',
         data  = mouse_500ug_1wk_qpcr_pk %>%
           rename(ionis_id = treatment) %>%
           filter(analyte != 'PK') %>%
           mutate(readout  = paste(tissue, analyte)) %>%
           mutate(ionis_id = as.character(ionis_id)) %>%
           select(ionis_id, readout, value)),

    list(label = 'Mouse 500 µg 1x 1wk', dataset_key = 'Mouse 500ug 1x 1wk FOB',
         data  = mouse_500ug_1wk_fob_qpcr %>%
           rename(ionis_id = isis_no) %>%
           mutate(ionis_id = as.character(ionis_id)) %>%
           select(ionis_id, readout, value)),

    list(label = 'Mouse 500 µg 3x 10wk', dataset_key = 'Mouse 500ug 3x 10wk',
         data  = mouse_500ug_3xicv_10wk %>%
           mutate(ionis_id = as.character(ionis_id)) %>%
           select(ionis_id, readout, value)),

    list(label = 'Mouse 500 µg AxD 4wk', dataset_key = 'Mouse 500ug AxD 4wk',
         data  = mouse_axd_4wk %>%
           filter(dose=='500 ug') %>%
           mutate(ionis_id = as.character(ionis_id)) %>%
           select(ionis_id, readout, value)),

    list(label = 'Mouse 1000 µg AxD 4wk', dataset_key = 'Mouse 1000ug AxD 4wk',
         data  = mouse_axd_4wk %>%
           filter(dose=='1000 ug') %>%
           mutate(ionis_id = as.character(ionis_id)) %>%
           select(ionis_id, readout, value)),

    list(label = 'Rat 3 mg 1x 1wk', dataset_key = 'Rat 3mg 1x 1wk',
         data  = rat_3mg_1wk_qpcr %>%
           rename(ionis_id = treatment) %>%
           filter(ionis_id != 'PBS') %>%
           mutate(ionis_id = as.character(ionis_id),
                  readout  = paste(tissue, target)) %>%
           select(ionis_id, readout, value)),

    list(label = 'Rat 1 mg 3x 4wk', dataset_key = 'Rat 1mg 3x 4wk',
         data  = rat_1mg_3xicv_4wk %>%
           mutate(ionis_id = as.character(ionis_id)) %>%
           select(ionis_id, readout, value)),

    list(label = 'Rat 3 mg 1x 8wk', dataset_key = 'Rat 3mg 1x 8wk',
         data  = rat_3mg_8wk %>%
           rename(ionis_id = isis_no) %>%
           mutate(ionis_id = as.character(ionis_id)) %>%
           select(ionis_id, readout, value))
  )

  # Count panels to set figure height
  n_panels = sum(sapply(s2_datasets, function(ds) length(unique(ds$data$readout))))
  n_cols   = 6
  n_rows   = ceiling(n_panels / n_cols)

  resx = 600
  png('display_items/figure-s2.png', width = 6.5 * resx, height = 9 * resx, res = resx)
  layout(matrix(seq_len(n_cols * n_rows), nrow = n_rows, ncol = n_cols, byrow = T))
  panel = 1
  EXTLETTERS = c(LETTERS, paste0('A',LETTERS))

  for (ds in s2_datasets) {
    for (ro in sort(unique(ds$data$readout))) {
      df_ro = ds$data %>% filter(readout == ro, !is.na(value)) %>% select(ionis_id, value)
      meta  = build_meta(unique(df_ro$ionis_id))
      if (nrow(meta) == 0 || nrow(df_ro) == 0) next
      # Look up human-readable label; fall back to raw readout name if not found
      ro_label = s2_readout_labels %>%
        filter(dataset == ds$dataset_key, readout == ro) %>%
        pull(label)
      ro_label = if (length(ro_label) == 1) ro_label else ro
      xlims = plot_strip(df_ro, meta,
                 title = paste0(ds$label, '\n', ro_label),
                 ylab  = ro_label)
      mtext(EXTLETTERS[panel], side=3, cex=1, at=min(xlims), adj=1, line = 0.5); panel = panel + 1
    }
  }

  # Consolidated source data for all Figure S2 panels (long format)
  map_dfr(s2_datasets, function(ds) ds$data %>%
            mutate(value = as.character(value), dataset = ds$label, dataset_key = ds$dataset_key)) %>%
    filter(!is.na(value)) %>%
    inner_join(treatments %>% select(ionis_id, display_name), by='ionis_id') %>%
    left_join(s2_readout_labels, by=c('dataset_key'='dataset', 'readout'='readout')) %>%
    mutate(readout_label = coalesce(label, readout)) %>%
    select(dataset, display_name, readout = readout_label, value) -> fig_s2_source
  write_supp_table(fig_s2_source,
                   'Figure S2: tolerability readouts (FOB, qPCR, IHC) across rodent studies, all panels')

  silence_is_golden = dev.off()
}




## Figure 2.  #### 
tell_user('done.\nCreating Figure 2...')
if ('figure-2'=='figure-2') {
  
  resx=600
  png('display_items/figure-2.png',width=6.5*resx,height=4.5*resx,res=resx)

  # Row 1: A (UMAP), B (volcano: toxic vs PBS), C (DGE scores binary)
  # Row 2: D (DGE scores 23-0afa), E (gene length vs L2FC), F (base mean vs L2FC)
  layout_matrix = matrix(c(1,1,2,2,3,3,
                           4,4,5,5,6,6), nrow=2, byrow=T)
  layout(layout_matrix)
  par(mar=c(3,4,3,1))
  panel = 1

  ### A. UMAP ####
  umap_pad = 1.5
  umap_xlim = c(floor(min(umap_data$umap1)) - umap_pad, ceiling(max(umap_data$umap1)) + umap_pad)
  umap_ylim = c(floor(min(umap_data$umap2)) - umap_pad, ceiling(max(umap_data$umap2)) + umap_pad)
  plot(NA, xlim=umap_xlim, ylim=umap_ylim, xlab='', ylab='', axes=FALSE, frame.plot=FALSE)
  for (cond in unique(umap_data$condition)) {
    for (inoc in unique(umap_data$inoculum)) {
      sub = umap_data %>% filter(condition==cond, inoculum==inoc)
      if (nrow(sub) == 0) next
      pch_val = if (inoc == 'RML') 16 else 1
      points(sub$umap1, sub$umap2,
             col=condition_colors[cond], pch=pch_val, cex=1.2,
             lwd=if(inoc=='CBH') 2 else 1)
    }
  }
  axis(1, at=pretty(umap_xlim), lwd=0.5, cex.axis=0.8)
  axis(2, at=pretty(umap_ylim), las=1, lwd=0.5, cex.axis=0.8)
  mtext('UMAP 1', side=1, line=2, cex=0.8)
  mtext('UMAP 2', side=2, line=2, cex=0.8)
  legend('topright', bty='n', cex=0.7,
         legend = c('toxic','tolerated','PBS','RML','CBH'),
         col    = c(condition_colors[c('toxic','tolerated','PBS')], 'black','black'),
         pch    = c(16, 16, 16, 16, 1),
         pt.lwd = c(1, 1, 1, 1, 2))
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1
  
  write_supp_table(umap_data %>% select(sample_id, aso_name, condition, inoculum, umap1, umap2),
                   'Figure 2A: UMAP coordinates per RNA-seq sample')

  ### B. Volcano: toxic ASO vs. PBS (21-01ac, RML-infected) ####
  volcano_data %>%
    left_join(cell_type_enrichments %>% select(gene, celltype, celltype_color = color),
              by = 'gene') %>%
    mutate(celltype_color = replace_na(celltype_color, '#D7D7D7'),
           celltype       = replace_na(celltype, '0_none'),
           pval           = padj) %>%
    arrange(celltype) -> volcano_plot_data

  maxy=30
  alpha_threshold = 0.05
  bonf_threshold = 0.05 / nrow(volcano_plot_data)
  volcano(volcano_plot_data, clipwidth = 5, maxy = maxy,
          colorcolname = 'celltype_color', pvalcolname = 'pval',
          log2fccolname = 'log2fold_change',
          title = '')
  mtext(side=3, line=0.25, text='Toxic ASO vs. PBS\n(RML-infected)', cex=0.7)
  rect(xleft=-5, xright=5, ybottom=0, ytop=-log10(alpha_threshold), col='#FFFFFFAA', border = NA)
  abline(h=c(-log10(alpha_threshold)), lty=3)
  mtext(side=4, line=0.15, at=c(-log10(alpha_threshold)), las=2, text=c('adjusted P < 0.05'), cex=0.3)
  mtext(LETTERS[panel], side=3, cex=1, adj = -0.2, line = 0.5); panel = panel + 1

  volcano_plot_data %>%
    filter(celltype == '3_astrocyte', -log10(padj) > 10 & (-log10(padj) > 17 | log2fold_change < -1)) %>%
    arrange(pval) %>%
    head(15) -> astro_label
  text(clipdist(astro_label$log2fold_change, -5, 5), pmin(maxy, -log10(astro_label$pval)),
       labels = astro_label$gene, cex = 0.6, pos = 2, font=3)

  volcano_plot_data %>%
    distinct(celltype, celltype_color) %>%
    mutate(celltype_disp = gsub('[0-9]_', '', celltype)) %>%
    arrange(celltype_disp) -> vol_leg
  legend('topright', legend = vol_leg$celltype_disp, col = vol_leg$celltype_color,
         pch = 20, bty = 'n', cex = 0.55)

  write_supp_table(volcano_plot_data %>%
                     select(gene, log2fold_change, base_mean, pvalue, padj, celltype) %>%
                     arrange(desc(padj)),
                   'Figure 2B: toxic ASO vs. PBS differential expression (per gene)')
  
  write_supp_table(g_tox_genes, 'Toxicity gene signature (g_tox_genes)')
  
  ### C. DGE scores (24-0432, continuous) ####
  # Label ASOs by the canonical naming in treatments.tsv (the aso_name/label columns in the
  # study-24 metadata used an outdated, incorrect lettering scheme).
  # Include the PBS reference at G_tox = 0 (PBS is the contrast reference, so it has no
  # ASO-vs-PBS score row of its own; its score is 0 by definition).
  dge_scores_24 %>%
    filter(condition != 'modified') %>%
    left_join(treatments %>% select(ionis_id, display_name), by = 'ionis_id') %>%
    mutate(label = paste0(display_name, ' (', inoculum, ')')) %>%
    bind_rows(tibble(contrast_id = 'RML_PBS_vs_RML_PBS', condition = 'PBS',
                     inoculum = 'RML', gtox_score = 0, display_name = 'PBS', label = 'PBS')) %>%
    arrange(gtox_score) -> dge_main

  n_asos = nrow(dge_main)
  score_xticks = pretty(range(dge_main$gtox_score), n=5)
  score_xlim = range(score_xticks)
  par(mar=c(3,9,3,1))
  plot(NA, xlim=score_xlim, ylim=c(0.5, n_asos+0.5),
       xlab='', ylab='', axes=FALSE, frame.plot=FALSE)
  for (i in seq_len(n_asos)) {
    row = dge_main[i, ]
    pch_val = if (row$inoculum == 'RML') 16 else 1
    points(row$gtox_score, i,
           col=condition_colors[row$condition], pch=pch_val, cex=1.4,
           lwd=if(row$inoculum=='CBH') 2 else 1)
  }
  axis(1, at=score_xticks, lwd=0.5, cex.axis=0.8)
  axis(2, at=seq_len(n_asos), labels=dge_main$label,
       las=2, cex.axis=0.75, lwd=0, tck=0)
  mtext('tox score', side=1, line=2, cex=0.8)
  abline(v=0, lty=3, col='grey50')
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(dge_main %>% select(contrast_id, display_name, condition, inoculum, gtox_score),
                   'Figure 2C: G_tox DGE scores per ASO (study 24-0432)')

  # Continuous DGE score per contrast: sum of L2FC * reference sign over g_tox genes.
  # Mirrors compute_gtox_continuous() in dan_contrib/src/analysis.R.
  # No padj cutoff — every g_tox gene contributes; signal accumulates, noise cancels.
  # Computed here (formerly in Figure 2E); used by panel D below and by Figure 3 (barplot + volcanoes).
  dge_scores_all = contrasts %>%
    inner_join(g_tox_genes %>% select(gene, sign), by = 'gene') %>%
    group_by(dataset, group1, group2) %>%
    summarise(
      gtox_score = sum(log2fold_change * sign, na.rm = TRUE),
      n_gtox     = n(),
      .groups    = 'drop'
    )

  write_supp_table(dge_scores_all %>% filter(!grepl('ASO1-MOE|ASO1-OMe2', group1)),
                   'Figure 2C, 2D, 3E: DGE scores per contrast')

  ### D. DGE scores dot plot (23-0afa) ####
  mri_timeseries_tox %>%
    mutate(mri_label = if_else(`90 dpi` != '\u2014', 'toxic', 'tolerated')) %>%
    select(ionis_id, mri_label) -> mri_90dpi_labels

  dge_scores_all %>%
    filter(dataset == '23-0afa') %>%
    left_join(treatments %>% select(ionis_id, display_name),
              by = c('group1' = 'ionis_id')) %>%
    left_join(mri_90dpi_labels, by = c('group1' = 'ionis_id')) %>%
    mutate(mri_label = coalesce(mri_label, 'tolerated')) %>%
    left_join(toxtol_colors, by = c('mri_label' = 'label')) %>%
    arrange(gtox_score) -> dge_23

  n_asos_23   = nrow(dge_23)
  score_xlim  = range(pretty(dge_23$gtox_score))
  score_xats  = pretty(dge_23$gtox_score)

  par(mar = c(3, 6, 3, 1))
  plot(NA, xlim = score_xlim, ylim = c(0.5, n_asos_23 + 0.5),
       xlab = '', ylab = '', axes = FALSE, frame.plot = FALSE)
  abline(v = 0, lty = 3, col = 'grey50', lwd = 0.5)
  points(dge_23$gtox_score, seq_len(n_asos_23),
         col = dge_23$color, pch = 16, cex = 1.2)
  axis(side = 1, at = score_xats, tck = -0.04, labels = NA)
  axis(side = 1, at = score_xats, labels = score_xats,
       lwd = 0, line = -0.5, cex.axis = 0.6)
  axis(side = 2, at = seq_len(n_asos_23), labels = dge_23$display_name,
       las = 2, cex.axis = 0.8, lwd = 0, tck = 0)
  mtext(side = 1, line = 2, text = 'tox score', cex = 0.8)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  volcano_data %>%
    clean_names() %>%
    inner_join(lens %>% select(-tss, -aaa), by='gene') %>%
    arrange(desc(padj)) -> dge_properties

  col_down  = '#8856A7'
  col_up    = '#2CA25F'
  col_nonsig = '#AAAAAA'
  dge_properties = dge_properties %>%
    mutate(pt_col = case_when(
      padj >= 0.05               ~ col_nonsig,
      log2fold_change < 0        ~ col_down,
      TRUE                       ~ col_up
    ))

  length_xticks   = c(1e2, 1e3, 1e4, 1e5, 1e6)
  basemean_xticks = c(1, 10, 100, 1000, 10000)
  xats = rep(1:9, 6) * 10^(rep(0:5, each=9))

  ### E. Gene length vs L2FC ####
  ylims_lfc = range(pretty(dge_properties$log2fold_change))
  par(mar=c(3,4,3,1))
  plot(NA, NA,
       xlim = range(length_xticks), ylim = ylims_lfc,
       axes = FALSE, ann = FALSE, log = 'x')
  abline(h = 0, lwd = 0.5, col = '#888888')
  points(dge_properties$length, dge_properties$log2fold_change,
         pch = 20, cex = 0.25, col = dge_properties$pt_col)
  axis(side = 1, at = length_xticks, labels = NA, tck = -0.04)
  axis(side = 1, at = length_xticks, labels = scales::comma(length_xticks),
       lwd=0, cex.axis = 0.7, line=-0.5)
  axis(side = 1, at = xats, labels = NA, tck = -0.015)
  axis(side = 2, at = pretty(ylims_lfc), tck = -0.02, labels = NA)
  axis(side = 2, at = pretty(ylims_lfc), labels = pretty(ylims_lfc),
       lwd = 0, line = -0.5, las = 2, cex.axis = 0.6)
  mtext(side = 1, line = 1.75,   text = 'Gene length (bp)', cex = 0.8)
  mtext(side = 2, line = 1.75, text = expression(log[2]~fold~change), cex = 0.8)
  legend('topright', legend = c('down', 'up', 'n.s.'),
         col = c(col_down, col_up, col_nonsig), pch = 20,
         text.col = c(col_down, col_up, col_nonsig),
         bty = 'n', cex = 0.55, pt.cex = 0.8)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  ### F. Base mean vs L2FC ####
  par(mar=c(3,4,3,1))
  plot(NA, NA,
       xlim = range(basemean_xticks), ylim = ylims_lfc,
       axes = FALSE, ann = FALSE, log = 'x')
  abline(h = 0, lwd = 0.5, col = '#888888')
  points(dge_properties$base_mean, dge_properties$log2fold_change,
         pch = 20, cex = 0.25, col = dge_properties$pt_col)
  axis(side = 1, at = basemean_xticks, labels = NA,
       tck = -0.04)
  axis(side = 1, at = basemean_xticks, labels = scales::comma(basemean_xticks),
       lwd=0, cex.axis = 0.7, line=-0.5)
  axis(side = 1, at = xats, labels = NA, tck = -0.015)
  axis(side = 2, at = pretty(ylims_lfc), tck = -0.02, labels = NA)
  axis(side = 2, at = pretty(ylims_lfc), labels = pretty(ylims_lfc),
       lwd = 0, line = -0.5, las = 2, cex.axis = 0.6)
  mtext(side = 1, line = 1.75,   text = 'Base expression level (TPM)', cex = 0.8)
  mtext(side = 2, line = 1.75, text = expression(log[2]~fold~change), cex = 0.8)
  legend('topright', legend = c('down', 'up', 'n.s.'),
         col = c(col_down, col_up, col_nonsig), pch = 20,
         text.col = c(col_down, col_up, col_nonsig),
         bty = 'n', cex = 0.55, pt.cex = 0.8)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(dge_properties %>%
                     select(gene, gene_length_bp = length, base_mean, log2fold_change, padj) %>%
                     arrange(desc(padj)),
                   'Figure 2E,F: per-gene length, base expression, and toxic-signature fold change')

  silence_is_golden = dev.off()
}
### end Figure 2 ####


## Figure 3.  ####
tell_user('done.\nCreating Figure 3...')
if ('figure-3'=='figure-3') {

  resx=600
  png('display_items/figure-3.png',width=6.5*resx,height=5*resx,res=resx)

  # Row 1: A-D volcanoes of RML-infected vs. uninfected at 60, 89, 120, 149 dpi (21-01ac)
  # Row 2: E DGE score barplot across timepoints/tranches (21-01ac)
  layout_matrix = matrix(c(1,2,3,4,
                           5,5,5,5), nrow=2, byrow=T)
  layout(layout_matrix)
  panel = 1

  maxy            = 30
  alpha_threshold = 0.05

  ### A-D. Volcanoes: RML-infected vs. uninfected over disease course (21-01ac) ####
  volc_meta = tribble(
    ~dpi, ~group1,       ~group2,
    60,   'RML_60_NO',   'CBH_60_NO',
    89,   'RML_89_NO',   'CBH_89_NO',
    120,  'RML_120_PBS', 'CBH_120_PBS',
    149,  'RML_149_PBS', 'CBH_149_PBS'
  )

  for (i in seq_len(nrow(volc_meta))) {
    vm = volc_meta[i, ]
    contrasts %>%
      filter(dataset == '21-01ac', group1 == vm$group1, group2 == vm$group2) %>%
      left_join(cell_type_enrichments %>% select(gene, celltype, celltype_color = color),
                by = 'gene') %>%
      mutate(celltype_color = replace_na(celltype_color, '#D7D7D7'),
             celltype       = replace_na(celltype, '0_none'),
             pval           = padj) %>%
      arrange(celltype) -> rml_plot_data

    par(mar = c(3, 4, 3, 1))
    volcano(rml_plot_data, clipwidth = 5, maxy = maxy,
            colorcolname = 'celltype_color', pvalcolname = 'pval',
            log2fccolname = 'log2fold_change',
            title = '')
    mtext(side=3, line=0.25, text=paste0(vm$dpi, ' dpi'), cex=0.7)
    rect(xleft=-5, xright=5, ybottom=0, ytop=-log10(alpha_threshold), col='#FFFFFFAA', border = NA)
    abline(h=c(-log10(alpha_threshold)), lty=3)

    rml_plot_data %>%
      filter(celltype %in% c('3_astrocyte','2_microglia')) %>%
      filter(-log10(padj) > 2 &
               (rank(padj) < 3 | 
               (-log10(padj) > 15 & (-log10(padj) > 17 | log2fold_change > 1)))) %>%
      arrange(pval) %>%
      head(15) -> rml_astro_label
    if (nrow(rml_astro_label) > 0) {
      par(xpd=T)
      text(clipdist(rml_astro_label$log2fold_change, -5, 5), pmin(maxy, -log10(rml_astro_label$pval)),
           labels = rml_astro_label$gene, cex = 0.5, pos = 4, font = 3)
      par(xpd=F)
    }

    if (i == 1) {
      rml_plot_data %>%
        distinct(celltype, celltype_color) %>%
        mutate(celltype_disp = gsub('[0-9]_', '', celltype)) %>%
        arrange(celltype_disp) -> rml_vol_leg
      legend('topleft', legend = rml_vol_leg$celltype_disp, col = rml_vol_leg$celltype_color,
             pch = 20, bty = 'n', cex = 0.5)
    }
    mtext(LETTERS[panel], side=3, cex=1, adj = -0.2, line = 0.5); panel = panel + 1
  }

  write_supp_table(contrasts %>%
                     filter(dataset == '21-01ac') %>%
                     inner_join(volc_meta, by = c('group1', 'group2')) %>%
                     select(dpi, gene, log2fold_change, base_mean, pvalue, padj) %>%
                     arrange(desc(padj)),
                   'Figure 3A-D: RML-infected vs. uninfected differential expression by timepoint')

  ### E. DGE scores barplot (21-01ac, continuous) ####
  # dge_scores_all computed in Figure 2
  score_panel_meta = read_tsv('analytic/score_panel_meta.tsv', col_types=cols()) %>%
    filter(!is.na(x))

  dge_scores_all %>%
    inner_join(score_panel_meta, by = c('dataset', 'group1', 'group2')) %>%
    mutate(
      dpi = as.integer(gsub('^(?:CBH|RML)_(\\d+)_.*', '\\1', group1, perl = TRUE)),
      bar_color = case_when(
        color == 'gray'      ~ '#A7A7A7',
        color == 'tolerated' ~ toxtol_colors$color[toxtol_colors$label == 'tolerated'],
        color == 'toxic'     ~ toxtol_colors$color[toxtol_colors$label == 'toxic']
      )
    ) %>%
    mutate(tranche = gsub('PBS ','PBS\n',tranche)) -> score_plot_data

  write_supp_table(score_plot_data %>%
                     select(dataset, group1, group2, dpi, tranche, gtox_score, n_gtox),
                   'Figure 3E: G_tox DGE scores across the disease course (study 21-01ac)')

  # Tranche summary: x range and label per tranche
  tranche_meta = score_plot_data %>%
    filter(!is.na(tranche)) %>%
    group_by(tranche) %>%
    summarise(x_lo = min(x) - 0.4, x_hi = max(x) + 0.4, x_mid = mean(x), .groups = 'drop')

  xlims = c(0.5, max(score_plot_data$x) + 0.5)
  ylims = range(pretty(c(0, score_plot_data$gtox_score)))
  yats  = pretty(c(0, score_plot_data$gtox_score))

  par(mar=c(5,4,3,1))
  plot(NA, NA, xlim = xlims, ylim = ylims, axes = F, ann = F, xaxs = 'i', yaxs = 'i')
  abline(h = 0, lwd = 0.5, col = '#888888')

  # Primary x-axis: dpi labels at bar positions
  axis(side = 1, at = xlims, labels = NA, lwd.ticks = 0)
  mtext(side = 1, line=0.15, at = score_plot_data$x, text = score_plot_data$dpi, cex=0.6)
  mtext(side = 1, line=0.15, at = -0.5, text = 'dpi: ', cex=0.6)

  # Tranche brackets: upward-pointing ticks at bracket ends, label at midpoint
  for (i in seq_len(nrow(tranche_meta))) {
    tr = tranche_meta[i, ]
    axis(side = 1, at = c(tr$x_lo, tr$x_hi), labels = NA,
         tck = 0.025, lwd = 1, lwd.ticks = 1, line = 1.5)
    mtext(side = 1, at = tr$x_mid, text = tr$tranche,
          line = 1.2, cex = 0.5, las = 1, padj=1)
  }

  axis(side = 2, at = yats, tck = -0.02, labels = NA)
  axis(side = 2, at = yats, labels = yats, lwd = 0, line = -0.5, las = 2, cex.axis = 0.8)
  mtext(side = 2, line = 1.75, text = 'tox score', cex = 0.8)

  bar_w = 0.4
  solid_bars  = score_plot_data %>% filter(fill == 'solid')
  hollow_bars = score_plot_data %>% filter(fill == 'hollow')

  rect(xleft   = solid_bars$x - bar_w,
       xright  = solid_bars$x + bar_w,
       ybottom = 0,
       ytop    = solid_bars$gtox_score,
       col     = solid_bars$bar_color,
       border  = NA)

  rect(xleft   = hollow_bars$x - bar_w,
       xright  = hollow_bars$x + bar_w,
       ybottom = 0,
       ytop    = hollow_bars$gtox_score,
       col     = '#FFFFFF',
       border  = hollow_bars$bar_color)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  silence_is_golden = dev.off()
}
### end Figure 3 ####


## Figure 4.  ####
tell_user('done.\nCreating Figure 4...')
if ('figure-4'=='figure-4') {

  resx=600
  png('display_items/figure-4.png',width=6.5*resx,height=3.25*resx,res=resx)

  # Row 1: A (ASO 6 dose-response), B (divalent siRNA dose-response)
  layout(matrix(c(1,2), nrow=1))
  par(mar=c(3,4,3,1))
  panel = 1

  dp_dose_meta = tibble(dose=c(0,70,244),x=c(1,2,3))
  inoc_meta = tibble(dp_inoculum=c('RML','none'), drwwo_inoculum=c('RML','CBH'), disp_inoc=c('RML','uninfected'), color=c('#902323','#23CDBE'))
  dp_qpcr %>%
    separate(sample_group, sep=',', into = c('dose','inoculum'), fill='right') %>%
    mutate(inoculum = trimws(replace_na(inoculum,'none'))) %>%
    mutate(dose = replace_na(suppressWarnings(as.numeric(gsub(' ug','',dose))), 0)) %>%
    inner_join(dp_dose_meta, by='dose') %>%
    inner_join(inoc_meta, by=c('inoculum'='dp_inoculum')) %>%
    rename(residual_rna=mean) -> dp

  drwwo_raw = read_tsv('analytic/PRP230201_dose_response_qpcr.tsv', col_types=cols()) %>% clean_names() %>% rename(ionis_id = aso)
  drwwo_dose_meta = tibble(dose=c(0,100,300),x=c(1,2,3))
  drwwo_raw %>%
    inner_join(treatments %>% select(ionis_id, display_name), by='ionis_id') %>%
    inner_join(drwwo_dose_meta, by='dose') %>%
    inner_join(inoc_meta, by=c('inoculum'='drwwo_inoculum')) %>%
    rename(residual_rna = prnp_wholehemi_normed) %>%
    select(display_name, inoculum, disp_inoc, x, dose, color, residual_rna) -> drwwo

  ### A. ASO ####
  xlims = c(0.5, 3.5)
  ylims = c(0, 1.5)
  yats = 0:6/4
  ybigs = 0:3/2
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  mtext(side=1, at=drwwo_dose_meta$x, line=0.25, text=drwwo_dose_meta$dose)
  axis(side=2, at=yats, tck=-0.02, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, line=-0.5, las=2)
  mtext(side=1, line=1.6, text='active ASO 6 dose (µg)')
  mtext(side=2, line=2.7, text='residual Prnp mRNA')
  abline(h=1, lty=3)
  points(drwwo$x, drwwo$residual_rna, col=alpha(drwwo$color, ci_alpha), pch=19)
  drwwo %>%
    group_by(inoculum, x, dose, color) %>%
    summarize(.groups='keep',
              n = n(),
              mean = mean(residual_rna),
              sd   =   sd(residual_rna),
              l95 = lower(residual_rna),
              u95 = upper(residual_rna)) %>%
    ungroup() %>%
    arrange(x, dose, inoculum) -> drwwo_smry
  drwwo %>%
    group_by(x, dose) %>%
    summarize(.groups='keep',
              pval = t.test(residual_rna[disp_inoc=='RML'], residual_rna[disp_inoc=='uninfected'])$p.value,
              diff = mean(residual_rna[disp_inoc=='uninfected']) - mean(residual_rna[disp_inoc=='RML'])) %>%
    ungroup() -> drwwo_p
  barwidth = 0.3
  segments(x0=drwwo_smry$x - barwidth/2, x1=drwwo_smry$x + barwidth/2, y0=drwwo_smry$mean, lwd=2, col=drwwo_smry$color)
  arrows(x0=drwwo_smry$x, y0=drwwo_smry$l95, y1=drwwo_smry$u95, code=3, angle=90, length=0.05, lwd=1.5, col=drwwo_smry$color)
  mtext(side=3, at=drwwo_p$x[2:3], line=0, cex=0.7, text=paste0('P',format_p(drwwo_p$pval[2:3])))
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  drwwo_smry %>%
    inner_join(drwwo_p, by=c('x','dose')) -> drwwo_out
  write_supp_table(drwwo_out, 'Figure 4A: ASO 6 potency in RML 105 dpi vs. uninfected mice at 30 days post-dose.')

  ### B. divalent siRNA ####
  xlims = c(0.5, 3.5)
  ylims = c(0, 1.5)
  yats = 0:6/4
  ybigs = 0:3/2
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  mtext(side=1, at=dp_dose_meta$x, line=0.25, text=dp_dose_meta$dose)
  axis(side=2, at=yats, tck=-0.02, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, line=-0.5, las=2)
  mtext(side=1, line=1.6, text='siRNA 1682-s4 dose (µg)')
  mtext(side=2, line=2.7, text='residual Prnp mRNA')
  abline(h=1, lty=3)
  points(dp$x, dp$residual_rna, col=alpha(dp$color, ci_alpha), pch=19)
  dp %>%
    group_by(inoculum, x, dose, color) %>%
    summarize(.groups='keep',
              n = n(),
              mean = mean(residual_rna),
              sd   =   sd(residual_rna),
              l95 = lower(residual_rna),
              u95 = upper(residual_rna)) %>%
    ungroup() %>%
    arrange(x, dose, inoculum) -> dp_smry
  dp %>%
    group_by(x, dose) %>%
    summarize(.groups='keep',
              pval = t.test(residual_rna[disp_inoc=='RML'], residual_rna[disp_inoc=='uninfected'])$p.value,
              diff = mean(residual_rna[disp_inoc=='uninfected']) - mean(residual_rna[disp_inoc=='RML'])) %>%
    ungroup() -> dp_p
  barwidth = 0.3
  segments(x0=dp_smry$x - barwidth/2, x1=dp_smry$x + barwidth/2, y0=dp_smry$mean, lwd=2, col=dp_smry$color)
  arrows(x0=dp_smry$x, y0=dp_smry$l95, y1=dp_smry$u95, code=3, angle=90, length=0.05, lwd=1.5, col=dp_smry$color)
  mtext(side=3, at=dp_p$x[2:3], line=0, cex=0.7, text=paste0('P',format_p(dp_p$pval[2:3])))
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  dp_smry %>%
    inner_join(dp_p, by=c('x','dose')) -> dp_out
  write_supp_table(dp_out, 'Figure 4B: Divalent siRNA 1682-s4 potency in RML 105 dpi vs. uninfected mice at 30 days post-dose.')

  par(xpd=T)
  legend(x=-0.6, y=-0.1, inoc_meta$disp_inoc, col=inoc_meta$color, bty='n', pch=19, lwd=1, cex=0.7)
  par(xpd=F)

  silence_is_golden = dev.off()
}
### end Figure 4 ####



## Figure S5.  ####
tell_user('done.\nCreating Figure S5...')
if ('figure-s5'=='figure-s5') {

  resx=600
  png('display_items/figure-s5.png',width=6.5*resx,height=5*resx,res=resx)
  
  alpha_threshold = 0.05
  
  sn20_supp2 %>%
    select(GeneID, log2FC_4wpi_main, log2FC_8wpi_main, log2FC_12wpi_main, log2FC_14wpi_main, log2FC_16wpi_main, log2FC_18wpi_main, log2FC_20wpi_main, log2FC_term_main,
           Pval_4wpi_main, Pval_8wpi_main, Pval_12wpi_main, Pval_14wpi_main, Pval_16wpi_main, Pval_18wpi_main, Pval_20wpi_main, Pval_term_main) %>%
    pivot_longer(cols=-GeneID) %>%
    rename(gene = GeneID) %>%
    mutate(wpi = gsub('wpi','',gsub('(log2FC|Pval)_','',gsub('_main','',name)))) %>%
    mutate(dpi_label = case_when(wpi=='term' ~ 'terminal',
                                 TRUE ~ paste0(as.character(suppressWarnings(as.integer(wpi))*7), ' dpi'))) %>%
    mutate(dpi = case_when(wpi=='term' ~ 170,
                           TRUE ~ suppressWarnings(as.integer(wpi))*7)) %>%
    mutate(variable = tolower(gsub('_.*','',name))) %>%
    pivot_wider(id_cols=c(gene, wpi, dpi_label, dpi), names_from=variable, values_from=value) %>%
    group_by(wpi, dpi_label, dpi) %>%
    mutate(p_bonf = pmin(1, pval * n())) %>%
    ungroup() %>%
    mutate(p_color = case_when(pval >= alpha_threshold ~ '#C9C9C9',
                             pval < alpha_threshold & p_bonf >= alpha_threshold & log2fc < 0 ~ '#FFCCCC',
                             pval < alpha_threshold & p_bonf >= alpha_threshold & log2fc > 0 ~ '#CCCCFF',
                             p_bonf < alpha_threshold & log2fc < 0 ~ '#BB1111',
                             p_bonf < alpha_threshold & log2fc > 0 ~ '#1111BB')) %>%
    left_join(cell_type_enrichments %>% select(gene, celltype, celltype_color=color), by='gene') %>%
    mutate(celltype_color = replace_na(celltype_color, '#D7D7D7')) %>%
    mutate(celltype = replace_na(celltype, '0_none')) %>%
    arrange(dpi, celltype) -> sn20
  
  write_supp_table(sn20 %>%
                     select(gene, dpi, dpi_label, log2fc, pval, p_bonf, celltype) %>%
                     arrange(desc(p_bonf)),
                   'Figure S5A-L: Sorce & Nuvolone scrapie differential expression by timepoint')

  layout_matrix = matrix(c(1:12, rep(13,4)), nrow=4, byrow=T)
  layout(layout_matrix, heights=c(1,1,1,0.3))
  par(mar=c(3,4,3,1))
  panel = 1
  
  ### re-plot of Aguzzi's data ####
  for (each_dpi in unique(sn20$dpi)) {
    subs = sn20 %>% filter(dpi==each_dpi)
    volcano(subs, title=paste0(each_dpi,' dpi'), maxy=30, colorcolname='celltype_color')
    rect(xleft=-5, xright=5, ybottom=0, ytop=-log10(alpha_threshold), col='#FFFFFFDD', border = NA)
    rect(xleft=-5, xright=5, ybottom=-log10(alpha_threshold), ytop=-log10(alpha_threshold/nrow(subs)), col='#FFFFFFAA', border = NA)
    abline(h=c(-log10(alpha_threshold), -log10(alpha_threshold/nrow(subs))), lty=3)
    mtext(side=4, line=0.15, at=c(-log10(alpha_threshold), -log10(alpha_threshold/nrow(subs))), las=2, text=c('nominal','Bonferroni'), cex=0.3)
    mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1
  }
  
  
  
  sn20_dge_scores = sn20 %>%
    inner_join(g_tox_genes %>% select(gene, sign), by = 'gene') %>%
    group_by(wpi, dpi_label, dpi) %>%
    summarise(
      .groups    = 'keep',
      gtox_score = sum(log2fc * sign, na.rm = TRUE),
      n_gtox     = n()
    ) %>%
    ungroup()
  
  contrasts %>%
    filter(group1=='RML_60_NO'   & group2=='CBH_60_NO' |
           group1=='RML_89_NO'   & group2=='CBH_89_NO' |
           group1=='RML_120_PBS' & group2=='CBH_120_PBS' | 
           group1=='RML_149_PBS' & group2=='CBH_149_PBS') %>%
    mutate(our_dpi = case_when(
      group1=='RML_60_NO'   ~ 60,
      group1=='RML_89_NO'   ~ 89,
      group1=='RML_120_PBS' ~ 120,
      group1=='RML_149_PBS' ~ 149
    )) %>%
    mutate(aguzzi_dpi = case_when(
      group1=='RML_60_NO'   ~ 56,
      group1=='RML_89_NO'   ~ 84,
      group1=='RML_120_PBS' ~ 126,
      group1=='RML_149_PBS' ~ 140
    )) %>%
    inner_join(sn20, by=c('gene','aguzzi_dpi'='dpi')) %>%
    arrange(aguzzi_dpi, celltype) -> our_data_vs_aguzzi
  
  write_supp_table(our_data_vs_aguzzi %>%
                     select(our_dpi, aguzzi_dpi, gene,
                            this_study_log2fc = log2fold_change, sorce_nuvolone_log2fc = log2fc,
                            this_study_padj = padj, sorce_nuvolone_p_bonf = p_bonf, celltype) %>%
                     arrange(desc(this_study_padj)),
                   'Figure S5 (comparison): this study vs. Sorce & Nuvolone fold changes, matched timepoints')

  ### comparison of our data vs. Aguzzi's ####
  xlims = c(-5, 5)
  xats = -5:5
  xbigs = c(-5, 5)
  yats = -5:5
  ybigs = c(-5, 5)
  ylims = c(-5, 5)
  for (each_dpi in sort(unique(our_data_vs_aguzzi$aguzzi_dpi))) {
    subs = our_data_vs_aguzzi %>% filter(aguzzi_dpi==each_dpi)
    our_dpi = subs$our_dpi[1]
    subs$plot_color = alpha(subs$celltype_color,
                            case_when(subs$padj > 0.05 & subs$p_bonf > 0.05 ~ 0.05,
                                      TRUE ~ 1))
    plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
    points(x=subs$log2fold_change, y=subs$log2fc,
         col=subs$plot_color, pch=20, cex=1)
    signif = subs %>%
      filter(padj < 0.05 | p_bonf < 0.05)
    m = lm(log2fc ~ log2fold_change, data=signif)
    abline(m, lwd=0.5)
    pearsons_obj = cor.test(signif$log2fc, signif$log2fold_change)
    mtext(side=3, line=0.4, cex=0.4,text=paste0('r^2=',formatC(pearsons_obj$estimate,digits=2),' P',format_p(pearsons_obj$p.value)))
    axis(side=1, pos=0, at=xats,  tck=-0.02, labels=NA)
    axis(side=1, pos=0, at=xbigs, tck=-0.05, labels=NA)
    axis(side=1, pos=0, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
    axis(side=2, pos=0, at=yats,  tck=-0.02, labels=NA)
    axis(side=2, pos=0, at=ybigs, tck=-0.05, labels=NA)
    axis(side=2, pos=0, at=ybigs, labels=ybigs, lwd=0, line=-0.5, las=2)
    mtext(side=1, line=1, text=paste0('L2FC\n(this study, ',our_dpi,' dpi)'), cex=0.6)
    mtext(side=2, line=0.5, text=paste0('L2FC\n(Sorce & Nuvolone, ',each_dpi,' dpi)'), cex=0.6)
    mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1
  }
  
  ### legend panel (13, no letter assigned) ####
  par(mar=c(0,0,0,0))
  plot(NA, NA, xlim=0:1, ylim=0:1, axes=F, ann=F)
  our_data_vs_aguzzi %>% 
    distinct(celltype, celltype_color) %>%
    mutate(celltype_disp = gsub('[0-9]_','',celltype)) %>%
    arrange(celltype_disp) -> leg
  legend('center',leg$celltype_disp,col=leg$celltype_color,pch=20,horiz=T, title='enriched cell type',title.col='#000000',bty='n', cex=0.8)
  

  silence_is_golden = dev.off()
}
### end Figure S5 ####






## Figure S3. ####
tell_user('done.\nCreating Figure S3...')
if ('figure-s3'=='figure-s3') {

  ### Nucleolus gene-list enrichment in toxic signature ####
  ### I think this is the most biologically correct way to do the overrepresentation test
  
  nucleolus_enrichment = bind_rows(
    test_set_enrichment(fazal_nucleolus_genes_mouse, 'Fazal 2019 nucleolus', direction='up'),
    test_set_enrichment(fazal_nucleolus_genes_mouse, 'Fazal 2019 nucleolus', direction='down'),
    test_set_enrichment(engel_nucleolus_genes_mouse, 'Engel 2022 nucleolus', direction='up'),
    test_set_enrichment(engel_nucleolus_genes_mouse, 'Engel 2022 nucleolus', direction='down')
  )
  write_supp_table(nucleolus_enrichment,
                   'Figure S3A: Hypergeometric enrichment of nucleolus gene lists in the toxic transcriptional signature.')

  ### Paraspeckle gene list — placeholder for future analysis ####
  paraspeckle_enrichment = bind_rows(
    test_set_enrichment(nakaya_fus_genes, 'Nakaya 2018 FUS R495X', direction='up'),
    test_set_enrichment(nakaya_fus_genes, 'Nakaya 2018 FUS R495X' , direction='down')
  )

  write_supp_table(paraspeckle_enrichment,
                   'Figure S3A: hypergeometric enrichment of a paraspeckle (FUS) gene list in the toxic transcriptional signature.')

  ### Top 10 GSEA pathways across Hallmark / Reactome / GO:BP (24-0432, toxic vs PBS) ####
  # Rank metric: signed -log10(p), matching prior toxseq GSEA convention.
  gene_ranks_24 = volcano_data %>%
    filter(!is.na(pvalue), !is.na(gene), gene != '', base_mean >= 10) %>%
    group_by(gene) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(rank_metric = -log10(pvalue) * sign(log2fold_change)) %>%
    arrange(desc(rank_metric))
  gene_ranks_24 = setNames(gene_ranks_24$rank_metric, gene_ranks_24$gene)


  # NOTE: This uses human pathways (better annotated, more relevant?) mapped to mouse genes.
  # Change this to db_species = "MM" species = "Mus musculus" if you want to use mouse pathways and mouse genes,
  # or db_species = "MM" and species = "Homo sapiens" if you want to convert mouse pathways to human gene symbols
  msigdb_sets = suppressMessages(msigdbr(species = 'Mus musculus'))
  build_pathways = function(coll, subcoll = NULL) {
    x = msigdb_sets %>% filter(gs_collection == coll)
    if (!is.null(subcoll)) x = x %>% filter(gs_subcollection == subcoll)
    split(x$gene_symbol, x$gs_name)
  }
  run_gsea = function(pathways, label) {
    set.seed(42)
    suppressWarnings(
      fgsea(pathways = pathways, stats = gene_ranks_24,
            minSize = 15, maxSize = 500, nPermSimple = 10000)
    ) %>%
      as_tibble() %>%
      arrange(padj, pval) %>%
      mutate(collection = label)
  }

  bind_rows(
    run_gsea(build_pathways('H'),                 'Hallmark'),
    run_gsea(build_pathways('C2', 'CP:REACTOME'), 'Reactome'),
    run_gsea(build_pathways('C5', 'GO:BP'),       'GO:BP')
  ) %>%
    select(collection, pathway, NES, padj, pval, size) -> all_pathways
  
  all_pathways %>%
    arrange(padj, desc(abs(NES))) %>%
    head(10) -> top10_pathways

  write_supp_table(top10_pathways,
                   'Figure S3B: Top 10 GSEA pathways enriched in the toxic transcriptional signature (24-0432, toxic ASOs vs PBS) ranked by adjusted P value across MSigDB Hallmark, Reactome, and GO:BP collections.')

  ### Plot ####
  col_up   = '#2CA25F'
  col_down = '#8856A7'

  resx = 600
  png('display_items/figure-s3.png', width=6.5*resx, height=7*resx, res=resx)

  layout(matrix(1:2, nrow=2), heights=c(0.4, 1))
  panel = 1

  ### A. Nucleolus + paraspeckle enrichment (up direction only, one point per gene list) ####
  bind_rows(
    nucleolus_enrichment  %>% filter(direction == 'up'),
    paraspeckle_enrichment %>% filter(direction == 'up')
  ) %>%
    arrange(set_name) %>%
    mutate(y = row_number()) -> enrich_plot

  xlims = range(pretty(c(enrich_plot$log2_fold_enrich, 0)))
  ylims = c(0.5, nrow(enrich_plot) + 0.5)
  xats  = pretty(xlims)
  par(mar=c(3, 12, 3, 1))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, bty='n')
  abline(h=enrich_plot$y, col='gray95', lty=1)
  abline(v=0, lty=2, col='gray60')
  points(enrich_plot$log2_fold_enrich, enrich_plot$y,
         pch=ifelse(enrich_plot$pval < 0.05, 19, 21),
         col='#444444', bg='#FFFFFF', cex=1.6, lwd=1.5)
  axis(side=1, at=xats, tck=-0.04, labels=NA)
  axis(side=1, at=xats, labels=xats, lwd=0, line=-0.5, cex.axis=0.7)
  axis(side=2, at=enrich_plot$y, labels=enrich_plot$set_name,
       las=1, lwd=0, lwd.ticks=0, cex.axis=0.7, line=-0.5)
  mtext(side=1, line=1.8, text=expression(log[2]~fold~enrichment), cex=0.7)
  par(xpd=NA)
  legend(x=-3, y=0,
         legend=c('P < 0.05', 'P \u2265 0.05'), pch=c(19, 21),
         col='#444444', pt.bg='#FFFFFF', bty='n', cex=0.7, horiz=TRUE)
  par(xpd=FALSE)
  mtext(LETTERS[panel], side=3, cex=1, adj=0.0, line=0.5); panel = panel + 1

  ### B. Top 10 GSEA pathways ####
  pw = top10_pathways %>%
    mutate(disp_pathway = gsub('_', ' ',
                          gsub('^(HALLMARK_|REACTOME_|GOBP_)', '', pathway)),
           color = if_else(NES >= 0, col_up, col_down)) %>%
    arrange(NES) %>%
    mutate(y = row_number())

  xlims = range(pretty(c(pw$NES, 0)))
  ylims = c(0.5, nrow(pw) + 0.5)
  xats  = pretty(xlims)
  par(mar=c(3, 22, 3, 1))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, bty='n')
  abline(h=pw$y, col='gray95', lty=1)
  abline(v=0, lty=2, col='gray60')
  points(pw$NES, pw$y,
         pch=ifelse(pw$padj < 0.05, 19, 21),
         col=pw$color, bg=pw$color, cex=1.6, lwd=1.5)
  axis(side=1, at=xats, tck=-0.04, labels=NA)
  axis(side=1, at=xats, labels=xats, lwd=0, line=-0.5, cex.axis=0.7)
  axis(side=2, at=pw$y,
       labels=paste0('[', pw$collection, '] ', pw$disp_pathway),
       las=1, lwd=0, lwd.ticks=0, cex.axis=0.55, line=-0.5)
  mtext(side=1, line=1.8, text='NES', cex=0.7)
  par(xpd=NA)
  legend(x=-4, y=0,
         legend=c('padj < 0.05', 'padj \u2265 0.05'), pch=c(19, 21),
         col='#444444', pt.bg='#444444', bty='n', cex=0.7, horiz=TRUE)
  par(xpd=FALSE)
  mtext(LETTERS[panel], side=3, cex=1, adj=0.0, line=0.5); panel = panel + 1

  silence_is_golden = dev.off()
}
### end Figure S3 ####




## Figure S6. ####
tell_user('done.\nCreating Figure S6...')
if ('figure-s6'=='figure-s6') {

  resx = 600
  png('display_items/figure-s6.png', width=9.75*resx, height=3*resx, res=resx)

  layout(matrix(c(1,2,3,
                  4,4,4), nrow=2, byrow=T), heights=c(1,.25))
  par(mar=c(3,4,3,1))
  panel = 1

  sn20_col    = '#7B2D8B'  # purple for Sorce & Nuvolone
  d0432_col   = '#00AAAA'  # cyan for 24-0432
  d2101ac_col = '#3B1C08'  # espresso for 21-01ac

  # --- sn20 data: gene counts per timepoint for Bonferroni correction ---
  sn20_supp2 %>%
    select(GeneID, matches('^Pval_.*_main$')) %>%
    pivot_longer(-GeneID, names_to='col', values_to='val') %>%
    filter(!is.na(val)) %>%
    mutate(wpi = gsub('Pval_|_main', '', col)) %>%
    group_by(wpi) %>%
    summarise(n_genes = n(), .groups='drop') -> sn20_n_tp

  # --- sn20 data: L2FC and significance for all three genes ---
  sn20_supp2 %>%
    filter(GeneID %in% c('Lgals3', 'Lgals9', 'Lgals3bp')) %>%
    select(GeneID, matches('(log2FC|Pval|FDR)_.*_main')) %>%
    pivot_longer(-GeneID, names_to='col', values_to='val') %>%
    mutate(
      metric = case_when(
        grepl('^log2FC', col) ~ 'log2fc',
        grepl('^Pval',   col) ~ 'pval',
        grepl('^FDR',    col) ~ 'fdr'
      ),
      wpi = gsub('log2FC_|Pval_|FDR_|_main', '', col)
    ) %>%
    pivot_wider(id_cols=c(GeneID, wpi), names_from=metric, values_from=val) %>%
    left_join(sn20_n_tp, by='wpi') %>%
    mutate(
      dpi    = if_else(wpi == 'term', 170L, suppressWarnings(as.integer(gsub('wpi','',wpi)) * 7L)),
      padj   = fdr,
      p_bonf = pmin(1, pval * n_genes),
      source = 'sn20'
    ) %>%
    select(gene=GeneID, dpi, log2fc, padj, p_bonf, source) -> lgals_sn20

  # --- 24-0432: RML_PBS vs CBH_PBS at 127 dpi ---
  contrasts_0432 %>%
    filter(group1 == 'RML_PBS', group2 == 'CBH_PBS') %>%
    nrow() -> n_0432_rml_cbh

  contrasts_0432 %>%
    filter(group1 == 'RML_PBS', group2 == 'CBH_PBS',
           gene %in% c('Lgals3', 'Lgals9', 'Lgals3bp')) %>%
    mutate(
      dpi    = 127L,
      log2fc = log2fold_change,
      p_bonf = pmin(1, pvalue * n_0432_rml_cbh),
      source = '24-0432'
    ) %>%
    select(gene, dpi, log2fc, padj, p_bonf, source) -> lgals_0432

  # --- 21-01ac: RML vs CBH NO and PBS at 60, 89, 120, 149 dpi ---
  contrasts %>%
    filter(dataset == '21-01ac') %>%
    group_by(group1, group2) %>%
    summarise(n_contrast = n(), .groups='drop') -> n_2101ac_contrasts

  contrasts %>%
    filter(dataset == '21-01ac',
           group1 %in% c('RML_60_NO','RML_89_NO','RML_120_PBS','RML_149_PBS'),
           group2 %in% c('CBH_60_NO','CBH_89_NO','CBH_120_PBS','CBH_149_PBS'),
           gene == 'Lgals3bp') %>%
    left_join(n_2101ac_contrasts, by=c('group1','group2')) %>%
    mutate(
      dpi = case_when(
        group1 == 'RML_60_NO'   ~ 60L,
        group1 == 'RML_89_NO'   ~ 89L,
        group1 == 'RML_120_PBS' ~ 120L,
        group1 == 'RML_149_PBS' ~ 149L
      ),
      log2fc = log2fold_change,
      p_bonf = pmin(1, pvalue * n_contrast),
      source = '21-01ac'
    ) %>%
    select(gene, dpi, log2fc, padj, p_bonf, source) -> lgals_2101ac

  # --- Combined ---
  bind_rows(lgals_sn20, lgals_0432, lgals_2101ac) %>%
    mutate(
      sig    = padj < 0.05 | p_bonf < 0.05,
      pt_col = case_when(
        source == 'sn20'    ~ sn20_col,
        source == '24-0432' ~ d0432_col,
        source == '21-01ac' ~ d2101ac_col
      ),
      pt_bg  = if_else(sig, pt_col, 'white')
    ) -> lgals_plot

  xlims = c(0, 180)
  xats  = seq(0, 180, 30)
  ylims = range(pretty(lgals_plot$log2fc))

  write_supp_table(lgals_plot %>%
                     select(gene, dpi, source, log2fc, padj, p_bonf),
                   'Figure S6: galectin (Lgals3, Lgals9, Lgals3bp) fold changes across datasets and timepoints')

  ### A, B & C. Lgals3, Lgals9, Lgals3bp ####
  for (g in c('Lgals3', 'Lgals9', 'Lgals3bp')) {
    sub = lgals_plot %>% filter(gene == g) %>% arrange(source, dpi)
    par(mar=c(3,4,3,1))
    plot(NA, xlim=xlims, ylim=ylims, axes=FALSE, ann=FALSE, xaxs='i')
    abline(h=0, lty=3, col='#888888', lwd=0.5)
    for (src in unique(sub$source)) {
      s = sub %>% filter(source == src)
      lines(s$dpi, s$log2fc, col=s$pt_col[1], lwd=1)
    }
    points(sub$dpi, sub$log2fc,
           pch=21, bg=sub$pt_bg, col=sub$pt_col, cex=1.2, lwd=1)
    axis(side=1, at=xats, tck=-0.04, labels=NA)
    axis(side=1, at=xats, labels=xats, lwd=0, line=-0.5, cex.axis=0.7)
    axis(side=2, at=pretty(ylims), tck=-0.04, labels=NA)
    axis(side=2, at=pretty(ylims), labels=pretty(ylims),
         lwd=0, line=-0.5, las=2, cex.axis=0.7)
    mtext(side=1, line=1.8, text='dpi', cex=0.7)
    mtext(side=2, line=2.5, text=expression(log[2]~fold~change), cex=0.7)
    mtext(side=3, line=0.25, text=bquote(italic(.(g))), cex=0.9)
    
    mtext(LETTERS[panel], side=3, cex=1, adj=0.0, line=0.5); panel = panel + 1
  }
  
  par(mar=c(0,0,0,0))
  plot(NA, xlim=0:1, ylim=0:1, axes=FALSE, ann=FALSE, xaxs='i')
  legend(x=0.1,y=0.9, bty='n', cex=0.6, pt.cex=1.1,
         legend=c('Sorce & Nuvolone 2020', 'RML PBS (N=12) vs. CBH PBS (N=6) single timepoint',
                  'RML none/PBS (N=4) vs. CBH none/PBS (N=4) time series'),
         pch   = c(21, 21, 21),
         col   = c(sn20_col, d0432_col, d2101ac_col),
         pt.bg = c(sn20_col, d0432_col, d2101ac_col))
  legend(x=0.4, y=0.9, bty='n', cex=0.6, pt.cex=1.1,
         legend=c('significant', 'not significant'),
         pch   = c( 21, 21),
         col   = c( '#555555', '#555555'),
         pt.bg = c( '#555555', 'white'))
  
  silence_is_golden = dev.off()
}
### end Figure S6 ####





## Figure S4. ####
tell_user('done.\nCreating Figure S4...')
if ('figure-s4'=='figure-s4') {
  
  resx=600
  png('display_items/figure-s4.png',width=6.5*resx,height=4.5*resx,res=resx)

  layout_matrix = matrix(c(1,2,
                           3,3), nrow=2, byrow=T)
  layout(layout_matrix)
  
  panel = 1

  ### A. potency ####
  prp240918_potency %>%
    rename(cord = thoracic_cord) %>%
    pivot_longer(cols=c(cortex, cord, brainstem), names_to='tissue', values_to='value') %>%
    group_by(animal_id, ionis_id) %>%
    summarize(.groups='keep', residual_rna = mean(value)) %>%
    ungroup() %>%
    select(ionis_id, residual_rna) -> potency_240918

  prp250210_qpcr %>%
    rename(animal_id = animal) %>%
    pivot_longer(cols=c(cortex, cord, brainstem), names_to='tissue', values_to='value') %>%
    group_by(animal_id, ionis_id) %>%
    summarize(.groups='keep', residual_rna = mean(value)) %>%
    ungroup() %>%
    select(ionis_id, residual_rna) -> potency_250210

  potency_meta = tibble(ionis_id = c('ASO2-MsPA2','ASO2-MsPA1','ASO2-MOE','ASO2-OMe2','ASO2','PBS'), x = 1:6) %>%
    left_join(treatments %>% select(ionis_id, display_name), by='ionis_id') %>%
    mutate(display_name = coalesce(display_name, ionis_id))

  bind_rows(potency_240918, potency_250210) %>%
    inner_join(potency_meta, by='ionis_id') %>%
    left_join(treatments %>% select(ionis_id, color), by='ionis_id') %>%
    mutate(color = coalesce(color, '#999999')) -> potency

  par(mar=c(3,9,3,1))
  xlims = c(0, 1.5)
  ylims = c(0.5, 6.5)
  xats  = 0:6/4
  xbigs = 0:3/2
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xats,  tck=-0.02, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5)
  mtext(side=1, line=1.8, text='residual Prnp mRNA')
  axis(side=2, at=c(ylims[1], potency_meta$x), labels=c('', potency_meta$display_name), las=2, tck=0, cex.axis=0.7)
  abline(v=1, lty=3)
  abline(v=mean(potency$residual_rna[potency$ionis_id=='ASO2']), lty=2, col='#CCCCCC')
  points(potency$residual_rna, potency$x, col=potency$color, bg='white', pch=21)
  potency %>%
    group_by(ionis_id, x, color) %>%
    summarize(.groups='keep',
              n    = n(),
              mean = mean(residual_rna),
              sd   =   sd(residual_rna),
              l95  = lower(residual_rna),
              u95  = upper(residual_rna)) %>%
    ungroup() -> potency_smry
  barwidth = 0.3
  segments(x0=potency_smry$mean, y0=potency_smry$x - barwidth/2, y1=potency_smry$x + barwidth/2,
           lwd=2, col=potency_smry$color)
  arrows(y0=potency_smry$x, x0=potency_smry$l95, x1=potency_smry$u95,
         code=3, angle=90, length=0.05, lwd=1.5, col=potency_smry$color)
  pbs_rna = potency$residual_rna[potency$ionis_id == 'PBS']
  potency %>%
    filter(ionis_id != 'PBS') %>%
    group_by(ionis_id, x) %>%
    summarize(.groups='keep',
              pval = t.test(residual_rna, pbs_rna)$p.value) %>%
    ungroup() -> potency_p
  #mtext(side=4, at=potency_p$x, line=0.25, cex=0.7, las=2, text=paste0('P',format_p(potency_p$pval)))
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(potency %>% select(display_name, residual_prnp_mrna = residual_rna),
                   'Figure S4A: residual Prnp mRNA per animal (modified ASO 2 potency)')

  ### B. DGE scores (24-0432, ASOs in Fig S4A) ####
  fig3a_asos = potency_meta %>% filter(ionis_id != 'PBS') %>% pull(ionis_id)

  contrasts_0432 %>%
    filter(group1 %in% paste0('RML_', fig3a_asos) & group2 == 'RML_PBS') %>%
    mutate(ionis_id = gsub('^RML_', '', group1)) %>%
    inner_join(g_tox_genes %>% select(gene, sign), by = 'gene') %>%
    group_by(ionis_id) %>%
    summarise(
      gtox_score = sum(log2fold_change * sign, na.rm = TRUE),
      n_gtox     = n(),
      .groups    = 'drop'
    ) %>%
    left_join(treatments %>% select(ionis_id, display_name, color), by = 'ionis_id') %>%
    arrange(gtox_score) -> dge_fig3a

  n_asos_3a  = nrow(dge_fig3a)
  score_xlim = c(0,70)
  score_xats = seq(0,70,10)

  par(mar = c(3, 9, 3, 1))
  ylims = c(0.5, n_asos_3a + 0.5)
  plot(NA, xlim = score_xlim, ylim = ylims,
       xlab = '', ylab = '', axes = FALSE, frame.plot = FALSE, xaxs='i', yaxs='i')
  abline(v = 0, lty = 3, col = 'grey50', lwd = 0.5)
  points(dge_fig3a$gtox_score, seq_len(n_asos_3a),
         col = dge_fig3a$color, pch = 16, cex = 1.2)
  axis(side = 1, at = score_xats, tck = -0.04, labels = NA)
  axis(side = 1, at = score_xats, labels = score_xats,
       lwd = 0, line = -0.5, cex.axis = 0.6)
  axis(side = 2, at = ylims, lwd.ticks=0, labels=NA)
  axis(side = 2, at = seq_len(n_asos_3a), labels = dge_fig3a$display_name,
       las = 2, cex.axis = 0.7, lwd = 0, tck = 0)
  mtext(side = 1, line = 2, text = 'tox score', cex = 0.6)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  par(mar=c(3,4,3,1))

  write_supp_table(dge_fig3a %>% select(display_name, gtox_score, n_gtox),
                   'Figure S4B: tox scores for modified ASO 2 variants')

  ### C. toxmit survival ####
  toxmit_survival %>%
    mutate(dpi = as.integer(as.Date(dod) - as.Date(inoculation_date)),
           event = 1) %>%
    left_join(treatments %>% select(ionis_id, display_name, color), by='ionis_id') %>%
    mutate(color        = coalesce(color, '#999999'),
           display_name = coalesce(display_name, as.character(ionis_id))) -> toxmit_surv_data

  toxmit_surv_meta = toxmit_surv_data %>%
    distinct(ionis_id, disp2=display_name, color) %>%
    arrange(disp2)

  sf = survfit(Surv(dpi, event) ~ ionis_id, data=toxmit_surv_data)
  sf$color = toxmit_surv_meta$color[match(gsub('ionis_id=', '', names(sf$strata)), toxmit_surv_meta$ionis_id)]

  plot_sf(sf, xlims=c(0,260), xbreaks=seq(0,300,50), xminor=seq(0,300,10), legend_pos='bottomleft')
  legend('bottomleft', legend=toxmit_surv_meta$disp2, col=toxmit_surv_meta$color,
         text.col=toxmit_surv_meta$color, lwd=2, bty='n', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1, adj = 0.0, line = 0.5); panel = panel + 1

  write_supp_table(toxmit_surv_data %>% select(animal, display_name, survival_dpi = dpi, event),
                   'Figure S4C: modified ASO 2 variant individual animal survival')

  sd_fullmoe = survdiff(Surv(dpi, event) ~ ionis_id, data=toxmit_surv_data %>%
             filter(ionis_id %in% c('PBS','ASO2-MOE')))
  sd_fullmoe$pvalue
  
  sd_2ome = survdiff(Surv(dpi, event) ~ ionis_id, data=toxmit_surv_data %>% 
                          filter(ionis_id %in% c('PBS','ASO2-OMe2')))
  sd_2ome$pvalue

  silence_is_golden = dev.off()
}
### end Figure S4 ####






# SUPPLEMENT ####

tell_user('done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)

elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))
