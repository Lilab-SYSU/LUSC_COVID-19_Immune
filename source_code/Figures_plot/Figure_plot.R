library(plotrix)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(VennDiagram)
library(grDevices)


# Figure 1F
davideMatrix <- read_csv(file = 'DAVIDE_pathway.csv')
davideMatrix2 <- davideMatrix %>% 
  filter(Benjamini<0.05) %>% 
  mutate(`-log(P.adjust)` = -log10(Benjamini)) %>% 
  mutate(class = factor(class), detail = factor(detail)) %>% 
  mutate(detail = str_to_upper(detail))
p <- ggdotchart(davideMatrix2,
           y="score",
           x="detail", 
           color = "-log(P.adjust)",
           group = "class",
           dot.size = 10,
           rotate = T,
           add = 'segments',
           add.params = list(color = "lightgray", size = 2),
           label = "score",
           font.label = list(color = "white", size = 9, vjust = 0.5)
           #y.text.col = T
           ) +
  scale_color_continuous(high='#8B008B', low='#FF00FF') +
  theme(panel.border = element_rect(colour="black", fill=NA))
  #theme_bw(base_family="Arial")
ggsave(filename = '棒棒糖图.pdf', plot = p, width = 10, height = 10)


# Figure 2D
KEGG_pathway_gsea <- read_csv(file = 'KEGG_pathway_gsea_report_for_Exh_1597073331292.csv')
KEGG_pathway_gsea2 <- KEGG_pathway_gsea %>% 
  filter(`FDRq-val`<0.01) %>% 
  mutate(`-log(P.adjust)` = -log10(`FDRq-val`+0.0001)) %>% 
  mutate(xAxis = 'Immune Exhaustion vs. Rest classes') %>% 
  mutate(NAME = str_replace(NAME, pattern = 'KEGG_', replacement = '')) %>% 
  arrange(NES) %>% 
  mutate(NAME=factor(.$NAME, levels = c("PANCREATIC_CANCER", "MELANOMA", "GLIOMA", "CALCIUM_SIGNALING_PATHWAY", "ENDOCYTOSIS", "GLYCOSPHINGOLIPID_BIOSYNTHESIS_GANGLIO_SERIES", "TGF_BETA_SIGNALING_PATHWAY", "ACUTE_MYELOID_LEUKEMIA", "DORSO_VENTRAL_AXIS_FORMATION", "NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION", "FC_EPSILON_RI_SIGNALING_PATHWAY", "RENIN_ANGIOTENSIN_SYSTEM", "AXON_GUIDANCE", "GRAFT_VERSUS_HOST_DISEASE", "PRIMARY_IMMUNODEFICIENCY", "ASTHMA", "APOPTOSIS", "ALLOGRAFT_REJECTION", "PATHWAYS_IN_CANCER", "LYSOSOME", "COMPLEMENT_AND_COAGULATION_CASCADES", "MAPK_SIGNALING_PATHWAY", "TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS", "B_CELL_RECEPTOR_SIGNALING_PATHWAY", "ANTIGEN_PROCESSING_AND_PRESENTATION", "TYPE_I_DIABETES_MELLITUS", "ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC", "NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY", "REGULATION_OF_ACTIN_CYTOSKELETON", "VASCULAR_SMOOTH_MUSCLE_CONTRACTION", "INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION", "SMALL_CELL_LUNG_CANCER", "T_CELL_RECEPTOR_SIGNALING_PATHWAY", "NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "ECM_RECEPTOR_INTERACTION", "DILATED_CARDIOMYOPATHY", "AUTOIMMUNE_THYROID_DISEASE", "HEMATOPOIETIC_CELL_LINEAGE", "HYPERTROPHIC_CARDIOMYOPATHY_HCM", "LEISHMANIA_INFECTION", "VIRAL_MYOCARDITIS", "FOCAL_ADHESION", "CHEMOKINE_SIGNALING_PATHWAY", "LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION", "JAK_STAT_SIGNALING_PATHWAY", "CELL_ADHESION_MOLECULES_CAMS", "CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"), ordered=T))

p22 <- ggplot(KEGG_pathway_gsea2, aes(x=xAxis, y=NAME, colour=`-log(P.adjust)`, size=NES)) +
  geom_point() +
  #scale_size_area(max_size = 6) +
  scale_radius(range = c(3, 8.5)) +
  scale_color_continuous(high='#EE0000', low='#1E90FF') +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black", fill=NA, size = 1.3))
ggsave(filename = 'KEGG_sizeNES_colorP.pdf', plot = p22, width = 13, height = 15)


# Figure 3B
LUSC1_diffgenes <- as_tibble(read.table(file = 'LUSC1_diffgenes_filter_lgFC1_pvalue5_UP.xls', header = T, sep="\t"))
SARS_CoV_2_diffgenes <- as_tibble(read.table(file = 'SARS_CoV_2_diffgenes_filter_lgFC1_pvalue5_UP.xls', header = T, sep="\t"))

venn.plot <- venn.diagram(
  x = list(
    Exhausted_LUSC_UP = LUSC1_diffgenes$gene_name,
    `SARS-CoV-2_UP` = SARS_CoV_2_diffgenes$gene_name
  ),
  filename = NULL, 
  # imagetype = "pdf",
  lwd = 3,
  fill = c("darkorchid1", "cornflowerblue"),
  alpha = 0.6,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",  
  fontface = "bold",
  cat.col = c("darkorchid1", "cornflowerblue"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.06, 0.05),
  cat.pos = c(0, 20),
  col = 'white'
)


pdf(file = 'Venn2.pdf')
grid.draw(venn.plot)
dev.off()


# Figure 3E
exhauGSEA <- as_tibble(read.table("gsea_report_for_Exhausted_1598350006543.xls", header = T, sep="\t"))
sarsCov2GSEA <- as_tibble(read.table('gsea_report_for_SARS_CoV2_1598858151469.xls', header = T, sep="\t"))

exhauGSEA2 <- exhauGSEA %>%
  filter(NOM.p.val<0.01)
sarsCov2GSEA2 <- sarsCov2GSEA %>%
  filter(NOM.p.val<0.01)

venn.plot <- venn.diagram(
  x = list(
    Exhausted = exhauGSEA2$NAME,
    `SARS-CoV-2` = sarsCov2GSEA2$NAME
  ),
  filename = NULL, 
  #imagetype = "tiff",
  lwd = 3,
  fill = c("darkorchid1", "cornflowerblue"),
  alpha = 0.6,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",  
  fontface = "bold",
  cat.col = c("darkorchid1", "cornflowerblue"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.06, 0.05),
  cat.pos = c(-20, 20),
  col = 'white'
)

pdf(file = 'Venn1.pdf')
grid.draw(venn.plot)
dev.off()


# Figure 4C
meta <- read_csv(file = 'All_cells_meta_info.csv')

sum_each_sample_meta <- meta %>%
  group_by(orig.ident) %>% 
  mutate(sumOfEachSample = n()) %>% 
  select(orig.ident, sumOfEachSample) %>% 
  distinct()

prop_sample_group_celltype <- meta %>%
  group_by(group, orig.ident, Celltype) %>% 
  tally() %>% 
  left_join(sum_each_sample_meta, by = 'orig.ident') %>% 
  rowwise() %>% 
  mutate(prop = n/sumOfEachSample)

cellList <- prop_sample_group_celltype$Celltype %>% unique()

#for (cell in c("Dendritic cells","Epithelial cells","Exhausted T cells","Macrophages","Monocytes","Neutrophils","NK cells","Plasmacytoid dendritic cells","T cells","Plasma cells","B cells")){
for (cell in cellList){
  prop_celltype <- prop_sample_group_celltype %>% 
    filter(Celltype==cell) %>% 
    mutate(c = case_when(group == 'healthy control' ~ '#E64B35FF',
                         group == 'mild' ~ '#00A1D5FF',
                         group == 'severe' ~ '#00A087FF'))
  maxProp <- prop_celltype %>%
    pull(var=prop) %>%
    max()
  p <- ggboxplot(prop_celltype,
                 x='group',
                 y='prop',
                 color = 'black',
                 #palette = 'npg',
                 add = "jitter",
                 add.params = list(color=prop_celltype$c, size=2.5)) +
    stat_compare_means(comparisons = list(c("healthy control", "mild"), c("healthy control", "severe"), c("mild", "severe"))) +
    stat_compare_means(label.y = maxProp+maxProp/2) +
    theme(panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA)) 
  ggsave(filename = paste0(cell, '.pdf'), plot = p, width = 4, height = 4)
}


# Figure 4D
Macro_venn_dat <- readRDS('Macro_venn_dat.RDS')
MacroT_venn_dat <- readRDS('MacroT_venn_dat.RDS')

venn.plot <- venn.diagram(
  x = list(
    CD3D = Macro_venn_dat$CD3D,
    CD163 = Macro_venn_dat$CD163,
    LAG3 = Macro_venn_dat$LAG3
  ),
  filename = NULL, 
  #imagetype = "pdf",
  lwd = 3,
  fill = c("#E64B35FF", "#00A087FF", '#4DBBD5FF'),
  alpha = 0.8,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",  
  fontface = "bold",
  cat.col = c("#E64B35FF", "#00A087FF", '#4DBBD5FF'),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.06, 0.05, 0.05),
  cat.pos = c(-20, 20, 180)
)

pdf(file = 'Macro_venn_dat.pdf')
grid.draw(venn.plot)
dev.off()

venn.plot <- venn.diagram(
  x = list(
    CD3D = MacroT_venn_dat$CD3D,
    CD163 = MacroT_venn_dat$CD163,
    LAG3 = MacroT_venn_dat$LAG3
  ),
  filename = NULL, 
  #imagetype = "tiff",
  lwd = 3,
  fill = c("#E64B35FF", "#00A087FF", '#4DBBD5FF'),
  alpha = 0.8,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",  
  fontface = "bold",
  cat.col = c("#E64B35FF", "#00A087FF", '#4DBBD5FF'),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.06, 0.05, 0.05),
  cat.pos = c(-20, 20, 180)
)

pdf(file = 'MacroT_venn_dat.pdf')
grid.draw(venn.plot)
dev.off()


# Figure 4F
yOrdered <- read_csv(file = 'interaction.csv')

ExhTCD82macro %>% 
  select(interacting_pair, secreted) %>% 
  write_csv(file = 'ExhTCD82macro.csv')

ExhTCD82macro_pAndMean <- ExhTCD82macro %>% 
  pivot_longer(cols = `B|Macro`:`pDC|Macro`, names_to='cell', values_to='p') %>% 
  left_join(means_longer, by=c("interacting_pair"="interacting_pair", "cell"="cell")) %>% 
  mutate(p2 = p+0.0001)

rowSorted <- ExhTCD82macro_pAndMean %>% 
  pull(var = cell) %>% 
  unique() %>% 
  sort()

p_8 <- ggplot(ExhTCD82macro_pAndMean, aes(x=cell, y=interacting_pair)) +
  geom_point(aes(size=-log10(p2), color=log2(m))) +
  scale_color_gradient('Log2(Mean)', low = '#63B8FF', high = '#FF4040') +
  theme_bw() +
  scale_y_discrete(limits=rev(yOrdered$yOrder)) +
  scale_x_discrete(limits=c(rowSorted[5:17], rowSorted[1:4], rowSorted[18:25])) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
ggsave(filename='气泡图ExhTCD82macro_yOrdered.pdf', plot = p_8, width = 8.3, height = 10.5)


# Figure 4G
pval <- read_tsv(file = 'pvalues.txt')

epi2ExhTCD8 <- pval %>%
  select(interacting_pair, secreted, contains('Epi')) %>% 
  filter(!str_detect(interacting_pair, 'complex')) #%>% 
  #filter(`Epi|Exh T_CD8` < 0.05)

secreted_epi2ExhTCD8 <- epi2ExhTCD8 %>%
  filter(secreted == T) %>% 
  pull(var = interacting_pair)
nonSecreted_epi2ExhTCD8 <- epi2ExhTCD8 %>%
  filter(secreted == F) %>% 
  pull(var = interacting_pair)

epi2ExhTCD8 %>%
  select(interacting_pair, secreted) %>% 
  write_csv(file = 'epi2ExhTCD8.csv')

epi2ExhTCD8_pAndMean <- epi2ExhTCD8 %>%
  pivot_longer(cols = `B|Epi`:`pDC|Epi`, names_to='cell', values_to='p') %>% 
  left_join(means_longer, by=c("interacting_pair"="interacting_pair", "cell"="cell")) %>% 
  mutate(p2 = p+0.0001)

rowSorted <- epi2ExhTCD8_pAndMean %>% 
  pull(var = cell) %>% 
  unique() %>% 
  sort()

p_2 <- ggplot(epi2ExhTCD8_pAndMean, aes(x=cell, y=interacting_pair)) +
  geom_point(aes(size=-log10(p2), color=log2(m))) +
  scale_color_gradient('Log2(Mean)', low = '#63B8FF', high = '#FF4040') +
  theme_bw() +
  scale_y_discrete(limits=c(secreted_epi2ExhTCD8, nonSecreted_epi2ExhTCD8)) +
  scale_x_discrete(limits=c(rowSorted[3:16], rowSorted[1:2], rowSorted[17:27])) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
ggsave(filename='气泡图epi2ExhTCD8.pdf', plot = p_2, width = 9.2, height = 10)


# Figure 5B
BALF_Cell_types_annot_metadata_no_scale <- read_csv(file = 'BALF_Cell_types_annot_metadata_no_scale.csv')

p <- ggboxplot(BALF_Cell_types_annot_metadata_no_scale,
               x='new_celltype',
               y='cytokine_storms_no_scale1',
               color = 'black',
               outlier.shape = NA,
               fill = 'group',
               palette = c('#98DF8AFF', '#F7B6D2FF', '#D62728FF'),
               #add = "jitter",
               #alpha = 0.9,
               bxp.errorbar = T
               #add.params = list(fill=BALF_Cell_types_annot_metadata_no_scale$c)
               ) +
  #stat_compare_means(comparisons = list(c('Health', 'Moderate'), c('Health', 'Severe'), c('Moderate', 'Severe'))) +
  #stat_compare_means(label.y = 0.01) +
  theme(panel.border = element_rect(linetype = 'solid', size = 1.2, fill = NA)) 
ggsave(filename = 'BALF_Cell_types_annot_metadata_no_scale_1.pdf', plot = p, width = 20, height = 7)

p2 <- ggboxplot(BALF_Cell_types_annot_metadata_no_scale,
               x='new_celltype',
               y='exhaustion_no_scale1',
               color = 'black',
               outlier.shape = NA,
               fill = 'group',
               palette = c('#98DF8AFF', '#F7B6D2FF', '#D62728FF'),
               bxp.errorbar = T
               #outlier.shape = NA,
               #palette = 'npg',
               #add = "jitter",
               #alpha = 0.5,
               #shape = 'group'
               #add.params = list(color=prop_celltype$c, size=2.5)
) +
  #stat_compare_means(comparisons = list(c('Health', 'Moderate'), c('Health', 'Severe'), c('Moderate', 'Severe'))) +
  #stat_compare_means(label.y = 0.01) +
  theme(panel.border = element_rect(linetype = 'solid', size = 1.2, fill = NA)) 
ggsave(filename = 'BALF_Cell_types_annot_metadata_no_scale_2.pdf', plot = p2, width = 20, height = 7)


# Figure 5D,5E
macroMeta <- read_csv(file = 'Macro_meta_data.csv')

macroMeta2 <- macroMeta %>% 
  mutate(cell2 = paste0('g', cell)) %>% 
  select(orig.ident, group, cell2)

sumEachCell <- macroMeta2 %>%
  group_by(orig.ident) %>% 
  mutate(sumOfEachCell = n()) %>% 
  select(orig.ident, sumOfEachCell) %>% 
  distinct()

propCellType <- macroMeta2 %>%
  group_by(group, orig.ident, cell2) %>% 
  tally() %>% 
  left_join(sumEachCell, by = 'orig.ident') %>%
  rowwise() %>%
  mutate(prop = n/sumOfEachCell)

cellList <- propCellType$cell2 %>% 
  unique()

for (cell in cellList){
  propCellType2 <- propCellType %>% 
    filter(cell2==cell) %>% 
    mutate(c = case_when(group == 'Healthy' ~ '#E64B35FF',
                         group == 'Moderate' ~ '#00A1D5FF',
                         group == 'Severe' ~ '#00A087FF'))
  maxProp <- propCellType2 %>%
    pull(var=prop) %>%
    max()
  p <- ggboxplot(propCellType2,
                 x='group',
                 y='prop',
                 color = 'black',
                 #palette = 'npg',
                 add = "jitter",
                 add.params = list(color=propCellType2$c, size=2.5)) +
    stat_compare_means(comparisons = list(c("Healthy", "Moderate"), c("Healthy", "Severe"), c("Moderate", "Severe"))) +
    stat_compare_means(label.y = maxProp+maxProp/2) +
    theme(panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA)) 
  ggsave(filename = paste0(cell, '.pdf'), plot = p, width = 4, height = 4)
}


# Figure 6B
meta <- read_csv(file = 'PBMC_meta.csv')

sum_each_sample_meta <- meta %>%
  group_by(orig.ident) %>% 
  mutate(sumOfEachSample = n()) %>% 
  select(orig.ident, sumOfEachSample) %>% 
  distinct()

prop_sample_group_celltype <- meta %>% 
  group_by(group, orig.ident, Celltype) %>% 
  tally() %>% 
  left_join(sum_each_sample_meta, by = 'orig.ident') %>% 
  rowwise() %>% 
  mutate(prop = n/sumOfEachSample)
  
for (cell in c("B","ELEP","Macro_c1","Macro_c2","Macro_c3","Macro_c4","Macro_c5","Macro_c6","Mono","Neu","NK","pDC","Plasma","Platelets","Progenitors","T","T_CD4","T_CD8")){
  prop_celltype <- prop_sample_group_celltype %>% 
    filter(Celltype==cell) %>% 
    mutate(c = case_when(group == 'Health' ~ '#E64B35FF',
                         group == 'Moderate' ~ '#00A1D5FF',
                         group == 'Severe' ~ '#00A087FF'))
  # maxProp <- prop_celltype %>% 
  #   pull(var=prop) %>% 
  #   max()
  p <- ggboxplot(prop_celltype,
                 x='group',
                 y='prop',
                 color = 'black',
                 #palette = 'npg',
                 add = "jitter",
                 add.params = list(color=prop_celltype$c, size=2.5)) +
    stat_compare_means(comparisons = list(c('Health', 'Moderate'), c('Health', 'Severe'), c('Moderate', 'Severe'))) +
    stat_compare_means(label.y = 0.01) +
    theme(panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA)) 
  ggsave(filename = paste0(cell, '.pdf'), plot = p, width = 4, height = 4)
}
