.libPaths(c("/usr/local/lib64/R/library","/sc/arion/projects/ad-omics/ricardo/Rlib4"))
library(vcfR)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggpubr)

# FROM SVTK
#LOF "Gene(s) on which the SV is predicted to have a loss-of-function effect."
#DUP_LOF "Gene(s) on which the SV is predicted to have a loss-of-function effect via intragenic exonic duplication."
#COPY_GAIN "Gene(s) on which the SV is predicted to have a copy-gain effect."
#DUP_PARTIAL "Gene(s) which are partially overlapped by an SV's duplication, such that an unaltered copy is preserved."
#MSV_EXON_OVR "Gene(s) on which the multiallelic SV would be predicted to have a LOF, DUP_LOF, COPY_GAIN, or DUP_PARTIAL annotation if the SV were biallelic."
#INTRONIC "Gene(s) where the SV was found to lie entirely within an intron."
#INV_SPAN "Gene(s) which are entirely spanned by an SV's inversion."
#UTR "Gene(s) for which the SV is predicted to disrupt a UTR."
#NEAREST_TSS "Nearest transcription start site to intragenic variants."
#promoter "Nearest promoter."
#INTERGENIC "SV does not overlap coding sequence."
#NONCODING_SPAN "Classes of noncoding elements spanned by SV."
#NONCODING_BREAKPOINT "Classes of noncoding elements disrupted by SV breakpoint."

work_dir = "~/ad-omics/ricardo/MyRepo/snakeSV"

vcf_annot = read.vcfR(paste0(work_dir,"/results_study_case_2/merged_cohort/gt_merged.annot.vcf.gz"))

info = as.data.frame(cbind(as.data.frame(getFIX(vcf_annot)), extract_info_tidy(vcf_annot)))
annot_df = unique(info[,c("ID","LOF","DUP_LOF","COPY_GAIN","INV_SPAN","DUP_PARTIAL","INTRONIC","UTR","INTERGENIC")])
colnames(annot_df) = toupper(colnames(annot_df))
annot_df$LOF = ifelse(is.na(annot_df$LOF),FALSE,TRUE)
annot_df$DUP_LOF = ifelse(is.na(annot_df$DUP_LOF),FALSE,TRUE)
annot_df$COPY_GAIN = ifelse(is.na(annot_df$COPY_GAIN),FALSE,TRUE)
annot_df$INV_SPAN = ifelse(is.na(annot_df$INV_SPAN),FALSE,TRUE)
annot_df$DUP_PARTIAL = ifelse(is.na(annot_df$DUP_PARTIAL),FALSE,TRUE)
annot_df$INTRONIC = ifelse(is.na(annot_df$INTRONIC),FALSE,TRUE)
annot_df$UTR = ifelse(is.na(annot_df$UTR),FALSE,TRUE)

annot_df$INTRONIC[(annot_df$DUP_PARTIAL == T) & (rowSums(annot_df[,-1]) == 1)] <- T
annot_df$INTRONIC[annot_df$INTRONIC == T & annot_df$UTR == T ] <- F
annot_df$INTRONIC[annot_df$INTRONIC == T & annot_df$LOF == T ] <- F
annot_df$INTRONIC[annot_df$INTRONIC == T & annot_df$DUP_LOF == T ] <- F

annot_df$EXONIC = ifelse(annot_df$LOF == T & rowSums(annot_df[,-1]) == 1, T,F)
annot_df$EXONIC = ifelse(annot_df$DUP_LOF == T & rowSums(annot_df[,-1]) == 1, T,annot_df$EXONIC)

rownames(annot_df) = annot_df$ID

#a = annot_df[, c("ID","EXONIC","UTR","INTRONIC","INTERGENIC")]
#all(sort(rowSums(a[-1]))==1)

annot_df_region <- annot_df[, c("ID","EXONIC","UTR","INTRONIC","INTERGENIC")] %>% 
  pivot_longer(-ID, names_to = "REGION") %>%
  filter(value==T) %>% select(-value)

annot_df_consequence <- annot_df[, c("ID","LOF","DUP_LOF","COPY_GAIN","INV_SPAN","DUP_PARTIAL")] %>% 
  pivot_longer(-ID, names_to = "CONSEQUENCE") %>%
  filter(value==T) %>% select(-value)

# Custom annot
annot_df_custom = na.omit(unique(info[,c("ID","NONCODING_BREAKPOINT")]))
colnames(annot_df_custom) = c("ID","CUSTOM")

annot_df_full <- annot_df_region %>% left_join(annot_df_consequence) %>% left_join(annot_df_custom)
annot_df_full$REGION = factor(annot_df_full$REGION, levels = rev(c("EXONIC","UTR","INTRONIC","INTERGENIC")))
annot_df_full$CONSEQUENCE = factor(annot_df_full$CONSEQUENCE, levels = rev(c("LOF","DUP_LOF","COPY_GAIN","INV_SPAN","DUP_PARTIAL")))

pdf(paste0(work_dir,"/results/annot.pdf"), width = 5, height = 1.8)
ggplot(annot_df_full, aes(y=REGION, fill=CONSEQUENCE)) + 
  geom_bar() +
  scale_fill_nejm(na.value="gray") +
  theme_classic() +
  labs(x = "Number of SVs", y = "Region", fill = "Consequence")
dev.off()

png(paste0(work_dir,"/results/annot.png"), width = 5, height = 1.8, units = "in", res = 1200)
ggplot(annot_df_full, aes(y=REGION, fill=CONSEQUENCE)) + 
  geom_bar() +
  scale_fill_nejm(na.value="gray") +
  theme_classic() +
  labs(x = "Number of SVs", y = "Region", fill = "Consequence")
dev.off()

custom_report = info[,c("ID","SVTYPE","SVLEN","NEAREST_TSS","UTR","INTRONIC")] %>% right_join(annot_df_full[!is.na(annot_df_full$CUSTOM), ])
custom_report = custom_report %>% unite("GENE", NEAREST_TSS:INTRONIC, na.rm = TRUE, remove = FALSE) %>%
  select(-NEAREST_TSS, -INTRONIC, -UTR) 

break_custom_annot <- function(custom_report){
  col_fields <- unique(unlist(stringr::str_split(string = custom_report$CUSTOM, pattern = ",")))
  ret <- stringr::str_split(string = custom_report$CUSTOM, pattern = ",") %>%
    lapply(function(x) {
      vals <- unlist(lapply(x, function(z) !is.na(z[1])))
      names(vals) <- unlist(lapply(x, function(z) z[1]))
      unname(vals[col_fields])
    }) %>% 
    unlist
  if(is.null(ret)){
    ret <- matrix(nrow = 0, ncol = length(col_fields), byrow = TRUE)
  } else {
    ret <- matrix(ret, ncol = length(col_fields), byrow = TRUE)
  }
  ret <- as.data.frame(ret, stringsAsFactors = FALSE) %>%
    setNames(col_fields) %>%
    tibble::as_tibble()
  return(cbind(custom_report, ret) %>% select(-CUSTOM))
}

custom_report2 <- break_custom_annot(custom_report)
rownames(custom_report2) = custom_report2$ID

head(custom_report2)

library(ComplexHeatmap)
library(circlize)

custom_only = !is.na(custom_report2[,7:ncol(custom_report2)]) %>% t() %>%
  as.data.frame()
custom_only = custom_only+0
for (i in 1:nrow(custom_only)){
  custom_only[i,custom_only[i,]==1] <- i
}

mat = custom_only
colnames(mat) <- gsub("(.*):(.*)","\\1",colnames(mat))

col_fun = c("white",pal_nejm("default")(nrow(custom_only)))
names(col_fun) = as.character(0:nrow(custom_only))

svtype_col = pal_npg()(2)
names(svtype_col) = c("DEL","INS")
svlen_col = colorRamp2(c(min(custom_report2$SVLEN), max(custom_report2$SVLEN)/2, max(custom_report2$SVLEN)), 
                       rev(RColorBrewer::brewer.pal(n = 3, name ="RdBu"))) 
region_col = pal_npg()(5)[3:5]
names(region_col) = c("INTERGENIC","INTRONIC","UTR")
  
ha = HeatmapAnnotation(CONSEQUENCE = custom_report2$CONSEQUENCE,
                       REGION = custom_report2$REGION, 
                       SVLEN = custom_report2$SVLEN,
                       SVTYPE = custom_report2$SVTYPE,
                       col = list(SVTYPE = svtype_col,
                                  SVLEN = svlen_col,
                                  REGION = region_col),
                       border = T,
                       annotation_name_gp = list("cex" = 0.7))

pdf(paste0(work_dir,"/results/annotcustom.pdf"), width = 12, height = 2.8)
Heatmap(mat, col = col_fun, bottom_annotation = ha, column_names_gp = gpar(fontsize = 10),
        show_row_dend = F, show_column_dend = F, rect_gp = gpar(col = "white", lwd = 1))
dev.off()





###


annot_df_m = reshape2::melt(annot_df, id.vars = "ID")
annot_df_m_region = annot_df_m %>% filter(variable %in% c("INTRONIC","UTR","INTERGENIC","NONCODING_BREAKPOINT"))
colnames(annot_df_m_region) = c("ID","REGION","REGION_VALUE")

annot_df_m_consequence = annot_df_m %>% filter(variable %in% c("LOF","DUP_LOF","COPY_GAIN","DUP_PARTIAL"))
colnames(annot_df_m_consequence) = c("ID","CONSEQUENCE","REGION_CONSEQUENCE")

annot_df_m2 = annot_df_m_region %>% left_join(annot_df_m_consequence)

g1 <- ggplot(annot_df_m_consequence, aes(x = CONSEQUENCE, y = as.numeric(REGION_CONSEQUENCE), fill = CONSEQUENCE)) + 
  stat_summary(fun = sum, geom = "bar", position = "stack") +
  scale_fill_manual(values = pal_nejm("default")(8)[1:4]) +
  theme_classic() +
  labs(y = "Number of SVs")

g2 <- ggplot(annot_df_m_region, aes(x = REGION, y = as.numeric(REGION_VALUE), fill = REGION)) + 
  stat_summary(fun = sum, geom = "bar", position = "stack") +
  scale_fill_manual(values = pal_nejm("default")(8)[5:8]) +
  theme_classic() +
  labs(y = "Number of SVs")



pdf(paste0(work_dir,"/results/annot.pdf"), width = 6, height = 6)
ggarrange(g1,g2, ncol = 1)
dev.off()

png(paste0(work_dir,"/results/annot.png"), width = 12, height = 3, res = 600, units = "in")
ggarrange(g1,g2, ncol = 2, labels = c("a","b"))
dev.off()

info$
