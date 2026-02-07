################################################################################
################################# Figure S1 ####################################
################################################################################

library(openxlsx) # version 4.2.5.2

# Soil sample grouping information in field
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil sample abundance information in field
Field_otu_raw <- read.xlsx("Field_ASVs_raw_data.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
rownames(Field_otu_raw) <- Field_otu_raw$`#OTU.ID`
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]
Field_otu_raw[1:6, 1:6]
dim(Field_otu_raw)
Field_otu_raw = Field_otu_raw[rowSums(Field_otu_raw) > 0, ]
#View(as.data.frame(rowSums(Field_otu_raw)))

# loading sample grouping information in greenhouse exp.
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
colnames(Green_group)

# loading sample grouping information in greenhouse exp.
Green_otu_raw <- read.xlsx("Greenhouse_ASVs_raw_data.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
rownames(Green_otu_raw) = Green_otu_raw$ASV_ID
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]
dim(Green_otu_raw)
Green_otu_raw = Green_otu_raw[rowSums(Green_otu_raw) > 0, ]
Green_otu_raw[1:6, 1:6]
sum(Green_otu_raw)

# show venn plot
field_ASV = data.frame(ASV = rownames(Field_otu_raw), type = "Field")
green_ASV = data.frame(ASV = rownames(Green_otu_raw), type = "Green")

library(limma) # version 3.50.3
venn_data = rbind(field_ASV,green_ASV)
venn_data$abun = 1
venn_data <- as.data.frame(reshape2 ::acast(venn_data, formula = ASV ~ type , value.var = "abun", fill = 0))
head(venn_data)
colnames(venn_data) <- c( "Field survey", "Greenhouse experiment")
venn_data=venn_data[rowSums(venn_data)>0,]
v_venn_data=vennCounts(venn_data)
vennDiagram(v_venn_data, circle.col = c('#3A648C', '#B79C64'), lwd=4, cex=1, scale=F)

# common ASVs id
common_ASVs <- intersect(rownames(Field_otu_raw), rownames(Green_otu_raw))
length(common_ASVs)


####### 评估共有ASVs对温室实验组成变异贡献
Green_otu_rel <- apply(Green_otu_raw, 2, function(x) x / sum(x))
Green_otu_rel <- as.data.frame(Green_otu_rel)
colSums(Green_otu_rel)

Field_otu_rel <- apply(Field_otu_raw, 2, function(x) x / sum(x))
Field_otu_rel <- as.data.frame(Field_otu_rel)
colSums(Field_otu_rel)


# 分年份站点依次计算共有类群对样本组成的变异贡献
Year = unique(Field_group$Year)
Site = unique(Field_group$Site)
diff_BC_merge_all = NULL

for (i in Year) {
  for (ii in Site) {
    sub_group = subset(Field_group, Year == i & Site == ii)
    
    # Field survey
    Field_otu_rel_sub = Field_otu_rel[,sub_group$Sample_ID]
    #Field_otu_rel_sub = Field_otu_rel_sub[rowSums(Field_otu_rel_sub) > 0, ]
    
    ## Step 1: 计算全体类群组成变异 ----
    # 计算样本间平均绝对差（类似 Bray–Curtis，但无抽平）
    x <- apply(combn(ncol(Field_otu_rel_sub), 2), 2, function(x) sum(abs(Field_otu_rel_sub[,x[1]] - Field_otu_rel_sub[,x[2]]))/2 )
    x_names <- apply(combn(ncol(Field_otu_rel_sub), 2), 2, function(x) paste(colnames(Field_otu_rel_sub)[x], collapse=' - '))
    Field_BC_full <- data.frame(x, x_names)
    
    ## Step 2. 共有类群的样本间差异
    Field_common_sub <- Field_otu_rel_sub[common_ASVs, ]
    x <- apply(combn(ncol(Field_common_sub), 2), 2, function(x) sum(abs(Field_common_sub[, x[1]] - Field_common_sub[, x[2]]))/2)
    x_names <- apply(combn(ncol(Field_common_sub), 2), 2, function(x) paste(colnames(Field_common_sub)[x], collapse=' - '))
    Field_BC_common <- data.frame(x, x_names)
    
    # add site, year information
    diff_BC_Field <- left_join(Field_BC_full, Field_BC_common, by=c('x_names'))
    #head(diff_BC_Field)
    diff_BC_Field$diff_BC <- diff_BC_Field$x.y/diff_BC_Field$x.x
    #mean(diff_BC_Field$diff_BC)
    #dim(diff_BC_Field)
    diff_BC_Field$Years = i; diff_BC_Field$Site = ii; diff_BC_Field$Type = "Field"
    
    # add species information
    split_names <- strsplit(diff_BC_Field$x_names, " - ")
    diff_BC_Field$Sample_ID1 <- sapply(split_names, function(x) x[1])
    diff_BC_Field$Sample_ID2 <- sapply(split_names, function(x) x[2])
    
    # 直接添加物种信息
    diff_BC_Field$Species_1 <- sub_group$Species[match(diff_BC_Field$Sample_ID1, rownames(sub_group))]
    diff_BC_Field$Species_2 <- sub_group$Species[match(diff_BC_Field$Sample_ID2, rownames(sub_group))]
    
    diff_BC_Field$Species_pair = paste0(diff_BC_Field$Species_1, " - ", diff_BC_Field$Species_2)
    colnames(diff_BC_Field)[2] = "Sample_pair"
    
    
    
    # Greenhouse experiment
    
    ## Step 1: 计算全体类群组成变异 ----
    # 计算整体样本间平均绝对差
    Green_otu_rel_sub = Green_mean_matrix[, sub_group$Species]
    
    x <- apply(combn(ncol(Green_otu_rel_sub), 2), 2, function(x) sum(abs(Green_otu_rel_sub[,x[1]] - Green_otu_rel_sub[,x[2]]))/2 )
    x_names <- apply(combn(ncol(Green_otu_rel_sub), 2), 2, function(x) paste(colnames(Green_otu_rel_sub)[x], collapse=' - '))
    Green_BC_full <- data.frame(x, x_names)
    
    
    ## Step 2. 共有类群的样本间差异
    Green_common_matrix <- Green_otu_rel_sub[common_ASVs, ]
    x <- apply(combn(ncol(Green_common_matrix), 2), 2, function(x) sum(abs(Green_common_matrix[, x[1]] - Green_common_matrix[, x[2]]))/2)
    x_names <- apply(combn(ncol(Green_common_matrix), 2), 2, function(x) paste(colnames(Green_common_matrix)[x], collapse=' - '))
    Green_BC_common <- data.frame(x, x_names)
    
    # add site, year information
    diff_BC_Green <- left_join(Green_BC_full, Green_BC_common, by=c('x_names'))
    #head(diff_BC_Green)
    diff_BC_Green$diff_BC <- diff_BC_Green$x.y/diff_BC_Green$x.x
    #mean(diff_BC_Green$diff_BC)
    diff_BC_Green$Years = i; diff_BC_Green$Site = ii; diff_BC_Green$Type = "Greenhouse"
    
    # add species information
    split_names <- strsplit(diff_BC_Green$x_names, " - ")
    diff_BC_Green$Species_1 <- sapply(split_names, function(x) x[1])
    diff_BC_Green$Species_2 <- sapply(split_names, function(x) x[2])
    
    # 直接添加物种信息
    diff_BC_Green$Sample_ID1 <- sub_group$Sample_ID[match(diff_BC_Green$Species_1, sub_group$Species)]
    diff_BC_Green$Sample_ID2 <- sub_group$Sample_ID[match(diff_BC_Green$Species_2, sub_group$Species)]
    
    diff_BC_Green$Sample_pair = paste0(diff_BC_Green$Sample_ID1, " - ", diff_BC_Green$Sample_ID2)
    colnames(diff_BC_Green)[2] = "Species_pair"
    
    # merge dataset
    #colnames(diff_BC_Green) %in% colnames(diff_BC_Field)
    # reorder
    diff_BC_Green = diff_BC_Green[,colnames(diff_BC_Field)]
    
    diff_BC_merge = rbind(diff_BC_Field, diff_BC_Green)
    
    # sloop 
    diff_BC_merge_all = rbind(diff_BC_merge_all, diff_BC_merge)
    
  }
}

site_colors <- c("Guangzhou" = "#F6DD61", "Guilin" = "#94684E", "Changsha" = "#BF5B1D",
                 "Wuhan" = "#3E91B7", "Zhengzhou" = "#0E4879", "Tai'an" = "#CFBD9F")

diff_BC_merge_all$Site <- factor(diff_BC_merge_all$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
diff_BC_merge_all$Type <- factor(diff_BC_merge_all$Type, levels = c("Field","Greenhouse"))
diff_BC_merge_all$Years <- as.factor(diff_BC_merge_all$Years)

diff_BC_merge_all %>% 
  group_by(Type) %>% 
  summarise(
    mean_diff_BC = mean(diff_BC*100),
    sd_diff_BC = sd(diff_BC*100),
    se_diff_BC = sd(diff_BC*100) / sqrt(n()))


site_colors <- c("Guangzhou" = "#0E4879", "Guilin" = "#3E91B7", "Changsha" = "#999999",
                 "Wuhan" = "#CFBD9F", "Zhengzhou" = "#94684E", "Tai'an" = "#BF5B1D")

type_colors <- c("Field" = "#3A648C", "Greenhouse" = "#B79C64")

ggplot(diff_BC_merge_all, aes(x = Site, y = diff_BC*100, fill = Years, color = Years)) +
  geom_boxplot(color = "black", outlier.shape = 21, outlier.size = 0.8, width = 0.6,
               position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = Years, fill = Years), fun.y = mean, geom = "point", shape = 23, size=1.8, #fill = "white",
               color = "black", position = position_dodge(width = 0.7)) + 
  #scale_color_manual(values = c("Field" = "#1C3C63", "Greenhouse" = "#93C8C0")) +
  scale_fill_manual(values = c("2018" = alpha("#0A636D", 0.7), "2020" = alpha("#F8A318", 0.7), "2021" = alpha("#76322A", 0.7))) +
  scale_y_continuous(labels = scales::label_comma(accuracy = 1), limits = c(0,100), 
                     breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.2))) + 
  geom_hline(yintercept = 50, linetype = 2) +
  ggh4x::facet_grid2( Type ~. , #switch = "y", 
                      strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = site_colors),
                                                  background_y = ggh4x::elem_list_rect(fill = type_colors))) +
  theme_classic() +
  theme(legend.position = c(0.9,0.15), 
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
        panel.grid=element_blank(), 
        axis.line.y = element_line(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        axis.title.x = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black', size = 14),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(colour = 'black', size = 12, angle = 35, vjust = 1, hjust = 1),
        strip.background = element_rect(color=NA, size=0.5, linetype="solid"),
        # strip.placement = "outside",
        #panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 14, face = "bold")) +
  labs(y = 'Shared taxa contributions\nto fungal community dissimilarity (%)', x = NULL, tag = "b") +
  geom_segment(aes(x = 0, xend = 0, y = 1, yend = 3), color = "black")


