################################################################################
################################## Table S4 ####################################
################################################################################

# loading R packages
library(openxlsx) # version 4.2.5.2
library(vegan) # version 2.6-4
library(GUniFrac) # version 1.5

# Load the grouping metadata of soil samples
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Load the information of pathogenic and arbuscular mycorrhizal (AM) and saprophytic fungi 
ASV_tax_information <- read.xlsx("ASV_tax_information.xlsx", colNames = T)

# Data Transformation
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Year <- as.factor(Field_group$Year)
Field_group$Site <- factor(Field_group$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_group$Origin <- factor(Field_group$Origin, levels = c("Native","Exotic"))

# notes: I have completed the above work, so I directly load the completed file
fungi_raw_abun <- read.xlsx("Field_ASVs_raw_data.xlsx", sheet = "raw_ASVs", rowNames = T, colNames = T)
fungi_Flattening <- fungi_raw_abun[,Field_group$Sample_ID]
dim(fungi_Flattening)
#fungi_Flattening <- fungi_Flattening[rowSums(fungi_Flattening) > 0, ]
colSums(fungi_Flattening)

# Z-scored standardized before analyses
Field_group <- Field_group[colnames(fungi_Flattening), ] 
colnames(Field_group)
select_var_scale <- c("Soil_ph", "Wcont","Soil_N","Tave","Prec",
                      "Chol","SLA","LDMC","SRL","FRR","RMF")
Field_group_scale = Field_group
Field_group_scale[select_var_scale] = scale(Field_group_scale[select_var_scale])
Field_group_scale$Group = paste0(Field_group_scale$Site, "|", Field_group_scale$Year) 

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
fungi_relative <- decostand(fungi_Flattening, method = "total", MARGIN = 2)

################### community composition of overall fungi #####################
set.seed(1234)
Overall_fungi_perM <- with(Field_group_scale, GUniFrac::adonis3(t(fungi_relative) ~ Field_SR + Origin * (Tave + Prec + Soil_ph + Wcont + Soil_N), method = "bray",
                                                                data = Field_group_scale, permutations = 999, strata = Group))

Overall_perM_tab <- as.data.frame(Overall_fungi_perM$aov.tab)
Overall_perM_tab$R2 <- round(Overall_perM_tab$R2,3)
Overall_perM_tab$df <- paste0(Overall_perM_tab$Df, ",", Overall_perM_tab$Df[length(Overall_perM_tab$Df)-1])
Overall_perM_tab$`F` <- round(Overall_perM_tab$F.Model,2)
Overall_perM_tab$Predictors <- c("Fungal richness", "Origin", "Temperature", "Precipitation", "Soil pH", "Wcont", "Soil N", 
                                 "Origin × Temperature", "Origin × Precipitation", "Origin × Soil pH", "Origin × Wcont", 
                                 "Origin × Soil N", "Residuals", "Total")
Overall_perM_tab <- Overall_perM_tab[-c(14) ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Overall_perM_tab) <- NULL
print(Overall_perM_tab)
