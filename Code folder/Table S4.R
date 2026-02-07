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
Overall_fungi_perM <- with(Field_group_scale, GUniFrac::adonis3(t(fungi_relative) ~ Origin * (Tave + Prec + Soil_ph + Wcont + Soil_N), method = "bray", by = "margin",
                                                                data = Field_group_scale, permutations = 999, strata = Group))

Overall_perM_tab <- as.data.frame(Overall_fungi_perM$aov.tab)
Overall_perM_tab$R2 <- round(Overall_perM_tab$R2,3)
Overall_perM_tab$df <- paste0(Overall_perM_tab$Df, ",", Overall_perM_tab$Df[length(Overall_perM_tab$Df)-1])
Overall_perM_tab$`F` <- round(Overall_perM_tab$F.Model,2)
Overall_perM_tab$Predictors <- c("Origin", "Temperature", "Precipitation", "Soil pH", "Wcont", "Soil N", 
                                 "Origin × Temperature", "Origin × Precipitation", "Origin × Soil pH", "Origin × Wcont", 
                                 "Origin × Soil N", "Residuals", "Total")
Overall_perM_tab <- Overall_perM_tab[-c(12:13) ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Overall_perM_tab) <- NULL
print(Overall_perM_tab)


################## community composition of pathogenic fungi ###################
Path_field_hel_no <- as.data.frame(t(fungi_Flattening))[Field_group_scale$Sample_ID ,subset(ASV_tax_information, Guilds == "Plant pathogen")$ASV_ID]
Path_field_hel_no <- Path_field_hel_no[, colSums(Path_field_hel_no) > 0]
Path_field_rela <- as.data.frame(decostand(t(Path_field_hel_no), method = "total", MARGIN = 2))
# colSums(Path_field_rela)

set.seed(1234)
Path_fungi_perM <- with(Field_group_scale, GUniFrac::adonis3(t(Path_field_rela) ~ Origin * (Tave + Prec + Soil_ph + Wcont + Soil_N), method = "bray", by = "margin",
                                                             data = Field_group_scale, permutations = 999, strata = Group))

Path_perM_tab <- as.data.frame(Path_fungi_perM$aov.tab)
Path_perM_tab$R2 <- round(Path_perM_tab$R2,3)
Path_perM_tab$df <- paste0(Path_perM_tab$Df, ",", Path_perM_tab$Df[length(Path_perM_tab$Df)-1])
Path_perM_tab$`F` <- round(Path_perM_tab$F.Model,2)
Path_perM_tab$Predictors <- c("Origin", "Temperature", "Precipitation", "Soil pH", "Wcont", "Soil N", 
                                 "Origin × Temperature", "Origin × Precipitation", "Origin × Soil pH", "Origin × Wcont", 
                                 "Origin × Soil N", "Residuals", "Total")
Path_perM_tab <- Path_perM_tab[-c(12:13) ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Path_perM_tab) <- NULL
print(Path_perM_tab)


###################### community composition of AM fungi #######################
AM_field_hel_no <- as.data.frame(t(fungi_Flattening))[Field_group_scale$Sample_ID ,subset(ASV_tax_information, Guilds == "Arbuscular Mycorrhizal")$ASV_ID]
AM_field_hel_no <- AM_field_hel_no[, colSums(AM_field_hel_no) > 0]
dim(AM_field_hel_no)

AM_field_hel_no_add <- AM_field_hel_no
AM_field_hel_no_add[rowSums(AM_field_hel_no_add) == 0, ] <- 1e-6

AM_field_rela <- as.data.frame(decostand(t(AM_field_hel_no_add), method = "total", MARGIN = 2))
# colSums(AM_field_rela)

set.seed(1234)
AM_fungi_perM <- with(Field_group_scale, GUniFrac::adonis3(t(AM_field_rela) ~ Origin * (Tave + Prec + Soil_ph + Wcont + Soil_N), method = "bray", by = "margin",
                                                            data = Field_group_scale, permutations = 999, strata = Group))

AM_perM_tab <- as.data.frame(AM_fungi_perM$aov.tab)
AM_perM_tab$R2 <- round(AM_perM_tab$R2,3)
AM_perM_tab$df <- paste0(AM_perM_tab$Df, ",", AM_perM_tab$Df[length(AM_perM_tab$Df)-1])
AM_perM_tab$`F` <- round(AM_perM_tab$F.Model,2)
AM_perM_tab$Predictors <- c("Origin", "Temperature", "Precipitation", "Soil pH", "Wcont", "Soil N", 
                              "Origin × Temperature", "Origin × Precipitation", "Origin × Soil pH", "Origin × Wcont", 
                              "Origin × Soil N", "Residuals", "Total")
AM_perM_tab <- AM_perM_tab[-c(12:13) ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(AM_perM_tab) <- NULL
print(AM_perM_tab)


################# community composition of saprophytic fungi ###################
Sap_field_hel_no <- as.data.frame(t(fungi_Flattening))[Field_group_scale$Sample_ID ,subset(ASV_tax_information, Guilds == "Saprotroph")$ASV_ID]
Sap_field_hel_no <- Sap_field_hel_no[, colSums(Sap_field_hel_no) > 0]
Sap_field_rela <- as.data.frame(decostand(t(Sap_field_hel_no), method = "total", MARGIN = 2))
# colSums(Sap_field_rela)

set.seed(1234)
Sap_fungi_perM <- with(Field_group_scale, GUniFrac::adonis3(t(Sap_field_rela) ~ Origin * (Tave + Prec + Soil_ph + Wcont + Soil_N), method = "bray", by = "margin",
                                                            data = Field_group_scale, permutations = 999, strata = Group))

Sap_perM_tab <- as.data.frame(Sap_fungi_perM$aov.tab)
Sap_perM_tab$R2 <- round(Sap_perM_tab$R2,3)
Sap_perM_tab$df <- paste0(Sap_perM_tab$Df, ",", Sap_perM_tab$Df[length(Sap_perM_tab$Df)-1])
Sap_perM_tab$`F` <- round(Sap_perM_tab$F.Model,2)
Sap_perM_tab$Predictors <- c("Origin", "Temperature", "Precipitation", "Soil pH", "Wcont", "Soil N", 
                              "Origin × Temperature", "Origin × Precipitation", "Origin × Soil pH", "Origin × Wcont", 
                              "Origin × Soil N", "Residuals", "Total")
Sap_perM_tab <- Sap_perM_tab[-c(12:13) ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Sap_perM_tab) <- NULL
print(Sap_perM_tab)



