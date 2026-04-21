################################################################################
################################## Table S3 ####################################
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

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
fungi_relative <- decostand(fungi_Flattening, method = "total", MARGIN = 2)
# colSums(fungi_relative)
Bray_dist_rela <- vegdist(t(fungi_relative), method = 'bray')

################### community composition of overall fungi #####################
set.seed(1234)
Overall_fungi_perM <- GUniFrac::adonis3(t(fungi_relative) ~ Field_SR + Year + Site + Origin/Species + Year:Site + 
                                          Origin:Year + Origin:Site + Origin:Year:Site + 
                                          (Origin/Species):Year + (Origin/Species):Site, method = "bray", 
                                        data = Field_group_scale, permutations = 999)

Overall_perM_tab <- as.data.frame(Overall_fungi_perM$aov.tab)
Overall_perM_tab$R2 <- round(Overall_perM_tab$R2,3)
Overall_perM_tab$df <- paste0(Overall_perM_tab$Df, ",", Overall_perM_tab$Df[length(Overall_perM_tab$Df)-1])
Overall_perM_tab$`F` <- round(Overall_perM_tab$F.Model,2)
Overall_perM_tab$Predictors <- c("Fungal richness", "Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                  "Year × Site × Origin", "Year × Species", "Site × Species","Residuals","Total")
Overall_perM_tab <- Overall_perM_tab[c(1:8,10,11,9,12), c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Overall_perM_tab) <- NULL
print(Overall_perM_tab)


##################### community richness of overall fungi ######################
Overall_field_SR_mod <- lm(Field_SR ~ Year + Site + Origin/Species + Year:Site + 
                             Origin:Year + Origin:Site + Origin:Year:Site + 
                             (Origin/Species):Year + (Origin/Species):Site, data = Field_group)
shapiro.test(residuals(Overall_field_SR_mod))

Overall_SR_mod_anova <- as.data.frame(anova(Overall_field_SR_mod))
Overall_SR_mod_anova[, 2:4] <- round(Overall_SR_mod_anova[, 2:4], digits = 2)
Overall_SR_mod_anova$`Pr(>F)` <- round(Overall_SR_mod_anova$`Pr(>F)`, 3)
Overall_SR_mod_anova$df <- paste0(Overall_SR_mod_anova$Df, ",", Overall_SR_mod_anova$Df[length(Overall_SR_mod_anova$Df)])
Overall_SR_mod_anova$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                  "Year × Site × Origin", "Year × Species", "Site × Species","Residuals")
Overall_SR_mod_anova <- Overall_SR_mod_anova[c(1:7,9,10,8,11), c(7,6,2:5)] # Reorder
rownames(Overall_SR_mod_anova) <- NULL
print(Overall_SR_mod_anova) 
