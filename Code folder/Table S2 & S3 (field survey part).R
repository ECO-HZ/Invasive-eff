################################################################################
###################### Table S2 & S3 field survey part #########################
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
Overall_fungi_perM <- GUniFrac::adonis3(t(fungi_relative) ~ Year + Site + Origin/Species + Year:Site + 
                                          Origin:Year + Origin:Site + Origin:Year:Site + 
                                          (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                        data = Field_group_scale, permutations = 999)

Overall_perM_tab <- as.data.frame(Overall_fungi_perM$aov.tab)
Overall_perM_tab$R2 <- round(Overall_perM_tab$R2,3)
Overall_perM_tab$df <- paste0(Overall_perM_tab$Df, ",", Overall_perM_tab$Df[length(Overall_perM_tab$Df)-1])
Overall_perM_tab$`F` <- round(Overall_perM_tab$F.Model,2)
Overall_perM_tab$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                  "Year × Site × Origin", "Year × Species", "Site × Species","Residuals","Total")
Overall_perM_tab <- Overall_perM_tab[c(1:7,9,10,8), c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Overall_perM_tab) <- NULL
print(Overall_perM_tab)


##################### community richness of overall fungi ######################
field_overall_SR <- as.data.frame(specnumber(t(fungi_relative)))
colnames(field_overall_SR) <- "Overall_SR"
field_overall_SR$Sample_ID <- rownames(field_overall_SR)
Field_group <- Field_group %>% left_join(field_overall_SR)

Overall_field_SR_mod <- lm(Overall_SR ~ Year + Site + Origin/Species + Year:Site + 
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

################## community composition of pathogenic fungi ###################
Path_field_hel_no <- as.data.frame(t(fungi_Flattening))[ ,subset(ASV_tax_information, Guilds == "Plant pathogen")$ASV_ID]
Path_field_hel_no <- Path_field_hel_no[, colSums(Path_field_hel_no) > 0]
Path_field_rela <- as.data.frame(decostand(t(Path_field_hel_no), method = "total", MARGIN = 2))
# colSums(Path_field_rela)

set.seed(1234)
Path_fungi_perM <- GUniFrac::adonis3(t(Path_field_rela) ~ Year + Site + Origin/Species + Year:Site + 
                                       Origin:Year + Origin:Site + Origin:Year:Site + 
                                       (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                     data = Field_group_scale, permutations = 999)

Path_perM_tab <- as.data.frame(Path_fungi_perM$aov.tab)
Path_perM_tab$R2 <- round(Path_perM_tab$R2,3)
Path_perM_tab$df <- paste0(Path_perM_tab$Df, ",", Path_perM_tab$Df[length(Path_perM_tab$Df)-1])
Path_perM_tab$`F` <- round(Path_perM_tab$F.Model,2)
Path_perM_tab$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                 "Year × Site × Origin", "Year × Species", "Site × Species","Residuals","Total")
Path_perM_tab <- Path_perM_tab[c(1:7,9,10,8), c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Path_perM_tab) <- NULL
print(Path_perM_tab)

################### community richness of pathogenic fungi #####################
field_path_SR <- as.data.frame(specnumber(t(Path_field_rela)))
colnames(field_path_SR) <- "Path_SR"
field_path_SR$Sample_ID <- rownames(field_path_SR)
Field_group <- Field_group %>% left_join(field_path_SR)
rownames(Field_group)  <- Field_group$Sample_ID

Path_field_SR_mod <- lm(Path_SR ~ Year + Site + Origin/Species + Year:Site + 
                          Origin:Year + Origin:Site + Origin:Year:Site + 
                          (Origin/Species):Year + (Origin/Species):Site, data = Field_group)
shapiro.test(residuals(Path_field_SR_mod))

Path_SR_mod_anova <- as.data.frame(anova(Path_field_SR_mod))
Path_SR_mod_anova[, 2:4] <- round(Path_SR_mod_anova[, 2:4], digits = 2)
Path_SR_mod_anova$`Pr(>F)` <- round(Path_SR_mod_anova$`Pr(>F)`, 3)
Path_SR_mod_anova$df <- paste0(Path_SR_mod_anova$Df, ",", Path_SR_mod_anova$Df[length(Path_SR_mod_anova$Df)])
Path_SR_mod_anova$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                     "Year × Site × Origin", "Year × Species", "Site × Species","Residuals")
Path_SR_mod_anova <- Path_SR_mod_anova[c(1:7,9,10,8,11), c(7,6,2:5)] # Reorder
rownames(Path_SR_mod_anova) <- NULL
print(Path_SR_mod_anova) 


###################### community composition of AM fungi #######################
AM_field_hel_no <- as.data.frame(t(fungi_Flattening))[ ,subset(ASV_tax_information, Guilds == "Arbuscular Mycorrhizal")$ASV_ID]
AM_field_hel_no <- AM_field_hel_no[, colSums(AM_field_hel_no) > 0]
dim(AM_field_hel_no)

AM_field_hel_no_add <- AM_field_hel_no
AM_field_hel_no_add[rowSums(AM_field_hel_no_add) == 0, ] <- 1e-6

AM_field_rela <- as.data.frame(decostand(t(AM_field_hel_no_add), method = "total", MARGIN = 2))
# colSums(AM_field_rela)

set.seed(1234)
AM_fungi_perM <- GUniFrac::adonis3(t(AM_field_rela) ~ Year + Site + Origin/Species + Year:Site + 
                                       Origin:Year + Origin:Site + Origin:Year:Site + 
                                       (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                     data = Field_group_scale, permutations = 999)

AM_perM_tab <- as.data.frame(AM_fungi_perM$aov.tab)
AM_perM_tab$R2 <- round(AM_perM_tab$R2,3)
AM_perM_tab$df <- paste0(AM_perM_tab$Df, ",", AM_perM_tab$Df[length(AM_perM_tab$Df)-1])
AM_perM_tab$`F` <- round(AM_perM_tab$F.Model,2)
AM_perM_tab$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                              "Year × Site × Origin", "Year × Species", "Site × Species","Residuals","Total")
AM_perM_tab <- AM_perM_tab[c(1:7,9,10,8), c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(AM_perM_tab) <- NULL
print(AM_perM_tab)

######################## community richness of AM fungi ########################
field_AM_SR <- as.data.frame(specnumber((AM_field_hel_no)))
colnames(field_AM_SR) <- "AM_SR"
field_AM_SR$Sample_ID <- rownames(field_AM_SR)
Field_group <- Field_group %>% left_join(field_AM_SR)
rownames(Field_group)  <- Field_group$Sample_ID

AM_field_SR_mod <- lm(AM_SR ~ Year + Site + Origin/Species + Year:Site + 
                          Origin:Year + Origin:Site + Origin:Year:Site + 
                          (Origin/Species):Year + (Origin/Species):Site, data = Field_group)
shapiro.test(residuals(AM_field_SR_mod))

AM_SR_mod_anova <- as.data.frame(anova(AM_field_SR_mod))
AM_SR_mod_anova[, 2:4] <- round(AM_SR_mod_anova[, 2:4], digits = 2)
AM_SR_mod_anova$`Pr(>F)` <- round(AM_SR_mod_anova$`Pr(>F)`, 3)
AM_SR_mod_anova$df <- paste0(AM_SR_mod_anova$Df, ",", AM_SR_mod_anova$Df[length(AM_SR_mod_anova$Df)])
AM_SR_mod_anova$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                  "Year × Site × Origin", "Year × Species", "Site × Species","Residuals")
AM_SR_mod_anova <- AM_SR_mod_anova[c(1:7,9,10,8,11), c(7,6,2:5)] # Reorder
rownames(AM_SR_mod_anova) <- NULL
print(AM_SR_mod_anova) 


################# community composition of saprophytic fungi ###################
unique(ASV_tax_information$Guilds)
Sap_field_hel_no <- as.data.frame(t(fungi_Flattening))[ ,subset(ASV_tax_information, Guilds == "Saprotroph")$ASV_ID]
Sap_field_hel_no <- Sap_field_hel_no[, colSums(Sap_field_hel_no) > 0]
Sap_field_rela <- as.data.frame(decostand(t(Sap_field_hel_no), method = "total", MARGIN = 2))
# colSums(Sap_field_rela)

set.seed(1234)
Sap_fungi_perM <- GUniFrac::adonis3(t(Sap_field_rela) ~ Year + Site + Origin/Species + Year:Site + 
                                       Origin:Year + Origin:Site + Origin:Year:Site + 
                                       (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                     data = Field_group_scale, permutations = 999)

Sap_perM_tab <- as.data.frame(Sap_fungi_perM$aov.tab)
Sap_perM_tab$R2 <- round(Sap_perM_tab$R2,3)
Sap_perM_tab$df <- paste0(Sap_perM_tab$Df, ",", Sap_perM_tab$Df[length(Sap_perM_tab$Df)-1])
Sap_perM_tab$`F` <- round(Sap_perM_tab$F.Model,2)
Sap_perM_tab$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                              "Year × Site × Origin", "Year × Species", "Site × Species","Residuals","Total")
Sap_perM_tab <- Sap_perM_tab[c(1:7,9,10,8), c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Sap_perM_tab) <- NULL
print(Sap_perM_tab)

################## community richness of saprophytic fungi #####################
field_Sap_SR <- as.data.frame(specnumber(t(Sap_field_rela)))
colnames(field_Sap_SR) <- "Sap_SR"
field_Sap_SR$Sample_ID <- rownames(field_Sap_SR)
Field_group <- Field_group %>% left_join(field_Sap_SR)
rownames(Field_group)  <- Field_group$Sample_ID

Sap_field_SR_mod <- lm(Sap_SR ~ Year + Site + Origin/Species + Year:Site + 
                          Origin:Year + Origin:Site + Origin:Year:Site + 
                          (Origin/Species):Year + (Origin/Species):Site, data = Field_group)
shapiro.test(residuals(Sap_field_SR_mod))

Sap_SR_mod_anova <- as.data.frame(anova(Sap_field_SR_mod))
Sap_SR_mod_anova[, 2:4] <- round(Sap_SR_mod_anova[, 2:4], digits = 2)
Sap_SR_mod_anova$`Pr(>F)` <- round(Sap_SR_mod_anova$`Pr(>F)`, 3)
Sap_SR_mod_anova$df <- paste0(Sap_SR_mod_anova$Df, ",", Sap_SR_mod_anova$Df[length(Sap_SR_mod_anova$Df)])
Sap_SR_mod_anova$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                  "Year × Site × Origin", "Year × Species", "Site × Species","Residuals")
Sap_SR_mod_anova <- Sap_SR_mod_anova[c(1:7,9,10,8,11), c(7,6,2:5)] # Reorder
rownames(Sap_SR_mod_anova) <- NULL
print(Sap_SR_mod_anova) 


