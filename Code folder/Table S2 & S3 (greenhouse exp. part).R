################################################################################
##################### Table S2 & S3 greenhouse exp. part #######################
################################################################################

# loading R packages
library(openxlsx) # version 4.2.5.2
library(vegan) # version 2.6-4
library(GUniFrac) # version 1.5

## load group data
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
# colnames(Green_group)

# loading database
green_otu <- read.xlsx("Greenhouse_ASVs_raw_data.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
rownames(green_otu) <- green_otu$ASV_ID
# green_otu[1:6,1:6]
green_otu <- green_otu[ ,Green_group$Sample_ID]
Green_group <- Green_group[colnames(green_otu), ]
#rownames(t(green_otu)) %in% rownames(Green_group)

## Tax INFORMATION
tax_default <- read.xlsx("Greenhouse_fungi_Flattening.xlsx", sheet = "ASV_tax", rowNames = F, colNames = T)[,1:3]
#head(tax_default)
rownames(tax_default)  <- tax_default$ASV_name
tax_default <- tax_default[colnames(t(green_otu)),]

OTU_tax <- tax_default %>% tidyr::separate(col = taxonomy_all, 
                                           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                           sep = ";\\s*", fill = "right")
head(OTU_tax[,1:6])
OTU_tax$GENUS <- sub("g__", "", OTU_tax$Genus)

# Loading FungalTraits database
FungalTraits <- read.xlsx("FungalTraits.xlsx", sheet = "FungalTraits", rowNames = F, colNames = T)
head(FungalTraits[,c(2:6,8)])
OTU_tax2 <- OTU_tax %>% left_join(FungalTraits[,c("GENUS","primary_lifestyle")], by = "GENUS")
head(OTU_tax2)

####  |Plant Pathogen|
Pathogen_id1 <- subset(OTU_tax2, primary_lifestyle == "plant_pathogen")$ASV_name
#### Fusarium
Pathogen_id2 <- subset(OTU_tax, Genus == "g__Fusarium")$ASV_name
### |Plant Pathogen| & Fusarium
Pathogens <- as.data.frame(unique(c(Pathogen_id1, Pathogen_id2)))
Pathogens$guild <- "Plant Pathogen"; colnames(Pathogens)[1] <- "ASV_name"

Green_pathogens = Pathogens

############################ Arbuscular Mycorrhizal ############################
####  Arbuscular Mycorrhizal
AMF_id1 <- subset(OTU_tax2, primary_lifestyle == "arbuscular_mycorrhizal")$ASV_name; length(AMF_id1)
#### Glomeromycota
AMF_id2 <- subset(OTU_tax, Phylum == "p__Glomeromycota")$ASV_name; length(AMF_id2)
### Arbuscular Mycorrhizal & Glomeromycota
AMF <- as.data.frame(unique(c(AMF_id1, AMF_id2)))
AMF$guild <- "Arbuscular Mycorrhizal"; colnames(AMF)[1] <- "ASV_name"

############################### saprophytic fungi ##############################
Saprotroph_type = c("soil_saprotroph","litter_saprotroph","unspecified_saprotroph","wood_saprotroph","nectar/tap_saprotroph","pollen_saprotroph")
Saprotroph <- as.data.frame(subset(OTU_tax2, primary_lifestyle %in% Saprotroph_type)$ASV_name)
Saprotroph$guild <- "Saprotroph"; colnames(Saprotroph)[1] <- "ASV_name"


##################### community richness of overall fungi ######################
Green_relative <- decostand(green_otu, method = "total", MARGIN = 2)
colSums(Green_relative)
Bray_dist_green_rela <- vegdist(t(Green_relative), method = 'bray')

set.seed(1234)
Overall_fungi_perM = GUniFrac::adonis3(t(Green_relative) ~ Origin/Species, method = "bray", by = "margin", 
                                    data = Green_group, permutations = 999)

Overall_perM_tab <- as.data.frame(Overall_fungi_perM$aov.tab)
Overall_perM_tab$R2 <- round(Overall_perM_tab$R2,3)
Overall_perM_tab$df <- paste0(Overall_perM_tab$Df, ",", Overall_perM_tab$Df[length(Overall_perM_tab$Df)-1])
Overall_perM_tab$`F` <- round(Overall_perM_tab$F.Model,2)
Overall_perM_tab$Predictors <- c("Origin", "Species","Residuals","Total")
Overall_perM_tab <- Overall_perM_tab[1:2 ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Overall_perM_tab) <- NULL
print(Overall_perM_tab)

##################### community richness of overall fungi ######################
green_overall_SR <- as.data.frame(specnumber(t(green_otu)))
colnames(green_overall_SR) <- "Overall_SR"
green_overall_SR$Sample_ID <- rownames(green_overall_SR)
Green_group = Green_group %>% left_join(green_overall_SR)

Overall_green_SR_mod <- lm(Overall_SR ~ Origin/Species, data = Green_group)
shapiro.test(residuals(Overall_green_SR_mod))

Overall_SR_mod_anova <- as.data.frame(anova(Overall_green_SR_mod))
Overall_SR_mod_anova[, 2:4] <- round(Overall_SR_mod_anova[, 2:4], digits = 2)
Overall_SR_mod_anova$`Pr(>F)` <- round(Overall_SR_mod_anova$`Pr(>F)`, 3)
Overall_SR_mod_anova$df <- paste0(Overall_SR_mod_anova$Df, ",", Overall_SR_mod_anova$Df[length(Overall_SR_mod_anova$Df)])
Overall_SR_mod_anova$Predictors <- c("Origin", "Species", "Residuals")
Overall_SR_mod_anova <- Overall_SR_mod_anova[, c(7,6,2:5)] # Reorder
rownames(Overall_SR_mod_anova) <- NULL
print(Overall_SR_mod_anova) 
Overall_SR_mod_anova$explan_var <- round((Overall_SR_mod_anova$`Sum Sq`/sum(Overall_SR_mod_anova$`Sum Sq`))*100, 3)


################## community composition of pathogenic fungi ###################
Path_Green_hel_no <- as.data.frame(t(green_otu))[ ,Pathogens$ASV_name]
Path_Green_rela <- as.data.frame(decostand(t(Path_Green_hel_no), method = "total", MARGIN = 2))
# colSums(Path_Green_rela)

set.seed(1234)
Path_fungi_perM = GUniFrac::adonis3(t(Path_Green_rela) ~ Origin/Species, method = "bray", by = "margin", 
                                       data = Green_group, permutations = 999)

Path_perM_tab <- as.data.frame(Path_fungi_perM$aov.tab)
Path_perM_tab$R2 <- round(Path_perM_tab$R2,3)
Path_perM_tab$df <- paste0(Path_perM_tab$Df, ",", Path_perM_tab$Df[length(Path_perM_tab$Df)-1])
Path_perM_tab$`F` <- round(Path_perM_tab$F.Model,2)
Path_perM_tab$Predictors <- c("Origin", "Species","Residuals","Total")
Path_perM_tab <- Path_perM_tab[1:2 ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Path_perM_tab) <- NULL
print(Path_perM_tab)

################### community richness of pathogenic fungi #####################
green_Path_SR <- as.data.frame(specnumber((Path_Green_hel_no)))
colnames(green_Path_SR) <- "Path_SR"
green_Path_SR$Sample_ID <- rownames(green_Path_SR)
Green_group = Green_group %>% left_join(green_Path_SR)

Path_green_SR_mod <- lm(Path_SR ~ Origin/Species, data = Green_group)
shapiro.test(residuals(Path_green_SR_mod))

Path_SR_mod_anova <- as.data.frame(anova(Path_green_SR_mod))
Path_SR_mod_anova[, 2:4] <- round(Path_SR_mod_anova[, 2:4], digits = 2)
Path_SR_mod_anova$`Pr(>F)` <- round(Path_SR_mod_anova$`Pr(>F)`, 3)
Path_SR_mod_anova$df <- paste0(Path_SR_mod_anova$Df, ",", Path_SR_mod_anova$Df[length(Path_SR_mod_anova$Df)])
Path_SR_mod_anova$Predictors <- c("Origin", "Species", "Residuals")
Path_SR_mod_anova <- Path_SR_mod_anova[, c(7,6,2:5)] # Reorder
rownames(Path_SR_mod_anova) <- NULL
Path_SR_mod_anova$explan_var <- round((Path_SR_mod_anova$`Sum Sq`/sum(Path_SR_mod_anova$`Sum Sq`))*100, 3)
print(Path_SR_mod_anova) 


###################### community composition of AM fungi #######################
AM_Green_hel_no <- as.data.frame(t(green_otu))[ ,AMF$ASV_name]
AM_Green_hel_no_add  <- AM_Green_hel_no
AM_Green_hel_no_add[rowSums(AM_Green_hel_no_add) == 0, ] <- 1e-6

AM_Green_rela <- as.data.frame(decostand(t(AM_Green_hel_no_add), method = "total", MARGIN = 2))
# colSums(AM_Green_rela)

set.seed(1234)
AM_fungi_perM = GUniFrac::adonis3(t(AM_Green_rela) ~ Origin/Species, method = "bray", by = "margin", 
                                    data = Green_group, permutations = 999)

AM_perM_tab <- as.data.frame(AM_fungi_perM$aov.tab)
AM_perM_tab$R2 <- round(AM_perM_tab$R2,3)
AM_perM_tab$df <- paste0(AM_perM_tab$Df, ",", AM_perM_tab$Df[length(AM_perM_tab$Df)-1])
AM_perM_tab$`F` <- round(AM_perM_tab$F.Model,2)
AM_perM_tab$Predictors <- c("Origin", "Species","Residuals","Total")
AM_perM_tab <- AM_perM_tab[1:2 ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(AM_perM_tab) <- NULL
print(AM_perM_tab)

######################## community richness of AM fungi ########################
green_AM_SR <- as.data.frame(specnumber((AM_Green_hel_no)))
colnames(green_AM_SR) <- "AM_SR"
green_AM_SR$Sample_ID <- rownames(green_AM_SR)
Green_group = Green_group %>% left_join(green_AM_SR)

AM_green_SR_mod <- lm(AM_SR ~ Origin/Species, data = Green_group)
shapiro.test(residuals(AM_green_SR_mod))

AM_SR_mod_anova <- as.data.frame(anova(AM_green_SR_mod))
AM_SR_mod_anova[, 2:4] <- round(AM_SR_mod_anova[, 2:4], digits = 2)
AM_SR_mod_anova$`Pr(>F)` <- round(AM_SR_mod_anova$`Pr(>F)`, 3)
AM_SR_mod_anova$df <- paste0(AM_SR_mod_anova$Df, ",", AM_SR_mod_anova$Df[length(AM_SR_mod_anova$Df)])
AM_SR_mod_anova$Predictors <- c("Origin", "Species", "Residuals")
AM_SR_mod_anova <- AM_SR_mod_anova[, c(7,6,2:5)] # Reorder
rownames(AM_SR_mod_anova) <- NULL
AM_SR_mod_anova$explan_var <- round((AM_SR_mod_anova$`Sum Sq`/sum(AM_SR_mod_anova$`Sum Sq`))*100, 3)
print(AM_SR_mod_anova) 


################# community composition of saprophytic fungi ###################
Sap_Green_hel_no <- as.data.frame(t(green_otu))[ ,Saprotroph$ASV_name]
Sap_Green_rela <- as.data.frame(decostand(t(Sap_Green_hel_no), method = "total", MARGIN = 2))
# colSums(Sap_Green_rela)

set.seed(1234)
Sap_fungi_perM = GUniFrac::adonis3(t(Sap_Green_rela) ~ Origin/Species, method = "bray", by = "margin", 
                                    data = Green_group, permutations = 999)

Sap_perM_tab <- as.data.frame(Sap_fungi_perM$aov.tab)
Sap_perM_tab$R2 <- round(Sap_perM_tab$R2,3)
Sap_perM_tab$df <- paste0(Sap_perM_tab$Df, ",", Sap_perM_tab$Df[length(Sap_perM_tab$Df)-1])
Sap_perM_tab$`F` <- round(Sap_perM_tab$F.Model,2)
Sap_perM_tab$Predictors <- c("Origin", "Species","Residuals","Total")
Sap_perM_tab <- Sap_perM_tab[1:2 ,c("Predictors", "df" ,"F", "R2", "Pr(>F)")] # Reorder
rownames(Sap_perM_tab) <- NULL
print(Sap_perM_tab)

################## community richness of saprophytic fungi #####################
green_Sap_SR <- as.data.frame(specnumber((Sap_Green_hel_no)))
colnames(green_Sap_SR) <- "Sap_SR"
green_Sap_SR$Sample_ID <- rownames(green_Sap_SR)
Green_group = Green_group %>% left_join(green_Sap_SR)

Sap_green_SR_mod <- lm(Sap_SR ~ Origin/Species, data = Green_group)
shapiro.test(residuals(Sap_green_SR_mod))

Sap_SR_mod_anova <- as.data.frame(anova(Sap_green_SR_mod))
Sap_SR_mod_anova[, 2:4] <- round(Sap_SR_mod_anova[, 2:4], digits = 2)
Sap_SR_mod_anova$`Pr(>F)` <- round(Sap_SR_mod_anova$`Pr(>F)`, 3)
Sap_SR_mod_anova$df <- paste0(Sap_SR_mod_anova$Df, ",", Sap_SR_mod_anova$Df[length(Sap_SR_mod_anova$Df)])
Sap_SR_mod_anova$Predictors <- c("Origin", "Species", "Residuals")
Sap_SR_mod_anova <- Sap_SR_mod_anova[, c(7,6,2:5)] # Reorder
rownames(Sap_SR_mod_anova) <- NULL
Sap_SR_mod_anova$explan_var <- round((Sap_SR_mod_anova$`Sum Sq`/sum(Sap_SR_mod_anova$`Sum Sq`))*100, 3)
print(Sap_SR_mod_anova) 


