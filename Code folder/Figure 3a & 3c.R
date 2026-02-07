################################################################################
############################## Figure 3a & 3c ##################################
################################################################################

# loading R packages
library(openxlsx) # version 4.2.5.2
library(dplyr) # version 1.1.1
library(ggplot2) # version 3.5.2
library(vegan) # version 2.6-4
library(ggtext) # version 0.1.2
library(phytools) # version 2.1-1
library(funrar) # version 1.5.0
library(patchwork) # version 1.3.1

# Custom style
mytheme = theme(
  legend.position = "none",
  panel.grid=element_blank(), 
  strip.background = element_rect(color="black", fill="white", size=0.5, linetype="solid"),
  strip.text.x = element_text(size = 9, color = "black"), # face = "bold.italic"
  legend.title = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 10),
  legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
  axis.ticks = element_line(color='black'),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(colour='black', size=14),
  axis.title.y = element_text(colour='black', size=14),
  axis.text = element_text(colour='black',size=12),
  plot.title = element_textbox(
    size = 14, color = "black", fill = "#E6E5E5",
    box.color = "black",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
    halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
  plot.tag = element_text(size = 16, face = "bold")) 


####################### Loading field survey database ##########################
# Soil sample grouping information in field survey
Field_group = read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

# Data Transformation
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Year <- as.factor(Field_group$Year)
Field_group$Site <- factor(Field_group$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_group$Origin <- factor(Field_group$Origin, levels = c("Native","Exotic"))

# Loading samples-abundance table in field survey
Field_raw_abun <- read.xlsx("Field_ASVs_raw_data.xlsx", sheet = "raw_ASVs", rowNames = T, colNames = T)
Field_raw_abun <- Field_raw_abun[,Field_group$Sample_ID]
# dim(Field_raw_abun)
Field_raw_abun <- Field_raw_abun[rowSums(Field_raw_abun) > 0, ]
# colSums(Field_raw_abun)

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
Field_fungi_relative <- decostand(Field_raw_abun, method = "total", MARGIN = 2)
# colSums(Field_fungi_relative)
Field_Bray_dist_rela <- vegdist(t(Field_fungi_relative), method = 'bray')

################### Loading greenhouse experiment database #####################
# Soil sample grouping information in greenhouse exp.
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)

# Soil samples-abundance table in greenhouse exp.
Green_raw_abun <- read.xlsx("Greenhouse_ASVs_raw_data.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_raw_abun <- Green_raw_abun[ ,Green_group$Sample_ID]
# Green_raw_abun[1:6,1:6]
# colSums(Green_raw_abun)

Green_fungi_relative <- decostand(Green_raw_abun, method = "total", MARGIN = 2)
# colSums(Green_fungi_relative)
Green_Bray_dist_rela <- vegdist(t(Green_fungi_relative), method = 'bray')
# save(Green_Bray_dist_rela, file = "Green_Bray_dist_rela.RData")

# note:
# To match the rhizosphere fungal data with the corresponding species from the 
# field survey, we calculated the mean pairwise dissimilarity values among 
# species (averaged across three replicate samples).
Green_dist <- as.matrix(Green_Bray_dist_rela)
Green_dist_data <- reshape2::melt(Green_dist, varnames = c("Sample_ID_A", "Sample_ID_B"),
                                  value.name = "dist", na.rm = T)
colnames(Green_dist_data)[1] <- "Sample_ID"
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
colnames(Green_dist_data)[c(1,2)] <- c("Sample_ID2","Sample_ID")
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
Green_mean_dist <- Rmisc::summarySE(Green_dist_data, measurevar = c("dist"), groupvars = c("Species.x", "Species.y"))
bc_dist_green_mean <- reshape2::dcast(Green_mean_dist, Species.x  ~ Species.y , value.var = "dist")
rownames(bc_dist_green_mean) <- bc_dist_green_mean$Species.x
bc_dist_green_mean <- bc_dist_green_mean[,-1]
diag(bc_dist_green_mean) <- 0
Green_Bray_dist_mean <- bc_dist_green_mean

# Loop
# Input Bray-Curtis distance matrix
Field_input_dist = Field_Bray_dist_rela
Green_input_dist = Green_Bray_dist_mean

# create null database for save data
Field_BC_data_all = NULL
Green_BC_data_all = NULL

Year = unique(Field_group$Year)
Site = unique(Field_group$Site)

for(i in Year){
  for(ii in Site) {
    
    ## Field survey
    select_group <- subset(Field_group, Year == i & Site == ii)
    Native_sample <- subset(select_group, Origin == "Native")$Sample_ID
    Exotic_sample <- subset(select_group, Origin == "Exotic")$Sample_ID
    
    ## select sp latin names
    Native_sp <- subset(select_group, Origin == "Native")$Species
    Exotic_sp <- subset(select_group, Origin == "Exotic")$Species
    
    ## Native
    Field_Native_matrix <- as.matrix(Field_input_dist)[Native_sample, Native_sample]
    dim(Field_Native_matrix)
    Field_lower_tri_mask <- lower.tri(Field_Native_matrix, diag = FALSE)
    Field_Native_matrix_lower <- Field_Native_matrix
    Field_Native_matrix_lower[!Field_lower_tri_mask] <- NA
    Field_Native_BC_long <- reshape2::melt(Field_Native_matrix_lower, na.rm = TRUE, 
                                           varnames = c("Sample_ID1", "Sample_ID2"), value.name = "Field_dist")
    Field_Native_BC_long$Origin <- "Native"
    
    ## Exotic
    Field_Exotic_matrix <- as.matrix(Field_input_dist)[Exotic_sample, Exotic_sample]
    dim(Field_Exotic_matrix)
    Field_lower_tri_mask <- lower.tri(Field_Exotic_matrix, diag = FALSE)
    Field_Exotic_matrix_lower <- Field_Exotic_matrix
    Field_Exotic_matrix_lower[!Field_lower_tri_mask] <- NA
    Field_Exotic_BC_long <- reshape2::melt(Field_Exotic_matrix_lower, na.rm = TRUE, 
                                           varnames = c("Sample_ID1", "Sample_ID2"), value.name = "Field_dist")
    Field_Exotic_BC_long$Origin <- "Exotic"
    
    ## merge data sets
    Field_Pairwise_BC_data <- rbind(Field_Native_BC_long, Field_Exotic_BC_long)
    
    ## adding species Latin names
    Sample_ID1 = select_group[,c("Sample_ID", "Species")]; colnames(Sample_ID1) = c("Sample_ID1", "Species1")
    Sample_ID2 = select_group[,c("Sample_ID", "Species")]; colnames(Sample_ID2) = c("Sample_ID2", "Species2")
    Field_Pairwise_BC_data = Field_Pairwise_BC_data %>% left_join(Sample_ID1) %>% left_join(Sample_ID2)
    
    Field_Pairwise_BC_data$Year <- i; Field_Pairwise_BC_data$Site <- ii
    
    ## merge data sets
    Field_BC_data_all = rbind(Field_BC_data_all, Field_Pairwise_BC_data)
    
    ########################## Greenhouse experiment ###########################
    ## Native
    Green_Native_matrix <- as.matrix(Green_input_dist)[Native_sp, Native_sp]
    dim(Green_Native_matrix)
    Green_lower_tri_mask <- lower.tri(Green_Native_matrix, diag = FALSE)
    Green_Native_matrix_lower <- Green_Native_matrix
    Green_Native_matrix_lower[!Green_lower_tri_mask] <- NA
    Green_Native_BC_long <- reshape2::melt(Green_Native_matrix_lower, na.rm = TRUE, 
                                           varnames = c("Species1", "Species2"), value.name = "Green_dist")
    Green_Native_BC_long$Origin <- "Native"
    
    ## Exotic
    Green_Exotic_matrix <- as.matrix(Green_input_dist)[Exotic_sp, Exotic_sp]
    dim(Green_Exotic_matrix)
    Green_lower_tri_mask <- lower.tri(Green_Exotic_matrix, diag = FALSE)
    Green_Exotic_matrix_lower <- Green_Exotic_matrix
    Green_Exotic_matrix_lower[!Green_lower_tri_mask] <- NA
    Green_Exotic_BC_long <- reshape2::melt(Green_Exotic_matrix_lower, na.rm = TRUE, 
                                           varnames = c("Species1", "Species2"), value.name = "Green_dist")
    Green_Exotic_BC_long$Origin <- "Exotic"
    
    ## merge data
    Green_Pairwise_BC_data <- rbind(Green_Native_BC_long, Green_Exotic_BC_long)
    Green_Pairwise_BC_data$Year <- i; Green_Pairwise_BC_data$Site <- ii
    
    ## Merge data sets
    Green_BC_data_all = rbind(Green_BC_data_all, Green_Pairwise_BC_data)
    
  }
}

# Compositional dissimilarity in rhizosphere fungal communities among natives or 
# aliens co-occurred at the same site and time in the field and greenhouse exp.
with_bc_data <- Field_BC_data_all %>% left_join(Green_BC_data_all)

# Add temperature information
with_bc_data <- with_bc_data %>% left_join(unique(Field_group[,c("Site", "Year", "Tave")]))

# Estimated the environmental effects
with_bc_data$env_Effect_on_fungi <- log(with_bc_data$Field_dist/with_bc_data$Green_dist)

mod = lm(env_Effect_on_fungi ~ Origin * Year * Site, data = with_bc_data)
anova(mod)

t.test(with_bc_data$env_Effect_on_fungi, mu = 0)
t.test(subset(with_bc_data, Origin == "Native")$ env_Effect_on_fungi,
       subset(with_bc_data, Origin == "Exotic")$ env_Effect_on_fungi)

# Loading traits database
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Loading plant phylogenetic tree
plant_tree <- read.newick("IQ_tree_plant_2025.newick")
plant_tree_max <- cophenetic(plant_tree)

traits_var = c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF")
Year = unique(Field_group$Year)
Site = unique(Field_group$Site)

# Loop
Field_traits_dist_all = NULL

for(i in Year){
  for(ii in Site) {
    
    ## Field survey
    select_group <- subset(Field_group, Year == i & Site == ii)
    Native_sample <- subset(select_group, Origin == "Native")$Sample_ID
    Exotic_sample <- subset(select_group, Origin == "Exotic")$Sample_ID
    
    ## select sp latin names
    Native_sp <- subset(select_group, Origin == "Native")$Species
    Exotic_sp <- subset(select_group, Origin == "Exotic")$Species
    
    # six functional traits
    traits_max = Field_group[,c("Sample_ID", "Chol", "SLA", "LDMC", "SRL", "FRR", "RMF")]
    traits_max$Sample_ID = NULL
    traits_max_dist = compute_dist_matrix(traits_max, metric = "euclidean", scale = T, center = T)
    ## Native
    N_Native_matrix <- as.matrix(traits_max_dist)[Native_sample, Native_sample]
    N_lower_tri_mask <- lower.tri(N_Native_matrix, diag = FALSE)
    N_Native_matrix_lower <- N_Native_matrix
    N_Native_matrix_lower[!N_lower_tri_mask] <- NA
    all_Native_BC_long <- reshape2::melt(N_Native_matrix_lower, na.rm = TRUE, 
                                         varnames = c("Sample_ID1", "Sample_ID2"), value.name = "all_traits")
    all_Native_BC_long$Origin <- "Native"
    
    ## Exotic
    N_Exotic_matrix <- as.matrix(traits_max_dist)[Exotic_sample, Exotic_sample]
    N_lower_tri_mask <- lower.tri(N_Exotic_matrix, diag = FALSE)
    N_Exotic_matrix_lower <- N_Exotic_matrix
    N_Exotic_matrix_lower[!N_lower_tri_mask] <- NA
    all_Exotic_BC_long <- reshape2::melt(N_Exotic_matrix_lower, na.rm = TRUE, 
                                         varnames = c("Sample_ID1", "Sample_ID2"), value.name = "all_traits")
    all_Exotic_BC_long$Origin <- "Exotic"
    
    ## merge data
    all_Pairwise_BC_data <- rbind(all_Native_BC_long, all_Exotic_BC_long) 
    
    # plant phylogenetic distance
    ## Native
    N_Native_matrix <- as.matrix(plant_tree_max)[Native_sp, Native_sp]
    N_lower_tri_mask <- lower.tri(N_Native_matrix, diag = FALSE)
    N_Native_matrix_lower <- N_Native_matrix
    N_Native_matrix_lower[!N_lower_tri_mask] <- NA
    phy_Native_BC_long <- reshape2::melt(N_Native_matrix_lower, na.rm = TRUE, 
                                         varnames = c("Species_ID1", "Species_ID2"), value.name = "Phylo")
    phy_Native_BC_long$Origin <- "Native"
    
    ## Exotic
    N_Exotic_matrix <- as.matrix(plant_tree_max)[Exotic_sp, Exotic_sp]
    N_lower_tri_mask <- lower.tri(N_Exotic_matrix, diag = FALSE)
    N_Exotic_matrix_lower <- N_Exotic_matrix
    N_Exotic_matrix_lower[!N_lower_tri_mask] <- NA
    phy_Exotic_BC_long <- reshape2::melt(N_Exotic_matrix_lower, na.rm = TRUE, 
                                         varnames = c("Species_ID1", "Species_ID2"), value.name = "Phylo")
    phy_Exotic_BC_long$Origin <- "Exotic"
    
    ## merge data
    phy_Pairwise_BC_data <- rbind(phy_Native_BC_long, phy_Exotic_BC_long) 
    
    sample_group1 <- subset(Field_group, Year == i & Site == ii)[,c("Sample_ID", "Species")]
    colnames(sample_group1) = c("Sample_ID1", "Species_ID1")
    
    sample_group2 <- subset(Field_group, Year == i & Site == ii)[,c("Sample_ID", "Species")]
    colnames(sample_group2) = c("Sample_ID2", "Species_ID2")
    
    phy_Pairwise_BC_data = phy_Pairwise_BC_data %>% left_join(sample_group1) %>% left_join(sample_group2)
    
    ## Merge data sets
    Field_traits_dist = all_Pairwise_BC_data %>% 
      left_join(phy_Pairwise_BC_data)
    
    Field_traits_dist$Year <- i; Field_traits_dist$Site <- ii
    
    ##
    Field_traits_dist_all = rbind(Field_traits_dist_all, Field_traits_dist)
  }
}

head(Field_traits_dist_add)

# plant richness per site per year
plant_SR = Field_group %>% group_by(Origin, Tave) %>%
  summarise(plant_sr = n())

# merge all information
with_bc_data$Year = as.factor(with_bc_data$Year)
Field_traits_dist_add$Year = as.factor(Field_traits_dist_add$Year)
with_bc_data_add <- with_bc_data %>% left_join(Field_traits_dist_add)

# averaging the data per site per year
with_bc_data_add_mean = with_bc_data_add %>%
  group_by(Origin, Tave) %>% 
  summarise(mean_env = mean(env_Effect_on_fungi, na.rm = TRUE),
            n = sum(!is.na(env_Effect_on_fungi)),                     
            se_env = sd(env_Effect_on_fungi, na.rm = TRUE) / sqrt(n), 
            ci_lower = mean_env - 1.96 * se_env,                     
            ci_upper = mean_env + 1.96 * se_env,                       
            mean_traits = mean(all_traits),
            mean_phylo = mean(Phylo)) %>%
  left_join(plant_SR) %>% as.data.frame()
# str(with_bc_data_add_mean)

# rename
with_bc_data_add_mean$Origin[with_bc_data_add_mean$Origin == "Exotic"] <- "Alien"
with_bc_data_add_mean$Origin = factor(with_bc_data_add_mean$Origin, levels = c("Native","Alien"))
with_bc_data_add_mean$plant_sr <- as.numeric(with_bc_data_add_mean$plant_sr)

# 
with_bc_data_add_mean2 = with_bc_data_add_mean
with_bc_data_add_mean2 = with_bc_data_add_mean2 %>% left_join(unique(Field_group[,c("Tave", "Site_pool")]))
mod0 <- lm(mean_env ~ Origin*Site_pool, data = with_bc_data_add_mean2)
anova(mod0)

# 计算站点水平下真菌richness 均值
fungal_richness_site = Field_group %>% group_by(Tave, Origin) %>% 
  summarise(site_SR = mean(Field_SR))
fungal_richness_site$Origin[fungal_richness_site$Origin == "Exotic"] <- "Alien"
with_bc_data_add_mean2 = with_bc_data_add_mean2 %>% left_join(fungal_richness_site)

mod0 <- lm(mean_env ~ Origin*site_SR, data = with_bc_data_add_mean2)
anova(mod0)

# relationship between plant richness and environmental effects
mod1 <- lm(mean_env ~ Origin*plant_sr, data = with_bc_data_add_mean)
anova(mod1)

# relationship between pairwise distance of plant functional traits and environmental effects
mod2 <- lm(mean_env ~ Origin*mean_traits, data = with_bc_data_add_mean)
anova(mod2)

# relationship between pairwise distance of plant phylogeny and environmental effects
mod3 <- lm(mean_env ~ Origin*mean_phylo, data = with_bc_data_add_mean)
anova(mod3)

################################# Figure 3a ####################################
ggplot(data = with_bc_data_add_mean, aes(x = Tave, y = mean_env, fill = Origin, color = Origin, group = Origin)) + 
  geom_smooth(method = "lm", se = F, aes(linetype = Origin), size = 0.5) +
  geom_errorbar(aes(ymin = mean_env - se_env*1.96,
                    ymax = mean_env + se_env*1.96), width = 0, linewidth = 0.5) +
  geom_point(size = 3, pch = 21, color = "black") + 
  ggpmisc::stat_poly_eq(aes(color = Origin, group = Origin, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, label.y.npc = "top", rr.digits = 3) + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85,0.15),
        panel.grid=element_blank(), 
        legend.background = element_rect(fill = NA),
        axis.title = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        plot.tag = element_text(size = 14, face = "bold")) + 
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_linetype_manual(values = c(1,2)) + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.1))) + 
  labs(x = "Annual mean temperature of sampling site (℃)",
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Pairwise~dissimilarity["estimated in the field"], 
                                       Pairwise~dissimilarity["estimated in greenhouse"]) ~ ")")), tag = "a",
       title = NULL) -> Figure_3a_right; Figure_3a_right

Rmisc::summarySE(with_bc_data_add, groupvars = "Origin", measurevar = "env_Effect_on_fungi")


with_bc_data_add_mean_all = with_bc_data_add %>%
  group_by(Origin) %>% 
  summarise(mean_env = mean(env_Effect_on_fungi, na.rm = TRUE),
            n = sum(!is.na(env_Effect_on_fungi)),                     
            se_env = sd(env_Effect_on_fungi, na.rm = TRUE) / sqrt(n), 
            ci_lower = mean_env - 1.96 * se_env,                     
            ci_upper = mean_env + 1.96 * se_env) 

with_bc_data_add_mean_all$Origin[with_bc_data_add_mean_all$Origin == "Exotic"] <- "Alien"
with_bc_data_add_mean_all$Origin = factor(with_bc_data_add_mean_all$Origin, levels = c("Native","Alien"))

ggplot(data = with_bc_data_add_mean_all, 
       aes(x = Origin, y = mean_env, fill = Origin)) +
  geom_point(size = 3.5, color = "black", pch = 21) + 
  geom_errorbar(aes(x = Origin, 
                    ymax = mean_env + se_env*1.96, 
                    ymin = mean_env - se_env*1.96),
                width=0,alpha = 1, color = "black")+
  theme_bw() + mytheme + 
  #geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) + 
  theme(legend.position = "none",
        panel.grid=element_blank(), 
        legend.background = element_rect(fill = NA),
        axis.title = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        plot.tag = element_text(size = 14, face = "bold")) + 
  labs(x = NULL,
       y = NULL,
       #y = bquote(atop("Environmental effects", 
       #                Ln ~ "(" ~ frac(Pairwise~dissimilarity["estimated in the field"], 
       #                                 Pairwise~dissimilarity["estimated in greenhouse"]) ~ ")")),
       title = NULL) + 
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) -> Figure_3a_left; Figure_3a_left


Figure_3a_right + annotation_custom(grob=ggplotGrob(Figure_3a_left), ymin = -0.45, ymax = -0.05, xmin=18, xmax=21.5) -> Figure_3a; Figure_3a


################################# Figure 3c ####################################
# difference in environmental effects between native and alien species (natives minus aliens)
with_bc_data_mean_nat <- subset(with_bc_data_add_mean, Origin == "Native")
colnames(with_bc_data_mean_nat)[-c(1:2)] <- paste0("nat_",colnames(with_bc_data_mean_nat)[-c(1:2)])
with_bc_data_mean_exo <- subset(with_bc_data_add_mean, Origin == "Alien")
colnames(with_bc_data_mean_exo)[-c(1:2)] <- paste0("exo_",colnames(with_bc_data_mean_exo)[-c(1:2)])

with_bc_data_mean_merg = with_bc_data_mean_nat[,-1] %>% left_join(with_bc_data_mean_exo[,-1])
with_bc_data_mean_merg$diff = with_bc_data_mean_merg$nat_mean_env - with_bc_data_mean_merg$exo_mean_env
with_bc_data_mean_merg$diff_se = sqrt(with_bc_data_mean_merg$nat_se_env^2 + with_bc_data_mean_merg$exo_se_env^2)

cor.test(with_bc_data_mean_merg$Tave, with_bc_data_mean_merg$diff, method = "spearman")

ggplot(data = with_bc_data_mean_merg, aes(x = Tave, y = diff)) + 
  geom_smooth(data = subset(with_bc_data_mean_merg, Tave != 14.70000), 
              mapping = aes(x = Tave, y = diff),
              method = "lm", se = F, formula = y ~ x, 
              linetype = 1, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = diff - diff_se*1.96,
                    ymax = diff + diff_se*1.96), width = 0, linewidth = 0.5) +
  geom_point(size = 3, pch = 21, color = "black", fill = "grey") + 
  geom_point(data = subset(with_bc_data_mean_merg, Tave == 14.70000), mapping = aes(x = Tave, y = diff),
             size = 3, pch = 21, color = "black", fill = "black") +
  ggpmisc::stat_poly_eq(data = subset(with_bc_data_mean_merg, Tave != 14.70000), 
                        mapping = aes(x = Tave, y = diff, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, rr.digits = 3,
                        label.y.npc = 0.98, label.x.npc = 0.95) + 
  annotate("text", x = 17.5, y = 0.53, label = "Excluding Tai'an in 2021:", size = 4) + 
  ggpmisc::stat_poly_eq(data = with_bc_data_mean_merg, 
                        mapping = aes(x = Tave, y = diff, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, rr.digits = 3,
                        label.y.npc = 0.92, label.x.npc = 0.95) + 
  annotate("text", x = 18.5, y = 0.479, label = "All sites:", size = 4) + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.80,0.25),
        panel.grid=element_blank(), 
        legend.background = element_rect(fill = NA),
        axis.title = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        plot.tag = element_text(size = 14, face = "bold")) + 
  scale_linetype_manual(values = c(1,2)) + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  labs(x = "Annual mean temperature of sampling site (℃)",
       y = "Difference in environmental effect\n(Natives - Aliens)", tag = "c",
       title = NULL) -> Figure_3c; Figure_3c


# 11.26 x 5.11
Figure_3a/Figure_3c

