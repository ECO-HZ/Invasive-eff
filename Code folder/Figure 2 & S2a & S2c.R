################################################################################
########################### Figure 2 & S2a & S2c ###############################
################################################################################

# loading R packages
library(openxlsx) # version 4.2.5.2
library(ggplot2) # version 3.5.2
library(vegan) # version 2.6-4
library(ggridges) # version 0.5.6
library(ggh4x) # version 0.2.8
library(ggtext) # version 0.1.2
library(dplyr) # version 1.1.1
library(ggplotify) # version 0.1.2
library(aplot) # version 0.2.2

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

site_colors <- c("Guangzhou" = "#0E4879", "Guilin" = "#3E91B7", "Changsha" = "#999999",
                 "Wuhan" = "#CFBD9F", "Zhengzhou" = "#94684E", "Tai'an" = "#BF5B1D")

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
Field_t_test_result_BC = NULL
Field_BC_data_all = NULL

Green_t_test_result_BC = NULL
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
    
    ## merge data
    Field_Pairwise_BC_data <- rbind(Field_Native_BC_long, Field_Exotic_BC_long)
    
    ## adding species Latin names
    Sample_ID1 = select_group[,c("Sample_ID", "Species")]; colnames(Sample_ID1) = c("Sample_ID1", "Species1")
    Sample_ID2 = select_group[,c("Sample_ID", "Species")]; colnames(Sample_ID2) = c("Sample_ID2", "Species2")
    Field_Pairwise_BC_data = Field_Pairwise_BC_data %>% left_join(Sample_ID1) %>% left_join(Sample_ID2)
    
    Field_Pairwise_BC_data$Year <- i; Field_Pairwise_BC_data$Site <- ii
    
    ## Student's t-test 
    t_test_mod <- t.test(Field_dist ~ Origin, data = Field_Pairwise_BC_data)
    Field_t_test_result <- data.frame(Type = "Field", Year = i, Site = ii,
                                      t_value = round(t_test_mod$statistic, 2),
                                      p_value = round(t_test_mod$p.value, 3))
    
    ## Merge data sets
    Field_t_test_result_BC = rbind(Field_t_test_result_BC, Field_t_test_result)
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
    
    ## Student's t-test 
    t_test_mod <- t.test(Green_dist ~ Origin, data = Green_Pairwise_BC_data)
    Greem_t_test_result <- data.frame(Type = "Green", Year = i, Site = ii,
                                      t_value = round(t_test_mod$statistic, 2),
                                      p_value = round(t_test_mod$p.value, 3))
    
    ## Merge data sets
    Green_t_test_result_BC = rbind(Green_t_test_result_BC, Greem_t_test_result)
    Green_BC_data_all = rbind(Green_BC_data_all, Green_Pairwise_BC_data)
    
  }
}

head(Field_BC_data_all)
head(Green_BC_data_all)

#pair_bc_data_add = Field_BC_data_all %>% left_join(Green_BC_data_all)

################################## Figure 2a ###################################
Field_BC_data_all$Site = factor(Field_BC_data_all$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

# rename
Field_BC_data_all$Origin[Field_BC_data_all$Origin == "Exotic"] <- "Alien"
Field_BC_data_all$Origin = factor(Field_BC_data_all$Origin, levels = c("Native","Alien"))
Field_BC_data_all$Year = factor(Field_BC_data_all$Year, levels = rev(c("2018", "2020", "2021")))

# reorder by origin and site
Field_pair_BC_sorted <- Field_BC_data_all %>% arrange(Site, Year, Origin)

# Student’s t-tests for total database in field survey
t.test(Field_pair_BC_sorted$Field_dist ~ Field_pair_BC_sorted$Origin)
(0.7951722 - 0.7759438)/0.7759438 # decreased 2.5% 

# plot data for total database
ggplot(Field_pair_BC_sorted, aes(x = Field_dist, y = "Overall", fill = Origin, color = Origin)) +
  geom_density_ridges(rel_min_height = 0.01, scale = 0.3,alpha = 0.5, linewidth = 0.7,
                      quantile_lines = TRUE, quantile_fun = mean) +
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01),limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25),expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0.36, 0))) +
  theme_minimal() +
  theme(plot.title = element_textbox(size = 14, color = "black", fill = "grey90",
                                     box.color = "grey50", padding = margin(5, 5, 5, 5), 
                                     margin = margin(b = 0),halign = 0.5, width = grid::unit(1, "npc")),
        legend.position = "none", legend.title = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 10),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_line(color = 'black'),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_blank(),
        strip.background = element_rect(color=NA, size=0.5, linetype="solid"),
        strip.placement = "outside",
        strip.text.y = element_text(size = 12, colour = "black"),
        panel.spacing = unit(0, "lines")) +
  labs(y = NULL,x = NULL, title = "Field survey") -> p1; p1


# plot data per site per year
print(subset(Field_t_test_result_BC, p_value <= 0.05))

ggplot(Field_pair_BC_sorted, aes(x = Field_dist, y = Year, fill = Origin, color = Origin)) +
  geom_density_ridges(rel_min_height = 0.01, scale = 0.9, alpha = 0.5, linewidth = 0.7,
                      quantile_lines = TRUE, quantile_fun = mean) +
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01), limits = c(0,1), 
                     breaks = seq(0, 1, by = 0.25), expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0.36, 0))) +
  ggh4x::facet_grid2(Site ~ ., #switch = "y", 
                     strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = site_colors))) +
  theme_minimal() +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        axis.line.x = element_line(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        axis.title.x = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12),
        strip.background = element_rect(color=NA, size=0.5, linetype="solid"),
        strip.placement = "outside",
        strip.text.y = element_text(size = 12, colour = "black"),
        panel.spacing = unit(0, "lines")) +
  labs(x = 'Pairwise Bray–Curtis dissimilarities', y = NULL) +
  geom_segment(aes(x = 0, xend = 0, y = 1, yend = 3), color = "black") -> p2; p2

## 5.83 x 11.00
mian_Fig_2_left <- p2 %>% insert_top(p1,height = 0.1) %>% as.ggplot() 
mian_Fig_2_left # 


############################## Greenhouse experiment ###########################
Green_BC_data_all$Site <- factor(Green_BC_data_all$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

# rename
Green_BC_data_all$Origin[Green_BC_data_all$Origin == "Exotic"] <- "Alien"
Green_BC_data_all$Origin <- factor(Green_BC_data_all$Origin, levels = c("Native","Alien"))
Green_BC_data_all$Year <- factor(Green_BC_data_all$Year, levels = rev(c("2018", "2020", "2021")))

# reorder by origin and site
Green_pair_BC_sorted <- Green_BC_data_all %>% arrange(Site, Year, Origin)

## Student’s t-tests for total database in greenhouse experiment
Green_fungi_relative <- decostand(Green_raw_abun, method = "total", MARGIN = 2)
# colSums(Green_fungi_relative)
Green_Bray_dist_rela <- vegdist(t(Green_fungi_relative), method = 'bray')

native_sample <- subset(Green_group, Origin == "Native")$Sample_ID
exotic_sample <- subset(Green_group, Origin == "Exotic")$Sample_ID

# native
Green_Bray_nat <- as.matrix(Green_Bray_dist_rela)[native_sample, native_sample]
# dim(Green_Bray_nat)
Green_lower_tri_mask <- lower.tri(Green_Bray_nat, diag = FALSE)
Green_Native_matrix_lower <- Green_Bray_nat
Green_Native_matrix_lower[!Green_lower_tri_mask] <- NA
Green_Native_BC_long <- reshape2::melt(Green_Native_matrix_lower, na.rm = TRUE, 
                                       varnames = c("Sample_ID1", "Sample_ID2"), value.name = "Green_dist")
Green_Native_BC_long$Origin <- "Native"

# Alien
Green_Bray_exo <- as.matrix(Green_Bray_dist_rela)[exotic_sample, exotic_sample]
# dim(Green_Bray_exo)
Green_lower_tri_mask <- lower.tri(Green_Bray_exo, diag = FALSE)
Green_Exotic_matrix_lower <- Green_Bray_exo
Green_Exotic_matrix_lower[!Green_lower_tri_mask] <- NA
Green_Exotic_BC_long <- reshape2::melt(Green_Exotic_matrix_lower, na.rm = TRUE, 
                                       varnames = c("Sample_ID1", "Sample_ID2"), value.name = "Green_dist")
Green_Exotic_BC_long$Origin <- "Alien"

# merge data
Green_Pairwise_BC_data <- rbind(Green_Native_BC_long, Green_Exotic_BC_long)
Green_Pairwise_BC_data$Origin <- factor(Green_Pairwise_BC_data$Origin, levels = c("Native","Alien"))
t.test(Green_dist ~ Origin,Green_Pairwise_BC_data) # t = 0.4151, p = 0.6781

# plot data for total database
ggplot(Green_Pairwise_BC_data, aes(x = Green_dist, y = "Overall", fill = Origin, color = Origin)) +
  geom_density_ridges(rel_min_height = 0.01, scale = 0.3,alpha = 0.5, linewidth = 0.7,
                      quantile_lines = TRUE, quantile_fun = mean) +
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01),limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25),expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0.36, 0))) +
  theme_minimal() +
  theme(plot.title = element_textbox(size = 14, color = "black", fill = "grey90",
                                     box.color = "grey50", padding = margin(5, 5, 5, 5), 
                                     margin = margin(b = 0),halign = 0.5, width = grid::unit(1, "npc")),
        legend.position = "none", legend.title = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 10),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_line(color = 'black'),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_blank(),
        strip.background = element_rect(color=NA, size=0.5, linetype="solid"),
        strip.placement = "outside",
        strip.text.y = element_text(size = 12, colour = "black"),
        panel.spacing = unit(0, "lines")) +
  labs(y = NULL,x = NULL, title = "Greenhouse experiment") -> p3; p3

# plot data per site per year
print(subset(Green_t_test_result_BC, p_value <= 0.05))

t.test(Green_dist ~ Origin, subset(Green_pair_BC_sorted, Site == "Tai'an" & Year == "2018"))

ggplot(Green_pair_BC_sorted, aes(x = Green_dist, y = Year, fill = Origin, color = Origin)) +
  geom_density_ridges(rel_min_height = 0.01, scale = 0.9, alpha = 0.5, linewidth = 0.7,
                      quantile_lines = TRUE, quantile_fun = mean) +
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01), limits = c(0,1), 
                     breaks = seq(0, 1, by = 0.25), expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0.36, 0))) +
  ggh4x::facet_grid2(Site ~ ., #switch = "y", 
                     strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = site_colors))) +
  theme_minimal() +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        axis.line.x = element_line(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        axis.title.x = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 12),
        strip.background = element_rect(color=NA, size=0.5, linetype="solid"),
        strip.placement = "outside",
        strip.text.y = element_text(size = 12, colour = "black"),
        panel.spacing = unit(0, "lines")) +
  labs(x = 'Pairwise Bray–Curtis dissimilarities', y = NULL) -> p4; p4


## 5.83 x 11.00
mian_Fig_2_right <- p4 %>% insert_top(p3,height = 0.1) %>% as.ggplot() 
mian_Fig_2_right 

mian_Fig_2_left|mian_Fig_2_right -> mian_Fig_2; mian_Fig_2


############################### Figure S2a & 2c #################################
# The relationships between the mean pairwise fungal compositional dissimilarity 
# among co-occurring natives (b) or aliens (c) estimated in the field survey and 
# those estimated in the greenhouse experiment. 

pair_BC_sorted_merge = Field_pair_BC_sorted %>% left_join(Green_pair_BC_sorted)

# native database
pair_BC_sorted_merge_nat = subset(pair_BC_sorted_merge, Origin == "Native")

t.test(pair_BC_sorted_merge_nat$Field_dist, pair_BC_sorted_merge_nat$Green_dist)
(0.7951722 - 0.6984232)/0.6984232 # increased 13.9%

# Student t.test
Year = unique(pair_BC_sorted_merge_nat$Year)
Site = unique(pair_BC_sorted_merge_nat$Site)
all_t_test_result = NULL

for (i in Year) {
  for (ii in Site) {
    sub_group = subset(pair_BC_sorted_merge_nat, Year == i & Site == ii)
    
    t.test(sub_group$Field_dist, sub_group$Green_dist)
    ## Student's t-test 
    t_test_mod <- t.test(sub_group$Field_dist, sub_group$Green_dist, data = sub_group)
    t_test_result <- data.frame(Type = "Field vs. Greenhouse", Year = i, Site = ii,
                                t_value = round(t_test_mod$statistic, 2),
                                p_value = round(t_test_mod$p.value, 3))
    
    ## Merge data sets
    all_t_test_result = rbind(all_t_test_result, t_test_result)
  }
}

all_t_test_result$sig = ifelse(all_t_test_result$p_value <= 0.05, 1, 0)
print(subset(all_t_test_result, sig == 1))

colnames(β_Bray_Native)
β_Bray_Native = pair_BC_sorted_merge_nat %>%
  group_by(Year, Site) %>%
  summarise(
    Field_mean = mean(Field_dist, na.rm = TRUE),
    Field_sd = sd(Field_dist, na.rm = TRUE),
    Field_se = Field_sd / sqrt(n()),
    #
    Green_mean = mean(Green_dist, na.rm = TRUE),
    Green_sd = sd(Green_dist, na.rm = TRUE),
    Green_se = Green_sd / sqrt(n()),
    .groups = "drop") %>%
  left_join(all_t_test_result)

fit <- deming::deming(Field_mean ~ Green_mean, ystd = Field_sd, xstd = Green_sd, data = β_Bray_Native)
print(fit)

β_Bray_Native$Site <- factor(β_Bray_Native$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  #geom_abline(intercept = 0, slope = 1, color = "#8B0000", linetype = 1, size = 1) +
  geom_abline(intercept = -0.1287923, slope = 1.3337283, linetype = 2, linewidth = 0.5) +
  geom_errorbar(data = β_Bray_Native, mapping = aes(x = Green_mean, ymax = Field_mean + 1.96*Field_se, ymin = Field_mean - 1.96*Field_se, color = Site), width = 0, size = 0.5) +
  geom_errorbarh(data = β_Bray_Native, mapping = aes(y = Field_mean, xmax = Green_mean + 1.96*Green_se, xmin = Green_mean - 1.96*Green_se, color = Site), height = 0, size = 0.5) +
  geom_point(β_Bray_Native, mapping = aes(x = Green_mean,y = Field_mean, shape = Year, fill = Site, color = Site), size = 3, color = "black") +
  # add signal symbol
  #geom_point(subset(β_Bray_Native, sig == "0"), mapping = aes(x = Green_mean,y = Field_mean, shape = Year, fill = Site, color = Site), fill = "white",size = 2.5) +
  scale_color_manual(values = site_colors)+
  scale_fill_manual(values = site_colors)+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.54,0.88)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.58,0.92)) + 
  #annotate("segment", x = 0.56, xend = 0.58, y = 0.90, yend = 0.90, size = 1, color = "#8B0000") +
  #annotate("text", x = 0.60, y = 0.901, label = "1:1 line", size = 4) + 
  #annotate("segment", x = 0.56, xend = 0.58, y = 0.88, yend = 0.88, size = 1, color = "black") + 
  #annotate("text", x = 0.605, y = 0.881, label = "Model fit", size = 4) + 
  theme_bw() + mytheme + 
  #theme(legend.position = c(0.12, 0.60),
  #      legend.title = element_text(size = 10, color = "black")) + 
  labs(x = NULL, 
       y = "Mean pair-wise Bray–Curtis dissimilarities\nestimated in the field",
       tag = "a", title = c("Among co-occurring native species")) -> Figure_S2a; Figure_S2a


################################ Figure 2c #####################################
# Alien database
pair_BC_sorted_merge_exo = subset(pair_BC_sorted_merge, Origin == "Alien")

t.test(pair_BC_sorted_merge_exo$Field_dist, pair_BC_sorted_merge_exo$Green_dist)
(0.7759438 - 0.7281435 )/0.7281435 # increased 6.6%

# Student t.test
Year = unique(pair_BC_sorted_merge_exo$Year)
Site = unique(pair_BC_sorted_merge_exo$Site)
all_t_test_result = NULL

for (i in Year) {
  for (ii in Site) {
    sub_group = subset(pair_BC_sorted_merge_exo, Year == i & Site == ii)
    
    t.test(sub_group$Field_dist, sub_group$Green_dist)
    ## Student's t-test 
    t_test_mod <- t.test(sub_group$Field_dist, sub_group$Green_dist, data = sub_group)
    t_test_result <- data.frame(Type = "Field vs. Greenhouse", Year = i, Site = ii,
                                t_value = round(t_test_mod$statistic, 2),
                                p_value = round(t_test_mod$p.value, 3))
    
    ## Merge data sets
    all_t_test_result = rbind(all_t_test_result, t_test_result)
  }
}

all_t_test_result$sig = ifelse(all_t_test_result$p_value <= 0.05, 1, 0)
print(subset(all_t_test_result, sig == 1))


colnames(pair_BC_sorted_merge_exo)
β_Bray_Exotic = pair_BC_sorted_merge_exo %>%
  group_by(Year, Site) %>%
  summarise(
    Field_mean = mean(Field_dist, na.rm = TRUE),
    Field_sd = sd(Field_dist, na.rm = TRUE),
    Field_se = Field_sd / sqrt(n()),
    #
    Green_mean = mean(Green_dist, na.rm = TRUE),
    Green_sd = sd(Green_dist, na.rm = TRUE),
    Green_se = Green_sd / sqrt(n()),
    .groups = "drop") %>%
  left_join(all_t_test_result)

fit <- deming::deming(Field_mean ~ Green_mean, ystd = Field_sd, xstd = Green_sd, data = β_Bray_Exotic)
print(fit)

β_Bray_Exotic$Site <- factor(β_Bray_Exotic$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  #geom_abline(intercept = 0, slope = 1, color = "#96393D", linetype = 1, size = 1)+
  geom_abline(intercept = -3.030977, slope = 5.517322, linetype = 2, linewidth = 0.5)+
  geom_errorbar(data = β_Bray_Exotic, mapping = aes(x = Green_mean, ymax = Field_mean + 1.96*Field_se, ymin = Field_mean - 1.96*Field_se, color = Site), width = 0, size = 0.5) +
  geom_errorbarh(data = β_Bray_Exotic, mapping = aes(y = Field_mean, xmax = Green_mean + 1.96*Green_se, xmin = Green_mean - 1.96*Green_se, color = Site), height = 0, size = 0.5) +
  geom_point(β_Bray_Exotic, mapping = aes(x = Green_mean,y = Field_mean, shape = Year, fill = Site, color = Site), size = 3, color = "black") +
  # add signal symbol
  #geom_point(subset(β_Bray_Exotic, sig == "0"), mapping = aes(x = Green_mean,y = Field_mean, shape = Year, fill = Site, color = Site), fill = "white",size = 2.5) +
  scale_color_manual(values = site_colors)+
  scale_fill_manual(values = site_colors)+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.57,0.83)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.54,0.92)) + 
  theme_bw() + mytheme + 
  #theme(legend.position = c(0.12, 0.60),
  #      legend.title = element_text(size = 10, color = "black")) + 
  labs(x = "Mean pair-wise Bray–Curtis\ndissimilarities estimated in the greenhouse experiment", 
       y = "", title = c("Among co-occurring alien species"), tag = "c") -> Figure_S2c; Figure_S2c


