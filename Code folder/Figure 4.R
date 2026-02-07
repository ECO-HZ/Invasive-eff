################################################################################
################################# Figure 4 #####################################
################################################################################

# loading R packages
library(openxlsx) # version 4.2.5.2
library(dplyr) # version 1.1.1
library(ggplot2) # version 3.5.2
library(vegan) # version 2.6-4
library(ggtext) # version 0.1.2
library(phytools) # version 2.1-1
library(funrar) # version 1.5.0
library(ggpubr) # version 0.6.0
library(patchwork) # version 1.3.1
library(emmeans) # version 1.10.6
library(iCAMP) # version 1.5.12
library(NST) # version 3.1.10 # need to be NST >=3.0.3

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

# Custom function
mean_ci_196 <- function(x) {
  x <- x[!is.na(x)]
  m  <- mean(x)
  se <- sd(x) / sqrt(length(x))
  c(y = m, ymin = m - 1.96 * se, ymax = m + 1.96 * se)
}

# common ASVS
common_ASVs = intersect(rownames(Field_raw_abun), rownames(Green_raw_abun))
length(common_ASVs)

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
#Field_group$Origin <- factor(Field_group$Origin, levels = c("Native","Exotic"))

# Loading samples-abundance table in field survey
Field_raw_abun <- read.xlsx("Field_ASVs_raw_data.xlsx", sheet = "raw_ASVs", rowNames = T, colNames = T)
rownames(Field_raw_abun) = Field_raw_abun$`#OTU.ID`; Field_raw_abun$`#OTU.ID` = NULL    
Field_raw_abun <- Field_raw_abun[,Field_group$Sample_ID]
# dim(Field_raw_abun)
Field_raw_abun <- Field_raw_abun[rowSums(Field_raw_abun) > 0, ]
# colSums(Field_raw_abun)

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
Field_fungi_relative <- decostand(Field_raw_abun, method = "total", MARGIN = 2)
# colSums(Field_fungi_relative)
Field_Bray_dist_rela <- vegdist(t(Field_fungi_relative), method = 'bray')

# common_ASV DIST IN FIELD
Field_Bray_dist_rela <- vegdist(t(Field_fungi_relative[common_ASVs, ]), method = 'bray')


################### Loading greenhouse experiment database #####################
# Soil sample grouping information in greenhouse exp.
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)

# Soil samples-abundance table in greenhouse exp.
Green_raw_abun <- read.xlsx("Greenhouse_ASVs_raw_data.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
rownames(Green_raw_abun) = Green_raw_abun$ASV_ID; Green_raw_abun$ASV_ID = NULL
Green_raw_abun <- Green_raw_abun[ ,Green_group$Sample_ID]
# Green_raw_abun[1:6,1:6]
# colSums(Green_raw_abun)

Green_fungi_relative <- decostand(Green_raw_abun, method = "total", MARGIN = 2)
# colSums(Green_fungi_relative)
Green_Bray_dist_rela <- vegdist(t(Green_fungi_relative), method = 'bray')
# save(Green_Bray_dist_rela, file = "Green_Bray_dist_rela.RData")

# common_ASV DIST IN GREENHOUSE
Green_Bray_dist_rela <- vegdist(t(Green_fungi_relative[common_ASVs, ]), method = 'bray')


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

# Loading plant phylogenetic tree
plant_tree <- read.newick("IQ_tree_plant_2025.newick")
plant_tree_max <- cophenetic(plant_tree)

# Loop
# Input Bray-Curtis distance matrix
Field_input_dist = Field_Bray_dist_rela
Green_input_dist = Green_Bray_dist_mean

# create null database for save data
Field_BC_data_all = NULL

Year = unique(Field_group$Year)
Site = unique(Field_group$Site)

for(i in Year){
  for(ii in Site) {
    
    ## Field survey
    select_group <- subset(Field_group, Year == i & Site == ii)
    group_sample <- select_group$Sample_ID
    
    ## select sp latin names
    group_sp <- select_group$Species
    
    ## all species per site per year
    Field_all_matrix <- as.matrix(Field_input_dist)[group_sample, group_sample]
    dim(Field_all_matrix)
    Field_lower_tri_mask <- lower.tri(Field_all_matrix, diag = FALSE)
    Field_all_matrix_lower <- Field_all_matrix
    Field_all_matrix_lower[!Field_lower_tri_mask] <- NA
    Field_all_BC_long <- reshape2::melt(Field_all_matrix_lower, na.rm = TRUE, 
                                        varnames = c("Sample_ID1", "Sample_ID2"), value.name = "Field_dist")
    Field_all_BC_long$Group <- "Field"
    
    ## adding species Latin names
    Sample_ID1 = select_group[,c("Sample_ID", "Species")]; colnames(Sample_ID1) = c("Sample_ID1", "Species1")
    Sample_ID2 = select_group[,c("Sample_ID", "Species")]; colnames(Sample_ID2) = c("Sample_ID2", "Species2")
    Field_all_BC_long = Field_all_BC_long %>% left_join(Sample_ID1) %>% left_join(Sample_ID2)
    
    Field_all_BC_long$Year <- i; Field_all_BC_long$Site <- ii
    
    # functional traits distance
    traits_max = select_group[,c("Sample_ID", "Chol", "SLA", "LDMC", "SRL", "FRR", "RMF")]
    traits_max$Sample_ID = NULL
    traits_max_dist = compute_dist_matrix(traits_max, metric = "euclidean", scale = T, center = T)
    dim(traits_max_dist)
    
    ## all species
    traits_all_matrix <- as.matrix(traits_max_dist)
    traits_lower_tri_mask <- lower.tri(traits_all_matrix, diag = FALSE)
    traits_all_matrix_lower <- traits_all_matrix
    traits_all_matrix_lower[!traits_lower_tri_mask] <- NA
    all_traits_BC_long <- reshape2::melt(traits_all_matrix_lower, na.rm = TRUE, 
                                         varnames = c("Sample_ID1", "Sample_ID2"), value.name = "all_traits")
    all_traits_BC_long$Group <- "Field"
    
    # plant phylogenetic distance
    phylo_all_matrix <- as.matrix(plant_tree_max)[group_sp, group_sp]
    phylo_lower_tri_mask <- lower.tri(phylo_all_matrix, diag = FALSE)
    phylo_all_matrix_lower <- phylo_all_matrix
    phylo_all_matrix_lower[!phylo_lower_tri_mask] <- NA
    all_phylo_BC_long <- reshape2::melt(phylo_all_matrix_lower, na.rm = TRUE, 
                                        varnames = c("Species_ID1", "Species_ID2"), value.name = "Phylo")
    all_phylo_BC_long$Group <- "Field"
    
    sample_group1 <- subset(Field_group, Year == i & Site == ii)[,c("Sample_ID", "Species")]
    colnames(sample_group1) = c("Sample_ID1", "Species_ID1")
    
    sample_group2 <- subset(Field_group, Year == i & Site == ii)[,c("Sample_ID", "Species")]
    colnames(sample_group2) = c("Sample_ID2", "Species_ID2")
    
    all_phylo_BC_long = all_phylo_BC_long %>% left_join(sample_group1) %>% left_join(sample_group2)
    
    # merge all data sets
    Field_all_BC_long_add <- Field_all_BC_long %>% left_join(all_traits_BC_long) %>% left_join(all_phylo_BC_long)
    
    Field_BC_data_all = rbind(Field_BC_data_all, Field_all_BC_long_add)
    
  }
}

cor.test(Field_BC_data_all$Field_dist, Field_BC_data_all$all_traits)
cor.test(Field_BC_data_all$Field_dist, Field_BC_data_all$Phylo)
plot(Field_BC_data_all$Phylo, Field_BC_data_all$Field_dist)

################################# Figure 4a ####################################
ggplot(data = Field_BC_data_all, aes(x = all_traits, y = Field_dist)) + 
  geom_point(size = 3, pch = 21, color = "black", fill = "grey") + 
  geom_smooth(method = "lm", se = F, formula = y ~ x, linetype = 2,linewidth=0.5, color = "black") +
  ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, label.y.npc = "bottom", rr.digits = 3) + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = NA)) + 
  labs(x = "Pairwise plant functional distance",
       y = "Pairwise Bray–Curtis dissimilarities\nestimated in the field", tag = "a") -> p1; p1

################################# Figure 4b ####################################
ggplot(data = with_bc_data_add, aes(x = Phylo, y = Field_dist)) + 
  geom_point(size = 3, pch = 21, color = "black", fill = "grey") + 
  geom_smooth(method = "lm", se = F, formula = y ~ x, linetype = 2,linewidth=0.5, color = "black") +
  ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, label.y.npc = "bottom", rr.digits = 3) + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = NA)) +
  labs(x = "Pairwise plant phylogenetic distance",
       y = "Pairwise Bray–Curtis dissimilarities\nestimated in the field", tag = "b") -> p2; p2


################################# Figure 4c ####################################
Functinal_taxa_group <- Field_group
Functinal_taxa_group$Origin[Functinal_taxa_group$Origin == "Exotic"] <- "Alien"
Functinal_taxa_group$Origin <- factor(Functinal_taxa_group$Origin, levels = c("Native", "Alien"))

t.test(Functinal_taxa_group$Rela_generalist*100 ~ Functinal_taxa_group$Origin)

ggplot(data = Functinal_taxa_group, 
       aes(x = Origin, y = (Rela_generalist*100), fill = Origin)) +
  stat_summary(fun = mean, geom = "point", size = 3.5, color = "black", pch = 21) +
  stat_summary(fun.data = mean_ci_196, geom = "errorbar", 
               width = 0, color = "black") +
  theme_bw() + mytheme +
  labs(x = "", y = "Relative abundance (%)", tag = "c", title = "Generalist") +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) + 
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) -> p3_a; p3_a


t.test(Functinal_taxa_group$Rela_specialist*100 ~ Functinal_taxa_group$Origin)

ggplot(data = Functinal_taxa_group, 
       aes(x = Origin, y = (Rela_specialist*100), fill = Origin)) +
  stat_summary(fun = mean, geom = "point", size = 3.5, color = "black", pch = 21) +
  stat_summary(fun.data = mean_ci_196, geom = "errorbar", 
               width = 0, color = "black") +
  theme_bw() + mytheme +
  labs(x = "", y = "Relative abundance (%)", title = "Specialist") +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) + 
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) -> p3_b; p3_b

# combination plot
(p1|p2|p3_a|p3_b) + plot_layout(widths = c(0.5,0.5,0.2,0.2)) -> Figure_4_top
# ggsave("Figure_4_top.pdf", plot = Figure_4_top, width = 13.32, height = 4.85, units = "in", dpi = 300)


############################# Figure 4d & 4e & 4f ##############################
# We quantified the complexity of rhizosphere fungal interaction networks using 
# the community cohesion metric, which is an abundance-weighted and null 
# model–corrected pairwise taxon correlations, following the methods by 
# Herren CM, McMahon KD. Cohesion: a method for quantifying the connectivity 
# of microbial communities. ISME 11, 2426-2438 (2017).

## definition function
# Find the number of zeros in the data
zero <- function(vec){  
  num.zero <- length(which(vec == 0))  
  return(num.zero)
}

# Create a function that averages only negative values in the data
neg.mean <- function(vector){  
  neg.vals <- vector[which(vector < 0)]  
  n.mean <- mean(neg.vals)  
  if(length(neg.vals) == 0) n.mean <- 0  
  return(n.mean)
}

# Create a function that averages only positive values in the data
pos.mean <- function(vector){  
  pos.vals <- vector[which(vector > 0)]  
  p.mean <- mean(pos.vals)  
  if(length(pos.vals) == 0) p.mean <- 0  
  return(p.mean)
}


# define parameters
# The higher the iter value, the longer the script runs
iter <- 100

# Decide whether to use taxonomic group/column shuffling (tax.shuffle=T) or 
# row shuffling algorithm (tax.shuffle=F)
tax.shuffle <- T

# Enter your own correlation table options
use.custom.cors <- F

# input data (36 data sets (native or alien per site per year))
sub_field_group = subset(Field_group, Site == "Wuhan" & Year == "2020" & Origin == "Native")
otu_filtered = fungi_Flattening[,sub_field_group$Sample_ID]
otu_filtered = as.data.frame(t(otu_filtered[rowSums(otu_filtered) > 0, ]))
otu_filtered[1:3,1:3]
dim(otu_filtered)

b <- otu_filtered
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]
dim(c)

# 保存原始矩阵中每个样本的总个体数。如果数据为相对丰度，则此值为1，但如果矩阵c为计数数据，则不是。
rowsums.orig <- rowSums(c)

# 只保留至少存在于 3 个样本中的 ASV
d <- c[, colSums(c > 0) >= 3]

# 移除因 ASV 被筛掉而没有任何丰度的样本
d <- d[rowSums(d) > 0, ]
dim(d)
# View(d)
# 创建相对丰度矩阵。
rel.d <- d / rowsums.orig
# 可选，检查剔除分类群后保留的群落比例
#hist(rowSums(rel.d))

# 创建观察到的相关性矩阵
cor.mat.true <- cor(rel.d)
dim(cor.mat.true)


# 创建保存初始分类群的中位otu-otu相关性的向量
med.tax.cors <- vector()
# 运行此循环以获取期望的成对相关性（空模型）
# 如果选择输入自定义相关性矩阵，则绕过空模型
if(use.custom.cors == F) {
  if(tax.shuffle == T) { 
    pb <- txtProgressBar(min = 0, max = dim(rel.d)[2], style = 3)  # 初始化进度条
    for(which.taxon in 1:dim(rel.d)[2]){        
      # 创建向量以保存每个单独分类群的每次排列的相关性    
      ## perm.cor.vec.mat代表排列相关性向量矩阵    
      perm.cor.vec.mat <- vector()        
      for(i in 1:iter){      
        # 创建与rel.d相同维度的空矩阵      
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])      
        rownames(perm.rel.d) <- rownames(rel.d)      
        colnames(perm.rel.d) <- colnames(rel.d)            
        # 对每个分类群      
        for(j in 1:dim(rel.d)[2]){         
          # 用排列的分类群向量替换原始分类群向量        
          perm.rel.d[, j ] <- sample(rel.d[ ,j ])       
        }            
        # 不随机化焦点列      
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]            
        # 计算排列矩阵的相关性矩阵      
        cor.mat.null <- cor(perm.rel.d)            
        # 对每次迭代，保存焦点分类群与其他分类群之间的空矩阵相关性向量      
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])          
      }    
      # 保存焦点分类群与所有其他分类群之间的中位相关性    
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      setTxtProgressBar(pb, which.taxon)  # 更新进度条
    }
    close(pb)  # 关闭进度条
  } else {  
    pb <- txtProgressBar(min = 0, max = dim(rel.d)[2], style = 3)  # 初始化进度条
    for(which.taxon in 1:dim(rel.d)[2]){       
      # 创建向量以保存每个单独分类群的每次排列的相关性    
      ## perm.cor.vec.mat代表排列相关性向量矩阵    
      perm.cor.vec.mat <- vector()        
      for(i in 1:iter){      
        # 创建用于洗牌丰度的重复矩阵      
        perm.rel.d <- rel.d             
        # 对每个分类群      
        for(j in 1:dim(rel.d)[1]){         
          which.replace <- which(rel.d[j, ] > 0 )         
          # 如果焦点分类群大于零，则从替换向量中去掉，以保持焦点丰度不变        
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]                
          # 用一个向量替换原始分类群向量，该向量中大于0的值已被随机排列        
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal])       
        }
        # 计算排列矩阵的相关性矩阵      
        cor.mat.null <- cor(perm.rel.d)            
        # 对每次迭代，保存焦点分类群与其他分类群之间的空矩阵相关性向量      
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])          
      }    
      # 保存焦点分类群与所有其他分类群之间的中位相关性    
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))        
      # 对大型数据集来说，这可以帮助了解这个循环将运行多久    
      setTxtProgressBar(pb, which.taxon)  # 更新进度条
    } 
    close(pb)  # 关闭进度条
  }
}

# 保存观察与预期的相关性差。如果use.custom.cors = TRUE，则使用自定义相关性
if(use.custom.cors == T) {  
  obs.exp.cors.mat <- custom.cor.mat.sub} else{    
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors}
diag(obs.exp.cors.mat) <- 0


#### 
#### 生成所需的连接度和凝聚度向量
# 通过对观察到的正负相关性求平均来计算连接度
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
length(connectedness.pos)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)
length(connectedness.neg)

# 通过将相对丰度数据集与相关的连接度相乘来计算凝聚度
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

print(data.frame(cohesion.pos = cohesion.pos, cohesion.neg = cohesion.neg))

#### 将向量组合成一个列表
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", "Negative Cohesion", "Positive Cohesion")
print(output)

# Iteratively run the native and alien plant datasets for each year and each site, 
# sequentially retaining the calculated cohesion index data.

# Lading saved cohesion index data
cohesion_data <- read.xlsx("cohesion_data.xlsx", colNames = T, rowNames = F)

# add other group information
Field_group2 = Field_group %>% left_join(cohesion_data)
Field_group2$total.cohesion = Field_group2$cohesion.pos + abs(Field_group2$cohesion.neg)
Field_group2$Origin[Field_group2$Origin == "Exotic"] <- "Alien"
Field_group2$Origin <- factor(Field_group2$Origin, levels = c("Native", "Alien"))

################################# Figure 4d ####################################
t.test(Field_group2$total.cohesion ~ Field_group2$Origin)

ggplot(data = Field_group2, aes(x = Origin, y = total.cohesion, fill = Origin)) +
  geom_boxplot(outliers = T, width = 0.35, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=2, color = "black", aes(fill = Origin)) + 
  theme_bw() + mytheme + 
  #stat_compare_means(method = "t.test", 
  #                   comparisons = list(c("Native", "Alien")), 
  #                   label = "p.format", hide.ns = TRUE) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) + 
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "", y = "Tatal cohesion", tag = "d") -> Figure_4d_left; Figure_4d_left

# Slope 
mod1 <- lm(scale(total.cohesion) ~ Origin*scale(Tave), data = Field_group2)
anova(mod1)
mod1_emtrends <- emtrends(mod1, pairwise ~ Origin, var = "Tave")
test(mod1_emtrends, adjust = "BH")

ggplot(data = Field_group2, aes(x = Tave, y = total.cohesion, color = Origin, group = Origin)) + 
  geom_point(size = 3, pch = 21, aes(fill = Origin), color = "black") + 
  geom_smooth(method = "lm", se = F, formula = y ~ x, aes(linetype = Origin),linewidth=0.5) +
  ggpmisc::stat_poly_eq(aes(color = Origin, group = Origin, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, label.y.npc = "top", rr.digits = 3) + 
  scale_linetype_manual(values = c(1,1)) + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85, 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = NA)) + 
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  labs(x = expression("Annual mean temperature ("*degree*C*")"),
       y = NULL) -> Figure_4d_right; Figure_4d_right

# Extract the y-axis range and scale of right plot
range_y <- ggplot_build(Figure_4d_right)$layout$panel_params[[1]]$y.range
breaks_y <- ggplot_build(Figure_4d_right)$layout$panel_params[[1]]$y$breaks

# Set to the same y-axis scale
Figure_4d_left <- Figure_4d_left + scale_y_continuous(limits = range_y, breaks = breaks_y)
Figure_4d_right <- Figure_4d_right + scale_y_continuous(limits = range_y, breaks = breaks_y)

# joint plot
Figure_4d_left + Figure_4d_right + plot_layout(widths = c(0.25, 0.7)) -> Figure_4d; Figure_4d


################################# Figure 4e ####################################
t.test(Field_group2$cohesion.pos ~ Field_group2$Origin)

ggplot(data = Field_group2, aes(x = Origin, y = cohesion.pos, fill = Origin)) +
  geom_boxplot(outliers = T, width = 0.35, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=2, color = "black", aes(fill = Origin)) + 
  theme_bw() + mytheme + 
  #stat_compare_means(method = "t.test", 
  #                   comparisons = list(c("Native", "Alien")), 
  #                   label = "p.format", hide.ns = TRUE) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) + 
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "", y = "Positive cohesion", tag = "e") -> Figure_4e_left; Figure_4e_left

# Slope 
mod2 <- lm(scale(cohesion.pos) ~ Origin*scale(Tave), data = Field_group2)
anova(mod2)
mod2_emtrends <- emtrends(mod2, pairwise ~ Origin, var = "Tave")
test(mod2_emtrends, adjust = "BH")

ggplot(data = Field_group2, aes(x = Tave, y = cohesion.pos, color = Origin, group = Origin)) + 
  geom_point(size = 3, pch = 21, aes(fill = Origin), color = "black") + 
  geom_smooth(method = "lm", se = F, formula = y ~ x, aes(linetype = Origin),linewidth=0.5) +
  ggpmisc::stat_poly_eq(aes(color = Origin, group = Origin, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, label.y.npc = "top", rr.digits = 3) + 
  scale_linetype_manual(values = c(1,1)) + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85, 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = NA)) + 
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  labs(x = expression("Annual mean temperature ("*degree*C*")"),
       y = NULL) -> Figure_4e_right; Figure_4e_right

# Extract the y-axis range and scale of right plot
range_y <- ggplot_build(Figure_4e_right)$layout$panel_params[[1]]$y.range
breaks_y <- ggplot_build(Figure_4e_right)$layout$panel_params[[1]]$y$breaks

# Set to the same y-axis scale
Figure_4e_left <- Figure_4e_left + scale_y_continuous(limits = range_y, breaks = breaks_y)
Figure_4e_right <- Figure_4e_right + scale_y_continuous(limits = range_y, breaks = breaks_y)

# joint plot
Figure_4e_left + Figure_4e_right + plot_layout(widths = c(0.25, 0.7)) -> Figure_4e; Figure_4e


################################# Figure 4f ####################################
t.test(Field_group2$cohesion.neg ~ Field_group2$Origin)

ggplot(data = Field_group2, aes(x = Origin, y = cohesion.neg, fill = Origin)) +
  geom_boxplot(outliers = T, width = 0.35, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=2, color = "black", aes(fill = Origin)) + 
  theme_bw() + mytheme + 
  #stat_compare_means(method = "t.test", 
  #                   comparisons = list(c("Native", "Alien")), 
  #                   label = "p.format", hide.ns = TRUE) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) + 
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "", y = "Negative cohesion", tag = "f") -> Figure_4f_left; Figure_4f_left

# Slope 
mod3 <- lm(scale(cohesion.neg) ~ Origin*scale(Tave), data = Field_group2)
anova(mod3)
mod3_emtrends <- emtrends(mod3, pairwise ~ Origin, var = "Tave")
test(mod3_emtrends, adjust = "BH")

ggplot(data = Field_group2, aes(x = Tave, y = cohesion.neg, color = Origin, group = Origin)) + 
  geom_point(size = 3, pch = 21, aes(fill = Origin), color = "black") + 
  geom_smooth(method = "lm", se = F, formula = y ~ x, aes(linetype = Origin),linewidth=0.5) +
  ggpmisc::stat_poly_eq(aes(color = Origin, group = Origin, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, label.y.npc = "top", rr.digits = 3) + 
  scale_linetype_manual(values = c(1,1)) + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85, 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = NA)) + 
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  labs(x = expression("Annual mean temperature ("*degree*C*")"),
       y = NULL) -> Figure_4f_right; Figure_4f_right

# Extract the y-axis range and scale of right plot
range_y <- ggplot_build(Figure_4f_right)$layout$panel_params[[1]]$y.range
breaks_y <- ggplot_build(Figure_4f_right)$layout$panel_params[[1]]$y$breaks

# Set to the same y-axis scale
Figure_4f_left <- Figure_4f_left + scale_y_continuous(limits = range_y, breaks = breaks_y)
Figure_4f_right <- Figure_4f_right + scale_y_continuous(limits = range_y, breaks = breaks_y)

# joint plot
Figure_4f_left + Figure_4f_right + plot_layout(widths = c(0.25, 0.7)) -> Figure_4f; Figure_4f


################################# Figure 4g ####################################
# taxonomic NST
dist.method="bray"

# Sets the number of threads used by parallel operations
nworker=8

# randomization time for null model analysis. 
rand.time=1000

# Prefix name of output file
prefix="Field_NST"

#tnst_j = tNST(comm=t(Field_raw_abun), group=Group, dist.method=dist.method, 
#            abundance.weighted=TRUE, rand=rand.time, output.rand=TRUE, nworker=nworker, 
#            LB=FALSE, null.model="PF", between.group=TRUE, SES=TRUE, RC=TRUE)
#saveRDS(tnst_j, "tnst_j.RData")

load("tnst_j.RData")
DATA_nst = tnst_j$index.pair.grp

DATA_nst <- DATA_nst %>%
  tidyr::separate(group, into = c("Site", "Year"), sep = "_", remove = FALSE)

# conversion factor variable
DATA_nst$Site <- factor(DATA_nst$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
DATA_nst$Year <- as.factor(DATA_nst$Year)

# add information of temperature 
DATA_nst = DATA_nst %>% left_join(unique(Field_group[,c("Site", "Year", "Latitude", "Tave")]))

colnames(DATA_nst)

# add information on species origin grouping
name1_ID = Field_group[,c("Sample_ID", "Origin")]; colnames(name1_ID) = c("name1","Origin1")
name2_ID = Field_group[,c("Sample_ID", "Origin")]; colnames(name2_ID) = c("name2","Origin2")

DATA_nst <- DATA_nst %>% left_join(name1_ID) %>% left_join(name2_ID)
DATA_nst$Orign_within = paste0(DATA_nst$Origin1,"|",DATA_nst$Origin2)
unique(DATA_nst$Orign_within)

# among native or alien species datasets
DATA_nst_filter <- subset(DATA_nst, Orign_within == "Native|Native"|Orign_within == "Exotic|Exotic")
DATA_nst_filter$Orign_within[DATA_nst_filter$Orign_within == "Native|Native"] <- "Native"
DATA_nst_filter$Orign_within[DATA_nst_filter$Orign_within == "Exotic|Exotic"] <- "Alien"
DATA_nst_filter$Orign_within <- factor(DATA_nst_filter$Orign_within, levels = c("Native", "Alien"))

ggplot(data = DATA_nst_filter, aes(x = Orign_within, y = MST.ij.bray, fill = Orign_within)) +
  geom_boxplot(outliers = T, width = 0.35, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=2, color = "black", aes(fill = Orign_within)) + 
  theme_bw() + mytheme + 
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("Native", "Alien")), 
                     label = "p.format",hide.ns = TRUE) + 
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) + 
  labs(x = "", y = "Modified stochasticity ratio (%)", tag = "g") -> Figure_4g_left; Figure_4g_left


mod4 = lm(scale(MST.ij.bray) ~ Orign_within*scale(Tave), data = DATA_nst_filter)
anova(mod4)
mod4_emtrends <- emtrends(mod4, pairwise ~ Orign_within, var = "Tave")
test(mod4_emtrends, adjust = "BH")

ggplot(data = DATA_nst_filter, aes(x = Tave, y = MST.ij.bray, color = Orign_within)) + 
  geom_point(size = 3, pch = 21, aes(fill = Orign_within), color = "black") + 
  geom_smooth(method = "lm", se = F, formula = y ~ x, aes(linetype = Orign_within),linewidth=0.5) +
  ggpmisc::stat_poly_eq(aes(color = Orign_within, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, label.y.npc = "top", rr.digits = 3) + 
  scale_linetype_manual(values = c(1,2)) + 
  theme_bw() + mytheme + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = NA)) + 
  scale_color_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  scale_fill_manual(values = c("Native" = "#356A5D", "Alien" = "#ECD45D")) +
  labs(x = expression("Annual mean temperature ("*degree*C*")"),
       y = NULL) -> Figure_4g_right; Figure_4g_right

# 
range_y <- ggplot_build(Figure_4g_right)$layout$panel_params[[1]]$y.range
breaks_y <- ggplot_build(Figure_4g_right)$layout$panel_params[[1]]$y$breaks

# 
Figure_4g_left <- Figure_4g_left + scale_y_continuous(limits = range_y, breaks = breaks_y)
Figure_4g_right <- Figure_4g_right + scale_y_continuous(limits = range_y, breaks = breaks_y)

# 
Figure_4g_left + Figure_4g_right + plot_layout(widths = c(0.25, 0.75)) -> Figure_4g; Figure_4g



