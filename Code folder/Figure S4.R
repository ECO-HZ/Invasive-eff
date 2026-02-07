################################################################################
################################## Figure S4 ###################################
################################################################################

# Loading R packages
library(openxlsx) # version 4.2.5.2
library(vegan) # version 2.6-4

############################## (Field survey part) #############################
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil sample abundance information
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]

par(las = 1, mar = c(5, 5, 4, 2))
Field_min_depth <- min(rowSums(t(Field_otu_raw)))
rarecurve(as.data.frame(t(Field_otu_raw)), step = 100, 
          col = "black", cex = 0.6, label = FALSE, xlab = "", ylab = "", main = "")
abline(v = Field_min_depth, col = "black", lty = 2, lwd = 2)
title(main = "Field survey", xlab = "Number of sequences", ylab = "Number of ASVs")


###########################  (Greenhouse exp. part) ############################
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)

Green_otu_raw <- read.xlsx("Greenhouse_data_row_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]
dim(Green_otu_raw)

par(las = 1, mar = c(5, 5, 4, 2))
Green_min_depth <- min(rowSums(t(Green_otu_raw)))
rarecurve(as.data.frame(t(Green_otu_raw)), step = 100, 
          col = "black", cex = 0.6, label = FALSE, xlab = "", ylab = "", main = "")
abline(v = Green_min_depth, col = "black", lty = 2, lwd = 2)
title(main = "Greenhouse experiment", xlab = "Number of sequences", ylab = "Number of ASVs")

################################################################################
########## Calculation of asymptotic exponential Shannon diversity #############
################################################################################

# Loading the R packages
library(iNEXT)

################################# Field survey ################################# 
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil samples-abundance table
Field_otu_raw <- read.xlsx("Field_ASVs_raw_data.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]
Field_otu_raw[1:6, 1:6]
#colSums(Field_otu_raw)

# Calculation of asymptotic exponential Shannon diversity
results <- list()
for(i in 1:385) {
  results[[i]] <- ChaoShannon(as.data.frame(Field_otu_raw)[i], datatype = "abundance")
  print(paste0(i, ":", colnames(Field_otu_raw)[i]))
}

Field_chao_shannon <- do.call(rbind, results)
mod_Field <- lm(Estimator ~ Observed, data = Field_chao_shannon)
summary(mod_Field)
cor.test(Field_chao_shannon$Observed, Field_chao_shannon$Estimator)

############################# Greenhouse experiment ############################
# Soil sample grouping information in greenhouse exp.
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)

# Soil samples-abundance table in greenhouse exp.
Green_otu_raw <- read.xlsx("Greenhouse_ASVs_raw_data.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]
Green_otu_raw[1:6,1:6]
#colSums(Green_otu_raw)

# Calculation of asymptotic exponential Shannon diversity
results <- list()
for(i in 1:162) {
  results[[i]] <- ChaoShannon(as.data.frame(Green_otu_raw)[i], datatype = "abundance")
  print(paste0(i, ":", colnames(Green_otu_raw)[i]))
}

Green_chao_shannon <- do.call(rbind, results)
mod_Green <- lm(Estimator ~ Observed, data = Green_chao_shannon)
summary(mod_Green)
cor.test(Green_chao_shannon$Observed, Green_chao_shannon$Estimator)

# Therefore, we used the raw sequencing abundance data for all analysis.

