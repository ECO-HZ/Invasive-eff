################################################################################
################################# Figure S3 ####################################
################################################################################
library(ggplot2) # version 3.5.2
library(openxlsx) # version 4.2.5.2
library(emmeans) # version 1.10.6

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
    size = 14, color = "black", fill = "grey90",
    box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
    halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
  plot.tag = element_text(size = 16, face = "bold")) 

# Loading the grouping metadata of soil samples
Field_group = read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

# Data Transformation
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Year <- as.factor(Field_group$Year)
Field_group$Site <- factor(Field_group$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_group$Origin <- factor(Field_group$Origin, levels = c("Native","Exotic"))

# 创建偏移量映射
offset_mapping <- c("2018" = -0.2, "2020" = 0, "2021" = 0.2)

## Soil total nitrogen content
data_Soil_N <- Field_group[,c("Origin", "Year" ,"Site", "Latitude", "Soil_N")]

mod4 = lm(Soil_N ~ Year * Latitude, data = data_Soil_N)
anova(mod4)
shapiro.test(residuals(mod4))


mod4_emtrends <- emtrends(mod4, pairwise ~ Year, var = "Latitude")
test(mod4_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Soil_N, mapping = aes(x = Latitude + offset_mapping[as.character(Year)], 
                                        y = Soil_N, fill = Year, group = Year),
             position = position_dodge(width = 0), size = 2.5, pch = 21) +  
  geom_smooth(data_Soil_N, mapping = aes(x = Latitude, y = Soil_N, color = Year),
              method = "lm", formula = y ~ x, se = F, size = 0.8, linetype = 1, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Soil_N, mapping = aes(x = Latitude, y = Soil_N, color = Year, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) +
  geom_smooth(data = data_Soil_N, mapping = aes(x = Latitude, y = Soil_N), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.5) + 
  ggpmisc::stat_poly_eq(data = data_Soil_N,  mapping = aes(x = Latitude, y = Soil_N, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) + 
  scale_fill_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  labs(x = NULL, y = "Soil total nitrogen content (%, sqrt)", tag = "a") -> p1; p1


## Soil water content
data_Wcont <- Field_group[,c("Origin", "Year" ,"Site", "Latitude", "Wcont")]

mod5 = lm(Wcont ~ Origin * Year * Latitude, data = data_Wcont)
shapiro.test(residuals(mod5))
mod5_emtrends <- emtrends(mod5, pairwise ~ Year, var = "Latitude")
test(mod5_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Wcont, mapping = aes(x = Latitude + offset_mapping[as.character(Year)], 
                                       y = Wcont, fill = Year, group = Year),
             position = position_dodge(width = 0), size = 2.5, pch = 21) +  
  #geom_errorbar(data_Wcont_mean, mapping = aes(x = Latitude, y = mean_Wcont, ymin = mean_Wcont - se_Wcont, ymax = mean_Wcont + se_Wcont, color = Year, group = Year), 
  #              width = 0.5, size = 0.5, position = position_dodge(width = 0)) +
  geom_smooth(data_Wcont, mapping = aes(x = Latitude, y = Wcont, color = Year, linetype = Year),
              method = "lm", formula = y ~ x, se = F, size = 0.8) + 
  ggpmisc::stat_poly_eq(data_Wcont, mapping = aes(x = Latitude, y = Wcont, color = Year, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) +
  ## 整体数据集
  geom_smooth(data = data_Wcont, mapping = aes(x = Latitude, y = Wcont), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.5) + 
  ggpmisc::stat_poly_eq(data = data_Wcont,  mapping = aes(x = Latitude, y = Wcont, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) + 
  scale_fill_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) +
  scale_linetype_manual(values = c(1,1,2)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  #scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  #scale_x_continuous(breaks=c(23.1,25.2,27.9,30.5,34.6,36.2)) + 
  labs(x = NULL, y = "Soil water content (g g-1, sqrt)", tag = "b") -> p2; p2

## Soil pH
data_Soil_ph <- Field_group[,c("Origin", "Year" ,"Site", "Latitude", "Soil_ph")]

mod3 = lm(Soil_ph ~ Year * Latitude, data = data_Soil_ph)
anova(mod3)
mod3_emtrends <- emtrends(mod3, pairwise ~ Year, var = "Latitude")
test(mod3_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Soil_ph, mapping = aes(x = Latitude + offset_mapping[as.character(Year)], 
                                         y = Soil_ph, fill = Year, group = Year),
             position = position_dodge(width = 0), size = 2.5, pch = 21) +  
  geom_smooth(data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph, color = Year, linetype = Year),
              method = "lm", formula = y ~ x, se = F, size = 0.8, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph, color = Year, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) +
  geom_smooth(data = data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.5) + 
  ggpmisc::stat_poly_eq(data = data_Soil_ph,  mapping = aes(x = Latitude, y = Soil_ph, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) + 
  scale_fill_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) +
  scale_linetype_manual(values = c(2,1,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = NULL, y = "Soil pH", tag = "c") -> p3; p3


## annual average temperature 
data_Tave <- unique(Field_group[,c("Year" ,"Site", "Latitude", "Tave")])

mod1 = lm(Tave ~ Year * Latitude, data_Tave)
anova(mod1)
mod1_emtrends <- emtrends(mod1, pairwise ~ Year, var = "Latitude")
test(mod1_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Tave, mapping = aes(x = Latitude + offset_mapping[as.character(Year)]
                                      , y = Tave, fill = Year),
             position = position_dodge(width = 0), size = 2.5, pch = 21) + 
  geom_smooth(data_Tave, mapping = aes(x = Latitude, y = Tave, color = Year),
              method = "lm", formula = y ~ x, se = F, size = 0.8, linetype = 1, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Tave, mapping = aes(x = Latitude, y = Tave, color = Year, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) + 
  geom_smooth(data = data_Tave, mapping = aes(x = Latitude, y = Tave), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.5) + 
  ggpmisc::stat_poly_eq(data = data_Tave,  mapping = aes(x = Latitude, y = Tave, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85,0.85)) + 
  scale_color_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) + 
  scale_fill_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  labs(x = NULL, y = expression("Annual mean temperature ( " * degree * "C)"), tag = "d") -> p4; p4


## annual precipitation
data_Prec <- unique(Field_group[,c("Year" ,"Site", "Latitude", "Prec")])
# anova(lm(Prec ~ Latitude, subset(data_Prec, Year == "2020")))
mod2 = lm(Prec ~ Year * Latitude, data_Prec)
anova(mod2)
mod2_emtrends <- emtrends(mod2, pairwise ~ Year, var = "Latitude")
test(mod2_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Prec, mapping = aes(x = Latitude + offset_mapping[as.character(Year)], 
                                      y = Prec, fill = Year),
             position = position_dodge(width = 0), size = 2.5, pch = 21) + 
  geom_smooth(data_Prec, mapping = aes(x = Latitude, y = Prec, color = Year, linetype = Year),
              method = "lm", formula = y ~ x, se = F, size = 0.8, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Prec, mapping = aes(x = Latitude, y = Prec, color = Year, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) + 
  geom_smooth(data = data_Prec, mapping = aes(x = Latitude, y = Prec), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.5) + 
  ggpmisc::stat_poly_eq(data = data_Prec,  mapping = aes(x = Latitude, y = Prec, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme +
  scale_color_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) + 
  scale_fill_manual(values = c("2018" = "#549299", "2020" = "#FABF5E","2021" = "#9F706A")) +
  scale_linetype_manual(values = c(1,2,2)) + 
  #scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  #scale_x_continuous(breaks=c(23.1,25.2,27.9,30.5,34.6,36.2)) + 
  labs(x = NULL, y = expression("Annual precipitation (mm)"), tag = "e") -> p5; p5


(p1/p4)|(p2/p5)|(p3/p5)


