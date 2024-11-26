## Date 2024-11-19
## This is the code for visualizing data obtained from a behavioral assay.
## Worms are conditioned with either the vehicle (V) or the odorant (O) and
## tested for dispersal (no choice assay) and chemotaxis (V vs O choice assay).
## N ≤ 8 arenas with > 100 adult animals (n).

## set working directory
setwd("/Users/samihatasnim/Desktop/Codes_Figures") ## set working directory

## install and load packages for data frame manipulation and plotting
install.packages(tidyverse)
install.packages(readxl)
install.packages(ggpubr)
install.packages(effsize)
install.packages(svglite)
install.packages(entropy)
install.packages(gridExtra)
install.packages(knitr)
install.packages(ggplot2)
library(tidyverse)
library(readxl)
library(ggpubr)
library(effsize)
library(svglite)
library(entropy)
library(gridExtra)
library(knitr)
library(ggplot2)

#custom color palette (paste as appropriate below)
eth_vsbutanone= "#d55e00"
but_vsbutanone= "#8d78a3"
eth_vsbenzaldehyde="#c35090"
eth_vsnonanone= "#0072b2"

## set up data frame(s) and input relevant information
date = "20210125"
strain = "N2"
vehicle = "eth"
odorant = "but"
orientation = "horizontal"
person = "Samiha"
errorbars= "#000000"
V_dispcolor= "#939598"
O_dispcolor= "#bcbec0"
V_color= "#d55e00"
O_color= "#8d78a3"

df_LV <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/20241117 EtOH_disp.xlsx")) ## df_LV is dispersal counts post vehicle pre-exposure
df_LO <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/20241117 Odor_disp.xlsx")) ## df_LO is dispersal counts post odorant pre-exposure
df_CV <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/20241117 EtOH_chem.xlsx")) ## df_CV is chemotaxis counts post vehicle pre-exposure
df_CO <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/20241117 Odor_chem.xlsx")) ## df_CO is chemotaxis counts post odorant pre-exposure

## define starting theme; modified later during plot generation

my_theme <- theme(panel.background = element_blank(),
                  strip.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_line(linetype = 1, size = 0.75, lineend="butt"),
                  axis.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  axis.ticks = element_line(size = 0.75), axis.ticks.length = unit(.2, "cm"),
                  axis.text = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  plot.title = element_text(family = "Helvetica", colour = "black", size = 10, hjust = 0, vjust = 0),
                  legend.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10, hjust = 0, vjust = 0.5),
                  legend.text = element_text(size = 10, hjust = 0, vjust = 0.5),
                  strip.text = element_text(size = 10))

## processing df_LV (data processing for final plot(s))
df_LV <- df_LV[,colnames(df_LV)[colnames(df_LV) !="Quadrant"]] ## gets rid of first column
df_LV = df_LV[colSums(df_LV) > 100] ## picks columns/lanes where n > 100
Ndf_LV = ncol(df_LV) ## readout of N (number of columns/lanes picked)
df_LV_23ex <- as.data.frame(subset(df_LV[2:3,])) ## extracts data from rows Q2 and Q3 for processing
df_LV_23sum <-as.data.frame(t(colSums(df_LV_23ex[1:ncol(df_LV)]))) ## ## calculates sum of Q2 and Q3 for each column/lane and adds it to a dataframe
df_LV_calc <- rbind(as.data.frame(df_LV[1,]), as.data.frame(df_LV_23sum[1,]), as.data.frame(df_LV[4,])) ## dataframe on which statistical calculations will be done
row.names(df_LV_calc)<- c("Q1", "Q23", "Q4") ## renaming rows 
df_LV_calc_proportion <- as.data.frame(t(t(df_LV_calc)/colSums(df_LV_calc))) ## calculating proportions of worms in quadrants
df_LV_calc_mean <- rowMeans(df_LV_calc_proportion) ## calculating mean proportions
processed_df_LV <- cbind(df_LV_calc_proportion, 
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LV_calc_proportion))$Q1)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LV_calc_proportion))$Q23)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LV_calc_proportion))$Q4)$conf.int)))),
                         df_LV_calc_mean,
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LV_calc_proportion))$Q1, as.data.frame(t(df_LV_calc_proportion))$Q23)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LV_calc_proportion))$Q23, as.data.frame(t(df_LV_calc_proportion))$Q4)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LV_calc_proportion))$Q1, as.data.frame(t(df_LV_calc_proportion))$Q4)$p.value))))) ## binding everything together
a <- Ndf_LV
colnames(processed_df_LV) <- c(paste0("VO", 1:a), "lower_df_LV", "upper_df_LV", "mean_df_LV", "pval_df_LV") ## renaming columns
df_LV_calc_proportion_index_calc <- as.data.frame(t(df_LV_calc_proportion)) ## transposing dataframe for calculating dispersal
chemotaxis_LV <- as.data.frame(((df_LV_calc_proportion_index_calc[3] - df_LV_calc_proportion_index_calc[1]) / (df_LV_calc_proportion_index_calc[3] + df_LV_calc_proportion_index_calc[1]))) ## calculating chemotaxis indices for each column/lane

## processing df_LO (data processing for final plot(s))

df_LO <- df_LO[,colnames(df_LO)[colnames(df_LO) !="Quadrant"]] ## gets rid of first column
df_LO = df_LO[colSums(df_LO) > 100] ## picks columns/lanes where n > 100
Ndf_LO = ncol(df_LO) ## readout of N (number of columns/lanes picked)
df_LO_23ex <- as.data.frame(subset(df_LO[2:3,])) ## extracts data from rows Q2 and Q3 for processing
df_LO_23sum <-as.data.frame(t(colSums(df_LO_23ex[1:ncol(df_LO)]))) ## ## calculates sum of Q2 and Q3 for each column/lane and adds it to a dataframe
df_LO_calc <- rbind(as.data.frame(df_LO[1,]), as.data.frame(df_LO_23sum[1,]), as.data.frame(df_LO[4,])) ## dataframe on which statistical calculations will be done
row.names(df_LO_calc)<- c("Q1", "Q23", "Q4") ## renaming rows 
df_LO_calc_proportion <- as.data.frame(t(t(df_LO_calc)/colSums(df_LO_calc))) ## calculating proportions of worms in quadrants
df_LO_calc_mean <- rowMeans(df_LO_calc_proportion) ## calculating mean proportions
processed_df_LO <- cbind(df_LO_calc_proportion, 
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LO_calc_proportion))$Q1)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LO_calc_proportion))$Q23)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LO_calc_proportion))$Q4)$conf.int)))),
                         df_LO_calc_mean,
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LO_calc_proportion))$Q1, as.data.frame(t(df_LO_calc_proportion))$Q23)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LO_calc_proportion))$Q23, as.data.frame(t(df_LO_calc_proportion))$Q4)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_LO_calc_proportion))$Q1, as.data.frame(t(df_LO_calc_proportion))$Q4)$p.value))))) ## binding everything together
b <- Ndf_LO
colnames(processed_df_LO) <- c(paste0("VO", 1:b), "lower_df_LO", "upper_df_LO", "mean_df_LO", "pval_df_LO") ## renaming columns
df_LO_calc_proportion_index_calc <- as.data.frame(t(df_LO_calc_proportion)) ## transposing dataframe for calculating dispersal
chemotaxis_LO <- as.data.frame(((df_LO_calc_proportion_index_calc[3] - df_LO_calc_proportion_index_calc[1]) / (df_LO_calc_proportion_index_calc[3] + df_LO_calc_proportion_index_calc[1]))) ## calculating chemotaxis indices for each column/lane

## processing to plot dispersal counts as mean proportions in quadrants

to_plot_dispersal_spread <- as.data.frame(cbind(processed_df_LV$mean_df_LV, processed_df_LV$lower_df_LV, processed_df_LV$upper_df_LV, processed_df_LV$pval_df_LV,
                                                processed_df_LO$mean_df_LO, processed_df_LO$lower_df_LO, processed_df_LO$upper_df_LO, processed_df_LO$pval_df_LO))
colnames(to_plot_dispersal_spread) <- c("mean_df_LV", "lower_df_LV", "upper_df_LV", "ttest_df_LV",
                                        "mean_df_LO", "lower_df_LO", "upper_df_LO", "ttest_df_LO")
to_plot_dispersal_spread <- to_plot_dispersal_spread %>%  mutate_all(as.numeric) 
to_plot_dispersal_spread$Quadrant <- c("1", "2+3", "4")
condition <- c(paste(vehicle), paste(vehicle, "+", odorant))

## plot dispersal spread

plot_dispersal_spread <- ggplot(to_plot_dispersal_spread, aes(x=Quadrant, group = 1)) + 
  geom_point(aes(y=mean_df_LV, colour=paste(vehicle), size=4)) +
  geom_line(aes(y=mean_df_LV, colour=paste(vehicle), size=2)) +
  geom_errorbar(aes(ymin=lower_df_LV, ymax=upper_df_LV), size=1, width=.25,  alpha = 1.0, color = c(errorbars)) +
  geom_point(aes(y=mean_df_LO, colour=paste(odorant), size=4)) +
  geom_line(aes(y=mean_df_LO, colour=paste(odorant), size=2)) +
  geom_errorbar(aes(ymin=lower_df_LO, ymax=upper_df_LO), size=1, width=.25,  alpha = 1.0, color = c(errorbars)) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, .5, 1)) +
  my_theme + 
  theme(legend.position = c("top"), 
        legend.justification = c("left", "top"), 
        legend.box.just = "right", 
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1.0, "cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA)) +
  guides(size = "none") +
  guides(colour = guide_legend(override.aes = list(size=6)), hjust = 0, vjust = 0) + 
  scale_color_manual(values = c(V_dispcolor, O_dispcolor)) +
  labs(title = paste("dispersal of", strain), x = "quadrant", y = "mean proportion \nof worms", color = "cond:", size = "")
print(plot_dispersal_spread)

## processing to plot mean dispersal
create_df_LV <- function(df_LV, n) {
  return(df_LV[, n])} ## function to create a new dataframe with the first n columns
for (i in 1:ncol(df_LV)) {
  new_df_LV <- create_df_LV(df_LV, i)
  assign(paste0("df_LV", i), new_df_LV)} ## loop to create dataframes with 1 to N columns
calculate_entropy_LV <- function(df_LV) {
  entropy.empirical(df_LV, unit = c("log2"))} ## function to calculate entropy of a dataframe
entropy_values_LV <- c()
for (i in 1:ncol(df_LV)) {
  df_LV_name <- paste0("df_LV", i)
  entropy_value_LV <- calculate_entropy_LV(get(df_LV_name))
  entropy_values_LV <- c(entropy_values_LV, entropy_value_LV)} ## loop to calculate entropy for each dataframe
dispersal_LV <- print(entropy_values_LV) ## print the entropy values

create_df_LO <- function(df_LO, n) {
  return(df_LO[, n])} ## function to create a new dataframe with the first n columns
for (i in 1:ncol(df_LO)) {
  new_df_LO <- create_df_LO(df_LO, i)
  assign(paste0("df_LO", i), new_df_LO)} ## loop to create dataframes with 1 to N columns
calculate_entropy_LO <- function(df_LO) {
  entropy.empirical(df_LO, unit = c("log2"))} ## function to calculate entropy of a dataframe
entropy_values_LO <- c()
for (i in 1:ncol(df_LO)) {
  df_LO_name <- paste0("df_LO", i)
  entropy_value_LO <- calculate_entropy_LO(get(df_LO_name))
  entropy_values_LO <- c(entropy_values_LO, entropy_value_LO)} ## loop to calculate entropy for each dataframe
dispersal_LO <- print(entropy_values_LO) ## print the entropy values

dispersal_CI_df_LV <- as.data.frame(t(as.data.frame(t.test(dispersal_LV)$conf.int))) ## calculating confidence interval of mean dispersal
dispersal_CI_df_LO <- as.data.frame(t(as.data.frame(t.test(dispersal_LO)$conf.int))) ## calculating confidence interval of mean dispersal

mean_dispersal_LV <- mean(dispersal_LV)
dispersal_V <- as.data.frame(mean(dispersal_LV))
colnames(dispersal_V) <- c("mean")

mean_dispersal_LO <- mean(dispersal_LO)
dispersal_O <- as.data.frame(mean(dispersal_LO))
colnames(dispersal_O) <- c("mean")
dispersal_error <- rbind(as.data.frame(dispersal_CI_df_LV),
                         as.data.frame(dispersal_CI_df_LO))

dispersal <- rbind(as.data.frame(dispersal_V), as.data.frame(dispersal_O))
condition <- c(paste(vehicle), paste(vehicle, "+", odorant))
dispersal_plot <- cbind(as.data.frame(dispersal), as.data.frame(condition), as.data.frame(dispersal_error))

## plot dispersal bar

plot_dispersal_bar <- 
  ggplot(dispersal_plot, aes(x=fct_inorder(condition), y=mean)) +
  geom_bar(stat = "identity", fill = c(V_dispcolor, O_dispcolor)) +
  geom_errorbar(aes(ymin=V1, ymax=V2), size=1, width=.25, color = c(errorbars)) +
  scale_y_continuous(limits=c(0, 2), breaks=c(0, 1, 2)) +
  my_theme +
  theme(plot.title = element_text(family = "Helvetica", colour = "#000000", size = 10)) +
  theme(legend.position = c("top"), legend.justification = c("centre", "top"), legend.box.just = "right", legend.margin = margin(0, 0, 0, 0)) + guides(size = "none") + guides(colour = guide_legend(override.aes = list(size=1))) + 
  labs(title = paste("dispersal of", strain), 
       x = "condition", y = "dispersal")
print(plot_dispersal_bar)

## processing df_CV (data processing for final plot(s))
df_CV <- df_CV[,colnames(df_CV)[colnames(df_CV) !="Quadrant"]] ## gets rid of first column
df_CV = df_CV[colSums(df_CV) > 100] ## picks columns/lanes where n > 100
Ndf_CV = ncol(df_CV) ## readout of N (number of columns/lanes picked)
df_CV_23ex <- as.data.frame(subset(df_CV[2:3,])) ## extracts data from rows Q2 and Q3 for processing
df_CV_23sum <-as.data.frame(t(colSums(df_CV_23ex[1:ncol(df_CV)]))) ## ## calculates sum of Q2 and Q3 for each column/lane and adds it to a dataframe
df_CV_calc <- rbind(as.data.frame(df_CV[1,]), as.data.frame(df_CV_23sum[1,]), as.data.frame(df_CV[4,])) ## dataframe on which statistical calculations will be done
row.names(df_CV_calc)<- c("Q1", "Q23", "Q4") ## renaming rows 
df_CV_calc_proportion <- as.data.frame(t(t(df_CV_calc)/colSums(df_CV_calc))) ## calculating proportions of worms in quadrants
df_CV_calc_mean <- rowMeans(df_CV_calc_proportion) ## calculating mean proportions
processed_df_CV <- cbind(df_CV_calc_proportion, 
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CV_calc_proportion))$Q1)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CV_calc_proportion))$Q23)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CV_calc_proportion))$Q4)$conf.int)))),
                         df_CV_calc_mean,
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CV_calc_proportion))$Q1, as.data.frame(t(df_CV_calc_proportion))$Q23)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CV_calc_proportion))$Q23, as.data.frame(t(df_CV_calc_proportion))$Q4)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CV_calc_proportion))$Q1, as.data.frame(t(df_CV_calc_proportion))$Q4)$p.value))))) ## binding everything together
c <- Ndf_CV
colnames(processed_df_CV) <- c(paste0("VO", 1:c), "lower_df_CV", "upper_df_CV", "mean_df_CV", "pval_df_CV") ## renaming columns
df_CV_calc_proportion_index_calc <- as.data.frame(t(df_CV_calc_proportion)) ## transposing dataframe for calculating chemotaxis index
chemotaxis_df_CV <- as.data.frame(((df_CV_calc_proportion_index_calc[3] - df_CV_calc_proportion_index_calc[1]) / (df_CV_calc_proportion_index_calc[3] + df_CV_calc_proportion_index_calc[1]))) ## calculating chemotaxis indices for each column/lane
chemotaxis_mean_df_CV <- colMeans(chemotaxis_df_CV) ## calculating mean chemotaxis index
chemotaxis_CI_df_CV <- as.data.frame(t(as.data.frame(t.test(chemotaxis_df_CV$Q4)$conf.int))) ## calculating confidence interval of mean chemotaxis index

## processing df_CO (data processing for final plot(s))
df_CO <- df_CO[,colnames(df_CO)[colnames(df_CO) !="Quadrant"]] ## gets rid of first column
df_CO = df_CO[colSums(df_CO) > 100] ## picks columns/lanes where n > 100
Ndf_CO = ncol(df_CO) ## readout of N (number of columns/lanes picked)
df_CO_23ex <- as.data.frame(subset(df_CO[2:3,])) ## extracts data from rows Q2 and Q3 for processing
df_CO_23sum <-as.data.frame(t(colSums(df_CO_23ex[1:ncol(df_CO)]))) ## ## calculates sum of Q2 and Q3 for each column/lane and adds it to a dataframe
df_CO_calc <- rbind(as.data.frame(df_CO[1,]), as.data.frame(df_CO_23sum[1,]), as.data.frame(df_CO[4,])) ## dataframe on which statistical calculations will be done
row.names(df_CO_calc)<- c("Q1", "Q23", "Q4") ## renaming rows 
df_CO_calc_proportion <- as.data.frame(t(t(df_CO_calc)/colSums(df_CO_calc))) ## calculating proportions of worms in quadrants
df_CO_calc_mean <- rowMeans(df_CO_calc_proportion) ## calculating mean proportions
processed_df_CO <- cbind(df_CO_calc_proportion, 
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CO_calc_proportion))$Q1)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CO_calc_proportion))$Q23)$conf.int))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CO_calc_proportion))$Q4)$conf.int)))),
                         df_CO_calc_mean,
                         rbind(as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CO_calc_proportion))$Q1, as.data.frame(t(df_CO_calc_proportion))$Q23)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CO_calc_proportion))$Q23, as.data.frame(t(df_CO_calc_proportion))$Q4)$p.value))),
                               as.data.frame(t(as.data.frame(t.test(as.data.frame(t(df_CO_calc_proportion))$Q1, as.data.frame(t(df_CO_calc_proportion))$Q4)$p.value))))) ## binding everything together
d <- Ndf_CO
colnames(processed_df_CO) <- c(paste0("VO", 1:d), "lower_df_CO", "upper_df_CO", "mean_df_CO", "pval_df_CO") ## renaming columns
df_CO_calc_proportion_index_calc <- as.data.frame(t(df_CO_calc_proportion)) ## transposing dataframe for calculating chemotaxis index
chemotaxis_df_CO <- as.data.frame(((df_CO_calc_proportion_index_calc[3] - df_CO_calc_proportion_index_calc[1]) / (df_CO_calc_proportion_index_calc[3] + df_CO_calc_proportion_index_calc[1]))) ## calculating chemotaxis indices for each column/lane
chemotaxis_mean_df_CO <- colMeans(chemotaxis_df_CO) ## calculating mean chemotaxis index
chemotaxis_CI_df_CO <- as.data.frame(t(as.data.frame(t.test(chemotaxis_df_CO$Q4)$conf.int))) ## calculating confidence interval of mean chemotaxis index

to_plot_chemotaxis_spread <- as.data.frame(cbind(processed_df_CV$mean_df_CV, processed_df_CV$lower_df_CV, processed_df_CV$upper_df_CV, processed_df_CV$pval_df_CV,
                                                 processed_df_CO$mean_df_CO, processed_df_CO$lower_df_CO, processed_df_CO$upper_df_CO, processed_df_CO$pval_df_CO))
colnames(to_plot_chemotaxis_spread) <- c("mean_df_CV", "lower_df_CV", "upper_df_CV", "ttest_df_CV",
                                         "mean_df_CO", "lower_df_CO", "upper_df_CO", "ttest_df_CO")
to_plot_chemotaxis_spread <- to_plot_chemotaxis_spread %>%  mutate_all(as.numeric) 
to_plot_chemotaxis_spread$Quadrant <- c("1", "2+3", "4")

## processing chemotaxis barplot dataframe(s)
chemotaxis_V <- chemotaxis_mean_df_CV
chemotaxis_O <- chemotaxis_mean_df_CO
chemotaxis_error <- rbind(as.data.frame(chemotaxis_CI_df_CV),
                          as.data.frame(chemotaxis_CI_df_CO))
chemotaxis <- c(chemotaxis_V, chemotaxis_O)
condition <- c(paste(vehicle), paste(vehicle, "+", odorant))
chemotaxis_plot <- cbind(as.data.frame(chemotaxis), as.data.frame(condition), as.data.frame(chemotaxis_error))

# effect sizes
chemotaxis_vs_dispersal_V_ES <- round(cohen.d(chemotaxis_df_CV$Q4, chemotaxis_LV$Q4)$estimate, digits = 2)
chemotaxis_vs_dispersal_O_ES <- round(cohen.d(chemotaxis_df_CO$Q4, chemotaxis_LO$Q4)$estimate, digits = 2)
chemotaxis_OvsV_ES <- round(cohen.d(chemotaxis_df_CO$Q4, chemotaxis_df_CV$Q4)$estimate, digits = 2)

## plot chemotaxis spread
plot_chemotaxis_spread <- ggplot(to_plot_chemotaxis_spread, aes(x=Quadrant, group = 1)) + 
  geom_point(aes(y=mean_df_CV, colour=paste(vehicle), size=4)) +
  geom_line(aes(y=mean_df_CV, colour=paste(vehicle), size=2)) +
  geom_errorbar(aes(ymin=lower_df_CV, ymax=upper_df_CV), size=1, width=.25,  alpha = 1.0, color = c(errorbars)) +
  geom_point(aes(y=mean_df_CO, colour=paste(vehicle,"+", odorant), size=4)) +
  geom_line(aes(y=mean_df_CO, colour=paste(vehicle,"+", odorant), size=2)) + 
  geom_errorbar(aes(ymin=lower_df_CO, ymax=upper_df_CO), size=1, width=.25,  alpha = 1.0, color = c(errorbars)) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, .5, 1)) +
  my_theme + 
  theme(legend.position = c("top"), 
        legend.justification = c("left", "top"), 
        legend.box.just = "right", 
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1.0, "cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA)) +
  guides(size = "none") +
  guides(colour = guide_legend(override.aes = list(size=6)), hjust = 0, vjust = 0) + 
  scale_color_manual(values = c(V_color, O_color)) +
  labs(title = paste("chemotaxis spread of", strain), x = "quadrant", y = "mean proportion \nof worms", color = "cond:", size = "") +
  annotate("polygon", x = c(1, 3, 3), y = c(0.9, 0.9, 1.0), alpha=0.4, fill=V_dispcolor)
print(plot_chemotaxis_spread)

## plot chemotaxis bar
plot_chemotaxis_bar <- 
  ggplot(chemotaxis_plot, aes(x=fct_inorder(condition), y=chemotaxis)) +
  geom_bar(stat = "identity", fill = c(V_color, O_color)) +
  geom_errorbar(aes(ymin=V1, ymax=V2), size=1, width=.25, color = c(errorbars)) +
  scale_y_continuous(limits=c(-1, 1), breaks=c(-1, -.5, 0, .5, 1)) +
  my_theme +
  theme(plot.title = element_text(family = "Helvetica", colour = "#000000", size = 10)) +
  theme(legend.position = c("top"), legend.justification = c("centre", "top"), legend.box.just = "right", legend.margin = margin(0, 0, 0, 0)) + guides(size = "none") + guides(colour = guide_legend(override.aes = list(size=1))) + 
  labs(title = paste("chemotaxis of", strain),
       x = "condition", y = "chemotaxis index")
print(plot_chemotaxis_bar)


## plot 
figure <- ggarrange(plot_dispersal_spread, plot_dispersal_bar,
                    plot_chemotaxis_spread, plot_chemotaxis_bar,
                    ncol = 2, nrow = 2)
print(figure)

## useful data tables [this section doesn't need to be in the final plot but can be useful for a quick look at the data]
## comment out as necessary
Date = date
Strain = strain
Control = vehicle
Odorant = odorant
Researcher = person
experiment <- data.frame(c(Date, Researcher, Strain, Control, Odorant))
colnames(experiment) <- c("info")                         
rownames(experiment) <- c("Date",
                          "Researcher",
                          "Strain",
                          "Control",
                          "Odorant")
table_experiment <- tableGrob(experiment)
all_spread_data <- as.data.frame(rbind(cbind(round(processed_df_LV$mean_df_LV, digits = 2),
                                             round(processed_df_LV$lower_df_LV, digits = 2),
                                             round(processed_df_LV$upper_df_LV, digits = 2),
                                             round(processed_df_LV$pval_df_LV, digits = 2)),
                                       cbind(round(processed_df_LO$mean_df_LO, digits = 2),
                                             round(processed_df_LO$lower_df_LO, digits = 2),
                                             round(processed_df_LO$upper_df_LO, digits = 2),
                                             round(processed_df_LO$pval_df_LO, digits = 2)),
                                       cbind(round(processed_df_CV$mean_df_CV, digits = 2),
                                             round(processed_df_CV$lower_df_CV, digits = 2),
                                             round(processed_df_CV$upper_df_CV, digits = 2),
                                             round(processed_df_CV$pval_df_CV, digits = 2)),
                                       cbind(round(processed_df_CO$mean_df_CO, digits = 2),
                                             round(processed_df_CO$lower_df_CO, digits = 2),
                                             round(processed_df_CO$upper_df_CO, digits = 2),
                                             round(processed_df_CO$pval_df_CO, digits = 2))))
colnames(all_spread_data) <- c("µ_proportion", "lower_CI", "upper_CI", "ttest")
row.names(all_spread_data)<- c("disp V Q1", "dis V Q23", "disp V Q4", 
                               "disp O Q1", "disp O Q23", "disp O Q4",
                               "chem V Q1", "chem V Q23", "chem V Q4", 
                               "chem O Q1", "chem O Q23", "chem O Q4")
table_spread <- tableGrob(all_spread_data)
all_bar_data <- cbind(rbind(round(mean_dispersal_LV, digits =2),
                            round(mean_dispersal_LO, digits =2),
                            round(chemotaxis_V, digits =2),
                            round(chemotaxis_O, digits =2)
),
rbind(round(dispersal_CI_df_LV, digits =2),
      round(dispersal_CI_df_LO, digits =2),
      round(chemotaxis_CI_df_CV, digits =2),
      round(chemotaxis_CI_df_CO, digits =2)
),
rbind (Ndf_LV, Ndf_LO, Ndf_CV, Ndf_CO))
colnames(all_bar_data) <- c("µ", "CI_low", "CI_up", "N")
row.names(all_bar_data)<- c("disp V", 
                            "disp O",
                            "chem V", 
                            "chem O")
table_bar <- tableGrob(all_bar_data)
all_ES_data <- rbind (chemotaxis_OvsV_ES, chemotaxis_vs_dispersal_V_ES, chemotaxis_vs_dispersal_O_ES)
colnames(all_ES_data) <- c("effect size")
rownames(all_ES_data) <- c("CO_vs_CV",
                           "CV_vs_LV",
                           "CO_vs_LO")
table_ES <- tableGrob(all_ES_data)
blank_plot <- ggplot() + theme_void()
tables_3 <- grid.arrange(table_experiment, table_ES, blank_plot, table_bar, ncol = 2)
print(tables_3)
tables_4 <- grid.arrange(tables_3, table_spread, ncol = 2)
print(tables_4)

## final final plot 

final_plot <- ggarrange(tables_4, figure,
                        ncol = 1, nrow = 2)
print(final_plot)

# ## export as svg
# ggsave(
#   filename = paste(date, person, strain, vehicle, odorant, ".svg"),
#   plot = last_plot(),
#   device = "svg",
#   path = "/Users/samihatasnim/Desktop/Codes_Figures/Plots/PAPER_01",
#   scale = 1,
#   width = 10,
#   height = 10,
#   units = c("in"),
#   dpi = 1000)