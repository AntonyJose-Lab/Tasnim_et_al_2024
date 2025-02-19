## this is the code as of 20250212 for visualizing data obtained from behavioral assays
## worms are pre-exposed with either the vehicle (V) or the odorant (O)
## tested for dispersal (no choice assay) and chemotaxis (V vs O choice assay)
## n ≥ 100 (N ≤ 8 according to n)

## set working directory
setwd("/Users/samihatasnim/Desktop/Codes_Figures") 

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
preX_eth_gradient_butanone = "#d55e00"
preX_but_gradient_butanone = "#faab6b"
preX_eth_gradient_benzaldehyde = "#0072b2"
preX_eth_gradient_nonanone = "#c35090"
dispcolor = "#bcbec0"
errorbars= "#000000"

## set up data frame(s) and input relevant information
date = "20230523"
strain = "N2"
choice_disp = "none"
choice_chem = "eth vs but"
vehicle = "eth"
odorant = "but"
orientation = "horizontal"
person = "Samiha"
errorbars= "#000000"
dispcolor= "#bcbec0"
O_color= "#d55e00"

dataframe1 = "N2 EtOH dispersal"
dataframe2 = "N2 but dispersal"
dataframe3 = "N2 EtOH chemotaxis n<100"
dataframe4 = "N2 but chemotaxis"

df_LV <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/202502 EtOH_disp.xlsx")) ## df_LV is dispersal counts post vehicle pre-exposure
df_LO <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/202502 Odor_disp.xlsx")) ## df_LO is dispersal counts post odorant pre-exposure
df_CV <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/202502 EtOH_chem.xlsx")) ## df_CV is chemotaxis counts post vehicle pre-exposure
df_CO <- data.frame(read_excel ("/Users/samihatasnim/Desktop/Codes_Figures/Data/202502 Odor_chem.xlsx")) ## df_CO is chemotaxis counts post odorant pre-exposure

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
df_LV_calc_proportion_index_calc <- as.data.frame(t(df_LV_calc_proportion)) ## transposing dataframe for R to calculate dispersal index

chemotaxis_LV <- as.data.frame(((df_LV_calc_proportion_index_calc[3] - df_LV_calc_proportion_index_calc[1]) / (df_LV_calc_proportion_index_calc[3] + df_LV_calc_proportion_index_calc[1]))) ## calculating chemotaxis indices for each column/lane
chemotaxis_mean_LV <- colMeans(chemotaxis_LV) ## calculating chemotaxis indices for each column/lane
chemotaxis_error_LV <- as.data.frame(t(as.data.frame(t.test(chemotaxis_LV$Q4)$conf.int)))


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
df_LO_calc_proportion_index_calc <- as.data.frame(t(df_LO_calc_proportion)) ## transposing dataframe for R to calculate locomotion index

chemotaxis_LO <- as.data.frame(((df_LO_calc_proportion_index_calc[3] - df_LO_calc_proportion_index_calc[1]) / (df_LO_calc_proportion_index_calc[3] + df_LO_calc_proportion_index_calc[1]))) ## calculating chemotaxis indices for each column/lane
chemotaxis_mean_LO <- colMeans(chemotaxis_LO) ## calculating chemotaxis indices for each column/lane
chemotaxis_error_LO <- as.data.frame(t(as.data.frame(t.test(chemotaxis_LO$Q4)$conf.int)))

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
df_CV_calc_proportion_index_calc <- as.data.frame(t(df_CV_calc_proportion)) ## transposing dataframe for R to calculate chemotaxis index

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
df_CO_calc_proportion_index_calc <- as.data.frame(t(df_CO_calc_proportion)) ## transposing dataframe for R to calculate chemotaxis index

chemotaxis_df_CO <- as.data.frame(((df_CO_calc_proportion_index_calc[3] - df_CO_calc_proportion_index_calc[1]) / (df_CO_calc_proportion_index_calc[3] + df_CO_calc_proportion_index_calc[1]))) ## calculating chemotaxis indices for each column/lane
chemotaxis_mean_df_CO <- colMeans(chemotaxis_df_CO) ## calculating mean chemotaxis index
chemotaxis_CI_df_CO <- as.data.frame(t(as.data.frame(t.test(chemotaxis_df_CO$Q4)$conf.int))) ## calculating confidence interval of mean chemotaxis index

## processing to plot spread as mean proportions in quadrants

to_plot_V_spread <- as.data.frame(cbind(processed_df_LV$mean_df_LV, processed_df_LV$lower_df_LV, processed_df_LV$upper_df_LV,
                                        processed_df_CV$mean_df_CV, processed_df_CV$lower_df_CV, processed_df_CV$upper_df_CV))
colnames(to_plot_V_spread) <- c("mean_df_LV", "lower_df_LV", "upper_df_LV",
                                "mean_df_CV", "lower_df_CV", "upper_df_CV")
to_plot_V_spread <- to_plot_V_spread %>%  mutate_all(as.numeric)
to_plot_V_spread$Quadrant <- c("1", "2+3", "4")

to_plot_O_spread <- as.data.frame(cbind(processed_df_LO$mean_df_LO, processed_df_LO$lower_df_LO, processed_df_LO$upper_df_LO,
                                        processed_df_CO$mean_df_CO, processed_df_CO$lower_df_CO, processed_df_CO$upper_df_CO))
colnames(to_plot_O_spread) <- c("mean_df_LO", "lower_df_LO", "upper_df_LO",
                                "mean_df_CO", "lower_df_CO", "upper_df_CO")
to_plot_O_spread <- to_plot_O_spread %>%  mutate_all(as.numeric)
to_plot_O_spread$Quadrant <- c("1", "2+3", "4")

Nchem <- c(Ndf_CV, Ndf_LV, Ndf_CO, Ndf_LO)

## plot dispersal w/vehicle pre-exposure
plot_vehicle_spread <- ggplot(to_plot_V_spread, aes(x=Quadrant, group = 1)) + 
  geom_point(aes(y=mean_df_LV, colour=paste("none"), size=4)) +
  geom_line(aes(y=mean_df_LV, colour=paste("none"), size=2)) +
  geom_errorbar(aes(ymin=pmax(lower_df_LV, 0), ymax=pmin(upper_df_LV,1)), size=1, width=.25,  alpha = 1.0, color = c(dispcolor)) +
  geom_point(aes(y=mean_df_CV, colour=paste(odorant), size=4)) +
  geom_line(aes(y=mean_df_CV, colour=paste(odorant), size=2)) +
  geom_errorbar(aes(ymin=pmax(lower_df_CV, 0), ymax=pmin(upper_df_CV, 1)), size=1, width=.25,  alpha = 1.0, color = c(O_color)) +
  scale_y_continuous(limits=c(0, 1.2), breaks=c(0, 1)) +
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
  scale_color_manual(values = c(O_color, dispcolor)) +
  labs(title = paste("spread of", strain, vehicle, "pre-exposure"), x = "quadrant", y = "µPw", color = "choice:", size = "") +
  annotate("polygon", x = c(1, 3, 3), y = c(1.1, 1.1, 1.2), alpha=1, fill=O_color)
print(plot_vehicle_spread)

## plot dispersal w/odorant pre-exposure
plot_odorant_spread <- ggplot(to_plot_O_spread, aes(x=Quadrant, group = 1)) + 
  geom_point(aes(y=mean_df_LO, colour=paste("none"), size=4)) +
  geom_line(aes(y=mean_df_LO, colour=paste("none"), size=2)) +
  geom_errorbar(aes(ymin=pmax(lower_df_LO, 0), ymax=pmin(upper_df_LO,1)), size=1, width=.25,  alpha = 1.0, color = c(dispcolor)) +
  geom_point(aes(y=mean_df_CO, colour=paste(odorant), size=4)) +
  geom_line(aes(y=mean_df_CO, colour=paste(odorant), size=2)) +
  geom_errorbar(aes(ymin=pmax(lower_df_CO, 0), ymax=pmin(upper_df_CO, 1)), size=1, width=.25,  alpha = 1.0, color = c(O_color)) +
  scale_y_continuous(limits=c(0, 1.2), breaks=c(0, 1)) +
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
  scale_color_manual(values = c(O_color, dispcolor)) +
  labs(title = paste("spread of", strain, odorant, "pre-exposure"), x = "quadrant", y = "µPw", color = "choice:", size = "") +
  annotate("polygon", x = c(1, 3, 3), y = c(1.1, 1.1, 1.2), alpha=1, fill=O_color)
print(plot_odorant_spread)

## processing to plot mean dispersal
calculate_entropy <- function(x) {entropy.empirical(x, unit = c("log2"))} # function to calculate Shannon entropy

dispersal_df_CV <- data.frame(Column = names(df_CV), Entropy = sapply(df_CV, calculate_entropy)) # calculate entropy for each column and create a new dataframe
dispersal_CI_df_CV <- t.test(dispersal_df_CV$Entropy)$conf.int ## calculating confidence interval of mean dispersal
dispersal_CV <- as.data.frame(mean(dispersal_df_CV$Entropy))
colnames(dispersal_CV) <- c("mean")

dispersal_df_LV <- data.frame(Column = names(df_LV), Entropy = sapply(df_LV, calculate_entropy)) # calculate entropy for each column and create a new dataframe
dispersal_CI_df_LV <- t.test(dispersal_df_LV$Entropy)$conf.int ## calculating confidence interval of mean dispersal

dispersal_LV <- as.data.frame(mean(dispersal_df_LV$Entropy))
colnames(dispersal_LV) <- c("mean")

dispersal_df_CO <- data.frame(Column = names(df_CO), Entropy = sapply(df_CO, calculate_entropy)) # calculate entropy for each column and create a new dataframe
dispersal_CI_df_CO <- t.test(dispersal_df_CO$Entropy)$conf.int ## calculating confidence interval of mean dispersal
dispersal_CO <- as.data.frame(mean(dispersal_df_CO$Entropy))
colnames(dispersal_CO) <- c("mean")

dispersal_df_LO <- data.frame(Column = names(df_LO), Entropy = sapply(df_LO, calculate_entropy)) # calculate entropy for each column and create a new dataframe
dispersal_CI_df_LO <- t.test(dispersal_df_LO$Entropy)$conf.int ## calculating confidence interval of mean dispersal
dispersal_LO <- as.data.frame(mean(dispersal_df_LO$Entropy))
colnames(dispersal_LO) <- c("mean")

dispersal <- rbind(as.data.frame(dispersal_CV),
                   as.data.frame(dispersal_LV), 
                   as.data.frame(dispersal_CO),
                   as.data.frame(dispersal_LO))
condition_dispersal <- c(paste(vehicle, "chem "), paste(vehicle, "disp"), paste(odorant, "chem "), paste(odorant, "disp"))

dispersal_error <- rbind(dispersal_CI_df_CV,
                         dispersal_CI_df_LV,
                         dispersal_CI_df_CO,
                         dispersal_CI_df_LO)

dispersal_plot <- cbind(as.data.frame(dispersal), as.data.frame(condition_dispersal), as.data.frame(dispersal_error))

Ndisp <- c(Ndf_CV, Ndf_LV, Ndf_CO, Ndf_LO)

V_dispersal_stat <-signif(t.test(as.data.frame(dispersal_df_LV$Entropy), as.data.frame(dispersal_df_CV$Entropy))$p.value, digits = 3)
V_dispersal_ES <- signif(cohen.d(dispersal_df_LV$Entropy, dispersal_df_CV$Entropy)$estimate, digits = 3)

O_dispersal_stat <-signif(t.test(as.data.frame(dispersal_df_LO$Entropy), as.data.frame(dispersal_df_CO$Entropy))$p.value, digits = 3)
O_dispersal_ES <- signif(cohen.d(dispersal_df_LO$Entropy, dispersal_df_CO$Entropy)$estimate, digits = 3)

## plot dispersal bar

plot_dispersal_bar <- 
  ggplot(dispersal_plot, aes(x=fct_inorder(condition_dispersal), y=mean)) +
  geom_bar(stat = "identity", fill = c(dispcolor, dispcolor, dispcolor, dispcolor), width = 0.25) +
  geom_errorbar(aes(ymin=pmax(V1, 0), ymax=pmin(V2, 2)), size=1, width=.1, color = c(dispcolor)) +
  
  geom_text(aes(x = -Inf, y = 3, label = paste0("N(n>100)=")), vjust = 1, color = "red") +
  geom_text(aes(label = paste0(Ndisp)), y = 3, vjust = 0.5, color = "red") +
  
  geom_text(aes(x = -Inf, y = 2.5, label = paste0("pval", "\n", "ES")), vjust = 1, color = "red") +
  geom_text(aes(x = mean(c(1,2)), y = 2.5, label = paste0(V_dispersal_stat, "\n", V_dispersal_ES)), vjust = 1, color = "red") +
  geom_text(aes(x = mean(c(3,4)), y = 2.5, label = paste0(O_dispersal_stat, "\n", O_dispersal_ES)), vjust = 1, color = "red") +
  
  coord_cartesian(clip = "off") +
  
  scale_y_continuous(limits=c(0, 3), breaks=c(0, 2)) +
  my_theme +
  theme(plot.title = element_text(family = "Helvetica", colour = "#000000", size = 10)) +
  theme(legend.position = c("top"), legend.justification = c("centre", "top"), legend.box.just = "right", legend.margin = margin(0, 0, 0, 0)) + guides(size = "none") + guides(colour = guide_legend(override.aes = list(size=1))) + 
  labs(title = paste("dispersal of", strain), 
       x = "test", y = "dispersal")
print(plot_dispersal_bar)

## processing chemotaxis barplot dataframe(s)

chemotaxis_V <- chemotaxis_mean_df_CV
chemotaxis_O <- chemotaxis_mean_df_CO
chemotaxis_error <- rbind(as.data.frame(chemotaxis_CI_df_CV),
                          as.data.frame(chemotaxis_error_LV),
                          as.data.frame(chemotaxis_CI_df_CO),
                          as.data.frame(chemotaxis_error_LV))
chemotaxis <- c(chemotaxis_V, chemotaxis_mean_LV, chemotaxis_O, chemotaxis_mean_LV)
condition_chemotaxis <- c(paste(vehicle, "chem"), paste(vehicle, "disp"), paste(odorant, "chem"), paste(odorant, "disp"))


rratio_df_LV <- as.data.frame(((df_LV_calc_proportion_index_calc[3] + df_LV_calc_proportion_index_calc[1]) / (df_LV_calc_proportion_index_calc[3] + df_LV_calc_proportion_index_calc[2] + df_LV_calc_proportion_index_calc[1]))) ## calculating response ratio for each LOlumn/lane
rratio_mean_df_LV <- colMeans(rratio_df_LV) ## calculating mean response ratio
rratio_CI_df_LV <- as.data.frame(t(as.data.frame(t.test(rratio_df_LV$Q4)$conf.int))) ## calculating LOnfidence interval of mean chemotaxis index
rratio_df_LO <- as.data.frame(((df_LO_calc_proportion_index_calc[3] + df_LO_calc_proportion_index_calc[1]) / (df_LO_calc_proportion_index_calc[3] + df_LO_calc_proportion_index_calc[2] + df_LO_calc_proportion_index_calc[1]))) ## calculating response ratio for each LOlumn/lane
rratio_CI_df_LO <- as.data.frame(t(as.data.frame(t.test(rratio_df_LO$Q4)$conf.int))) ## calculating LOnfidence interval of mean chemotaxis index
rratio_mean_df_LO <- colMeans(rratio_df_LO) ## calculating mean response ratio
rratio_df_CV <- as.data.frame(((df_CV_calc_proportion_index_calc[3] + df_CV_calc_proportion_index_calc[1]) / (df_CV_calc_proportion_index_calc[3] + df_CV_calc_proportion_index_calc[2] + df_CV_calc_proportion_index_calc[1]))) ## calculating response ratio for each column/lane
rratio_mean_df_CV <- colMeans(rratio_df_CV) ## calculating mean response ratio
rratio_CI_df_CV <- as.data.frame(t(as.data.frame(t.test(rratio_df_CV$Q4)$conf.int))) ## calculating confidence interval of mean chemotaxis index
rratio_df_CO <- as.data.frame(((df_CO_calc_proportion_index_calc[3] + df_CO_calc_proportion_index_calc[1]) / (df_CO_calc_proportion_index_calc[3] + df_CO_calc_proportion_index_calc[2] + df_CO_calc_proportion_index_calc[1]))) ## calculating response ratio for each column/lane
rratio_CI_df_CO <- as.data.frame(t(as.data.frame(t.test(rratio_df_CO$Q4)$conf.int))) ## calculating confidence interval of mean chemotaxis index
rratio_mean_df_CO <- colMeans(rratio_df_CO) ## calculating mean response ratio

response_ratio <- c(signif(rratio_mean_df_CV, digits = 2), signif(rratio_mean_df_LV, digits = 2), signif(rratio_mean_df_CO, digits = 2), signif(rratio_mean_df_LO, digits = 2))
rratio_error <- rbind(as.data.frame(rratio_CI_df_CV),
                      as.data.frame(rratio_CI_df_LV),
                      as.data.frame(rratio_CI_df_CO),
                      as.data.frame(rratio_CI_df_LO))

to_plot_chemotaxis_response_bars <- cbind(as.data.frame(chemotaxis), as.data.frame(condition_chemotaxis), as.data.frame(chemotaxis_error))

colnames(to_plot_chemotaxis_response_bars) <- c("chemotaxis", "condition",
                                                "lower_chem", "upper_chem")

# other indices and ratios
chemotaxis_vs_dispersal_V_signif <- signif((t.test(chemotaxis_df_CV$Q4, chemotaxis_LV$Q4))$p.value, digits = 3)
chemotaxis_vs_dispersal_O_signif <- signif((t.test(chemotaxis_df_CO$Q4, chemotaxis_LO$Q4))$p.value, digits = 3)

chemotaxis_vs_dispersal_V_ES <- signif(cohen.d(chemotaxis_df_CV$Q4, chemotaxis_LV$Q4)$estimate, digits = 3)
chemotaxis_vs_dispersal_O_ES <- signif(cohen.d(chemotaxis_df_CO$Q4, chemotaxis_LO$Q4)$estimate, digits = 3)

significance_dispersal_chemotaxis <- c(chemotaxis_vs_dispersal_V_signif, chemotaxis_vs_dispersal_O_signif)
effectsize_dispersal_chemotaxis <- c(chemotaxis_vs_dispersal_V_ES, chemotaxis_vs_dispersal_O_ES)

chemotaxis_OvsV_signif <- signif(t.test(chemotaxis_df_CO$Q4, chemotaxis_df_CV$Q4)$p.value, digits = 3)
chemotaxis_OvsV_ES <- signif(cohen.d(chemotaxis_df_CO$Q4, chemotaxis_df_CV$Q4)$estimate, digits = 3)

signif_table <- c(chemotaxis_OvsV_signif, significance_dispersal_chemotaxis)

nc_df_LV <- as.data.frame(as.data.frame(t(df_LV_calc))[3] + as.data.frame(t(df_LV_calc))[1]) ## calculating responding number for each LOlumn/lane
nc_mean_df_LV <- signif(colMeans(nc_df_LV), digits = 3) ## calculating mean responding number
nc_df_LO <- as.data.frame(as.data.frame(t(df_LO_calc))[3] + as.data.frame(t(df_LO_calc))[1]) ## calculating responding number for each LOlumn/lane
nc_mean_df_LO <- signif(colMeans(nc_df_LO), digits = 3) ## calculating mean responding number
nc_df_CV <- as.data.frame(as.data.frame(t(df_CV_calc))[3] + as.data.frame(t(df_CV_calc))[1]) ## calculating responding number for each column/lane
nc_mean_df_CV <- signif(colMeans(nc_df_CV), digits = 3) ## calculating mean responding number
nc_df_CO <- as.data.frame(as.data.frame(t(df_CO_calc))[3] + as.data.frame(t(df_CO_calc))[1]) ## calculating responding number for each column/lane
nc_mean_df_CO <- signif(colMeans(nc_df_CO), digits = 3) ## calculating mean responding number

nc <- c(nc_mean_df_CV, nc_mean_df_LV, nc_mean_df_CO, nc_mean_df_LO)

## plot chemotaxis bar
plot_chemotaxis_bar <- 
  ggplot(to_plot_chemotaxis_response_bars, aes(x=fct_inorder(condition))) +
  geom_bar(aes(y = chemotaxis), stat = "identity", fill = c(O_color, dispcolor, O_color, dispcolor), width = 0.25) +
  
  geom_text(aes(x = -Inf, y = 6, label = paste0("signif", "\n", "ES")), vjust = 1, color = "red") +
  geom_text(aes(x = mean(c(1,3)), y = 6, label = paste0(chemotaxis_OvsV_signif, "\n", chemotaxis_OvsV_ES)), vjust = 1, color = "red") +
  
  geom_text(aes(x = -Inf, y = 5, label = paste0("signif", "\n", "ES")), vjust = 1, color = "red") +
  geom_text(aes(x = mean(c(1,2)), y = 5, label = paste0(chemotaxis_vs_dispersal_V_signif, "\n", chemotaxis_vs_dispersal_V_ES)), vjust = 1, color = "red") +
  geom_text(aes(x = mean(c(3,4)), y = 5, label = paste0(chemotaxis_vs_dispersal_O_signif, "\n", chemotaxis_vs_dispersal_O_ES)), vjust = 1, color = "red") +
  
  geom_text(x = -Inf, y = 4, aes(label = paste0("N(n>100)=", "\n", "nchem=", "\n", "Dispersal=", "\n", "RRatio=", "\n", "CI=")), vjust = 1, color = "red") +
  geom_text(y = 4, aes(label = paste0(Nchem, "\n", nc, "\n", signif(dispersal_plot$mean, digits = 2), "\n", response_ratio, "\n", signif(to_plot_chemotaxis_response_bars$chemotaxis, digits = 2))), vjust = 1, color = "red") +
  
  
  coord_cartesian(clip = "off") +
  
  geom_errorbar(aes(ymin=pmax(lower_chem, -1) , ymax=pmin(upper_chem, 1)), size=1, width=.1, color = c(O_color, dispcolor, O_color, dispcolor)) +
  scale_y_continuous(limits=c(-1, 6), breaks=c(-1, 0, 1)) +
  my_theme +
  theme(plot.title = element_text(family = "Helvetica", colour = "#000000", size = 10)) +
  theme(legend.position = c("top"), legend.justification = c("centre", "top"), legend.box.just = "right", legend.margin = margin(0, 0, 0, 0)) + guides(size = "none") + guides(colour = guide_legend(override.aes = list(size=1))) + 
  labs(title = paste("chemotaxis of", strain),
       x = "test", y = "chemotaxis index")
print(plot_chemotaxis_bar)

## plot 
figure1 <- ggarrange(plot_vehicle_spread, plot_odorant_spread,
                     ncol = 2, nrow = 1)
figure <- ggarrange(figure1,
                    plot_dispersal_bar,
                    plot_chemotaxis_bar,
                    ncol = 1, nrow = 3,
                    heights = c(1.5, 2, 2))
print(figure)


dispersal_range <- cbind(dispersal_plot$mean, (((dispersal_plot$V2)-(dispersal_plot$V1))/2) )
colnames(dispersal_range) <- c("µ", "±")
row.names(dispersal_range) <- c("dispersal_CI_df_CV",
                                "dispersal_CI_df_LV",
                                "dispersal_CI_df_CO",
                                "dispersal_CI_df_LO")
rratio_range <- cbind(rbind(rratio_mean_df_CV,
                            rratio_mean_df_LV,
                            rratio_mean_df_CO,
                            rratio_mean_df_LO),
                      (((rratio_error$V2)-(rratio_error$V1))/2) )
colnames(rratio_range) <- c("µ", "±")



## data table(s) for this plot [this portion doesn't need to be in the final plot but can be useful for a quick look at the data]
## comment out when necessary
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

all_spread_data <- as.data.frame(rbind(cbind(signif(processed_df_LV$mean_df_LV, digits = 4),
                                             signif(processed_df_LV$lower_df_LV, digits = 4),
                                             signif(processed_df_LV$upper_df_LV, digits = 4),
                                             signif(processed_df_LV$pval_df_LV, digits = 4)),
                                       cbind(signif(processed_df_LO$mean_df_LO, digits = 4),
                                             signif(processed_df_LO$lower_df_LO, digits = 4),
                                             signif(processed_df_LO$upper_df_LO, digits = 4),
                                             signif(processed_df_LO$pval_df_LO, digits = 4)),
                                       cbind(signif(processed_df_CV$mean_df_CV, digits = 4),
                                             signif(processed_df_CV$lower_df_CV, digits = 4),
                                             signif(processed_df_CV$upper_df_CV, digits = 4),
                                             signif(processed_df_CV$pval_df_CV, digits = 4)),
                                       cbind(signif(processed_df_CO$mean_df_CO, digits = 4),
                                             signif(processed_df_CO$lower_df_CO, digits = 4),
                                             signif(processed_df_CO$upper_df_CO, digits = 4),
                                             signif(processed_df_CO$pval_df_CO, digits = 4))))
colnames(all_spread_data) <- c("µ", "lower_CI", "upper_CI", "pvalue")
row.names(all_spread_data) <- c(paste(dataframe1, "Q1"), paste(dataframe1, "Q23"), paste(dataframe1, "Q4"),
                               paste(dataframe2, "Q1"), paste(dataframe2, "Q23"), paste(dataframe2, "Q4"),
                               paste(dataframe3, "Q1"), paste(dataframe3, "Q23"), paste(dataframe3, "Q4"), 
                               paste(dataframe4, "Q1"), paste(dataframe4, "Q23"), paste(dataframe4, "Q4"))
table_spread <- tableGrob(all_spread_data)

all_bar_data <- cbind(rbind(signif(dispersal_LV, digits =4),
                            signif(dispersal_LO, digits =4),
                            signif(dispersal_CV, digits =4),
                            signif(dispersal_CO, digits =4),
                            
                            signif(chemotaxis_mean_LV, digits =4),
                            signif(chemotaxis_mean_LO, digits =4),
                            signif(chemotaxis_V, digits =4),
                            signif(chemotaxis_O, digits =4)),
                      
                      rbind(signif(dispersal_CI_df_LV, digits =4),
                            signif(dispersal_CI_df_LO, digits =4),
                            signif(dispersal_CI_df_CV, digits =4),
                            signif(dispersal_CI_df_CO, digits =4),
                            
                            signif(chemotaxis_error_LV, digits =4),
                            signif(chemotaxis_error_LO, digits =4),
                            signif(chemotaxis_CI_df_CV, digits =4),
                            signif(chemotaxis_CI_df_CO, digits =4)),
                      rbind(Ndf_LV, Ndf_LO, Ndf_CV, Ndf_CO))
colnames(all_bar_data) <- c("µ", "CI_low", "CI_up", "N")
row.names(all_bar_data)<- c(paste("dispersal", dataframe1),
                            paste("dispersal", dataframe2),
                            paste("dispersal", dataframe3),
                            paste("dispersal", dataframe4),
                            paste("chemotaxis", dataframe1),
                            paste("chemotaxis", dataframe2),
                            paste("chemotaxis", dataframe3),
                            paste("chemotaxis", dataframe4))
table_bar <- tableGrob(all_bar_data)

tables_3 <- grid.arrange(table_experiment, table_spread, table_bar, ncol = 1)
print(tables_3)

## final final plot 

final_plot <- ggarrange(tables_3, figure,
                        ncol = 1,
                        heights = c(1, 1))
print(final_plot)

# ## export as svg
# timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
# ggsave(
#   filename = paste(timestamp, date, person, strain, vehicle, odorant, ".svg"),
#   plot = last_plot(),
#   device = "svg",
#   path = "/Users/samihatasnim/Desktop/Codes_Figures/Plots/PAPER_01",
#   scale = 1,
#   width = 6.5,
#   height = 20,
#   units = c("in"),
#   dpi = 1000)
