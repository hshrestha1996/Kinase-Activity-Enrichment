#setwd("../")
getwd()
library(KSEAapp)

database <- read.csv("Kinase_Substrate_Dataset.csv")
phosphodata_LBD <- read.csv("KSEA_LBDvsCtl_p0.05_2sd.csv")
phosphodata_AD <- read.csv("KSEA_ADvsCtl_p0.05_2sd.csv")

KSEA.Complete(KSData = database, PX = phosphodata_LBD, NetworKIN=TRUE, NetworKIN.cutoff=1, m.cutoff=2, p.cutoff=0.05)
KSEA.Complete(KSData = database, PX = phosphodata_AD, NetworKIN=TRUE, NetworKIN.cutoff=1, m.cutoff=2, p.cutoff=0.1)


KSEA.Scores.1 <- read.csv("KSEA Kinase Scores_lbd.csv")
KSEA.Scores.2<- read.csv("KSEA Kinase ScoresAD.csv")

KSEA.Heatmap(score.list=list(KSEA.Scores.1, KSEA.Scores.2), 
             sample.labels=c("Kinases activity in LBD", "Kinases activity in AD"), 
             stats="p.value", m.cutoff=2, p.cutoff=0.05, sample.cluster=TRUE)

##Enrichment plot for kinases passing 0.05 and subunit threshold##
KSEA_lbd <- KSEA.Scores.1[KSEA.Scores.1$m>1,]
KSEA_lbd <- KSEA_lbd[KSEA_lbd$p.value<0.05,]
KSEA_lbd <- KSEA_lbd[order(KSEA_lbd$z.score,decreasing = TRUE),]


KSEA_ad <- KSEA.Scores.2[KSEA.Scores.2$m>1,]
KSEA_ad <- KSEA_ad[KSEA_ad$p.value<0.05,]
KSEA_ad <- KSEA_ad[order(KSEA_ad$z.score,decreasing = TRUE),]


df <- KSEA_ad
df$`-log10(p-value)` <- -log10(df$p.value)

library(ggplot2)
library(RColorBrewer)
library(scales)

# Define a reversed color palette
palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))


bwr_palette <- colorRampPalette(c("blue", "white", "red"))
n_colors <- 100
palette <- bwr_palette(n_colors)

# Normalize the color mapping to the range -6 to 6
norm_range <- c(-3, 3)

# Ensure Kinase.Gene is a factor with levels in the order of df
df$Kinase.Gene <- factor(df$Kinase.Gene, levels = df$Kinase.Gene)
ggplot(df, aes(x = `-log10(p-value)`, y = Kinase.Gene, fill = Enrichment)) +
  geom_bar(stat = "identity", orientation = "y", color = "black", size = 0.7) +  # Add border around bars
  scale_fill_gradientn(colors = palette, limits = norm_range, name = "Enrichment score") +
  labs(x = "-log10(p-value)", y = "Kinases", title = "Kinase Activity Enrichment") +
  theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90,color = "black"),
    axis.text.y = element_text(size = 12, angle = 0, hjust = 0,color = "black"),  # No rotation for y-axis labels
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines for cleaner look
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot area
    axis.ticks = element_line(color = "black", size = 0.5)
  ) +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(axis.text.y = element_text(size = 12, angle = 0, hjust = 0, vjust = 0.5))
ggsave("enriched_kinase_activity_ADvsCtl.png", dpi = 300, width = 4, height = 5)
dev.off()
