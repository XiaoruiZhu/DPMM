
load(file = "data/ALL_D.Rdata")

Test1 <- ALLD$Results
DPMM_Data <- ALLD$Data
y <- DPMM_Data[, c("X1", "X2")]
Sim_M <-  floor(exp(c(0:8, 8.5, 9.5, 10)))
Test1$NTables

# Draw sequential scatter plot --------------------------------------------

library(ggplot2); require("ggrepel")
clusters <- as.data.frame(DPMM_Data, row.names = T)
clusters[,2] <- as.factor(clusters[,2])
centers <- as.data.frame(DPMM_Data[,5:6], row.names = T)
colnames(centers) <- c("dim1","dim2")
p3 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID), size=0.001) + 
  geom_point()
p3 + geom_point(data=centers, aes(x=dim1,y=dim2), size=2, 
                inherit.aes = FALSE) + 
  ggtitle("True data with clusters") +
  labs(x = "X1") + labs(y = "X2")

Pred_center2 <- as.data.frame(Test1$All_theta[,,7], row.names = T)
colnames(Pred_center2) <- c("dim1","dim2")
Table_ID_Est <- as.factor(Pred_center2[,1])
# Rename all levels
levels(Table_ID_Est) <- c(1:length(unique(Table_ID_Est)))

# DPMM_ 2 -----------------------------------------------------------------

library(extrafont)
loadfonts(device = "pdf")

setwd("Example")
require(ggplot2)
png(file="2D_gibbs%02d.png", width=700, height=600)

Test1$NTables

for (i in c(4:dim(Test1$NTables)[1])) {
  # i=4
  Pred_center2 <- as.data.frame(Test1$All_theta[,,i], row.names = T)
  colnames(Pred_center2) <- c("dim1","dim2")
  Table_ID_Est <- as.factor(Pred_center2[,1])
  # Rename all levels
  levels(Table_ID_Est) <- c(1:length(unique(Table_ID_Est)))
  p3 <- ggplot(data=clusters, aes(x=X1, y=X2)) +
    geom_point(aes(color=Table_ID_Est), alpha=0.5) + 
    stat_ellipse(aes(color=Table_ID_Est), type = "norm", linetype = 1) +
    geom_point(data=Pred_center2, 
               aes(x=dim1,y=dim2, shape=Table_ID_Est), color="black", size = 4) +
    geom_point(data=Pred_center2, 
               aes(x=dim1,y=dim2, color=Table_ID_Est, shape=Table_ID_Est), size = 2) + 
    ggtitle(paste("Clusters of Gibbs sampler when M=", Sim_M[i])) +
    labs(x = "X1") + labs(y = "X2") +
    geom_point(data=centers, aes(x=dim1,y=dim2), shape=21, size=2, 
               inherit.aes = FALSE) +
    geom_label_repel(data=unique(centers),
                     aes(x=dim1,y=dim2, label="TrueCenter"), fill="green3",
                     fontface = 'bold', color = 'white',
                     box.padding = 5, point.padding = 0.2,
                     segment.color = 'grey50') +
    xlim(-1.2, 4)+ ylim(-4.2,3)
  print(p3)
}
dev.off()
# A little bit different in Mac OS, only use "convert... " is enough, 
# don't need to use "magick"
system("convert -delay 200 2D_gibbs*.png 2D_Clusters.gif")

