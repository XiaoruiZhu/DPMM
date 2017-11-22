# center inference ----------------------------------------------
# generate data
DPMM_Data <- DPMM_CRP(1000, .3)
plot(
  table( DPMM_Data$Table_ID )
  ,xlab="Table Index", ylab="Frequency", 
  main = "Chinese Restaurant Process"
)
y <- DPMM_Data$X
load(file = "data/Inf_centers.Rdata")
require(ggplot2)
# fit
Test1 <- Covg_M(1, y, alpha = 0.3)
Test1$NTables; unique(DPMM_Data$Center); hist(DPMM_Data$Table_ID)
# table(DPMM_Data$Center); table(Test1$End_Tables)
CombD2 <- data.frame(x=c(DPMM_Data$Center, Test1$End_Tables), 
                     Label=c(rep("TRUE", 1000), rep("Gibbs", 1000)))
# Data generation
library(extrafont)
# loadfonts(device = "win")
ggplot(CombD2, aes(x = x, color=Label, fill=Label)) +   # basic graphical object
  geom_histogram(stat = "bin", binwidth = 0.1) +
  ggtitle("Gibbs Centers vs. True Centers (M=1)") +
  xlim(-6, 6) + ylim(0, 1000) +
  theme(text=element_text(size=16, family="Comic Sans MS"))

# DPMM ---------------------------------------------------------------------
# https://www.r-graph-gallery.com/166-basic-animated-graph-with-imagemagick/
Draw_n <- 30; alpha = 0.3
load(file = "data/Inf_centers.Rdata")
y <- DPMM_Data$X

# dir.create("examples")
setwd("Example")
require(ggplot2)
png(file="gibbs%02d.png", width=600, height=600)

for (i in c(1:Draw_n)) {
  # plot.new()
  Test1 <- Covg_M(1000, y, alpha = alpha)
  CombD2 <- data.frame(y=c(DPMM_Data$Center, Test1$End_Tables), 
                       Label=c(rep("TRUE", 100), rep("Gibbs", 100)))
  # Data generation
  p2 <- ggplot(CombD2, aes(x = y, color=Label, fill=Label)) +  
    geom_histogram(stat = "bin", binwidth = 0.1) +
    xlim(-6, 4) + ylim(0, 100)
  print(p2)
}
dev.off()
system("magick convert -delay 100 *.png example_1.gif")
# to not leave the directory with the single jpeg files
# I remove them.
file.remove(list.files(pattern=".png"))
file.remove(list.files(pattern=".gif"))


# DPMM --------------------------------------------------------------------
Draw_n <- c(1:10, 30, 50, 150, 300, 500, 1000, 3000); alpha = 0.3
load(file = "data/Inf_centers.Rdata")
y <- DPMM_Data$X
library(extrafont)
loadfonts(device = "win")

# dir.create("examples")
setwd("Example")
require(ggplot2)
png(file="gibbs%02d.png", width=600, height=600)

for (i in c(1:length(Draw_n))) {
  # plot.new()
  Test1 <- Covg_M(Draw_n[i], y, alpha = 0.3)
  # table(DPMM_Data$Center); table(Test1$End_Tables)
  CombD2 <- data.frame(y=c(DPMM_Data$Center, Test1$End_Tables), 
                       Label=c(rep("TRUE", 1000), rep("Gibbs", 1000)))
  # Data generation
  tit <- paste("Gibbs Centers vs. True Centers : M=", Draw_n[i])
  p2 <- ggplot(CombD2, aes(x = y, color=Label, fill=Label)) +   # basic graphical object
    geom_histogram(stat = "bin", binwidth = 0.1) +
    ggtitle(tit) +
    xlim(-6, 4) + ylim(0, 1000) + 
    theme(text=element_text(size=16, family="Comic Sans MS"))
  
  print(p2)
}
dev.off()
system("magick convert -delay 100 *.png example_1.gif")
# to not leave the directory with the single jpeg files
# I remove them.
file.remove(list.files(pattern=".png"))
file.remove(list.files(pattern=".gif"))


# DPMM_ 2 -----------------------------------------------------------------

Draw_n <- c(1:10, 30, 50, 150, 300, 500, 1000, 3000); alpha = 0.3
load(file = "data/Inf_centers.Rdata")
y <- DPMM_Data$X
library(extrafont)
loadfonts(device = "win")

# dir.create("examples")
setwd("Example")
require(ggplot2)
png(file="gibbs_Diff%02d.png", width=600, height=600)

for (i in c(1:length(Draw_n))) {
  # plot.new()
  Test1 <- Covg_M(Draw_n[i], y, alpha = 0.3)
  # table(DPMM_Data$Center); table(Test1$End_Tables)
  CombD2 <- data.frame(y=c(DPMM_Data$Center, Test1$End_Tables), 
                       Label=c(rep("TRUE", 1000), rep("Gibbs", 1000)))
  # Data generation
  tit <- paste("Gibbs Centers vs. True Centers : M=", Draw_n[i])
  p2 <- ggplot(CombD2, aes(x = y, color=Label, fill=Label)) +   # basic graphical object
    geom_histogram(stat = "bin", binwidth = 0.1) +
    ggtitle(tit) +
    xlim(-6, 4) + ylim(0, 1000) + 
    theme(text=element_text(size=16, family="Comic Sans MS"))
  
  print(p2)
}
dev.off()
system("magick convert -delay 100 gibbs_Diff*.png example_2.gif")
