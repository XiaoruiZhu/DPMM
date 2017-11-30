

# Data generating  --------------------------------------------------------


#' DPMM_CRP
#' Dirichlet Process Mixture Model generating by Dirichlet process (CRP generating method)
#'
#' @param num_customers number of customers in total
#' @param theta concentration parameter in Dirichlet process
#'
#' @return
#' @export
#'
#' @examples
DPMM_CRP <- function(num_customers, theta, sig2=1, sig2_0=3) {
  table_ID <- c(1); table_Customers <- c(1)
  Center_TRUE <-  rnorm(1, 0, sig2_0) # runif(1, -10, 10) ;
  Center <- c(Center_TRUE)
  next.table <- 2; X <- rnorm(1, Center_TRUE[1], sig2)
  for (i in 2:num_customers) {
    Draw_Rand <- runif(1, 0, 1)
    if (Draw_Rand > (1 - (theta) / (i - 1 + theta)) ) {
      # Add a new ball color/ table
      table_ID <- c(table_ID, next.table)
      table_Customers <- c(table_Customers, 1)
      Center_TRUE[next.table] <- rnorm(1, 0, sig2_0) # runif(1, -10, 10) 
      X[i] <- rnorm(1, Center_TRUE[next.table], sig2)
      Center[i] <- Center_TRUE[next.table]
      next.table <- next.table+1
    } else {
      # Pick out a ball from the urn, and add back a
      # ball of the same color.
      # The nth customer chooses the first unoccupied table with probability
      # α/(n−1+α), and an occupied table with probability c/(n−1+α), where c is the
      # number of people sitting at that table.
      # select.table <- table_ID[sample(1:length(table_ID), 1)]
      select_table <- min(which(cumsum(table_Customers/(i-1+theta))>=Draw_Rand))
      
      table_Customers[select_table] <- table_Customers[select_table] + 1
      X[i] <- rnorm(1, Center_TRUE[select_table], sig2)
      Center[i] <- Center_TRUE[select_table]
      table_ID <- c(table_ID, select_table)
    }
  }
  DPMM <- data.frame(Table_ID=table_ID, X=X, Center = Center)
  return(DPMM)
}

# Test
Test1 <- DPMM_CRP(2000, 2)
(table(Test1$Table_ID))
# table(Test1$Center)

2*log(2000)

plot(
  table( DPMM_CRP(num_customers = 1000, .4)$Table_ID )
  ,xlab="Table Index", ylab="Frequency", 
  main = "Chinese Restaurant Process"
)
# Test

#' DPMM_PitmanYor
#' Dirichlet Process Mixture Model generating by Pitman-Yor process (CRP generating method)
#'
#' @param num_customers number of customers in total
#' @param theta concentration parameter in Dirichlet process
#' @param alpha is discount parameter in Pitman-Yor process
#'
#' @return
#' @export
#'
#' @examples
PitmanYorPMM_CRP <- function(num_customers, theta, alpha) {
  table_ID <- c(1); table_Customers <- c(1)
  Center_TRUE <- rnorm(1, 0, 2)
  next.table <- 2; X <- rnorm(1, Center_TRUE[1], 1)
  for (i in 2:num_customers) {
    if (runif(1,0,1) < (theta+length(table_ID) * alpha) / (i - 1 + theta)) {
      # if (runif(1,0,1) < (theta+length(table) * beta) / (theta + i)) {
      # Add a new ball color/ table
      table_ID <- c(table_ID, next.table)
      table_Customers <- c(table_Customers, 1)
      Center_TRUE[table_ID] <- rnorm(1, 0, 2)
      X[i] <- rnorm(1, Center_TRUE[table_ID], 1)
      next.table <- next.table+1
    } else {
      # Pick out a ball from the urn, and add back a
      # ball of the same color.
      # The nth customer chooses the first unoccupied table with probability
      # α/(n−1+α), and an occupied table with probability c/(n−1+α), where c is the
      # number of people sitting at that table.
      # select.table <- table_ID[sample(1:length(table_ID), 1)]
      
      Draw_Rand <- runif(1, 0, 1)
      select_table <- min(which(cumsum(table_Customers/(i-1))>=Draw_Rand))
      
      table_Customers[select_table] <- table_Customers[select_table] + 1
      X[i] <- rnorm(1, Center_TRUE[select_table], 1) # could be change based on base measure
      table_ID <- c(table_ID, select_table)
    }
  }
  DPMM <- data.frame(Table_ID=table_ID, X=X)
  return(DPMM)
}

# Test
(table(PitmanYorPMM_CRP(200, 2, 0.1)$Table_ID))
2*log(200)

plot(
  table( PitmanYorPMM_CRP(num_customers = 1000, .5, 0.1)$Table_ID )
  ,xlab="Table Index", ylab="Frequency",main = "Chinese Restaurant Process"
)
# Test


# Inference Algorithm -----------------------------------------------------


# Gibbs Sampling: Algorithm 1 ---------------------------------------------

#' Covg_M 
#' is a function to show the convergence of Algorithm 1 or 3 in Neal, R. M. (2000). 
#' Markov chain sampling methods for Dirichlet process mixture models. 
#' Journal of computational and graphical statistics, 9(2), 249-265.
#'
#' @param Ser_M 
#' @param y 
#' @param alpha concentration parameter of Dirichlet
#' @param Initial either 1 when all initial are random, or 2 when all initials are 0.1 
#'
#' @return
#' @export
#'
#' @examples
Covg_M <- function(Ser_M, y, alpha=2, Initial=1) {
  n_y <- length(y)
  Num_Table <- rep(NA, length(Ser_M))
  # y <- Tables$X
  # alpha <- 2
  if (Initial==1) {
    ini_theta <- runif(n_y, -4, 4) # #rep(0.1, n_y) #
  } else if (Initial==2) {
    ini_theta <-  rep(0.1, n_y) #runif(n_y, -4, 4)
  }
  q <- F_yiTi <- matrix(NA, nrow = n_y, ncol = n_y)
  Int <- rep(NA, n_y)
  theta <- select_theta <- b <- r <- c()
  for (m in 1:max(Ser_M)) {
    for (i in 1:n_y) {
      F_yiTi[i, ] <- 1/sqrt(2*pi) * exp(-1/2 * (y[i] - ini_theta)^2)
      Int[i] <- 1/sqrt(6 * pi) * exp(-1/6 * (y[i])^2)
      b[i] <- 1/(sum(F_yiTi[i,-c(i)]) + alpha * Int[i])
      q[i, ] <- b[i] * F_yiTi[i, ]
      r[i] <- b[i] * alpha * Int[i]
      # cumuProb[i, ] <- cumulateP(c(q[i,-c(i)], r[i]), n_y)
      # select_theta[i] <- max(which(cumuProb[i, ] <= runif(1,0,1)))
      Draw_Rand <- runif(1, 0, 1)
      select_theta[i] <- min(which(cumsum( c(q[i,-c(i)], r[i]) ) >= Draw_Rand))
      
      if (select_theta[i] > (n_y-1)) {
        theta[i] <- rnorm(1, 2/3 * y[i], sd = sqrt(2/3))
      } else {
        theta[i] <- (ini_theta[-c(i)])[select_theta[i]]
      }
    }
    ini_theta <- theta
    Ser_ID <- which(Ser_M==m)
    if (length(Ser_ID)!=0) {
      Num_Table[Ser_ID] <- length(unique(theta))
    }
  }
  list(NTables=data.frame(M=Ser_M, Tables=Num_Table), End_Tables=theta)
}

# 
DPMM_Data <- DPMM_CRP(100, .4, sig2 = 1, sig2_0 = 2)
plot(
  table( DPMM_Data$Table_ID )
  ,xlab="Table Index", ylab="Frequency", 
  main = "Chinese Restaurant Process"
)
load(file = "data/Covg_M.Rdata")
DPMM_Data <- All$Data

y <- DPMM_Data$X

# Draw graph to show convergency of algorithm
Sim_M <-  floor(exp(c(0:5)))
Repp_1 <- replicate(100, 
                    {Test1 <- Covg_M(Sim_M[1], y)
                    Test1$NTables[,2]})
Repp_2 <- replicate(100, 
                    {Test1 <- Covg_M(Sim_M[2], y)
                    Test1$NTables[,2]})
Repp_3 <- replicate(100, 
                    {Test1 <- Covg_M(Sim_M[3], y)
                    Test1$NTables[,2]})
Repp_4 <- replicate(100, 
                    {Test1 <- Covg_M(Sim_M[4], y)
                    Test1$NTables[,2]})
Repp_5 <- replicate(100, 
                    {Test1 <- Covg_M(Sim_M[5], y)
                    Test1$NTables[,2]})
Repp_6 <- replicate(100, 
                    {Test1 <- Covg_M(Sim_M[6], y)
                    Test1$NTables[,2]})
Repp_try <- replicate(10, 
                    {Test1 <- Covg_M(Sim_M[7], y)
                    Test1$NTables[,2]})
Repp_try
par(mfrow=c(2,3))
hist(Repp_1, breaks = 20, main = "M=2", 
     xlim = c(min(dim(table(DPMM_Data$Table_ID))-1, min(Repp_1)), max(Repp_1))); 
abline(v= dim(table(DPMM_Data$Table_ID)), col="red")
hist(Repp_2, breaks = 20, main = "M=7", 
     xlim = c(min(dim(table(DPMM_Data$Table_ID))-1, min(Repp_2)), max(Repp_2))); 
abline(v= dim(table(DPMM_Data$Table_ID)), col="red")
hist(Repp_3, breaks = 20, main = "M=20", 
     xlim = c(min(dim(table(DPMM_Data$Table_ID))-1, min(Repp_3)), max(Repp_3))); 
abline(v= dim(table(DPMM_Data$Table_ID)), col="red")
hist(Repp_4, breaks = 20, main = "M=54", 
     xlim = c(min(dim(table(DPMM_Data$Table_ID))-1, min(Repp_4)), max(Repp_4))); 
abline(v= dim(table(DPMM_Data$Table_ID)), col="red")
hist(Repp_5, breaks = 20, main = "M=148", 
     xlim = c(min(dim(table(DPMM_Data$Table_ID))-1, min(Repp_5)), max(Repp_5))); 
abline(v= dim(table(DPMM_Data$Table_ID)), col="red")
hist(Repp_6, breaks = 20, main = "M=403", 
     xlim = c(min(dim(table(DPMM_Data$Table_ID))-1, min(Repp_6)), max(Repp_6))); 
abline(v= dim(table(DPMM_Data$Table_ID)), col="red")


# Convergency of Algorithm 1 ----------------------------------------------

DPMM_Data <- DPMM_CRP(1000, .3, sig2 = 1, sig2_0 = 2)
plot(
  table( DPMM_Data$Table_ID )
  ,xlab="Table Index", ylab="Frequency", 
  main = "Chinese Restaurant Process"
)
hist(DPMM_Data$X)
y <- DPMM_Data$X

# load(file = "data/OneSimu_Diff.Rdata")
# DPMM_Data <- OneSimu_Diff$Data
DPMM_Data$Table_ID <- as.factor(DPMM_Data$Table_ID)
library(gridExtra)

plot1 <- ggplot(DPMM_Data, aes(x = X)) +   # basic graphical object
  geom_histogram(stat = "bin", binwidth = 0.09) +
  ggtitle("Observed data")
print(plot1)

plot2 <- ggplot(DPMM_Data, aes(x = X, color=Table_ID, fill=Table_ID)) +   # basic graphical object
  geom_histogram(stat = "bin", binwidth = 0.09) +
  ggtitle("Mixture Model")
print(plot2)

# Sim ---------------------------------------------------------------------


Sim_M <-  floor(exp(c(0:6))) # floor(exp(c(1:2))) 
Sim_N <- 20
NT_i <- c(); Avg_NT <- c() # average number of tables

for (j in 1:length(Sim_M)) {
  for (i in 1:Sim_N) {
    NT_i[i] <- Covg_M(Sim_M[j], y, alpha=0.3, Initial = 1)$NTables$Tables
  }
  Avg_NT[j] <- mean(NT_i)
}

Result_Covg <- data.frame(Sim_M=Sim_M, Avg_NT=Avg_NT)
All <- list(Result_Covg=Result_Covg, Data=DPMM_Data)
save(All, file = "data/Covg_M1.Rdata")
# 
Sim_M <-  floor(exp(c(0:6))) # floor(exp(c(1:2))) 
Sim_N <- 20
NT_i <- c(); Avg_NT <- c() # average number of tables

for (j in 1:length(Sim_M)) {
  for (i in 1:Sim_N) {
    NT_i[i] <- Covg_M(Sim_M[j], y, alpha=0.3, Initial = 2)$NTables$Tables
  }
  Avg_NT[j] <- mean(NT_i)
}

Result_Covg <- data.frame(Sim_M=Sim_M, Avg_NT=Avg_NT)
All2 <- list(Result_Covg=Result_Covg, Data=DPMM_Data)
save(All2, file = "data/Covg_M2.Rdata")


# Another simulation ------------------------------------------------------


Sim_M <-  floor(exp(c(0:6))) # floor(exp(c(1:2))) 
Sim_N <- 1
NT_i <- c(); Avg_NT <- c() # average number of tables

for (j in 1:length(Sim_M)) {
  for (i in 1:Sim_N) {
    NT_i[i] <- Covg_M(Sim_M[j], y, alpha=0.3, Initial = 1)$NTables$Tables
  }
  Avg_NT[j] <- mean(NT_i)
}

Result_Covg <- data.frame(Sim_M=Sim_M, Avg_NT=Avg_NT)
OneSimu_Diff <- list(Result_Covg=Result_Covg, Data=DPMM_Data)
save(OneSimu_Diff, file = "data/OneSimu_Diff.Rdata")
# 
Sim_M <-  floor(exp(c(0:6))) # floor(exp(c(1:2))) 
Sim_N <- 1
NT_i <- c(); Avg_NT <- c() # average number of tables

for (j in 1:length(Sim_M)) {
  for (i in 1:Sim_N) {
    NT_i[i] <- Covg_M(Sim_M[j], y, alpha=0.3, Initial = 2)$NTables$Tables
  }
  Avg_NT[j] <- mean(NT_i)
}

Result_Covg <- data.frame(Sim_M=Sim_M, Avg_NT=Avg_NT)
OneSimu_Same <- list(Result_Covg=Result_Covg, Data=DPMM_Data)
save(OneSimu_Same, file = "data/OneSimu_Same.Rdata")

Covg_fig <- data.frame(rbind(OneSimu_Diff$Result_Covg, OneSimu_Same$Result_Covg), 
           Label = c(rep("AllDiff", 7), rep("AllSame", 7)))
Fig_D <- list(Covg_fig=Covg_fig, Theo_cluster=length(unique(DPMM_Data$Table_ID)))
save(Fig_D, file = "data/Fig_D.Rdata")

ggplot(data=Fig_D$Covg_fig, aes(x=Sim_M, y=Avg_NT, color=Label)) +
  geom_line() + geom_point() + 
  geom_hline(yintercept = Fig_D$Theo_cluster, color="blue") +
  labs(x = "Simulate Times of Gibbs Sampler") + labs(y = "Average number of clusters") +
  ggtitle("Algorithm convergency") +
  theme(plot.title = element_text(hjust = 0.5))


# Draw convergency --------------------------------------------------------


# DPMM --------------------------------------------------------------------
Draw_n <- c(1:10, 30, 50, 150, 300, 500, 1000, 2000); alpha = 0.3
y <- DPMM_Data$X
library(extrafont)
loadfonts(device = "pdf")

# dir.create("examples")
setwd("Example")
require(ggplot2)
png(file="M1_gibbs%02d.png", width=600, height=600)

for (i in c(1:length(Draw_n))) {
  # plot.new()
  Test1 <- Covg_M(Draw_n[i], y, alpha = 0.3, Initial = 1)
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
system("magick convert -delay 100 M1_gibbs*.png Try1.gif")

# to not leave the directory with the single jpeg files
# I remove them.
file.remove(list.files(pattern=".png"))
file.remove(list.files(pattern=".gif"))


# DPMM_ 2 -----------------------------------------------------------------

Draw_n <- c(1:10, 30, 50, 150, 300, 500, 1000, 2000); alpha = 0.3
y <- DPMM_Data$X
library(extrafont)
loadfonts(device = "pdf")

# dir.create("examples")
setwd("Example")
require(ggplot2)
png(file="M2_gibbs_Diff%02d.png", width=600, height=600)

for (i in c(1:length(Draw_n))) {
  # plot.new()
  Test1 <- Covg_M(Draw_n[i], y, alpha = 0.3, Initial = 2)
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
system("magick convert -delay 100 M2_gibbs_Diff*.png Try2.gif")


