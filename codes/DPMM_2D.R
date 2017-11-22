library(ggplot2);library(MASS)

#' DPMM_2D
#' the function to generate Chinese Restaurant Process with base measure of 2D Gaussian distribution
#'
#' @param num_customers the total number of customers for simulation
#' @param theta the concentration parameter in Dirichlet Process, or 
#'
#' @return
#' @export
#'
#' @examples
DPMM_2D <- function(num_customers, theta, sig2=0.1, sig2_0=2) {
  
  sigma <- matrix(c(sig2,0,0,sig2),nrow=2)
  mu <- matrix(NA, num_customers, 2)
  clusterList <- NULL
  orig_cluster <- round(theta*log(num_customers)+10)
  mudelta <- mvrnorm(n=orig_cluster, c(0,0), matrix(c(sig2_0,0,0,sig2_0),nrow = 2))
    # matrix(5*runif(orig_cluster*2, 0, 1) , 
    #                 nrow = orig_cluster, ncol = 2)
  
  table_ID <- c(1); table_Customers <- c(1)
  next_table <- 2
  X <- matrix(NA, nrow = num_customers, ncol = 2)
  clusterList <- matrix(NA, nrow = num_customers, ncol = 6)
  X[1,] <- mvrnorm(1, mudelta[table_ID,], sigma)
  clusterList[1,] <- c(table_Customers, table_ID, X[1,], mudelta[table_ID, ])
  for (i in 2:num_customers) {
    if (runif(1,0,1) < (theta) / (i - 1 + theta)) {
      # Add a new table
      table_ID <- c(table_ID, next_table)
      X[i,] <- mvrnorm(n=1,mudelta[next_table,],sigma)
      #put in the 2-D normal cluster using the cluster number as mean vector
      table_Customers <- c(table_Customers, 1)
      clusterList[i,] <- c(i, next_table, X[i,], mudelta[next_table,])
      next_table <- next_table+1
    } else {
      # The nth customer chooses the first unoccupied table with probability
      # theta/(n−1+theta), and an occupied table with probability c/(n−1+theta), where c is the
      # number of people sitting at that table.
      Draw_Rand <- runif(1, 0, 1)
      select_table <- min(which(cumsum(table_Customers/(i-1))>=Draw_Rand))
      
      table_Customers[select_table] <- table_Customers[select_table] + 1
      X[i,] <- mvrnorm(n=1, mudelta[select_table, ], sigma)
      clusterList[i,] <- c(i, select_table, X[i,], mudelta[select_table,])
      table_ID <- c(table_ID, select_table)
    }
  }
  DPMM <- clusterList
  colnames(DPMM) <- c("ID", "Table_ID", "X1", "X2", "Center1", "Center2")
    # data.frame(Table_ID=table_ID, X=X, Center = Center)
  return(DPMM)
  # list(clusterList=clusterList, mudelta=mudelta, table_ID=table_ID)
}


# Simple example to check function ----------------------------------------
n <- 1000; theta <- .4
DPMM01 <- DPMM_2D(n, theta, sig2 = .5, sig2_0=2)

# Draw graph to show clusters ---------------------------------------------
clusters <- as.data.frame(DPMM01, row.names = T)
clusters[,2] <- as.factor(clusters[,2])
centers <- as.data.frame(DPMM01[,5:6], row.names = T)
colnames(centers) <- c("dim1","dim2")
p <- ggplot(data=clusters, aes(x=X1, y=X2)) + 
  geom_point()
p + geom_point(data=centers,
             aes(x=dim1,y=dim2),inherit.aes = FALSE) +  ggtitle("Observed data")
p2 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID)) + 
  geom_point()
p2 + geom_point(data=centers,
               aes(x=dim1,y=dim2),inherit.aes = FALSE) + ggtitle("Mixture Model")

# double check power law of asymptotic behavior of the number of clusters 
# Proposition 10.4 in lecture notes.
max(as.numeric(clusters[,2]))
theta * log(n)
# seems good, but maybe we can change n to show its rate
# movie3d(play3d(), duration = 18,
#         dir = getwd(), convert = TRUE)
# 
# movie3d(spin3d(axis = c(0, 0, 0)), duration = 18,
#         dir = getwd(), convert = TRUE)

# Gibbs Sampling: Algorithm 1 ---------------------------------------------

#' Covg_M_2D
#' is a function to show the convergence of Algorithm 1 or 3 in Neal, R. M. (2000). 
#' Markov chain sampling methods for Dirichlet process mixture models. 
#' Journal of computational and graphical statistics, 9(2), 249-265.
#'
#' @param Ser_M 
#' @param y 
#' @param alpha concentration parameter of Dirichlet
#'
#' @return
#' @export
#'
#' @examples
Covg_M_2D <- function(Ser_M, y, alpha=2, sig2=0.01, sig2_0=2) {
  S1 <- sig2_0/(sig2_0 + sig2)
  n_y <- dim(y)[1]
  Num_Table <- rep(NA, length(Ser_M))
  All_theta <- array(NA, dim = c(n_y, 2, length(Ser_M)) )
  ini_theta <-  y #matrix(rep(0, 2*n_y), nrow = n_y, ncol = 2)  # matrix(runif(2*n_y, 0, 1), nrow = n_y, ncol = 2) # 
  theta <- matrix(NA, nrow = n_y, ncol = 2)
  q <- F_yiTi <- matrix(NA, nrow = n_y, ncol = n_y)
  Int <- rep(NA, n_y)
  select_cluster <- b <- r <- c()
  for (m in 1:max(Ser_M)) {
    for (i in 1:n_y) {
      F_yiTi[i, ] <- 1/(2*pi*sig2) * exp(-1/(2*sig2) * rowSums((y[i,] - ini_theta)^2))
      Int[i] <- 1/(2 * (sig2+sig2_0) * pi) * exp(-1/(2* (sig2+sig2_0)) * sum((y[i,])^2))
      b[i] <- 1/(sum(F_yiTi[i,-c(i)]) + alpha * Int[i])
      q[i, ] <- b[i] * F_yiTi[i, ]
      r[i] <- b[i] * alpha * Int[i]
      Draw_Rand <- runif(1, 0, 1)
      select_cluster[i] <- min(which(cumsum( c(q[i,-c(i)], r[i]) ) >= Draw_Rand))
      
      if (select_cluster[i] > (n_y-1)) {
        ini_theta[i, ] <- mvrnorm(n = 1, mu = sig2_0/(sig2+sig2_0) * y[i,], 
                              Sigma = matrix(c(sig2_0*sig2/(sig2+sig2_0),0,0,sig2_0*sig2/(sig2+sig2_0)), 2, 2, byrow = T))
      } else {
        ini_theta[i, ] <- (ini_theta[-c(i), ])[select_cluster[i], ]
      }
    }
    theta <- ini_theta
    Ser_ID <- which(Ser_M==m)
    if (length(Ser_ID)!=0) {
      All_theta[,,Ser_ID] <- theta
      Num_Table[Ser_ID] <- dim(unique(theta))[1]
    }
  }
  list(NTables=data.frame(M=Ser_M, Tables=Num_Table), All_theta=All_theta, 
       End_Tables=theta)
}

# 
sig2=.01; sig2_0=1
DPMM_Data <- DPMM_2D(num_customers = 1000, theta = 0.4, 
                     sig2 = sig2, sig2_0 = sig2_0)
plot(
  table( DPMM_Data[, "Table_ID"] )
  ,xlab="Table Index", ylab="Frequency", 
  main = "Chinese Restaurant Process"
)
# Draw graph to show clusters ---------------------------------------------
clusters <- as.data.frame(DPMM_Data, row.names = T)
clusters[,2] <- as.factor(clusters[,2])
centers <- as.data.frame(DPMM_Data[,5:6], row.names = T)
colnames(centers) <- c("dim1","dim2")
p <- ggplot(data=clusters, aes(x=clusters[,3], y=clusters[,4])) + 
  geom_point()
p + geom_point(data=centers,
               aes(x=dim1,y=dim2),inherit.aes = FALSE)
p2 <- ggplot(data=clusters, aes(x=clusters[,3], y=clusters[,4], color=clusters[,2])) + 
  geom_point()
p2 + geom_point(data=centers,
                aes(x=dim1,y=dim2),inherit.aes = FALSE)

# load(file = "data/Covg_M.Rdata")

y <- DPMM_Data[, c("X1", "X2")]
Sim_M <-  floor(exp(c(1:8, 9, 9.5, 10)))

Test1 <- Covg_M_2D(Ser_M = Sim_M[1:11],y = y, alpha = 0.4, 
                   sig2 = sig2, sig2_0 = sig2_0)

ALLD <- list(Data=DPMM_Data, Results=Test1)
save(ALLD, file = "data/ALL_D.Rdata")

load(file = "data/ALL_D.Rdata")
Test1 <- ALLD$Results
DPMM_Data <- ALLD$Data
Test1$NTables
dim(unique(Test1$End_Tables))[1]
unique(Test1$All_theta[,,6]); unique(Test1$All_theta[,,7]); 
unique(Test1$All_theta[,,8]); unique(Test1$All_theta[,,9]); 
unique(Test1$All_theta[,,10]);unique(Test1$All_theta[,,11])
unique(Test1$End_Tables)
unique(DPMM_Data[, c("Center1", "Center2")])

library(ggplot2)
clusters <- as.data.frame(DPMM_Data, row.names = T)
clusters[,2] <- as.factor(clusters[,2])
centers <- as.data.frame(DPMM_Data[,5:6], row.names = T)
colnames(centers) <- c("dim1","dim2")
p3 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID)) + 
  geom_point()
p3 + geom_point(data=centers, aes(x=dim1,y=dim2),
                inherit.aes = FALSE) + 
  ggtitle("True data with clusters") +
  labs(x = "X1") + labs(y = "X2")
  

Pred_center2 <- as.data.frame(Test1$All_theta[,,11], row.names = T)
colnames(Pred_center2) <- c("dim1","dim2")
Table_ID_Est <- as.factor(Pred_center2[,1])
# Rename all levels
levels(Table_ID_Est) <- c(1:length(unique(Table_ID_Est)))
p3 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID_Est)) + 
  geom_point() 
p3 + geom_point(data=Pred_center2, aes(x=dim1,y=dim2),
                inherit.aes = FALSE) + 
  ggtitle(paste("Clusters of Gibbs sampler when M=", Sim_M[11])) +
  labs(x = "X1") + labs(y = "X2")


Pred_center3 <- as.data.frame(Test1$All_theta[,,7], row.names = T)
colnames(Pred_center3) <- c("dim1","dim2")
p3 <- ggplot(data=clusters, aes(x=clusters[,3], y=clusters[,4], 
                                color=as.factor(Pred_center3[,1]))) + 
  geom_point()
p3 + geom_point(data=Pred_center3, aes(x=dim1,y=dim2),
                inherit.aes = FALSE)

Pred_center4 <- as.data.frame(Test1$All_theta[,,8], row.names = T)
colnames(Pred_center4) <- c("dim1","dim2")
p3 <- ggplot(data=clusters, aes(x=clusters[,3], y=clusters[,4], 
                                color=as.factor(Pred_center4[,1]))) + 
  geom_point()
p3 + geom_point(data=Pred_center4, aes(x=dim1,y=dim2),
                inherit.aes = FALSE)


Pred_center5 <- as.data.frame(Test1$All_theta[,,9], row.names = T)
colnames(Pred_center5) <- c("dim1","dim2")
p3 <- ggplot(data=clusters, aes(x=clusters[,3], y=clusters[,4], 
                                color=as.factor(Pred_center5[,1]))) + 
  geom_point()
p3 + geom_point(data=Pred_center5, aes(x=dim1,y=dim2),
                inherit.aes = FALSE)

Pred_center6 <- as.data.frame(Test1$All_theta[,,10], row.names = T)
colnames(Pred_center6) <- c("dim1","dim2")
p3 <- ggplot(data=clusters, aes(x=clusters[,3], y=clusters[,4], 
                                color=as.factor(Pred_center6[,1]))) + 
  geom_point()
p3 + geom_point(data=Pred_center6, aes(x=dim1,y=dim2),
                inherit.aes = FALSE)

Pred_center7 <- as.data.frame(Test1$All_theta[,,11], row.names = T)
colnames(Pred_center7) <- c("dim1","dim2")
p3 <- ggplot(data=clusters, aes(x=clusters[,3], y=clusters[,4], 
                                color=as.factor(Pred_center7[,1]))) + 
  geom_point()
p3 + geom_point(data=Pred_center7, aes(x=dim1,y=dim2),
                inherit.aes = FALSE)


# Draw graph to show convergency of algorithm
Sim_M <-  floor(exp(c(1:7)))
Repp_1 <- replicate(100, 
                    {Test1 <- Covg_M(Sim_M[1], y, alpha = 0.5)
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

Sim_M <-  floor(exp(c(1:7))) # floor(exp(c(1:2))) 
Sim_N <- 100
NT_i <- c(); Avg_NT <- c() # average number of tables

for (j in 1:length(Sim_M)) {
  for (i in 1:Sim_N) {
    NT_i[i] <- Covg_M(Sim_M[j], y)$NTables$Tables
  }
  Avg_NT[j] <- mean(NT_i)
}

Result_Covg <- data.frame(Sim_M=Sim_M, Avg_NT=Avg_NT)
All <- list(Result_Covg=Result_Covg, Data=DPMM_Data)