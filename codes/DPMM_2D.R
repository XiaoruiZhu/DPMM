library(ggplot2);library(MASS)

#' DPMM_2D
#' the function to generate DPMM with base measure of 2D Gaussian distribution
#'
#' @param num_customers the total number of customers for simulation
#' @param theta the concentration parameter in Dirichlet Process, or 
#' @param sig2 variances of y
#' @param sig2_0 variances of prior distribution of centers
#'
#' @return
#' @export
#'
#' @examples
DPMM_2D <- function(num_customers, theta, sig2=0.1, sig2_0=2) {
  
  sigma <- matrix(c(sig2,0,0,sig2),nrow=2, byrow = T)
  mu <- matrix(NA, num_customers, 2)
  clusterList <- NULL
  orig_cluster <- round(theta*log(num_customers)+10)
  mudelta <- mvrnorm(n=orig_cluster, c(0,0), matrix(c(sig2_0,0,0,sig2_0),nrow = 2, byrow = T))
    # matrix(5*runif(orig_cluster*2, 0, 1) , 
    #                 nrow = orig_cluster, ncol = 2)
  
  table_ID <- c(1); table_Customers <- c(1)
  next_table <- 2
  X <- matrix(NA, nrow = num_customers, ncol = 2)
  clusterList <- matrix(NA, nrow = num_customers, ncol = 6)
  X[1,] <- mvrnorm(1, mudelta[table_ID,], sigma)
  clusterList[1,] <- c(table_Customers, table_ID, X[1,], mudelta[table_ID, ])
  for (i in 2:num_customers) {
    Draw_Rand <- runif(1, 0, 1)
    if ( Draw_Rand > (1- (theta) / (i - 1 + theta)) ) {
      # Add a new table
      table_ID <- c(table_ID, next_table)
      X[i,] <- mvrnorm(n=1, mudelta[next_table,], sigma)
      #put in the 2-D normal cluster using the cluster number as mean vector
      table_Customers <- c(table_Customers, 1)
      clusterList[i,] <- c(i, next_table, X[i,], mudelta[next_table,])
      next_table <- next_table+1
    } else {
      # The nth customer chooses the first unoccupied table with probability
      # theta/(n−1+theta), and an occupied table with probability c/(n−1+theta), 
      # where c is the number of people sitting at that table.
      
      select_table <- min(which(cumsum(table_Customers/(i-1+theta))>=Draw_Rand))
      # Here is a typo, wrong probability for sitting at old table
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
DPMM01 <- DPMM_2D(n, theta, sig2 = .01, sig2_0=4)

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

# Gibbs Sampling: Algorithm 1 ---------------------------------------------

#' Covg_M_2D
#' function for Gibbs Sampler of 2D-DPMM, can be used to show the convergence of 
#' Algorithm 1 or 3 in Neal, R. M. (2000). 
#' Markov chain sampling methods for Dirichlet process mixture models. 
#' Journal of computational and graphical statistics, 9(2), 249-265.
#'
#' @param Ser_M iteration times series that want to be saved
#' @param y 2D observed data
#' @param alpha concentration parameter of Dirichlet
#' @param sig2 variances of y
#' @param sig2_0 variances of prior distribution of centers
#' @param Initial 1 for all different initials, 2 for all same initials, 3 for initials are same as y
#'
#' @return
#' @export
#'
#' @examples
Covg_M_2D <- function(Ser_M, y, alpha=2, sig2=0.01, sig2_0=2, Initial=1) {
  # S1 <- sig2_0/(sig2_0 + sig2)
  n_y <- dim(y)[1]
  Num_Table <- rep(NA, length(Ser_M))
  All_theta <- array(NA, dim = c(n_y, 2, length(Ser_M)) )
  if (Initial==1) {
    ini_theta <- matrix(runif(2*n_y, 0, 1), nrow = n_y, ncol = 2) # all different
  } else if (Initial==2) {
    ini_theta <- matrix(rep(0.1, 2*n_y), nrow = n_y, ncol = 2) # all same
  } else if (Initial==3) {
    ini_theta <- y
  }
  theta <- matrix(NA, nrow = n_y, ncol = 2)
  q <- F_yiTi <- matrix(NA, nrow = n_y, ncol = n_y)
  Int <- rep(NA, n_y)
  select_cluster <- b <- r <- c()
  for (m in 1:max(Ser_M)) {
    for (i in 1:n_y) {
      Temp_y <- matrix(rep(y[i,], n_y), nrow = n_y, ncol = 2, byrow = T)
      F_yiTi[i, ] <- 1/(2*pi*sig2) * exp(-1/(2*sig2) * rowSums((Temp_y - ini_theta)^2))
      Int[i] <- 1/(2 * (sig2+sig2_0) * pi) * exp(-1/(2* (sig2+sig2_0)) * sum((y[i,])^2))
      b[i] <- 1/(sum(F_yiTi[i,-c(i)]) + alpha * Int[i])
      q[i, ] <- b[i] * F_yiTi[i, ]
      r[i] <- b[i] * alpha * Int[i]
      Draw_Rand <- runif(1, 0, 1)
      select_cluster[i] <- min(which(cumsum( c(q[i,-c(i)], r[i]) ) >= Draw_Rand))
      
      if (select_cluster[i] > (n_y-1)) {
        theta[i, ] <- mvrnorm(n = 1, mu = sig2_0/(sig2+sig2_0) * y[i,], 
                              Sigma = matrix(c(sig2_0*sig2/(sig2+sig2_0),0,0,sig2_0*sig2/(sig2+sig2_0)), 2, 2, byrow = T))
      } else {
        theta[i, ] <- (ini_theta[-c(i), ])[select_cluster[i], ]
      }
    }
    ini_theta <- theta
    Ser_ID <- which(Ser_M==m)
    if (length(Ser_ID)!=0) {
      All_theta[,,Ser_ID] <- theta
      Num_Table[Ser_ID] <- dim(unique(theta))[1]
    }
  }
  list(NTables=data.frame(M=Ser_M, Tables=Num_Table), All_theta=All_theta, 
       End_Tables=theta)
}


# Simulation and show algorithm -------------------------------------------

sig2=0.05; sig2_0=10; true_theta=0.4
DPMM_Data <- DPMM_2D(num_customers = 500, theta = true_theta, 
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
p <- ggplot(data=clusters, aes(x=X1, y=X2)) + 
  geom_point()
p + geom_point(data=centers,
               aes(x=dim1,y=dim2),inherit.aes = FALSE)
p2 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID)) + 
  geom_point()
p2 + geom_point(data=centers,
                aes(x=dim1,y=dim2),inherit.aes = FALSE)

y <- DPMM_Data[, c("X1", "X2")]
Sim_M <-  floor(exp(c(0:8, 8.5, 9.5, 10)))

Test1 <- Covg_M_2D(Ser_M = Sim_M[1:9], y = y, alpha = true_theta, 
                   sig2 = sig2, sig2_0 = sig2_0, Initial = 2)

ALLD <- list(Data=DPMM_Data, Results=Test1)
save(ALLD, file = "data/ALL_D.Rdata")

load(file = "data/ALL_D.Rdata")
Test1 <- ALLD$Results
DPMM_Data <- ALLD$Data
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
  
Pred_center2 <- as.data.frame(Test1$All_theta[,,6], row.names = T)
colnames(Pred_center2) <- c("dim1","dim2")
Table_ID_Est <- as.factor(Pred_center2[,1])
# Rename all levels
levels(Table_ID_Est) <- c(1:length(unique(Table_ID_Est)))
p3 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID_Est)) + 
  geom_point() 
p3 + geom_point(data=Pred_center2, 
                aes(x=dim1,y=dim2, color=Table_ID_Est, shape=Table_ID_Est), size = 4) + 
  geom_point(data=Pred_center2, 
             aes(x=dim1,y=dim2, shape=Table_ID_Est), color="black", size = 2) + 
  ggtitle(paste("Clusters of Gibbs sampler when M=", Sim_M[6])) +
  labs(x = "X1") + labs(y = "X2") +
  geom_point(data=centers, aes(x=dim1,y=dim2), shape=21, size=2, 
             inherit.aes = FALSE) +
  geom_label_repel(data=unique(centers),
                   aes(x=dim1,y=dim2, label="TrueCenter"), fill="green3",
                   fontface = 'bold', color = 'white',
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = 'grey50') 

Pred_center2 <- as.data.frame(Test1$All_theta[,,7], row.names = T)
colnames(Pred_center2) <- c("dim1","dim2")
Table_ID_Est <- as.factor(Pred_center2[,1])
# Rename all levels
levels(Table_ID_Est) <- c(1:length(unique(Table_ID_Est)))
p3 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID_Est)) + 
  geom_point() 
p3 + geom_point(data=Pred_center2, 
             aes(x=dim1,y=dim2, color=Table_ID_Est, shape=Table_ID_Est), size = 4) + 
  geom_point(data=Pred_center2, 
             aes(x=dim1,y=dim2, shape=Table_ID_Est), color="black", size = 2) + 
  ggtitle(paste("Clusters of Gibbs sampler when M=", Sim_M[7])) +
  labs(x = "X1") + labs(y = "X2") +
  geom_point(data=centers, aes(x=dim1,y=dim2), shape=21, size=2, 
             inherit.aes = FALSE) +
  geom_label_repel(data=unique(centers),
                   aes(x=dim1,y=dim2, label="TrueCenter"), fill="green3",
    fontface = 'bold', color = 'white',
    box.padding = 0.35, point.padding = 0.5,
    segment.color = 'grey50') 

Pred_center2 <- as.data.frame(Test1$All_theta[,,8], row.names = T)
colnames(Pred_center2) <- c("dim1","dim2")
Table_ID_Est <- as.factor(Pred_center2[,1])
# Rename all levels
levels(Table_ID_Est) <- c(1:length(unique(Table_ID_Est)))
p3 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID_Est)) + 
  geom_point() 
p3 + geom_point(data=Pred_center2, 
                aes(x=dim1,y=dim2, color=Table_ID_Est, shape=Table_ID_Est), size = 4) + 
  geom_point(data=Pred_center2, 
             aes(x=dim1,y=dim2, shape=Table_ID_Est), color="black", size = 2) + 
  ggtitle(paste("Clusters of Gibbs sampler when M=", Sim_M[8])) +
  labs(x = "X1") + labs(y = "X2") +
  geom_point(data=centers, aes(x=dim1,y=dim2), shape=21, size=2, 
             inherit.aes = FALSE) +
  geom_label_repel(data=unique(centers),
                   aes(x=dim1,y=dim2, label="TrueCenter"), fill="green3",
                   fontface = 'bold', color = 'white',
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = 'grey50') 

Pred_center2 <- as.data.frame(Test1$All_theta[,,9], row.names = T)
colnames(Pred_center2) <- c("dim1","dim2")
Table_ID_Est <- as.factor(Pred_center2[,1])
# Rename all levels
levels(Table_ID_Est) <- c(1:length(unique(Table_ID_Est)))
p3 <- ggplot(data=clusters, aes(x=X1, y=X2, color=Table_ID_Est)) + 
  geom_point() 
p3 + geom_point(data=Pred_center2, 
                aes(x=dim1,y=dim2, color=Table_ID_Est, shape=Table_ID_Est), size = 4) + 
  geom_point(data=Pred_center2, 
             aes(x=dim1,y=dim2, shape=Table_ID_Est), color="black", size = 2) + 
  ggtitle(paste("Clusters of Gibbs sampler when M=", Sim_M[9])) +
  labs(x = "X1") + labs(y = "X2") +
  geom_point(data=centers, aes(x=dim1,y=dim2), shape=21, size=2, 
             inherit.aes = FALSE) +
  geom_label_repel(data=unique(centers),
                   aes(x=dim1,y=dim2, label="TrueCenter"), fill="green3",
                   fontface = 'bold', color = 'white',
                   box.padding = 0.9, point.padding = 0.9,
                   segment.color = 'grey50') 
