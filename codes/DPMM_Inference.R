

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
DPMM_CRP <- function(num_customers, theta) {
  table_ID <- c(1); table_Customers <- c(1)
  Center_TRUE <- rnorm(1, 0, 2); Center <- c(Center_TRUE)
  next.table <- 2; X <- rnorm(1, Center_TRUE[1], 1)
  for (i in 2:num_customers) {
    if (runif(1,0,1) < (theta) / (i - 1 + theta)) {
      # Add a new ball color/ table
      table_ID <- c(table_ID, next.table)
      table_Customers <- c(table_Customers, 1)
      Center_TRUE[next.table] <- rnorm(1, 0, 2)
      X[i] <- rnorm(1, Center_TRUE[next.table], 1)
      Center[i] <- Center_TRUE[next.table]
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
      X[i] <- rnorm(1, Center_TRUE[select_table], 1)
      Center[i] <- Center_TRUE[select_table]
      table_ID <- c(table_ID, select_table)
    }
  }
  DPMM <- data.frame(Table_ID=table_ID, X=X, Center = Center)
  return(DPMM)
}

# Test
Test1 <- DPMM_CRP(200, 2)
(table(Test1$Table_ID))
table(Test1$Center)

2*log(200)

plot(
  table( DPMM_CRP(num_customers = 1000, 100)$Table_ID )
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
  table( PitmanYorPMM_CRP(num_customers = 1000, 2, 0.5)$Table_ID )
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
#'
#' @return
#' @export
#'
#' @examples
Covg_M <- function(Ser_M, y) {
  n_y <- length(y)
  Num_Table <- rep(NA, length(Ser_M))
  # y <- Tables$X
  alpha <- 2
  ini_theta <- runif(n_y, 0, 1) #rep(0.1, n_y)
  q <- F_yiTi <- matrix(NA, nrow = n_y, ncol = n_y)
  cumuProb <- matrix(NA, nrow = n_y, ncol = n_y+1)
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
DPMM_Data <- DPMM_CRP(100, 2)
y <- DPMM_Data$X

# Draw graph to show convergency of algorithm
Sim_M <-  floor(exp(c(1:7)))
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

# Show parameters comparison
Test1 <- Covg_M(600, y)
Test1$NTables; table(DPMM_Data$Table_ID)

names(table(DPMM_Data$Center)); table(DPMM_Data$Center)
Data_Plot <- data.frame(T)
cbind(table( Test1$End_Tables) , names(table(DPMM_Data$Center)), table(DPMM_Data$Center))
plot(table( round(Test1$End_Tables, digits = 3)) )
plot(table(round(DPMM_Data$Center, digits = 3)))

cbind(DPMM_Data$Center, Test1$End_Tables)
# table(Test1$End_Tables); table(DPMM_Data$Center)

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
save(All, file = "data/Covg_M.Rdata")


# Gibbs Sampling: Algorithm 3 ---------------------------------------------
Ser_M <- c(1:50,100, 200)
Ser_M <- 10
DPMM_Data <- DPMM_CRP(1000, 2)
y <- DPMM_Data$X

n_y <- length(y)
Num_Table <- rep(NA, length(Ser_M))

alpha <- 2
ini_Table <- rep(1, n_y) ##c(1:n_y) 

Int <- rep(NA, n_y)
select_Table <- b <- r <- c()
Table <- matrix(NA, nrow = Ser_M, ncol = n_y)
for (m in 1:Ser_M) {
  for (i in 1:n_y) {
    # i <- +1
    n_ic <- table(ini_Table[-i])
    q <- Int_c <- matrix(NA, nrow = 1, ncol = dim(n_ic))
    cumuProb <- matrix(NA, nrow = 1, ncol = dim(n_ic)+2)
    # <- matrix(NA, nrow = n_y, ncol = dim(n_ic)+2)
    for (j in 1:dim(n_ic)) {
      sum_yj <- sum(y[ini_Table==n_ic[j]][-i])
      sum_yj2 <- sum((y[ini_Table==n_ic[j]][-i])^2)
      Int_c[1, j] <- n_ic[j]/(n_y-1+alpha) * (sqrt(3))^(n_ic[j]-1) / (sqrt(2*pi* (2*n_ic[j]+1)^3)) * 
        (2 * sum_yj - (2 * n_ic[j] -1) * y[i]) *
        exp((2 * n_ic[j] + 1)/4 * (2*(y[i] + sum_yj)/(2*n_ic[j]+1))^2 - 0.5 * (y[i])^2 - 1/3 * sum_yj2)
    }
    Int[i] <- alpha/(n_y-1+alpha) * 1/sqrt(6 * pi) * exp(-1/6 * (y[i])^2)
    b[i] <- 1/(sum(Int_c[1,]) + Int[i])
    q[1, ] <- b[i] * Int_c[1, ]
    r[i] <- b[i] *  Int[i]
    cumuProb[1, ] <- cumulateP(c(q[1,], r[i]), N = (dim(n_ic)+1))
    select_Table[i] <- max(which(cumuProb[1, ] <= runif(1,0,1)))
    if (select_Table[i] <= dim(n_ic)) {
      ini_Table[i] <- as.numeric(names(n_ic[select_Table[i]]))
    }  
    # if (select_Table[i] > (n_y-1)) {
    #   Table[i] <- rnorm(1, 2/3 * y[i], sd = sqrt(2/3))
    # } else {
    #   Table[i] <- (ini_Table[-c(i)])[select_Table[i]]
    # }
  }
  Table[m, ] <- ini_Table
  Ser_ID <- which(Ser_M==m)
  if (length(Ser_ID)!=0) {
    Num_Table[Ser_ID] <- length(unique(Table))
  }
}
unique(Table[Ser_M,])

