

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

