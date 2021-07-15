library(WriteXLS)
library(doParallel)
library(MisRepARMA)

### AR(1) simulation
alpha <- round(seq(0.1, 0.9, 0.1), 1); q <- round(seq(0.1, 0.9, 0.1), 1); 
w <- round(seq(0.1, 0.9, 0.1), 1); mu <- 500
res <- data.frame(expand.grid(alpha=alpha, q=q, w=w))

genEsts <- function(i)
{
  print(paste0("Simulation step ", i, " out of ", dim(res)[1]))
  x    <- as.numeric(arima.sim(model=list(ar=c(res$alpha[i])), n=1000)+mu)
  aux  <- rbinom(length(x), 1, res$w[i]) ### Randomly choose underreported observations
  y    <- ifelse(aux==0, x, res$q[i]*x)  ### Observed series
  ### Parameter estimation
  est_params <- try(fitMisRepARMA(y, 1e-8, 500, 1, 0, covars=NULL, "U"), silent=TRUE)
  j <- 1
  while (class(est_params)=="try-error" & j < 10)
  {
    print(paste0("try-error: ", j))
	  x    <- as.numeric(arima.sim(model=list(ar=c(res$alpha[i])), n=1000)+mu)
    aux  <- rbinom(length(x), 1, res$w[i]) ### Randomly choose underreported observations
    y    <- ifelse(aux==0, x, res$q[i]*x)  ### Observed series
	### Parameter estimation
	  est_params <- try(fitMisRepARMA(y, 1e-8, 500, 1, 0), silent=TRUE)
	  j <- j + 1
  }
  ### Keep the results
  alpha_hat    <- est_params$t0[1]
  alpha_hatLOW <- est_params$t0[1]-qnorm(0.975)*apply(est_params$t, 2, sd, na.rm=T)[1]
  alpha_hatUPP <- est_params$t0[1]+qnorm(0.975)*apply(est_params$t, 2, sd, na.rm=T)[1]
  q_hat        <- est_params$t0[4]
  q_hatLOW     <- est_params$t0[4]-qnorm(0.975)*apply(est_params$t, 2, sd, na.rm=T)[4]
  q_hatUPP     <- est_params$t0[4]+qnorm(0.975)*apply(est_params$t, 2, sd, na.rm=T)[4]
  w_hat        <- est_params$t0[5]
  w_hatLOW     <- est_params$t0[5]-qnorm(0.975)*apply(est_params$t, 2, sd, na.rm=T)[5]
  w_hatUPP     <- est_params$t0[5]+qnorm(0.975)*apply(est_params$t, 2, sd, na.rm=T)[5]  
  return(c(res$alpha[i], res$q[i], res$w[i], alpha_hat, alpha_hatLOW,
	         alpha_hatUPP, q_hat, q_hatLOW, q_hatUPP, w_hat, w_hatLOW, w_hatUPP))
}

nCores <- detectCores()
registerDoParallel(nCores)
dat.fin <- foreach(k=1:dim(res)[1], .combine=rbind) %dopar% genEsts(k)
colnames(dat.fin) <- c("alpha", "q", "w", "alpha_hat", "alpha_hatLOW",
                       "alpha_hatUPP", "q_hat", "q_hatLOW", "q_hatUPP",
                       "w_hat", "w_hatLOW", "w_hatUPP")

### Excel exportation
WriteXLS(as.data.frame(dat.fin), "../Results/sim_AR1_2021.xls")
