### HPV example
library(gdata)
library(WriteXLS)
library(forecast)
library(MisRepARMA)

hpv <- read.table("Data/HPV.csv", sep=";", header=TRUE)
hpv$X <- NULL
hpv <- hpv[, c(1, 2, 6)]

pob <- read.xls("Data/aec-245.xls")
pob$X <- NULL
pob <- pob[, c(1:2)]

hpv_inc <- merge(hpv, pob, by=c("Any"))

### The period 2010-2014 was analysed in SiM (Fernández-Fontelo et al., 2016)
hpv_inc <- hpv_inc[hpv_inc$Any>=2010 & hpv_inc$Any<=2014, ]

### Incidence x 100,000 persons
hpv_inc$incid <- hpv_inc$Girona.x/hpv_inc$Girona.y*100000
hpv_inc <- hpv_inc[order(hpv_inc$Any, hpv_inc$Setmana),]

### Figure 1
rho_obs <- acf(hpv_inc$incid, plot=FALSE)

### Plot
m1 <- lm(log(rho_obs$acf[2:10])~1+seq(1:9))
summary(m1)
  plot(log(as.numeric(rho_obs$acf[2:10])), ylab=expression(log(rho[k])), xlab="Lag", xlim=c(1,9),  xaxt='n')
  axis(side = 1, at=1:9)
  abline(m1)
  text(7, -4, expression(log(rho[k]) == -2.05-0.17*k))

### Analysis
set.seed(1234)
pr <- fitMisRepARMA(hpv_inc$incid, 1e-8, 500, 1, 0, covars=NULL, misReport="U")
summary(pr)

### Reconstruction of the hidden process
y_t <- reconstruct(pr)
mean(hpv_inc$incid)/mean(y_t)*100

### Check residuals
res1 <- arima(y_t, order=c(1, 0, 0))
res1$method <- "AR(1)"
  checkresiduals(res1, main="")

### Plot the registered and estimated series
hpv.ts <- ts(hpv_inc$incid, start=c(2010, 1), end=c(2014, 12), freq=52)
y_t.ts <- ts(y_t, start=c(2010, 1), end=c(2014, 12), freq=52)
comb.ts <- cbind(hpv.ts, y_t.ts)
  par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
    ts.plot(comb.ts, ylab="HPV incidence (x 100,000 individuals)", plot.type="single", gpars=list(col=c("black", "red")))
    legend("bottom", inset=c(0, -0.35), legend=c("Registered", "Estimated"), col=c("black", "red"), lty=1)