library(reshape2)
library(forecast)
library(MisRepARMA)
library(WriteXLS)

### Read the data
dades <- read.table("Data/COVID19_China.csv", header=T, sep=",")

COVID19 <- melt(dades, id.vars=c("Province.State", "Country.Region", "Population"))
colnames(COVID19) <- c("Province", "Country", "Population", "Date", "Cases")
rm(dades)

COVID19 <- COVID19[order(COVID19$Province), ]
COVID19$Date <- rep(seq(as.Date("2020-01-22", format="%Y-%m-%d"), as.Date("2020-02-26", format="%Y-%m-%d"), "days"), 31)
COVID19$Incidence <- COVID19$Cases/COVID19$Population*100000

### Keep only Heilongjiang region
COVIDHeilongjiang <- COVID19[COVID19$Province=="Heilongjiang", ]

### Estimation
set.seed(1234)
pr <- fitMisRepARMA(COVIDHeilongjiang$Incidence, 1e-8, 500, 0, 1, covars=NULL, misReport="U")
summary(pr)

### Reconstruction of the hidden process
y_t <- reconstruct(pr)
mean(COVIDHeilongjiang$Incidence)/mean(y_t)*100

### Check residuals
res1 <- arima(y_t, order=c(0, 0, 1))
res1$method <- "MA(1)"
  checkresiduals(res1, main="")

### Plot the registered and estimated series
inds <- seq(as.Date("2020-01-22"), as.Date("2020-02-26"), by = "day")
par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
  plot(COVIDHeilongjiang$Incidence~inds, xaxt = "n", type = "l", xlab="Time", ylab="COVID-19 incidence (x 100,000 individuals)")
  axis(1, inds, format(inds, "%Y-%m-%d"), cex.axis=0.7)
  lines(y_t~inds, col="red")
  legend("bottom", inset=c(0, -0.34), legend=c("Registered", "Estimated"), col=c("black", "red"), lty=1)
