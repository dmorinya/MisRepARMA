### Covid-19 Catalonia example
library(dplyr)
library(MisRepARMA)
library(forecast)
library(lubridate)

dades <- read.csv("Data/casos_sexe_municipi.csv", header=T, sep=";") # https://dadescovid.cat/static/csv/casos_sexe_municipi.zip
dades$TIPUSCASDATA <- as.Date(dades$TIPUSCASDATA, format="%d/%m/%Y")
dades <- dades[order(dades$TIPUSCASDATA), ]
covidCAT <- dades %>% group_by(TIPUSCASDATA) %>%
  summarise(Casos=sum(NUMCASOS))
covidCAT <- covidCAT[covidCAT$TIPUSCASDATA>="2021-05-16" & covidCAT$TIPUSCASDATA<="2021-06-20", ]
covidCAT$Incidence <- covidCAT$Casos/7716760*100000

### Covariates
t <- seq(1, length(covidCAT$Incidence), 1)
cov2 <- sin(2*pi*t/7)
cov3 <- cos(2*pi*t/7)
covars <- cbind(t, cov2, cov3)

### Estimation
set.seed(1234)
pr <- fitMisRepARMA(covidCAT$Incidence, 1e-8, 500, 2, 0, covars=covars, misReport="U")
summary(pr)

### Reconstruction of the hidden process
y_t <- reconstruct(pr)
mean(covidCAT$Incidence)/mean(y_t)*100

### Check residuals
res1 <- arima(y_t, order=c(2, 0, 0), xreg=covars)
res1$method <- "AR(2)"
  checkresiduals(res1, main="")

### Plot the registered and estimated series
inds <- seq(as.Date("2021-05-16"), as.Date("2021-06-20"), by = "day")
par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
  plot(covidCAT$Incidence~inds, xaxt = "n", type = "l", xlab="Time", ylab="COVID-19 incidence (x 100,000 individuals)",
       ylim=c(min(covidCAT$Incidence), max(y_t)+2))
  axis(1, inds, format(inds, "%Y-%m-%d"), cex.axis=0.7)
  lines(y_t~inds, col="red")
  legend("bottom", inset=c(0, -0.34), legend=c("Registered", "Estimated"), col=c("black", "red"), lty=1)
