### HPV example
library(gdata)
library(WriteXLS)
library(forecast)
library(MisRepARMA)
library(ggplot2)
library(ggfortify)

setwd("/home/dmorinya/Insync/dmorina@ub.edu/OneDrive Biz/Projectes/2018/0102018. BGSM/004UNDERREPORTEDHPV/Paper/Scientific Reports/Rev2")

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
hpv_inc <- hpv_inc[order(hpv_inc$Any, hpv_inc$Setmana), ]

### Figure 1
rho_obs <- acf(hpv_inc$incid, plot=FALSE)

### Plot
df <- data.frame(log(rho_obs$acf[2:10]), seq(1:9))
m1 <- lm(log(rho_obs$acf[2:10])~1+seq(1:9))

ggplot(df, aes(seq.1.9., log.rho_obs.acf.2.10..)) + geom_point(color="#de2d26") +
  geom_smooth(method="lm", aes(y=df$log.rho_obs.acf.2.10..,x=df$seq.1.9.), se=FALSE, color="black")+
  theme(legend.position = "none")+ylab(expression(hat(log(rho[k]))))+xlab("Lag")+scale_x_continuous(breaks=seq(1:9))+
  annotate(geom = 'text', x = 3, y = -4.4, 
                    label = "hat(log(rho[k])) == -2.05-0.17*k", 
                    parse = TRUE, color = 'black', size=5)+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),
                                                                            axis.text.x=element_text(size=12),
                                                                            axis.text.y=element_text(size=12))

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
comb.ts <- ts.union(hpv.ts, y_t.ts)
p <- autoplot(comb.ts, facets=FALSE)
p+scale_color_manual(labels = c("Actual", "Forecasted"),
                                     values=c("black", "#de2d26"))+aes(linetype = series)+
  ylab("HPV incidence (x 100,000 individuals)")+xlab("")+theme_bw()+theme(legend.position = "none")+
  theme(axis.title=element_text(size=14,face="bold"), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12))