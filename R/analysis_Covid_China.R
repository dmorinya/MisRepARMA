library(reshape2)
library(forecast)
library(MisRepARMA)
library(WriteXLS)
library(ggplot2)
library(ggfortify)
library(zoo)

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
  checkresiduals(res1, main="")+ theme(plot.title = element_text(hjust = 0.5))

### Plot the registered and estimated series
covid.ts <- with(COVIDHeilongjiang, zoo(Incidence, Date))
y_t.df <- data.frame(y_t, Date=COVIDHeilongjiang$Date)
y_t.ts <- with(y_t.df, zoo(y_t, Date))
comb.ts <- merge.zoo(covid.ts, y_t.ts)
p <- autoplot(comb.ts, facets=NULL)
p+scale_color_manual(labels = c("Actual", "Forecasted"),
                     values=c("black", "#de2d26"))+aes(linetype = Series)+
  ylab("COVID-19 incidence (x 100,000 individuals)")+theme_bw()+theme(legend.position = "none")+
  theme(axis.title=element_text(size=14,face="bold"), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12))+scale_x_date(date_breaks = "3 day", date_labels = "%Y-%m-%d")+
  xlab("")+theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
