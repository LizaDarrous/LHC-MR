### Generates simple plots for simulation
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(ggpubr)
library(extrafont)
library(SimDesign)
#font_import(pattern = "cmr10*")
#loadfonts()

#colnames = mLL - pX - pU - pY - h2X - h2Y - tX - tY - a - b - EggerEst - EggerErr - WMedianEst - WMedianErr - IVWEst - IVWErr - ModeEst - ModeErr - WModeEst - WModeErr
df=as.data.frame(read.csv("Standard_SIM.csv",header = FALSE))  ## result of 100 concatinated  result files from LHC_SIMarray.sh
df_fil=df[,-c(1,12,14,16,18,20)] # remove mLL and MR standard error columns, leave only estimated parameters.
colnames(df_fil)=c("pX","pU","pY","h2X","h2Y","tX", "tY", "a", "b","Egger", "WMedian","IVW","Mode","WMode")
## Separate parameters from causal effect estimates
dat <- tibble::tribble(
  ~Parameter, ~Type,
  "pX", " ",
  "pU", " ", 
  "pY",  " ",
  "h2X", " ",
  "h2Y", " ",
  "tX", " ",
  "tY", " ",
  "a", "causal estimate",
  "b",  " ",
  "Egger", "causal estimate",
  "WMedian", "causal estimate",
  "IVW", "causal estimate",
  "Mode" , "causal estimate",
  "WMode", "causal estimate"
)

# Make them ordered factors to appear on boxplots as such
dat$Parameter <- as.factor(dat$Parameter)
dat$Parameter <- forcats::fct_relevel(dat$Parameter, "pX", "pU", "pY", "h2X", "h2Y", "tX", "tY", "a", "b", "Egger", "WMedian", "IVW", "Mode" ,"WMode")
df2=cbind(dat, t(df_fil))
df3 = reshape2::melt(df2, id.vars = c("Parameter","Type"))

theme_set(theme_bw())

## Set up real/actual paramater values from data generation
tX = sqrt(0.3)*0.3  # sqrt(h2U)*qX
tY = sqrt(0.3)*0.2  # sqrt(h2U)*qY

data.segm = data.frame(x='a',xend='WMode',y=0.3,yend=0.3, Type="causal estimate")
data.points = data.frame(Parameter = factor(c("pX","pU","pY","h2X","h2Y","tX", "tY","b")),
           Value = c(0.1,0.05,0.15,0.25,0.2,tX,tY,0),
           Type=" ")

## Here is the plotting
ggplot(data = df3) + 
  geom_boxplot(aes(Parameter, value)) + 
  geom_point(data = data.points, 
             aes(x=Parameter,y=Value), inherit.aes = FALSE, color = "Blue", fill = 'Blue', shape=23, size=2 )+ 
  geom_segment(data=data.segm,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, color = 'Blue')+
  ylab("Estimated value") + xlab("Parameter") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(plot.title = element_text(size=20, face = "bold"), plot.subtitle = element_text(size=15), 
        axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size = 14)) +
  facet_wrap(~Type, scales='free')

### Bias, RMSE and Variance calculations
y = c(0.1,0.05,0.15,0.25,0.2,tX,tY,0.3,0,0.3,0.3,0.3,0.3,0.3)  # real/actual parameter values according to columns in df_fil
bias_df=t(as.data.frame(bias(df_fil,y)))
rmse_df=t(as.data.frame(RMSE(df_fil,y)))
var_df=t(as.data.frame(diag(var(df_fil))))

summ=rbind(" ",bias_df,rmse_df,var_df)
rownames(summ)=c("Standard","Bias","RMSE","Variance")
write.table("Standard Settings", "Standard_summary.csv", sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(summ, file = "Standard_summary.csv", sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE, append = TRUE)
