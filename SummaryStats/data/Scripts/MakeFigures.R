library(tidyverse)
library(data.table)
library(GGally)
library(mixtools)
library(extRemes)

df = fread("AllRes_parscale.csv")

df %>%
  filter(mLL<min(mLL+0.1)) -> parscale_filtered #Get all parameters 0.1 likelihood away from smallest mLL
main1= select(parscale_filtered, -c(2,3)) #remove iX iY, used for pair plot
#main1 %>% gather() %>% head()
main2=reshape2::melt(parscale_filtered, id.vars = NULL) #used for hist plot
#main2$variable <- factor(main2$variable, levels=c("pX","pU","pY","h2X","h2Y","tX","tY","b","rho"))

p_hist = ggplot(main2, aes(value)) +
  geom_histogram(bins = 50) +
  facet_wrap(. ~ variable, scales = c('free'), ncol = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 15),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 10)
  )
ggsave("Hist.pdf", plot = p_hist, device = "pdf", width = 297, height = 210, units ="mm")

p_pairs = ggpairs(main1,
                  axisLabels = "show", columnLabels = colnames(main1)) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 15),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12)
  )
ggsave("Pairs.pdf", plot = p_pairs, device = "pdf", width = 297, height = 210, units ="mm")

## Obtain bimodal mean values if possible
bimodal=list()
for (x in 1:dim(main1)[2]){
  #single=normalmixEM(df_fil$a, lambda=1, k=1, arbmean = TRUE)$loglik
  print(x)
  double=normalmixEM(main1[,x], k=2,fast=TRUE,maxit=10000) #assuming bimodal structure always
  bimodal[[x]]=c(double$lambda,double$mu,double$sigma, double$sigma^2)
}
names(bimodal)=colnames(main1)
test=as.data.frame(t(as.data.frame(bimodal)))
colnames(test)=c("lambda1","lambda2","mu1","mu2","sigma1","sigma2","var1","var2")
write.table(test, file="Allres_bimodal.csv", sep=",", row.names = TRUE, col.names=NA)

print(paste0(round(test['a','mu1'],2),"/",round(test['a','mu2'],2)))
print(paste0(round(test['b','mu1'],2),"/",round(test['b','mu2'],2)))
