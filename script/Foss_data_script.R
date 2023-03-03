#R RStudio 2022.07.2 Build 576
#R 4.1.2

library(tidyverse)
library(funModeling)
library(ggplot2)
library(ggpubr)
library(dunn.test)
library(mgcv)
library(car)
library(MASS)
library(glmmTMB)
library(DHARMa)
library(rstatix)
library(gridExtra)
library(cowplot)
library(ggeffects)


###################################
#### 1 - Load data and prepare ####
###################################

###Load main data set foss_final
foss <- read_csv("./data/foss_final.csv")
foss

###Add labels and set factors
foss$month <- factor(foss$month, ordered=TRUE, labels = c("September", "October", "November", "December", "January", "February"))
foss$year <- factor(foss$year, labels = c("2017/18", "2018/19", "2019/20"))
foss$lightperiod <- factor(foss$lightperiod, labels = c("Day", "Night"))
foss$photo <- factor(foss$photo, labels = c("Dawn","Day","Dusk","Night"))
foss$lvl_stage <- factor(foss$lvl_stage, labels = c("Rising","Falling","Stable \n (reference)","Stable \n (elevated)"))
foss$p_event <- factor(foss$p_event, labels = c("Operation \n one", "Operation \n two"))
foss$p_day <- factor(foss$p_day, labels = c("24h Pre-PO", "24h Post-PO"))
foss$barrier <- factor(foss$barrier, labels =c("Normal operation", "Pre-dawn barrier closure", "Post barrier trial", "Operational level"))
foss$barrier2 <- factor(foss$barrier2, labels =c("Normal operation", "Pre-dawn barrier closure", "Post barrier trial"))

####Prepare yearly data sets
fossy1 <- foss%>%filter(year=="2017/18")
fossy2 <- foss%>%filter(year=="2018/19")
fossy3 <- foss%>%filter(year=="2019/20")

#############################
#####Tests for normality#####
#############################
#Data expected to be non-normal due to multi-modial temporal variability and seasonal/yearly elements.
#Unlikely to achieve normality with transformation. Should not be considered here.

shapiro.test(foss$total)
qqnorm(foss$total, pch = 1, frame = FALSE)
qqline(foss$total, col = "steelblue", lwd = 2)
#fish count non-normally distributed

shapiro.test(fossy1$total)
qqnorm(fossy1$total, pch = 1, frame = FALSE)
qqline(fossy1$total, col = "steelblue", lwd = 2)
#fish count non-normally distributed

shapiro.test(fossy2$total)
qqnorm(fossy2$total, pch = 1, frame = FALSE)
qqline(fossy2$total, col = "steelblue", lwd = 2)
#fish count non-normally distributed

shapiro.test(fossy3$total)
qqnorm(fossy3$total, pch = 1, frame = FALSE)
qqline(fossy3$total, col = "steelblue", lwd = 2)
#fish count non-normally distributed

#####################################################
#### 2 - Data exploration and summary statistics ####
#####################################################

glimpse(foss)
status(foss)

#######Temporal fish count data#######

#total fish count by year
foss  %>%  group_by(year) %>%  summarise(n = sum(total))
ggplot(foss, aes(year,total))+
  geom_bar(stat="identity")+
  geom_text(aes(label=after_stat(y)),stat = 'summary', fun = sum)
#total fish count by year, light period
foss  %>%  group_by(year,lightperiod) %>%  summarise(n = sum(total), med = median(total), min = min(total), max = max(total), IQR = IQR(total))
ggplot(foss, aes(lightperiod,total))+
  geom_boxplot()
#total fish count by year, photoperiod
foss  %>%  group_by(year,photo) %>%  summarise(n = sum(total))
ggplot(foss, aes(photo,total))+
  geom_boxplot()
ggplot(foss, aes(photo,total))+
  geom_boxplot()+facet_wrap(~year)
#total fish count by hour
foss  %>%  group_by(time2) %>%  summarise(n = sum(total))
ggplot(foss, aes(as.factor(time),total))+
  geom_boxplot()
#total fish count by month
foss  %>%  group_by(year,month) %>%  summarise(n = sum(total), med = median(total), min = min(total), max = max(total), IQR = IQR(total))
ggplot(foss, aes(month,total))+
  geom_bar(stat="identity")+
  geom_text(aes(label=after_stat(y)),stat = 'summary', fun = sum)+facet_wrap(~year)
#total fish count by month, light period
foss  %>%  group_by(year,month, lightperiod) %>% summarise(n = sum(total), med = median(total), min = min(total), max = max(total), IQR = IQR(total))
#total fish count by month, river level grouped
foss  %>%  group_by(year,month, lvl_stage) %>% summarise(n = sum(total), med = median(total), min = min(total), max = max(total), IQR = IQR(total))
ggplot(foss, aes(as.factor(lvl_stage),total))+
  geom_boxplot()+facet_wrap(~year)


#######Post pump operation fish count data#######
#total fish count by day (pre, post operation)
foss  %>% filter(!is.na(p_day))%>% group_by(p_event,p_day) %>%  summarise(n = sum(p_total), med = median(p_total), min = min(p_total), max = max(p_total), IQR = IQR(p_total))
ggplot(foss%>%filter(!is.na(p_day)), aes(p_day,p_total))+
  geom_bar(stat="identity")+
  geom_text(aes(label=after_stat(y)),stat = 'summary', fun = sum)+facet_wrap(~p_event)


#######Floodgate operation fish count data#######
#total fish count by barrier operation (normal, pre-dawn, post-trial)
foss  %>% filter(!is.na(barrier2))%>% group_by(barrier2) %>%  summarise(n = sum(b_total2), med = median(b_total2), min = min(b_total2), max = max(b_total2), IQR = IQR(b_total2))
ggplot(foss%>%filter(!is.na(barrier2)), aes(barrier2,b_total2))+
  geom_boxplot()

#######Environmental data#######
#River level within observed data range - descriptive values presented in Table S1 represent full duration of study
foss  %>%  group_by(year,month) %>%  summarise(med = median(lvl), min = min(lvl), max = max(lvl), IQR = IQR(lvl))
ggplot(foss, aes(month,lvl))+
  geom_boxplot()+facet_wrap(~year)
#temperature within observed data range - descriptive values presented in Table S1 represent full duration of study
foss %>% filter(year!="2017/18")%>% group_by(year,month) %>%  summarise(med = median(temp), min = min(temp), max = max(temp), IQR = IQR(temp))
ggplot(foss%>%filter(year!="2017/18"), aes(month,temp))+
  geom_boxplot()+facet_wrap(~year)


##################################
#### 3 - Statistical analysis ####
##################################

### Data is heavily skewed and non-normal, supported by shapiro.test
### Use non-parametric testing throughout

#######Temporal fish count data#######

#Fish count differences between years
kruskal.test(foss$total~foss$year)
#Fish count differences between months
kruskal.test(foss$total~foss$month)
#Post-hoc dunn's test
dunn.test(x=foss$total, g=foss$month, method="bh")
#Fish count differences between day and night
wilcox.test(foss$total~foss$lightperiod)
#Fish count differences between photoperiods
kruskal.test(foss$total~foss$photo)
#Post-hoc dunn's test
dunn.test(x=foss$total, g=foss$photo, method="bh")

#######Post pump operation fish count data#######

#Fish count pre and post pump operation
wilcox.test(foss$p_total~foss$p_day)
#Individual operations
foss %>% filter(p_event=="Operation \n one") %>%
  wilcox.test(data=.,p_total~p_day)
foss %>% filter(p_event=="Operation \n two") %>%
  wilcox.test(data=.,p_total~p_day)

#######Floodgate operation fish count data#######

#Trial comp post
foss %>% filter(barrier2=="Post barrier trial"|barrier2=="Pre-dawn barrier closure") %>%
  wilcox.test(data=.,b_total2~barrier2)
#Trial comp hydro
foss %>% filter(barrier2=="Normal operation"|barrier2=="Pre-dawn barrier closure") %>%
  wilcox.test(data=.,b_total2~barrier2)
#Post comp hydro
foss %>%  filter(barrier2=="Normal operation"|barrier2=="Post barrier trial") %>%
  wilcox.test(data=.,b_total2~barrier2)

#######################################
#### 3.1 - Modelling - GAM and GLMM####
#######################################

####Generalised Additive Modelling####
#Non-linear time effect

###Construct GAMS to check for model fit before using geom_smooth
###Geom_smooth uses same function as mgcv GAM, 
###Therefore not concerned about plotting lines directly from model as no extra terms used

gam_y <- gam(data=foss, total ~ s(time, by=month), method = "REML")
gam.check(gam_y)
summary(gam_y)

gam_y1 <- gam(data=fossy1, total ~ s(time, by=month), method = "REML")
gam.check(gam_y1)
summary(gam_y1)

gam_y2 <- gam(data=fossy2, total ~ s(time, by=month), method = "REML")
gam.check(gam_y2)
summary(gam_y2)

gam_y3 <- gam(data=fossy3, total ~ s(time, by=month), method = "REML")
gam.check(gam_y3)
summary(gam_y3)


####Generalised Linear Mixed Modelling####
#Linear effects with environmental variables

#Start with GLM before considering adding random factor

#Data profiling

str(foss)

##4647 obs of 19 variables
##each row of total is a separate fish count
#year is a factor with three levels (2017/18, 2018/19, 2019/20)
#month is a factor with six levels
#light period is factor with two levels
#photo is a factor with four levels
#lvl_stage is a factor with four levels
#lvl is a continuous co variate
#temp is a continuous co variate

#check for missing data

colSums(is.na(foss))

#Year one temperature data not available, all other years have no missing data

#check balance of data across categorical variables
table(foss$lvl_stage)
#not well balanced on lvl_stage, however data are random
table(foss$month)
#not well balanced in November and January
table(foss$photo)
#biased towards night-time counts

#check within years

table(fossy1$lvl_stage) #no elevated
table(fossy1$month)
table(fossy1$photo)
table(fossy2$lvl_stage) #no elevated
table(fossy2$month) #well balanced
table(fossy2$photo)
table(fossy3$lvl_stage) #elevated
table(fossy3$month)
table(fossy3$photo)

#Check % of zeros in response variable
sum(foss$total == 0,
    na.rm = TRUE) * 100 / nrow(foss)
#20% of fish counts are 0, potential for zero-inflation
#if amount of observed zeros is larger than number of predicted, 
#then the model is under fitting zeros

##check for multi-collinearity in continuous variables

cor.test(foss$lvl, foss$temp, method="spearman", exact=FALSE)

plot(foss$temp, foss$lvl)

#Small degree of collinearity here, but is obscured by temporal information so unlikely to be relevant
#higher river levels occur at all temps

#check Variance Inflation Factor to determine collineartiy (VIF >3)

vif(glm(total ~ lvl + temp + lvl_stage + lightperiod + month + year,
        family = poisson,
        data = foss))

#confirms that temp shows collinearity with lvl but does not exceed 3, accepted

###Explore plots to determine relationships of interest

plot(foss$total, foss$lvl)# negative linear
plot(fossy2$temp, fossy2$total) #negative linear

plot(foss$total ~ foss$lvl_stage) #stable(reference) higher
plot(foss$total~ foss$photo) #dawn and dusk higher
plot(foss$total~ foss$year) # different variation between years
plot(foss$total~ foss$month) # different variation between months

### Now consider interactions 

#fish count, temperature and photo period
#fish count, temperature and month
#fish count, temperature and year
#fish count, lvl_stage and photo

ggplot(foss, aes(x=temp, y=total)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) +facet_grid(~photo)

#strength of interaction between temperature and night is slightly weaker, 
#possible interaction with light period and temperature

ggplot(fossy2, aes(x=temp, y=total)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) 

ggplot(fossy2, aes(x=lvl, y=total)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) 
#evidence of interaction, but expected as seasonal temperature decline

ggplot(foss, aes(x=temp, y=total)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) +facet_grid(~year)
#Stronger interaction in year two


##determining independence of the dependable variable is difficult as no way of knowing if fish counted belong to one shoal, or many, or belong to one species, or many
##Assume not independent

######Create independent models for each year#####

#year 1

Pois1 <- glm(total ~ lvl+ lvl_stage + month + photo,
             data = fossy1,
             family = poisson(link = log))
summary (Pois1)

#calculate psuedo R2 100 x (null deviance-residual deviance) / null deviance
#31% of variation fish count explained by this model

#year 2

Pois2 <- glm(total ~ temp+ lvl + lvl_stage + month + photo,
             data = fossy2,
             family = poisson(link = log))
summary (Pois2)

#52% of variation fish count explained by this model

#year 3

Pois3 <- glm(total ~ temp+ lvl + lvl_stage + month + photo,
             data = fossy3,
             family = poisson(link = log))
summary (Pois3)

#47% of variation fish count explained by this model

#check models for overdispersion

Pois1$deviance / Pois1$df.residual
Pois2$deviance / Pois2$df.residual
Pois3$deviance / Pois3$df.residual

###all years are overdispersed
###what is the source of the problem? 
#Wrong distribution
#zero inflation 
#Autocorrelation

# try changing distribution to negative binomial to treat the overdispersion
# First exploring with Y2 data as inclusive of temperature

model1 <-glm.nb(total~temp+lvl+lvl_stage+photo,data=fossy2)
summary(model1)

model1$deviance / model1$df.residual

#reduce dispersion to 1.2, but ignores the source of dependency associated with sample month, and potentially individual fish count observations (i.e. OLRE)
#Therefore, need to switch to a model that incorporates random effects -> GLMM

# Quick function to test for overdispersion as can not use pearson resid in GLMM
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#Add individual row ID and incorporate as random effect (OLRE)

foss$ID<-1:nrow(foss)
fossy1$ID<-1:nrow(fossy1)
fossy2$ID<-1:nrow(fossy2)
fossy3$ID<-1:nrow(fossy3)
model2 <-glmmTMB(total~temp+lvl+lvl_stage+photo+(1|ID),data=fossy2, family=poisson)
summary(model1)
overdisp_fun(model1)

# Adding an OLRE introduces underdispersion 

model3 <-glmmTMB(total~temp+lvl+lvl_stage+photo+(1|ID)+(1|month),data=fossy2, family=poisson)
summary(model3)
overdisp_fun(model3)

#adding month does improve the dispersion factor, but still under dispersed
#this suggests either zero inflation, or a negative bionomial distribution required 

###Year one models
#Check % of zeros in response variable
sum(fossy1$total == 0,
    na.rm = TRUE) * 100 / nrow(fossy1)
#13% zero

mod4_y1_nb <-glmmTMB(total~lvl+lvl_stage+photo+(1|month),data=fossy1)
summary(mod4_y1_nb)
overdisp_fun(mod4_y1_nb)

fittedModely1 <- mod4_y1_nb
simuout1 <- simulateResiduals(fittedModel = fittedModely1)
plot(simuout1)

testDispersion(simuout1)
testZeroInflation(simuout1) #zero inflation!!!

mod4_y1_zi <- glmmTMB(total~+lvl+lvl_stage+photo+(1|month), data=fossy1,ziformula=~1, family=nbinom2())
summary(mod4_y1_zi)

#Check goodness of fit via residual diagnostics (Zurr)

fittedModely1_zi <- mod4_y1_zi
simuout2 <- simulateResiduals(fittedModel = fittedModely1_zi)
plot(simuout2)

testDispersion(simuout2) # marginally underdispersed 0.97
testZeroInflation(simuout2) # 0.98 zero inflation treated 


######## mod4_y1_zi accepted as year one model ######


###Year two models
#Check % of zeros in response variable
sum(fossy2$total == 0,
    na.rm = TRUE) * 100 / nrow(fossy2)
#25% zero


mod4_y2_nb <-glmmTMB(total~temp+lvl+lvl_stage+photo+(1|month),data=fossy2)
summary(mod4_y2_nb)

fittedModely2 <- mod4_y2_nb
simuout3 <- simulateResiduals(fittedModel = fittedModely2)
plot(simuout3)

testDispersion(simuout3) 
testZeroInflation(simuout3) # zero inflation!

mod4_y2_zi <- glmmTMB(total~temp+lvl+lvl_stage+photo, data=fossy2,ziformula=~1, family=nbinom2())
summary(mod4_y2_zi)
###Random effect of month may interact oddly with ggpredict, consider removing for plot

fittedModely2_zi <- mod4_y2_zi
simuout4 <- simulateResiduals(fittedModel = fittedModely2_zi)
plot(simuout4)

testDispersion(simuout4) 
testZeroInflation(simuout4) # 1.2 some zero inflation present but much improved over previous model


######## mod4_y2_zi accepted as year two model ######

####Year three models
#Check % of zeros in response variable
sum(fossy3$total == 0,
    na.rm = TRUE) * 100 / nrow(fossy3)
#21% zero

mod4_y3_nb <-glmmTMB(total~temp+lvl+lvl_stage+photo+(1|month),data=fossy3)
summary(mod4_y3_nb)
overdisp_fun(mod4_y3_nb)

fittedModely3 <- mod4_y3_nb
simuout5 <- simulateResiduals(fittedModel = fittedModely3)
plot(simuout5)

testDispersion(simuout5) 
testZeroInflation(simuout5) # Zero inflation!

mod4_y3_zi <- glmmTMB(total~temp+lvl+lvl_stage+photo+(1|month), data=fossy3,ziformula=~1, family=nbinom2())
summary(mod4_y3_zi)


fittedModely3_zi <- mod4_y3_zi
simuout6 <- simulateResiduals(fittedModel = fittedModely3_zi)
plot(simuout6)

testDispersion(simuout6) # 0.89
testZeroInflation(simuout6) # 0.97 zero inflation treated 

######## mod4_y3_zi accepted as year three model ######

##################################################
#### 4 - Visualization and publication figures####
##################################################

#create colour vector

colours <- c("steelblue4", "red4", "palegreen4")

#Create theme

theme_JN <- function(base_size=10){ 
  theme_grey() %+replace%
    theme(
      axis.text = element_text(colour="black"),
      axis.title = element_text(colour="black"),
      axis.ticks = element_line(colour="black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    ) 
}

################################
#####Manuscript - Figure 2######
################################

fosstotal <- ggplot(transform(foss, year = factor(year,labels = c("Year one", "Year two", "Year three"))),
                    aes(x = time, y = total)) +
  geom_rect(aes(xmin=dplyr::lag(time)+1,xmax=(time)+1, ymin=-Inf,ymax=Inf,fill=photo))+
  geom_rect(data=data.frame(year='Year three', month='February'), inherit.aes=F, 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'white')+
  geom_rect(data=data.frame(year='Year one', month=c('September','October')), inherit.aes=F, 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = 'white')+
  scale_fill_manual(values=c('grey90','grey90', 'white', 'grey80')) +
  geom_jitter(size = 0.6, alpha=0.2)+
  geom_smooth(method="gam", colour="black")+
  labs(x = "Time of day (24h)", y = expression(Fish~count~("individuals·2m" ^2~h^-1)))+
  coord_cartesian(xlim = c(1,23))+
  scale_y_continuous(breaks = seq(0, 16, 4), limits=c(0, 16)) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),
                     labels=c("","01:00","","","","05:00","","","","09:00","","","", "13:00","", "","", "17:00", "", "", "", "21:00", "",""),
                     limits=c(-1,24),
                     expand=c(0,0))+
  theme_JN()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) + 
  rremove("legend")+
  facet_grid(factor(month, levels=c('September','October','November','December','January','February'))
             ~(factor(year, levels=c('Year one','Year two','Year three'))))
fosstotal

ggsave(filename="fosstotal.svg", plot=fosstotal, device = "svg",units="cm", width=14,height=20, dpi=600)

################################
#####Manuscript - Figure 3######
################################

###Create summary statistics table
stat.test1 <- foss %>%
  group_by(year) %>%
  wilcox_test(total ~ lvl_stage) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test1
#adding xy values
stat.test1 <- stat.test1 %>% add_xy_position(x = "lvl_stage")
#removing unwanted rows
stat.test1 <-stat.test1[-c(9,11),]
#adjusting y position manually
stat.test1["y.position"] <- c(22,24,22,22,24,22,22,24,22,22)
stat.test1[8, "xmax"] <- 4

box_lvl_stage <-  ggboxplot(
  foss, x = "lvl_stage", y = "total", fill = "lightperiod",
  facet.by = "year", scales = "free_x", outlier.shape = NA,
  width=0.5, bxp.errorbar=TRUE, bxp.errorbar.width=0.3) +
  scale_y_continuous(breaks = seq(0, 20, 4),limits=c(0, 26)) +
  labs(x = "D/S river level (Yorkshire Ouse) (mAOD)")+
  coord_cartesian(clip="off")+
  theme_JN()+
  theme(axis.title.y = element_blank(),
        legend.position = c(.07, 0.6),
        legend.title = element_blank(),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text(size=8),
        legend.background = element_rect(colour = 'black', fill = 'white'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x=unit(0, "lines")) +
  scale_fill_brewer(palette = "Greys")+
  stat_pvalue_manual(stat.test1) +
  ggpubr::rotate_x_text()

box_lvl_stage

#reorder factors for plot
foss$photo <- factor(foss$photo , levels=c("Dawn", "Dusk", "Day", "Night"))

###Create summary statistics table
stat.test2 <- foss %>%
  group_by(year) %>%
  wilcox_test(total ~ photo) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test2
#adding xy values
stat.test2 <- stat.test2 %>% add_xy_position(x = "photo")
#removing unwanted rows
stat.test2 <-stat.test2[-c(2,4,5,8,10,11,14,16,17),]
#adjusting y position manually
stat.test2["y.position"] <- c(22,24,22,22,24,22,22,24,22)


box_photo<- ggboxplot(
  foss, x = "photo", y = "total",
  facet.by = "year", outlier.shape = NA,
  width=0.5, bxp.errorbar=TRUE, bxp.errorbar.width=0.3,
  panel.labs = list(year = c("Year one", "Year two", "Year three"))) +
  scale_y_continuous(breaks = seq(0, 20, 4),limits=c(0, 26)) +
  labs(x = "Photoperiod")+
  coord_cartesian(clip="off")+
  theme_JN()+
  theme(axis.title.y = element_blank(),
        panel.spacing.x=unit(0, "lines"))+
  stat_pvalue_manual(stat.test2)
box_photo

combined_bind<- ggdraw() +
  draw_plot(box_photo, x = 0, y = 0.55, width = 1, height = 0.45) +
  draw_plot(box_lvl_stage, x =0, y=0, width=1, height=0.55)
combined_bind

finalbox <- annotate_figure(combined_bind, 
                            left=text_grob(expression (Fish~count~("individuals·2m" ^2~h^-1)),
                                           rot=90, hjust=0.4))
finalbox

ggsave(filename="finalbox.svg", plot=finalbox, device = "svg",units="cm", width=16,height=14, dpi=600)


################################
#####Manuscript - Figure 4######
################################


#Get predicted model fit line using ggpredict and store as dataframe
plot(ggpredict(mod4_y1_zi, terms = c("lvl[5.2:7.2, by=0.01]")))
plot(ggpredict(mod4_y2_zi, terms = c("lvl[5.2:7.2, by=0.01]")))
plot(ggpredict(mod4_y3_zi, terms = c("lvl[5.2:7.2, by=0.01]")))
plot(ggpredict(mod4_y2_zi, terms = c("temp[5:15, by=0.1]")))
plot(ggpredict(mod4_y3_zi, terms = c("temp[5:15, by=0.1]")))

modlvly1 <-ggpredict(mod4_y1_zi, terms = c("lvl[5.2:7.2, by=0.01]"))
#Convert column group from factor to numeric for row binding
#as.character first required to convert factor to character, and then to numeric
modlvly1$group = as.numeric(as.character(modlvly1$group))
modlvly2 <-ggpredict(mod4_y2_zi, terms = c("lvl[5.2:7.2, by=0.01]"))
#Rename group column to 2
modlvly2["group"] <- 2
modlvly3 <-ggpredict(mod4_y3_zi, terms = c("lvl[5.2:7.2, by=0.01]"))
#rename group column to 3
modlvly3["group"] <- 3
#bind data frames
foss_lvl_mod <- bind_rows(modlvly1, modlvly2, modlvly3)

modpredy2 <- ggpredict(mod4_y2_zi, terms = c("temp[5:15, by=0.1]"))
modpredy2["group"] <- 2
modpredy2$group = as.numeric(as.character(modpredy2$group))
modpredy3 <- ggpredict(mod4_y3_zi,terms = "temp[5:15, by=0.1]")
modpredy3["group"] <- 3
foss_temp_mod <- bind_rows(modpredy2, modpredy3)

foss_temp_mod<- foss_temp_mod %>% add_row(group = 1)

tempregress <-ggplot(foss_temp_mod, aes(x=x, y=predicted, color=as.factor(group),fill=as.factor(group)))+
  geom_line(lwd=1)+
  geom_ribbon(aes(x=x,ymin=conf.low, ymax=conf.high),alpha=0.3, colour="black",linetype=0)+
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  scale_y_continuous(breaks = seq(0, 16, 2) ,limits=c(0,16), expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,14,1), limits=c(5,14), expand=c(0,0)) +
  coord_cartesian(clip = "off") +
  labs(x = "River temperature \n (Yorkshire Ouse) (°C)")+
  theme_JN()+
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(5.5,10,5.5,5.5, "pt"))
tempregress

lvlregress <-ggplot(foss_lvl_mod, aes(x=x, y=predicted, color=as.factor(group),fill=as.factor(group)))+
  geom_line(lwd=1)+
  geom_ribbon(aes(x=x,ymin=conf.low, ymax=conf.high),alpha=0.3, colour="black",linetype=0)+
  scale_fill_manual(values=colours)+
  scale_colour_manual(values=colours)+
  scale_y_continuous(breaks = seq(0, 16, 2) ,limits=c(0,16), expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5.2,7.2,0.4), limits=c(5.2,7.2), expand=c(0,0)) +
  coord_cartesian(clip = "off") +
  labs(x = "D/S river level \n (Yorkshire Ouse) (mAOD)",
       y = expression (Fish~count~("individuals·2m" ^2~h^-1)))+
  theme_JN()+
  theme(legend.position = "none",
        plot.margin = margin(5.5,10,5.5,5.5, "pt"))
lvlregress

##
#combinedregress<- ggdraw() +
#  draw_plot(lvlregress, x = 0, y = 0, width = 0.53, height = 1) +
#  draw_plot(tempregress, x =0.53, y=0, width=0.47, height=1)+
#  draw_plot_label(label = c("a)", "b)", "CI = 95%", "CI = 95%"), 
#                  size = 10,
#                  x = c(0.10, 0.58,0.35, 0.8), 
#                  y = c(0.97, 0.97,0.97, 0.97))


combinedregress <-plot_grid(lvlregress, tempregress,
                            ncol = 2, nrow = 1, rel_widths = c(0.43,0.4), align = "h") 
combinedregress
ggsave(filename="combinedregress.svg", plot=combinedregress, device = "svg",units="cm", width=14,height=7, dpi=600)

################################
#####Manuscript - Figure 5######
################################

###Create summary statistics table
stat.pp <- foss %>%
  filter(!is.na(p_event)) %>%
  group_by(p_event) %>%
  wilcox_test(p_total ~ p_day) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.pp

#Quick plot to get x/y values for table
bxp <- ggplot(foss %>% filter(!is.na(p_event)), aes(x = p_event, fill = p_day, weight=p_total)) + 
  geom_bar(position="dodge", colour="black")
bxp
#adding xy values
stat.pp <- stat.pp %>% add_xy_position(x = "p_event")

#adjusting y position manually and adding KW tests to table 
stat.pp["y.position"] <- 240
stat.test2[c(2,5,8), "p.adj.signif"] <- "KW = <0.001"

#weight is used instead of y to allow geom_bar to count and create y
##inherit.aes added to stat call to prevent missing group from gg call 

barpp <- ggplot(foss %>% filter(!is.na(p_event)), aes(x = p_event)) + 
  geom_bar(aes(fill=p_day, weight=p_total),position="dodge", colour="black") +
  scale_y_continuous(breaks = seq(0, 200, 50), limits=c(0, 250), expand=c(0,0)) +
  labs(x = "Pumping operation", y = expression ("Total fish count"~("individuals·2m" ^-2)))+
  scale_color_grey() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.text.x=element_text(colour="black", size = 10),       
        legend.position = c(.3, 0.8),
        legend.title = element_blank(),
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text(size=8),
        legend.background = element_rect(colour = 'black', fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA)) +
  scale_fill_brewer(palette = "Greys") +
  labs(fill = "Day") +
  stat_pvalue_manual(stat.pp)
barpp

pp_lvl <-ggplot(foss,aes(x=p_date, y=p_lvl), inherit.aes = FALSE) +
  geom_line(color="steelblue4",lwd=1) +
  geom_hline(yintercept=7.6, linetype="dashed", color="black", lwd=0.5) +
  scale_y_continuous(breaks = seq(5, 8.5, 0.5),limits=c(5,9),expand=c(0,0), position = "right") +
  scale_x_continuous(breaks = c(25,	121,	217,	313,	409,	505,	601,	697),
                     label= c("02 Oct","06 Oct", "10 Oct","14 Oct", "18 Oct","22 Oct", "26 Oct","30 Oct"),
                     limits=c(1,744),expand=c(0,0)) +
  labs(x = "Date", y = "D/S river Level \n (Yorkshire Ouse) (mAOD)")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        panel.border = element_rect(colour = "black", fill=NA)) +
  ggpubr::rotate_x_text()
pp_lvl

# using cowplot is more effective, as can define relative width instead of specifying x and y 
# Also solves alignment issue with x axis 

pp_bind <-plot_grid(barpp, pp_lvl,
                    ncol = 2, nrow = 1, rel_widths = c(3.5, 6.5),align = "h")

################################
#####Manuscript - Figure 6######
################################


###summarise data, create ID for each row for continuous axis 

hist_sum <- foss %>% filter(barrier2=="Pre-dawn barrier closure")%>%
  group_by(b_time) %>% 
  summarise(
    med = median(b_total2))
hist_sum$barrier <-'Pre-dawn barrier closure'
hist_sum$x_ID<-1:nrow(hist_sum)

hist_sum2 <- foss %>% filter(barrier2=="Normal operation")%>%
  group_by(b_time) %>% 
  summarise(
    med = median(b_total2))
hist_sum2$barrier <-'Normal operation'
hist_sum2$x_ID<-1:nrow(hist_sum2)

hist_sum3 <- foss %>% filter(barrier2=="Post barrier trial")%>%
  group_by(b_time) %>% 
  summarise(
    med = median(b_total2))
hist_sum3$barrier <-'Post barrier trial'
hist_sum3$x_ID<-1:nrow(hist_sum3)

hist_sum_bind <- rbind(hist_sum, hist_sum2, hist_sum3)
hist_sum_bind$barrier <- as.factor(hist_sum_bind$barrier)

#######Create histogram

histbox <- ggplot(hist_sum_bind%>%filter(barrier=="Pre-dawn barrier closure"), aes(x=x_ID, y = med)) +
  geom_bar(stat="identity", width=1, alpha=0.2) +
  scale_y_continuous(breaks=seq(0, 18, 2),
                     limits=c(0, 18),
                     expand=c(0,0)) +
  scale_x_continuous(breaks = c(1, 5, 9, 13),
                     label= c("07:00","08:00","09:00","10:00"),
                     limits=c(0,14),
                     expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")  +
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") 
histbox


#####Create and extract only the boxplot, to overlay on histogram

boxplotbarrier <- ggplot(foss%>%filter(barrier2=="Pre-dawn barrier closure"), aes( y = b_total2)) +
  stat_boxplot(geom = "errorbar", width=1, position = position_dodge(0.5))+
  geom_boxplot(outlier.shape=NA,width=3,coef = 0, fill='white') +
  scale_y_continuous(breaks = seq(0, 18, 2),limits=c(0, 18), expand=c(0,0)) +
  theme_bw() + 
  xlim(-6, 6)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y=element_blank(),
        strip.text.y = element_blank())+
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.axis")+
  rremove("y.ticks")+
  rremove("y.text")+
  rremove("legend")
boxplotbarrier

aligned_plots1 <- align_plots(histbox, boxplotbarrier, align="hv", axis="tblr")
topbox1 <-ggdraw(aligned_plots1[[1]]) + draw_plot(aligned_plots1[[2]])
topbox1

barrierplot<- ggplot(foss%>%filter(barrier=="Pre-dawn barrier closure"), aes(x=b_timecode, y = b_total)) +
  annotate("rect",xmin=32,xmax=35, ymin=-Inf,ymax=Inf,fill="grey80")+
  annotate("rect",xmin=56,xmax=59, ymin=-Inf,ymax=Inf,fill="grey80")+
  annotate("rect",xmin=80,xmax=83, ymin=-Inf,ymax=Inf,fill="grey80")+
  annotate("rect",xmin=104,xmax=107, ymin=-Inf,ymax=Inf,fill="grey80")+
  annotate("rect",xmin=128,xmax=131, ymin=-Inf,ymax=Inf,fill="grey80")+
  geom_bar(stat="identity", width=1.1, colour=NA) +
  scale_y_continuous(breaks=seq(0, 18, 2),
                     limits=c(0, 18),
                     expand=c(0,0)) +
  scale_x_continuous(breaks = c(25,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,
                                92,96,100,104,108,112,116,120,124,128,132,136,140,144),
                     limits=c(20,148),
                     expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +   rremove("x.ticks")

barrierplot

names <- c("D/S YO level", "Operational level ")

barrierplot_hyd<- ggplot(foss%>%filter(barrier=="Pre-dawn barrier closure" | barrier=="Operational level"), 
                         aes(x=b_timecode, y=b_lvl, col=barrier, linetype=barrier, size=barrier)) +
  geom_line() +
  scale_color_manual(values=c("steelblue4", "black"),labels=names)+
  scale_linetype_manual(values=c(1,2),labels=names)+
  scale_size_manual(values=c(1, 0.5),labels=names)+
  scale_y_continuous(breaks=seq(5, 8, 0.5),
                     limits=c(5, 8),
                     expand=c(0,0), position = "right") +
  scale_x_continuous(breaks = c(25,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,
                                92,96,100,104,108,112,116,120,124,128,132,136,140,144,
                                limits=c(20,148)),
                     expand=c(0,0))+
  theme_half_open() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = c(.65, 0.7),
        legend.title = element_blank(),
        legend.key.height = unit(0.4,"cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text(size=8),
        legend.background = element_rect(colour = 'black', fill = 'white'))+
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") 

barrierplot_hyd

aligned_plots2 <- align_plots(barrierplot, barrierplot_hyd, align="hv", axis="tblr")
topbox2 <-ggdraw(aligned_plots2[[1]]) + draw_plot(aligned_plots2[[2]])
topbox2

top_bind<- ggdraw() +
  draw_plot(topbox1, x = 0, y = 0, width = 0.3, height = 1) +
  draw_plot(topbox2, x =0.3, y=0, width=0.7, height=1)
top_bind

histbox2 <- ggplot(hist_sum_bind%>%filter(barrier=="Post barrier trial"), aes(x=x_ID, y = med)) +
  geom_bar(stat="identity", width=1, alpha=0.2) +
  scale_y_continuous(breaks=seq(0, 18, 2),
                     limits=c(0, 18),
                     expand=c(0,0)) +
  scale_x_continuous(breaks = c(1, 5, 9, 13),
                     limits=c(0,14),
                     expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")  +
  rremove("x.ticks") 
histbox2

boxplotbarrier2 <- ggplot(foss%>%filter(barrier2=="Post barrier trial"), aes( y = b_total2)) +
  stat_boxplot(geom = "errorbar", width=1, position = position_dodge(0.5))+
  geom_boxplot(outlier.shape=NA,width=3,coef = 0, fill='white') +
  scale_y_continuous(breaks = seq(0, 18, 2),limits=c(0, 18), expand=c(0,0)) +
  theme_bw() + 
  xlim(-6, 6)+
  theme_half_open() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y=element_blank(),
        strip.text.y = element_blank())+
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.axis")+
  rremove("y.ticks")+
  rremove("y.text")+
  rremove("legend")
boxplotbarrier2

aligned_plots3 <- align_plots(histbox2, boxplotbarrier2, align="hv", axis="tblr")
middlebox1 <-ggdraw(aligned_plots3[[1]]) + draw_plot(aligned_plots3[[2]])
middlebox1

barrierplot2<- ggplot(foss%>%filter(barrier=="Post barrier trial"), aes(x=b_timecode, y = b_total)) +
  geom_bar(stat="identity", width=1.1,colour=NA,) +
  scale_y_continuous(breaks=seq(0, 18, 2),
                     limits=c(0, 18),
                     expand=c(0,0)) +
  scale_x_continuous(breaks = c(32,40,48,56,64,72,80,88,
                                96,104,112,120,128,136,144),
                     limits=c(20,148),
                     expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +   rremove("x.ticks")

barrierplot2

barrierplot_hyd2<- ggplot(foss%>%filter(barrier=="Post barrier trial" | barrier=="Operational level"), 
                          aes(x=b_timecode, y=b_lvl, col=barrier, linetype=barrier, size=barrier)) +
  geom_line() +
  scale_color_manual(values=c("steelblue4", "black"),labels=names)+
  scale_linetype_manual(values=c(1,2),labels=names)+
  scale_size_manual(values=c(1, 0.5),labels=names)+
  scale_y_continuous(breaks=seq(5, 8, 0.5),
                     limits=c(5, 8),
                     expand=c(0,0), position = "right") +
  scale_x_continuous(breaks = c(25,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,
                                92,96,100,104,108,112,116,120,124,128,132,136,140,144,
                                limits=c(20,148)),
                     expand=c(0,0))+
  theme_half_open() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks")

barrierplot_hyd2

aligned_plots4 <- align_plots(barrierplot2, barrierplot_hyd2, align="hv", axis="tblr")
middlebox2 <-ggdraw(aligned_plots4[[1]]) + draw_plot(aligned_plots4[[2]])

middle_bind<- ggdraw() +
  draw_plot(middlebox1, x = 0, y = 0, width = 0.3, height = 1) +
  draw_plot(middlebox2, x =0.3, y=0, width=0.7, height=1)
middle_bind


######### Bototm plot
histbox3 <- ggplot(hist_sum_bind%>%filter(barrier=="Normal operation"), aes(x=x_ID, y = med)) +
  geom_bar(stat="identity", width=1, alpha=0.2) +
  scale_y_continuous(breaks=seq(0, 18, 2),
                     limits=c(0, 18),
                     expand=c(0,0)) +
  scale_x_continuous(breaks = c(1, 5, 9, 13),
                     label= c("07:00","08:00","09:00","10:00"),
                     limits=c(0,14),
                     expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.text.x=element_text(colour="black", size = 10,  angle = 90,vjust = 0.5),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") 
histbox3

boxplotbarrier3 <- ggplot(foss%>%filter(barrier2=="Normal operation"), aes( y = b_total2)) +
  stat_boxplot(geom = "errorbar", width=1, position = position_dodge(0.5))+
  geom_boxplot(outlier.shape=NA,width=3,coef = 0, fill='white') +
  scale_y_continuous(breaks = seq(0, 18, 2),limits=c(0, 18), expand=c(0,0)) +
  theme_bw() + 
  xlim(-6, 6)+
  theme_half_open() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y=element_blank(),
        strip.text.y = element_blank())+
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.axis")+
  rremove("y.ticks")+
  rremove("y.text")+
  rremove("legend")
boxplotbarrier3

aligned_plots5 <- align_plots(histbox3, boxplotbarrier3, align="hv", axis="tblr")
bottombox1 <-ggdraw(aligned_plots5[[1]]) + draw_plot(aligned_plots5[[2]])
bottombox1

barrierplot3<- ggplot(foss%>%filter(barrier=="Normal operation"), aes(x=b_timecode, y = b_total)) +
  geom_bar(stat="identity", width=1.1,colour=NA,) +
  scale_y_continuous(breaks=seq(0, 18, 2),
                     limits=c(0, 18),
                     expand=c(0,0)) +
  scale_x_continuous(breaks = c(32,40,48,56,64,72,80,88,
                                96,104,112,120,128,136,144),
                     label= c("07:00","15:00","23:00",
                              "07:00","15:00","23:00",
                              "07:00","15:00","23:00",
                              "07:00","15:00","23:00",
                              "07:00","15:00","23:00"),
                     limits=c(20,148),
                     expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(colour="black", size = 10, angle = 90,vjust = 0.5),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")

barrierplot3

barrierplot_hyd3<- ggplot(foss%>%filter(barrier=="Normal operation" | barrier=="Operational level"), 
                          aes(x=b_timecode, y=b_lvl, col=barrier, linetype=barrier, size=barrier)) +
  geom_line() +
  scale_color_manual(values=c("steelblue4", "black"),labels=names)+
  scale_linetype_manual(values=c(1,2),labels=names)+
  scale_size_manual(values=c(1, 0.5),labels=names)+
  scale_y_continuous(breaks=seq(5, 8, 0.5),
                     limits=c(5, 8),
                     expand=c(0,0), position = "right") +
  scale_x_continuous(breaks = c(25,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,
                                92,96,100,104,108,112,116,120,124,128,132,136,140,144,
                                limits=c(20,148)),
                     expand=c(0,0))+
  theme_half_open() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks")

barrierplot_hyd3

aligned_plots6 <- align_plots(barrierplot3, barrierplot_hyd3, align="hv", axis="tblr")
bottombox2 <-ggdraw(aligned_plots6[[1]]) + draw_plot(aligned_plots6[[2]])
bottombox2

bottom_bind<- ggdraw() +
  draw_plot(bottombox1, x = 0, y = 0, width = 0.3, height = 1) +
  draw_plot(bottombox2, x =0.3, y=0, width=0.7, height=1)
bottom_bind

combinedbarrier <-ggarrange(top_bind, bottom_bind, ncol=1, nrow=2)
combinedbarrier

combined_bind<- ggdraw() +
  draw_plot(top_bind, x = 0, y = 0.69, width = 1, height = 0.31) +
  draw_plot(middle_bind, x =0, y=0.38, width=1, height=0.31)+
  draw_plot(bottom_bind, x =0, y=0, width=1, height=0.38)
combined_bind

finalbarrier <- annotate_figure(combined_bind, 
                                left=text_grob(expression("Fish count" ~ ("individuals·2m" ^-2)),
                                               rot=90, hjust=0.4),
                                right=text_grob(expression("D/S river level (Yorkshire Ouse) (mAOD)"),
                                                rot=90, hjust=0.4),
                                bottom=text_grob("Time of day (24h)", hjust=0.5))
finalbarrier
finalbarrier2 <- finalbarrier + draw_plot_label(label = c("a)", "b)", "c)", "d)"), size = 12,
                                                x = c(0.09, 0.34,0.09, 0.34), y = c(0.98, 0.98,0.54, 0.54))
finalbarrier2

ggsave(filename="finalbarrier.svg", plot=finalbarrier, device = "svg",units="cm", width=16,height=12, dpi=600)
