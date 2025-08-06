
library(readr)
library(lmerTest)
library(sjPlot)
library(simr)
library(car)
library(tidyverse) # data wrangling and visualization
library(sjPlot)    # to visualizing mixed-effects models
library(effects)   # to visualizing mixed-effects models
library(lme4)      # "golden standard" for mixed-effects modelling in R (no p-values)
library(lmerTest)  # p-values for MEMs based on the Satterthwaite approximation
library(report)    # mainly for an "report" function
library(emmeans)   # post-hoc analysis
library(knitr)     # beautifying tables
library(sjstats)   # ICC - intraclass-correlation coefficient
library(caret)     # ML, model comparison & utility functions
library(multcomp)
library(ggeffects)
library(parameters)



df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df =  subset(df,df$AI != 'NaN')
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))

## AI
# pshare
m_pshare_1 <- lmer(pshare ~  AI + site + (1|subj), data=df, REML = FALSE)
tab_model(m_pshare_1)
summary(m_pshare_1)
anova(m_pshare_1)

m_pshare_2 <- lmer(pshare ~ AI + age + iq + srs+ site + (1|subj), data=df, REML = FALSE) # winning model
tab_model(m_pshare_2)
summary(m_pshare_2)
anova(m_pshare_2)

m_pshare_3 <- lmer(pshare ~ AI + poly(age,degree = 2,raw = TRUE) + iq + (1|subj), data=df, REML = FALSE)
tab_model(m_pshare_3)
summary(m_pshare_3)
anova(m_pshare_3)

m_pshare_4 <- lmer(pshare ~ AI * (iq + age + srs) + site + (1|subj), data=df, REML = FALSE)
tab_model(m_pshare_4)
summary(m_pshare_4)
anova(m_pshare_4)

anova(m_pshare_2,m_pshare_4)



# ashare
m_ashare_1 <- lmer(ashare ~ AI + site + (1|subj), data=df, REML = FALSE)
m_ashare_1 <- lmer(ashare ~ AI + (1|subj), data=df, REML = FALSE)
tab_model(m_ashare_1)
summary(m_ashare_1)
anova(m_ashare_1)

m_ashare_2 <- lmer(ashare ~ age + iq + srs+AI + site +(1|subj), data=df, REML = FALSE)
m_ashare_2 <- lmer(ashare ~ age + iq + srs+AI +(1|subj), data=df, REML = FALSE)
tab_model(m_ashare_2)
summary(m_ashare_2)
anova(m_ashare_2)

m_ashare_3 <- lmer(ashare ~ AI * (iq + age + srs) + site + (1|subj), data=df, REML = FALSE) ## winning model
tab_model(m_ashare_3)
summary(m_ashare_3)
anova(m_ashare_3)

m_ashare_4 <- lmer(ashare ~ iq + poly(age,degree = 2,raw = TRUE)*AI + (1|subj), data=df, REML = FALSE)
tab_model(m_ashare_4)
summary(m_ashare_4)
anova(m_ashare_4)

anova(m_ashare_1,m_ashare_2,m_ashare_3,m_ashare_4)

df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df =  subset(df,df$AI != 'NaN')
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
# difft4t
m_dshare_1 <- lmer(nochoice ~ AI + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1)
summary(m_dshare_1)
anova(m_dshare_1)

m_dshare_2 <- lmer(keepa ~ AI + iq + age + srs + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_2)
summary(m_dshare_2)
anova(m_dshare_2)

m_dshare_3 <- lmer(nochoicea ~ AI*(iq + age + srs) + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_3)
summary(m_dshare_3)
anova(m_dshare_3)

anova(m_dshare_2,m_dshare_3)

##social
# pshare
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])

df =  subset(df,df$AI != 'NaN')
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))

m_pshare_1 <- lmer(pshare ~ social + site + (1|subj), data=df, REML = FALSE) 
tab_model(m_pshare_1)
summary(m_pshare_1)
anova(m_pshare_1)

m_pshare_2 <- lmer(pshare ~ social + age + iq + srs + site + (1|subj), data=df, REML = FALSE) ##winning model
tab_model(m_pshare_2)
summary(m_pshare_2)
anova(m_pshare_2)

m_pshare_3 <- lmer(pshare ~ social * (iq + age+srs) + site + (1|subj), data=df, REML = FALSE)
tab_model(m_pshare_3)
summary(m_pshare_3)
anova(m_pshare_3)

anova(m_pshare_1,m_pshare_2,m_pshare_3)

df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df =  subset(df,df$AI != 'NaN')
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))


#ashare 
m_ashare_s1 <- lmer(ashare ~ social + site + (1|subj), data=df, REML = FALSE) 
tab_model(m_ashare_s1)
summary(m_ashare_s1)
anova(m_ashare_s1)

m_ashare_s2 <- lmer(ashare ~ social + age+iq+srs+site+(1|subj), data=df, REML = FALSE)
tab_model(m_ashare_s2)
summary(m_ashare_s2)
anova(m_ashare_s2)

m_ashare_s3 <- lmer(ashare ~ social * (age + iq + srs) + site +(1|subj), data=df, REML = FALSE)
tab_model(m_ashare_s3)
summary(m_ashare_s3)
anova(m_ashare_s3)

anova(m_ashare_s2,m_ashare_s1)

summary(glht(m_pshare_2, linfct = mcp(AI="Tukey")))
g <- ggpredict(m_pshare_2, terms = "AI")
plot(g)

df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df =  subset(df,df$AI != 'NaN')
df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))

# difft4t
m_dshare_1 <- lmer(keepp ~ AI + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1)
summary(m_dshare_1)
anova(m_dshare_1)
m_dshare_1a <- lmer(keepp ~ AI + site + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1a)
summary(m_dshare_1a)
anova(m_dshare_1a)
m_dshare_1b <- lmer(keepa ~ social + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1b)
summary(m_dshare_1b)
anova(m_dshare_1b)
m_dshare_1c <- lmer(keepa ~ social + site + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1c)
summary(m_dshare_1c)
anova(m_dshare_1c)

m_dshare_2 <- lmer(nochoice ~ social + iq + age + srs + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_2)
summary(m_dshare_2)
anova(m_dshare_2)

m_dshare_3 <- lmer(nochoice ~ social*(iq + age + srs) + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_3)
summary(m_dshare_3)
anova(m_dshare_3)

m_dshare_1 <- lmer(nochoicea ~ social + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1)
summary(m_dshare_1)
anova(m_dshare_1)

m_dshare_2 <- lmer(nochoicea ~ social + iq + age + srs + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_2)
summary(m_dshare_2)
anova(m_dshare_2)

m_dshare_3 <- lmer(nochoicea ~ social*(iq + age + srs) + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_3)
summary(m_dshare_3)
anova(m_dshare_3)
anova(m_dshare_2,m_dshare_3)

m_dshare_1 <- lmer(sharep ~ social + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1)
summary(m_dshare_1)
anova(m_dshare_1)

m_dshare_2 <- lmer(sharep ~ social + iq + age + srs + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_2)
summary(m_dshare_2)
anova(m_dshare_2)

m_dshare_3 <- lmer(sharep ~ social*(iq + age + srs) + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_3)
summary(m_dshare_3)
anova(m_dshare_3)

m_dshare_1 <- lmer(sharea ~ social + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_1)
summary(m_dshare_1)
anova(m_dshare_1)

m_dshare_2 <- lmer(sharea ~ social + iq + age + srs + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_2)
summary(m_dshare_2)
anova(m_dshare_2)

m_dshare_3 <- lmer(sharep ~ social*(iq + age + srs) + (1|subj), data=df, REML = FALSE)
tab_model(m_dshare_3)
summary(m_dshare_3)
anova(m_dshare_3)
#####################
#glmer
###AI
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
df[,'social'] <- factor(df[,'social'], levels = c('fix30','fix70','computer','interactive'))

df =  subset(df,df$AI ==0)
m_adap <- glmer(pshare ~ tit4tat*srs+
              (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_adap)
summary(m_adap)
anova(m_adap)


m_adap1 <- glmer(pshare ~ iq + AI*tit4tat + age + 
                  (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa")) 
tab_model(m_adap1)
summary(m_adap1)
anova(m_adap1)

m_adap2 <- glmer(pshare ~ iq + AI*tit4tat*age + 
                   (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa")) ## winning model
tab_model(m_adap2)
summary(m_adap2)
anova(m_adap2)

anova(m_adap,m_adap1,m_adap2)

###nonsocial vs fix70
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df[,'social'] <- factor(df[,'social'], levels = c('fix70','computer'))

df[,'age']<-scale(df[,'age'], center = TRUE, scale = TRUE)
df[,'iq']<-scale(df[,'iq'], center = TRUE, scale = TRUE)

m_adap_0 <- glmer(pshare ~ age + iq + social*tit4tat +
                  (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))## winning model
tab_model(m_adap_0)
summary(m_adap_0)
anova(m_adap_0)

m_adap_1 <- glmer(pshare ~ iq + social*tit4tat*age +
                    (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_adap_1)
summary(m_adap_1)
anova(m_adap_1)

anova(m_adap_0,m_adap_1)

###social vs nonsocial
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))

df[,'age']<-scale(df[,'age'], center = TRUE, scale = TRUE)
df[,'iq']<-scale(df[,'iq'], center = TRUE, scale = TRUE)

m_adap_2 <- glmer(pshare ~ age + iq + social*tit4tat +
                    (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_adap_2)
summary(m_adap_2)
anova(m_adap_2)

m_adap_3 <- glmer(pshare ~ iq + social*tit4tat*age + 
                    (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))## winning model
tab_model(m_adap_3)
summary(m_adap_3)
anova(m_adap_3)

m_adap_4 <- glmer(pshare ~ age + iq  + social*tit4tat+
                    (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_adap_4)
summary(m_adap_4)
anova(m_adap_4)

anova(m_adap_2,m_adap_3,m_adap_4)

########trial
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))


m_trial_1 <- glmer(ashare ~ trial*AI +
                    (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_trial_1)
summary(m_trial_1)
anova(m_trial_1)

m_trial2 <- glmer(pshare ~ trial*AI + age + iq + srs +
                     (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_trial2)
summary(m_trial2)
anova(m_trial2)

m_trial2 <- glmer(pshare ~ trial*AI*(age + iq + srs) +
                    (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_trial2)
summary(m_trial2)
anova(m_trial2)

df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))

df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))
m_trial_2 <- glmer(pshare ~ trial*social +age + iq + srs +
                     (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_trial_2)
summary(m_trial_2)
anova(m_trial_2)

###pshare after share
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/lme_data.csv")
df[,'pshare'] <- as.numeric(df[,'pshare'])
df[,'subj']<-factor(df[,'subj'])
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
df =  subset(df,df$tit4tat == 1)
df[,'social'] <- factor(df[,'social'], levels = c('fix70','interactive'))

m_adap <- glmer(pshare ~ age+iq + social +
                  (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m_adap)
summary(m_adap)
anova(m_adap)

df[,'social'] <- factor(df[,'social'], levels = c('fix30','fix70','computer','interactive','real'))
df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive','real'))
df[,'social'] <- factor(df[,'social'], levels = c('Fix','computer','interactive'))
df[,'social'] <- factor(df[,'social'], levels = c('interactive','real'))

df[,'social'] <- factor(df[,'social'], levels = c('fix70','computer'))


df =  subset(df,df$AI == 3)

m1 <- glmer(pshare ~ tit4tat +
              (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m1)
summary(m1)
anova(m1)




#############
#age invest
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/lme_data.csv")
df =  subset(df,df$AI != 'NaN')
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
df[,'game_n'] <- factor(df[,'game_n'], levels = c('3','4','1','2'))

cbPalette <- c("#0072B2", "#D55E00", "#E69F00", "#CC79A7")

ggplot(df,aes(x=trial,y=pshare,color=game_n))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=game_n),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  theme(axis.line.x = element_line(color = 'black',size = 1),
        axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  labs(fill="game_n", x="trial", y="pshare investor")+
  # scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-0.1,1.1))+
  coord_cartesian(ylim=c(-0.05,1.05),xlim=c(1,20))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))



####################
df<- read.csv("/Users/wl/Documents/lme_data.csv")
df[,'social'] <- factor(df[,'social'], levels = c('fix30','fix70','computer','interactive'))
df[,'subj']<-factor(df[,'subj'])
df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive','real'))

model1 <- lmer(pshare ~ social*AI+ (1 | subj), data = df)
tab_model(model1)
summary(model1)

model2 <- lmer(ashare ~ social + (1 | subj), data = df)
tab_model(model2)
summary(model2)
##################
df<- read.csv("/Users/wl/Documents/lme_data_tit4tat.csv")
df[,'subj']<-factor(df[,'subj'])
df[,'social']<-factor(df[,'social'])

m2 <- glmer(investor ~ prediction + (1 | subj), data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
tab_model(m2)
summary(m2)
#############
#age invest
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data_within.csv")
df =  subset(df,df$AI != 'NaN')
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
cbPalette <- c("#0072B2", "#D55E00", "#E69F00", "#CC79A7")

df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))
df =  subset(df,df$social != 'NaN')

df$age <- df$age / 12

ggplot(df,aes(x=iq,y=ashare,color=AI))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=AI),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
 # theme(axis.line.x = element_line(color = 'black',size = 1),
 #       axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  labs(fill="AI", x="IQ", y="ashare investor")+
  scale_x_continuous(expand = c(0,0), limits=c(60,150)) +
  scale_y_continuous(expand = c(0,0), limits=c(-0.1,5.1))+
  coord_cartesian(ylim=c(-0.2,5.2),xlim=c(60,150))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

##age invest
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data.csv")
df =  subset(df,df$AI != 'NaN')
df[,'AI'] <- factor(df[,'AI'], levels = c('0','1'))
cbPalette <- c("#0072B2", "#D55E00", "#E69F00", "#CC79A7")

df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))
df =  subset(df,df$social != 'NaN')
df$age <- df$age / 12


ggplot(df,aes(x=age,y=difft4t,color=social))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=social),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  # theme(axis.line.x = element_line(color = 'black',size = 1),
  #       axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  labs(fill="social", x="age", y="%tit4tat")+
  scale_x_continuous(expand = c(0,0), limits=c(8,18)) +
  scale_y_continuous(expand = c(0,0), limits=c(-0.1,1.1))+
  coord_cartesian(ylim=c(-0.05,1.05),xlim=c(8,18))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))


df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data.csv")

df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))
df =  subset(df,df$social != 'NaN')
cbPalette <- c("#0072B2", "#D55E00", "#E69F00", "#CC79A7")

ggplot(df,aes(x=srs,y=sharea,color=social))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=social),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  # theme(axis.line.x = element_line(color = 'black',size = 1),
  #       axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  labs(fill="social", x="srs", y="share mount")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

ggplot(df,aes(x=srs,y=sharea))+
  geom_point(size=5)+
  geom_smooth(method="lm",size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  # theme(axis.line.x = element_line(color = 'black',size = 1),
  #       axis.line.y = element_line(color = 'black',size = 1))+
  labs(x="srs", y="%tit4tat")+
  scale_colour_manual(values=cbPalette)+
  scale_y_continuous(expand = c(0,0), limits=c(-0.1,1.1))+
  coord_cartesian(ylim=c(-0.01,1.01))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

############################################
df<- read.csv("/Users/wl/Documents/DSN_Lab/Inperson_trustgame/m_lme_data_model.csv")
df[,'subj']<-factor(df[,'subj'])
df[,'social'] <- factor(df[,'social'], levels = c('computer','interactive'))
df =  subset(df,df$social == 'interactive')
#df =  subset(df,df$AI != 'adap')

## AI
# pshare
m_pshare_1 <- lmer(p2 ~  AI +(1|subj), data=df, REML = FALSE)
tab_model(m_pshare_1)
summary(m_pshare_1)
anova(m_pshare_1)

m_pshare_1 <- lmer(p1 ~  AI + age + iq + (1|subj), data=df, REML = FALSE)
tab_model(m_pshare_1)
summary(m_pshare_1)
anova(m_pshare_1)

m_pshare_2 <- lm(p1 ~ poly(age,degree = 2,raw = TRUE)+ iq + srs, data=df, REML = FALSE)
tab_model(m_pshare_2)
summary(m_pshare_2)
anova(m_pshare_2)

m_pshare_3 <- lm(p2 ~ poly(age,degree = 2,raw = TRUE)+ iq + srs, data=df, REML = FALSE)
tab_model(m_pshare_3)
summary(m_pshare_3)
anova(m_pshare_3)

m_pshare_4 <- lm(BIC ~ poly(age,degree = 2,raw = TRUE)+ iq + srs, data=df, REML = FALSE)
tab_model(m_pshare_4)
summary(m_pshare_4)
anova(m_pshare_4)

df$age <- df$age / 12

cbPalette <- c("#D55E00", "#E69F00", "#CC79A7")
ggplot(df,aes(x=age,y=p1,color=social))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=social),formula = y ~ x + I(x^2),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  # theme(axis.line.x = element_line(color = 'black',size = 1),
  #       axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  labs(fill="social", x="age", y="decay p")+
  coord_cartesian(ylim=c(-0.01,1.01))+
  
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

cbPalette <- c("#D55E00", "#E69F00", "#CC79A7")
ggplot(df,aes(x=srs,y=p2,color=social))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=social),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  # theme(axis.line.x = element_line(color = 'black',size = 1),
  #       axis.line.y = element_line(color = 'black',size = 1))+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  labs(fill="social", x="srs", y="decision noise")+
  scale_y_continuous(expand = c(0,0), limits=c(-1.1,1.1))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))