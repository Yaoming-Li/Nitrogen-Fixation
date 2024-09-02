envANF <- read.csv('Data_VPA.csv',header=T)

library(car)
envanf <- data.frame(scale(envANF))
library(performance)
mod<-lm(ANF ~ ., data=envanf)
summary(mod)     
vif_value <- vif(mod)
vif_value
#remove fators that vif>10
dataselect<- envanf[c('Shannon','MAP','MAT','NPP','SOC','TN','TP','SM','DOC','NH4','NO3','pH','Fe','Mo','MDS1','ANF')]
mod1 <- lm(ANF ~ .,data=dataselect)
library(MuMIn)
options(na.action = "na.fail")
dd<-dredge(mod1)
de<-model.avg(dd, subset = delta < 2)
summary(de)
#best model(Minimum AIC value)
cc <- summary(get.models(dd, 1)[[1]])
library(ggthemes)
library(ggplot2)
library(ggpubr)
dataset <- read.csv('Result.csv',header=T)
#Definition sequence
dataset$index <- factor(dataset$index, levels = c('Mo','Fe','NH4','pH','TN','TP','MDS1','MAP'))
a1<-ggplot(dataset, aes(Estimate, index))+
  geom_point(size=4)+
  geom_errorbarh(aes(xmax = Estimate +`SE`, xmin = Estimate -`SE`),size= 0.5,height = 0.2) +
  scale_x_continuous(limits= c(-1, 1))+ 
  geom_vline(aes(xintercept = 0),color="gray",linetype="dashed", size = 0.6) +
  geom_text(aes(label=sig,x=Estimate +`SE`+0.1))+
  xlab('Effect size ')+ 
  ylab(' ')+
  theme_few()+
  theme(axis.text.x = element_text(size = 24, color = "black"))+
  theme(axis.text.y = element_text(size = 24, color = "black"))+
  theme(title=element_text(size=24))+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
a1

library(rdacca.hp)
ANF<-dataselect$ANF
env<-dataselect[c('MAP','TN','TP','NH4','pH','Fe','Mo','MDS1')]
hp<-rdacca.hp(ANF,env,method = 'RDA', type = 'adjR2', scale = FALSE)
hp
plot(hp)
hp.ie <- data.frame(hp$Hier.part, check.names = FALSE)
hp.ie$env <- rownames(hp.ie)
hp.ie$env <- factor(hp.ie$env, levels = rev(hp.ie$env))
hp.ie$exp_var <- ''
hp.ie
p2 <- ggplot(hp.ie, aes(x = exp_var, y = Individual*100, fill = env)) + 
  geom_col(width = 0.6) +  # relative effect of estimates
  scale_fill_manual(values = c('#d83c51','#ef6e46','#fdab68','#fae18e',
                               '#e3f19c','#b0dda2','#66c3a8','#3188ba','#5761C8',
                               '#817f80','#a0c155','#f7e599','#f9e4ce')) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.y = element_line(color = 'gray30')) +
  labs(x = paste('adjR2 =', round(hp$Total_explained_variation, 3)), y = 'relative effect of estimates (%)', fill = '') +
  scale_x_discrete(position = 'top') +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

p2
#I.perc(%)analyzethe contribution of each x
#ouput the overall R2 for the model
set.seed(123)
permu_hp <- permu.hp(ANF, env, method = 'RDA', type = 'adjR2', permutations = 999)
permu_hp
