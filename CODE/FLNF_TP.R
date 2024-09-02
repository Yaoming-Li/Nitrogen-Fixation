################ scripts for FLNF data analysis ######################
###################     authored by Yaoming      ##################### 
###################         Feb-15 2024          ##################### 
sessionInfo()
# R version 4.3.1
rm(list = ls())
getwd()
setwd("K:\\Research\\�о���\\����\\��ظ�ԭnifH\\Soil_Nif\\R\\NifH")
dir.create('data_raw')
dir.create('tem')
dir.create('Figs')
dir.create('Results')

#import micro & soil physiochemical data
library(readxl)
library(vegan)
library(reshape2)
FM <- as.data.frame(read_excel("data_raw/micro.xlsx", sheet = 2, col_names = TRUE)) 
rownames(FM) <- FM[,1]
FM[,1] <- NULL
plot(colSums(FM[,-c(1:2)]))

Alpha_D <- as.data.frame(read_excel("data_raw/micro.xlsx", sheet = 1, col_names = TRUE)) 
rownames(Alpha_D) <- Alpha_D[,1]
Alpha_D[,1] <- NULL

Phylum <- as.data.frame(read_excel("data_raw/micro.xlsx", sheet = 4, col_names = TRUE))
rownames(Phylum) <- Phylum[,1]
Phylum[,1] <- NULL

Fun_Clus <- as.data.frame(read_excel("data_raw/micro.xlsx", sheet = 5, col_names = TRUE))
rownames(Fun_Clus) <- Fun_Clus[,1]
Fun_Clus[,1] <- NULL

idx <- as.data.frame(read_excel("data_raw/Env.xlsx", sheet = 1, col_names = TRUE))
rownames(idx) <- idx[,1]
idx[,1] <- NULL

ANF <- as.data.frame(read_excel("data_raw/Env.xlsx", sheet = 3, col_names = TRUE))
rownames(ANF) <- ANF[,1]
ANF[,1] <- NULL

save.image(file = "data_raw/Nif.Rdata")

plot(idx$B16SQ, idx$nifHQ)
plot(log(idx$B16SQ), log(idx$nifHQ))

SIP_ASV <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 1, col_names = TRUE))
rownames(SIP_ASV) <- SIP_ASV[,1]
SIP_ASV[,1] <- NULL

SIP_tax <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 2, col_names = TRUE))
rownames(SIP_tax) <- SIP_tax[,1]
SIP_tax[,1] <- NULL

SIP_DNACon <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 3, col_names = TRUE))
rownames(SIP_DNACon) <- SIP_DNACon[,1]
SIP_DNACon[,1] <- NULL

SIP_NifCon <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 4, col_names = TRUE))
rownames(SIP_NifCon) <- SIP_NifCon[,1]
SIP_NifCon[,1] <- NULL

SIP_Tr <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 6, col_names = TRUE))
rownames(SIP_Tr) <- SIP_Tr[,1]
SIP_Tr[,1] <- NULL

save.image(file = "data_raw/SIP.Rdata")
##### Env Permanova
load(file = "data_raw/Nif.Rdata")
idx$MAT <- idx$MAT+5
df <- t(decostand(t(idx[,3:21]),"total"))  # decostand should site in rows and Spe/Envs in columns
colSums(df)
dis <- vegdist(df, method = "bray", na.rm = T)
sum(rownames(df) == rownames(idx))
adon <- adonis2(dis ~ Gty, data=idx, permutations=999)
adon

anos <- with(idx, anosim(dis, Gty, distance = "bray", permutations = 999))
summary(anos)

#### Env factor NMDS plot
dis <- vegdist(df, method = "bray", na.rm = T)
sum(rownames(df) == rownames(idx))
mds <- metaMDS(as.matrix(dis),k=2,trymax=100)
mds$stress
mds_po <- as.data.frame(mds$points)
NMDS <- merge(mds_po, idx, by = "row.names")

## ## ## ## Env factor  visualization ## ## ## ## ## ## 
library(ggplot2)
library(ggsci)
library(ggConvexHull)
p <- ggplot(NMDS, aes(MDS1, MDS2)) +
  geom_point(aes(color = Gty, shape = Gty),size=3.5) +
  geom_convexhull(alpha=.1, aes(fill=Gty))+ 
  labs(x="NMDS1",
       y="NMDS2")+
  scale_color_aaas()+
  scale_fill_aaas()+
  theme_classic()+
  theme(axis.title = element_text(size = 20,colour = "black"))+
  theme(axis.text = element_text(colour = "black",size = 20,angle = 0,vjust = 0.5))+
  theme(legend.title =  element_text(size = 20,colour = "black"), 
        legend.text = element_text(size = 20,colour = "black"),
        strip.text =  element_text(size = 20,colour = "black"),
        text = element_text(family = "serif"))+
  annotate("text", x=-0.2, y=0.39, label= "Sress=0.12",size = 6)
p
## save results
Gname <- paste("Figs/Env_NMDS",".pdf",sep = "")
ggsave(Gname, p, width = 8, height = 5)


##### Micro Permanova
df <- t(decostand(t(FM[,-c(1:2)]),"total"))  # decostand should site in rows and Spe/Envs in columns
dis <- vegdist(t(df), method = "bray", na.rm = T)
sum(rownames(t(df)) == rownames(idx))
adon <- adonis2(dis ~ Gty, data=idx, permutations=999)
adon

anos <- with(idx, anosim(dis, Gty, distance = "bray", permutations = 999))
summary(anos)

#### microbial NMDS plot
library(vegan)
df <- t(decostand(t(FM[,-c(1:2)]),"total"))  # decostand should site in rows and Spe/Envs in columns
dis <- vegdist(t(df), method = "bray", na.rm = T)
sum(rownames(t(df)) == rownames(idx))
adonis(dis ~ Gty, data=idx, permutations=99)
mds <- metaMDS(as.matrix(dis),k=2,trymax=100)
mds$stress
mds_po <- as.data.frame(mds$points)
NMDS <- merge(mds_po, idx, by = "row.names")

## ## ## ## NMDS  visualization ## ## ## ## ## ## 
## http://127.0.0.1:27913/graphics/plot_zoom_png?width=2048&height=1096
library(ggplot2)
library(ggsci)
library(ggConvexHull)
p <- ggplot(NMDS, aes(MDS1, MDS2)) +
  geom_point(aes(color = Gty, shape = Gty),size=3.5) +
  geom_convexhull(alpha=.1, aes(fill=Gty))+ 
  labs(x="NMDS1",
     y="NMDS2")+
  scale_color_aaas()+
  scale_fill_aaas()+
  theme_classic()+
  theme(axis.title = element_text(size = 20,colour = "black"))+
  theme(axis.text = element_text(colour = "black",size = 20,angle = 0,vjust = 0.5))+
  theme(legend.title =  element_text(size = 20,colour = "black"), 
        legend.text = element_text(size = 20,colour = "black"),
        strip.text =  element_text(size = 20,colour = "black"),
        text = element_text(family = "serif"))+
  annotate("text", x=-0.35, y=0.39, label= "Sress=0.2",size = 6)
p
## save results
Gname <- paste("Figs/FLNF_NMDS",".pdf",sep = "")
ggsave(Gname, p, width = 8, height = 5)

## beta dispersal
library("ggpubr")
groups <- factor(idx$Gty)
mod <- betadisper(dis, groups)
summary(mod)
anova(mod)
TukeyHSD(mod)
boxplot(mod)
distan <- as.data.frame(mod$distances)
betadis <- merge(distan, idx, by = "row.names")
# p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution
shapiro.test(subset(betadis, Gty=="Alpine_meadow")$`mod$distances`)
ggdensity(subset(betadis, Gty=="Alpine_meadow")$`mod$distances`)

shapiro.test(log(subset(betadis, Gty=="Alpine_steppe")$`mod$distances`))
ggdensity(subset(betadis, Gty=="Alpine_steppe")$`mod$distances`)

# if data is normal distribution
t.test(betadis$`mod$distances` ~ betadis$Gty)
# else
wilcox.test(betadis$`mod$distances` ~ betadis$Gty)


#### alpha diversity of FLNF 
load(file = "data_raw/Nif.Rdata")
library(reshape2)
df <- merge(Alpha_D,idx[,c(3,24)], by="row.names")
# p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution
Re <- as.data.frame(matrix(nrow=4, ncol=5))
colnames(Re) <- c("alpha_D","p","Alpine_meadow_NM","Alpine_steppe_NM","test")
for (i in 2:5){
  st1 <- shapiro.test(subset(df, Gty=="Alpine_meadow")[,i])
  st2 <- shapiro.test(subset(df, Gty=="Alpine_steppe")[,i])
  if (st1$p.value>0.05&st2$p.value>0.05)  {
   Res <- t.test(df[,i] ~ df$Gty)  
   Re[i-1,5] <- "T-test"
  }else{
   Res <-  wilcox.test(df[,i] ~ df$Gty)
   Re[i-1,5] <- "wilcox-test"
  }
  Re[i-1,1:4] <- c(colnames(df)[i], Res$p.value, st1$p.value, st2$p.value)
}
write.table(Re, "clipboard")

divers <- merge(distan, Alpha_D, by = "row.names")
rownames(divers) <- divers[,1]
divers[,1] <- NULL
colnames(divers)[1] <- "beta_Disp"
divers_M <- merge(divers, idx[,c(3,24)], by = "row.names")
divers_M_F <- divers_M[,c(2,3,5,8)] 
df <- melt(divers_M_F, id="Gty")
##### Beta dispersal visualization
level_order <- c("observed_features", "faith_pd","beta_Disp")
p <- ggplot(df, aes(x=Gty, y=value)) + 
  geom_violin()+ theme_classic()+
  geom_boxplot(width=0.2, color="black", alpha=0.2)+
  stat_summary(fun = "mean", geom = "point", shape = 8, size = 3, color = "midnightblue")+
  facet_wrap(~ factor(variable, level = level_order), scales = "free", nrow = 1)+  # grouped by variable  free_x free_y
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0))+  # change the text size of facet
  #theme(strip.background = element_rect(fill="gray"))+ # color of facet 
  theme(axis.text.x=element_text(angle=60,hjust=0.5, vjust=0.5)) +
  theme(legend.position="none") +
  labs(x='', y= 'Distance to centroid') +
  theme(axis.title.x = element_text(size = rel(1))) +
  theme(axis.title.y = element_text(size = rel(1))) +
  theme(axis.text.x = element_text(colour = 'black', size = 10, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = 'black', size = 10, hjust = 0.5, vjust = 0.5))
p
## save results
Gname <- paste("Figs/Diversity",".pdf",sep = "")
ggsave(Gname, p, width = 8, height = 5)


### Activity of Diazotrophs
A_df <- ANF[,c(1,27)]
A_df <- A_df[A_df$ANF!=0, ]
st1 <- shapiro.test(subset(A_df, Gty=="AM")[,2])
st2 <- shapiro.test(subset(A_df, Gty=="AS")[,2])
wilcox.test(A_df[,2] ~ A_df$Gty)

p <- ggplot(A_df, aes(x=Gty, y=ANF)) + 
  geom_violin()+ theme_classic()+
  geom_boxplot(width=0.2, color="black", alpha=0.2)+
  stat_summary(fun = "mean", geom = "point", shape = 8, size = 3, color = "midnightblue")+
  #facet_wrap(~ factor(variable, level = level_order), scales = "free", nrow = 1)+  # grouped by variable  free_x free_y
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0))+  # change the text size of facet
  #theme(strip.background = element_rect(fill="gray"))+ # color of facet 
  theme(axis.text.x=element_text(angle=60,hjust=0.5, vjust=0.5)) +
  theme(legend.position="none") +
  labs(x='', y= 'N2 fixation rate kg N��ha-1��yr-1') +
  theme(axis.title.x = element_text(size = rel(1))) +
  theme(axis.title.y = element_text(size = rel(1))) +
  theme(axis.text.x = element_text(colour = 'black', size = 10, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = 'black', size = 10, hjust = 0.5, vjust = 0.5))
p
## save results
Gname <- paste("Figs/Activity_D",".pdf",sep = "")
ggsave(Gname, p, width =5, height = 3)




#### Quantification of FLNF 
library(reshape2)
library(ggplot2)
library(GGally)
load(file = "data_raw/Nif.Rdata")

df <- idx[,22:24]
df <- df[!is.na(df$nifHQ), ]
df$RA <- df$nifHQ/df$B16SQ*100

# p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution
Re <- as.data.frame(matrix(nrow=4, ncol=5))
colnames(Re) <- c("ID","p","Alpine_meadow_NM","Alpine_steppe_NM","test")
for (i in c(1:2,4)){
  st1 <- shapiro.test(subset(df, Gty=="Alpine_meadow")[,i])
  st2 <- shapiro.test(subset(df, Gty=="Alpine_steppe")[,i])
  if (st1$p.value>0.05&st2$p.value>0.05)  {
    Res <- t.test(df[,i] ~ df$Gty)  
    Re[i-1,5] <- "T-test"
  }else{
    Res <-  wilcox.test(df[,i] ~ df$Gty)
    Re[i,5] <- "wilcox-test"
  }
  Re[i,1:4] <- c(colnames(df)[i], Res$p.value, st1$p.value, st2$p.value)
}
write.table(Re, "clipboard")

##### Q-PCR visualization
level_order <- c("nifHQ", "B16SQ","RA")
library(reshape2)
library(ggplot2)
df[df$RA>0.1,]$RA <- NA
df2 <- melt(df, id = 'Gty')
df2[df2$variable!="RA",]$value <- log(df2[df2$variable!="RA",]$value)
p <- ggplot(df2, aes(x=Gty, y=value)) + 
  geom_violin()+ theme_classic()+
  geom_boxplot(width=0.2, color="black", alpha=0.2)+
  stat_summary(fun = "mean", geom = "point", shape = 8, size = 3, color = "midnightblue")+
  facet_wrap(~ factor(variable, level = level_order), scales = "free", nrow = 1)+  # grouped by variable  free_x free_y
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0))+  # change the text size of facet
  #theme(strip.background = element_rect(fill="gray"))+ # color of facet 
  theme(axis.text.x=element_text(angle=60,hjust=0.5, vjust=0.5)) +
  theme(legend.position="none") +
  labs(x='', y= 'Quantification') +
  theme(axis.title.x = element_text(size = rel(1))) +
  theme(axis.title.y = element_text(size = rel(1))) +
  theme(axis.text.x = element_text(colour = 'black', size = 10, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = 'black', size = 10, hjust = 0.5, vjust = 0.5))
p
## save results
Gname <- paste("Figs/Q-PCR",".pdf",sep = "")
ggsave(Gname, p, width = 10, height = 5)


##### Phylum and cluster composition #######
library(vegan)
load(file = "data_raw/Nif.Rdata")
df <- decostand(t(Phylum),"total") # decostand should site in rows and Spe/Envs in columns
rowSums(df)
df2 <- merge(df, idx[,23:24], by = "row.names")
## difference test between AM and AS
# p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution
Re <- as.data.frame(matrix(nrow=8, ncol=9))
colnames(Re) <- c("Phylum","p","Alpine_meadow_Mean","Alpine_meadow_SE","Alpine_steppe_Mean","Alpine_steppe_SE","Alpine_meadow_NM","Alpine_steppe_NM","test")
for (i in 2:9){
  st1 <- shapiro.test(subset(df2, Gty=="Alpine_meadow")[,i])
  st2 <- shapiro.test(subset(df2, Gty=="Alpine_steppe")[,i])
  if (st1$p.value>0.05&st2$p.value>0.05)  {
    Res <- t.test(df2[,i] ~ df2$Gty)  
    Re[i-1,9] <- "T-test"
  }else{
    Res <-  wilcox.test(df2[,i] ~ df2$Gty)
    Re[i-1,9] <- "wilcox-test"
  }
  m1 <- mean(subset(df2, Gty=="Alpine_meadow")[,i])
  SE1 <- sd(subset(df2, Gty=="Alpine_meadow")[,i])/sqrt(length((subset(df2, Gty=="Alpine_meadow")[,i])))
  m2 <- mean(subset(df2, Gty=="Alpine_steppe")[,i])
  SE2 <- sd(subset(df2, Gty=="Alpine_steppe")[,i])/sqrt(length((subset(df2, Gty=="Alpine_steppe")[,i])))
  Re[i-1,1:8] <- c(colnames(df2)[i], Res$p.value, m1, SE1, m2, SE2, st1$p.value, st2$p.value)
}
write.table(Re, "clipboard")

df3<- melt(df2[,-c(1,11)], id='Gty')
df4 <- aggregate(value ~ Gty+variable, data = df3, mean)
df4$value <- df4$value*100
ord <- c("Others",	"Verrucomicrobia",	"Firmicutes",	"Cyanobacteria","Other_Proteobacteria","Gammaproteobacteria",	"Betaproteobacteria",	"Deltaproteobacteria",	"Alphaproteobacteria")
p <- ggplot(df4, aes(x=Gty, y=value, fill = factor(variable,level = ord))) + 
  geom_col(aes(color=factor(variable,level = ord)))+
  theme(strip.text.x = element_text(size = 32, colour = "black", angle = 0))+  # change the text size of facet
  #theme(strip.background = element_rect(fill="gray"))+ # color of facet 
  theme(axis.text.x=element_text(angle=60,hjust=0.5, vjust=0.5)) +
  labs(x='', y= 'Relative abundance (%)') +
  theme(axis.title.x = element_text(size = rel(2))) +
  theme(axis.title.y = element_text(size = rel(2))) +
  theme(axis.text.x = element_text(colour = 'black', size = 20, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = 'black', size = 20, hjust = 0.5, vjust = 0.5))+
  scale_color_aaas()+
  scale_fill_aaas()
p
## save results
Gname <- paste("Figs/Phylum_Composition",".pdf",sep = "")
ggsave(Gname, p, width = 5, height = 6)

#### #### ####  functional cluster composition ######### #### 
Fun_Clus <- t(Fun_Clus)
Fun_Clus <- decostand(Fun_Clus,"total") # decostand should site in rows and Spe/Envs in columns
rowSums(Fun_Clus)

df2 <- merge(Fun_Clus, idx[,23:24], by = "row.names")

Re <- as.data.frame(matrix(nrow=10, ncol=9))
colnames(Re) <- c("Phylum","p","Alpine_meadow_Mean","Alpine_meadow_SE","Alpine_steppe_Mean","Alpine_steppe_SE","Alpine_meadow_NM","Alpine_steppe_NM","test")
for (i in 2:10){
  st1 <- shapiro.test(subset(df2, Gty=="Alpine_meadow")[,i])
  st2 <- shapiro.test(subset(df2, Gty=="Alpine_steppe")[,i])
  if (st1$p.value>0.05&st2$p.value>0.05)  {
    Res <- t.test(df2[,i] ~ df2$Gty)  
    Re[i-1,9] <- "T-test"
  }else{
    Res <-  wilcox.test(df2[,i] ~ df2$Gty)
    Re[i-1,9] <- "wilcox-test"
  }
  m1 <- mean(subset(df2, Gty=="Alpine_meadow")[,i])
  SE1 <- sd(subset(df2, Gty=="Alpine_meadow")[,i])/sqrt(length((subset(df2, Gty=="Alpine_meadow")[,i])))
  m2 <- mean(subset(df2, Gty=="Alpine_steppe")[,i])
  SE2 <- sd(subset(df2, Gty=="Alpine_steppe")[,i])/sqrt(length((subset(df2, Gty=="Alpine_steppe")[,i])))
  Re[i-1,1:8] <- c(colnames(df2)[i], Res$p.value, m1, SE1, m2, SE2, st1$p.value, st2$p.value)
}
write.table(Re, "clipboard")

df3<- melt(df2[,-c(1,12)], id='Gty')
df4 <- aggregate(value ~ Gty+variable, data = df3, mean)
df4$value <- df4$value*100
ord <- c("Others",	"1P",	"1E",	"1D",	"1",	"3E",	"1B",	"1A",	"1K",	"1J")
p <- ggplot(df4, aes(x=Gty, y=value, fill = factor(variable,level = ord))) + 
  geom_col(aes(color=factor(variable,level = ord)))+
  theme(strip.text.x = element_text(size = 32, colour = "black", angle = 0))+  # change the text size of facet
  #theme(strip.background = element_rect(fill="gray"))+ # color of facet 
  theme(axis.text.x=element_text(angle=60,hjust=0.5, vjust=0.5)) +
  labs(x='', y= 'Relative abundance (%)') +
  theme(axis.title.x = element_text(size = rel(2))) +
  theme(axis.title.y = element_text(size = rel(2))) +
  theme(axis.text.x = element_text(colour = 'black', size = 20, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = 'black', size = 20, hjust = 0.5, vjust = 0.5))+
  scale_color_aaas()+
  scale_fill_aaas()
p
## save results
Gname <- paste("Figs/Functi_Cluster_Composition",".pdf",sep = "")
ggsave(Gname, p, width = 5, height = 6)

#########   Distance decay #############
#########Function to calculate geospatial distance between two points (lat,long) ############
library(vegan)
library(readxl)
load(file = "data_raw/Nif.Rdata")
df <- t(decostand(t(FM[,-c(1:2)]),"total"))  # decostand should site in rows and Spe/Envs in columns
dis <- vegdist(t(df), method = "bray", na.rm = T)
dist <- as.matrix(dis)

dist_Med <- dist[rownames(dist)%in%rownames(idx[idx$Gty=="Alpine_meadow",]), colnames(dist)%in%rownames(idx[idx$Gty=="Alpine_meadow",])]
dist_Stp <- dist[rownames(dist)%in%rownames(idx[idx$Gty=="Alpine_steppe",]), colnames(dist)%in%rownames(idx[idx$Gty=="Alpine_steppe",])]

xy <- t(combn(colnames(dist_Med), 2))
re <- data.frame(xy, dist=dist_Med[xy])
write.table(re, 'tem/distance_Med.txt')

xy <- t(combn(colnames(dist_Stp), 2))
re <- data.frame(xy, dist=dist_Stp[xy])
write.table(re, 'tem/distance_Stp.txt')

library(geosphere)
library(reshape2)
library(readxl)
GD <- read_xlsx("tem/tem.xlsx", sheet = 1, col_names = TRUE)
idx$MID <- rownames(idx)
colnames(GD)[1] <- "MID"
GD1 <- merge(GD,idx[,c(1:2,25)], by = "MID") 
colnames(GD1)[1:2] <- c("ID1", "MID")
GD2 <- merge(GD1,idx[,c(1:2,25)], by = "MID") 

for (i in 1:7359){
  re <- distm (GD2[i,5:6], GD2[i,7:8], fun = distHaversine)
  GD2[i, 9] <- re
}

cor.test(x = subset(GD2, Gty=="Alpine_meadow")$V9, y = subset(GD2, Gty=="Alpine_meadow")$dist)
cor.test(x = subset(GD2, Gty=="Alpine_steppe")$V9, y = subset(GD2, Gty=="Alpine_steppe")$dist)


library(ggplot2)
library(ggsci)
GD2$V9 <- GD2$V9/1000 # trans m to Km
df <- subset(GD2, V9>0)
df$dist <- 1 - df$dist
DD <- ggplot(df, aes(x=V9, y=dist)) + 
  geom_point(alpha=1, size=2, aes(colour=Gty))+
  #facet_wrap(~ GT+Depth, scale = 'free_y', nrow = 2)+ ##free in one dimension ("free_x", "free_y")
  geom_smooth(method = "lm", formula = y ~ poly(x,1, raw = T), aes(group=Gty, colour=Gty),)+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0))+  # change the text size of facet
  theme(strip.background = element_rect(fill="gray"))+ # color of facet 
  theme(axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5)) +
  labs(y = expression('Community Similarity'), x = expression('Geographic Distance (Km)'))+
  #scale_x_continuous(breaks=seq(0, 3000, 500))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 0.9, size=10,color="black"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=10,color="black"))+
  theme(axis.title.x = element_text(size = rel(1)))+
  theme(axis.title.y = element_text(size = rel(1)))+
  scale_color_aaas()+
  scale_fill_aaas()
DD
#geom_hline(aes(yintercept=0.5),colour="#990000", linetype="dashed")  #add horizon line
#geom_vline(aes(xintercept=0),colour="#990000", linetype="dashed")
tablename <- paste("Figs/Distandecay",".pdf",sep = "")
ggsave(tablename, DD, width = 8, height = 6)


######## Environmental effects #######
load(file = "data_raw/Nif.Rdata")
############################## alpha-diversity #####################
############################  contributions of environmental factors to alpha diversity change ##########
# Ref: https://mp.weixin.qq.com/s/80QDYIFejPBW4R2BPaiSrg
# this codes used GLMM model, and also can use random forest to assess env factor contributions
library(readxl)
library(clusterSim)
library(vegan)
library(GGally)#相关分析�?
library(car)

Ad_O <- merge(Alpha_D[,1:3], idx, by = 'row.names')
rownames(Ad_O) <- Ad_O[,1]
Ad_O[,1] <- NULL
Ad <- Ad_O[,-c(25:27)]

### colinear test for observed_features
a <- lm(observed_features~., data=log(Ad[,-c(2:3)]+1))
shapiro.test(a$residuals) ## p > 0.05 when data id normally distributed
summary(a)
VF <- vif(a)  #remove fators that vif>10
if (sum(VF > 10)>=1){
  namord <- order(VF, decreasing = T)[1]+3
  tem <- Ad[,-c(2:3, namord)]
  }

### repeat this part until sum(VF > 10)>1 is FAULSE##
a <- lm(observed_features~., data=log(tem))
VF <- vif(a)  #remove fators that vif>20
if (sum(VF > 10)>=1){
  namord <- order(VF, decreasing = T)[1]+1
  tem <- tem[,-c(namord)]
}
VF
sum(VF > 10)>1
namord
####################################################

Env_observed_features <- names(VF)


### colinear test for Faith's PD
a <- lm(faith_pd~., data=log(Ad[,-c(1,3)]))
shapiro.test(a$residuals) ## p > 0.05 when data id normally distributed
summary(a)
VF <- vif(a)  #remove fators that vif>10
if (sum(VF > 10)>=1){
  namord <- order(VF, decreasing = T)[1]+3
  tem <- Ad[,-c(1,3,4,namord)]
}
### repeat this part until sum(VF > 10)>1 is FAULSE##
a <- lm(faith_pd~., data=log(tem))
VF <- vif(a)  #remove fators that vif>20
if (sum(VF > 10)>=1){
  namord <- order(VF, decreasing = T)[1]+1
  tem <- tem[,-c(namord)]
}
VF
namord
sum(VF > 10)>1
#####################################################
Env_faith_pd <- names(VF)

ggpairs(tem[,c(3:17)]) #check self cor factors r >0.7
library(lme4)#新版混合模型�?
library(lmerTest)#输出自变量p�?
library(effectsize)#输出标准化回归系数�?
tem2 <- merge(log(Ad+1), Ad_O[,26:27], by = 'row.names')
##   observed_features
Env_observed_features
fit1<-lmer(observed_features ~ MAP+MAT+NPP+OC+TN+TP+SM+DOC+NH4+NO3+pH+EC+Fe+Mo+(1|Gty), tem2) #新版混合模型
summary(fit1) #significant test for each x
effectsize(fit1)  ### effect size of each x

library(MuMIn)# ouput the overall R2 for the model
library(glmm.hp)# analyzethe contribution of each x
r.squaredGLMM(fit1) # overall R2 for the model
a1<-glmm.hp(fit1) # R2 for each x
a1
plot(a1)

## faith_pd
fit2<-lmer(faith_pd ~ MAP+MAT+NPP+OC+TN+TP+SM+DOC+NH4+NO3+pH+EC+Fe+Mo+(1|Gty), tem2) #新版混合模型 正态分布不用log转换
summary(fit2) #significant test for each x
effectsize(fit2)  ### effect size of each x

library(MuMIn)# ouput the overall R2 for the model
library(glmm.hp)# analyzethe contribution of each x
r.squaredGLMM(fit2) # overall R2 for the model
a2<-glmm.hp(fit2) # R2 for each x
a2
plot(a2)

##### environmental factors cor with activity
### colinear test for ANF
library(vegan)
a <- lm(ANF~., data=ANF[,-c(1)])
shapiro.test(a$residuals) ## p > 0.05 when data id normally distributed
summary(a)
VF <- vif(a)  #remove fators that vif>10
VF
if (sum(VF > 10)>1){
  namord <- order(VF, decreasing = T)[1]+1
  tem <- ANF[,-c(1, namord)]
}
### repeat this part until sum(VF > 10)>1 is FAULSE##
a <- lm(ANF~., data=tem)
VF <- vif(a)  #remove fators that vif>20
if (sum(VF > 10)>=1){
  namord <- order(VF, decreasing = T)[1]
  tem <- tem[,-c(namord)]
}
sum(VF > 10)>1
namord
VF
#####################################################
Env_ANF <- names(VF)

ggpairs(tem) #check self cor factors r >0.7
library(lme4)#新版混合模型�?
library(lmerTest)#输出自变量p�?
library(effectsize)#输出标准化回归系数�?
tem2 <- merge(tem, ANF[,1:2], by = 'row.names')
##   ANF
fit1<-lmer(ANF ~ Ele+MAP+MAT+SOC+TN+ TP+ SM+ DOC + NH4+NO3+pH+ EC+ Fe+ Mo+ observed_features+ nifHQ+(1|Gty), tem2) #新版混合模型
summary(fit1) #significant test for each x
effectsize(fit1)  ### effect size of each x

library(MuMIn)# ouput the overall R2 for the model
library(glmm.hp)# analyzethe contribution of each x
r.squaredGLMM(fit1) # overall R2 for the model
a<-glmm.hp(fit1) # R2 for each x
a
plot(a)


## NMDS_1
fit3<-lmer(MDS1 ~ MAP+MAT+NPP+OC+TN+TP+SM+DOC+NH4+NO3+pH+EC+Fe+Mo+(1|Gty), tem2) #新版混合模型 正态分布不用log转换
summary(fit3) #significant test for each x
effectsize(fit3)  ### effect size of each x

library(MuMIn)# ouput the overall R2 for the model
library(glmm.hp)# analyzethe contribution of each x
r.squaredGLMM(fit3) # overall R2 for the model
a3<-glmm.hp(fit3) # R2 for each x
a3
plot(a3)

## Nif Quantification
tem3 <- tem2[!is.na(tem2$nifHQ), ]
fit4<-lmer(nifHQ ~ MAP+MAT+NPP+OC+TN+TP+SM+DOC+NH4+NO3+pH+EC+Fe+Mo+(1|Gty), tem3) #新版混合模型 正态分布不用log转换
summary(fit4) #significant test for each x
effectsize(fit4)  ### effect size of each x

library(MuMIn)# ouput the overall R2 for the model
library(glmm.hp)# analyzethe contribution of each x
r.squaredGLMM(fit4) # overall R2 for the model
a4<-glmm.hp(fit4) # R2 for each x
a4
plot(a4)

############ plot contribution of env factors
## plot
library(ggplot2)
library(ggsci)
library(readxl)
library(reshape2)
CorLa <- as.data.frame(read_excel("tem/tem.xlsx", sheet = 2, col_names = TRUE)) # one column, named 'Group'
CorLa <- melt(CorLa, id = "ID")

RF <- as.data.frame(read_excel("tem/tem.xlsx", sheet = 3, col_names = TRUE)) # one column, named 'Group'
RF <- melt(RF, id = "ID")
  
level_od_y <- c('MAP','MAT','NPP','OC','TN','TP','SM','DOC','NH4','NO3','pH','EC','Fe','Mo')
level_od_x <- c("NifQ","SR","PD",	"MDS","ANF"	)
heatmap<- ggplot()+
  geom_tile(data = CorLa, aes(x = factor(variable,levels = level_od_x), y = factor(ID,levels = level_od_y), fill = value))+ #热图
  scale_fill_gradientn(colors = c('#2D6DB1', 'white', '#DC1623'),
                       limit = c(-2, 2))+
  geom_point(data = RF, aes(x = factor(variable,levels = level_od_x), y = factor(ID,levels = level_od_y), 
                            size = value), shape = 21)+ #圆圈大小
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, color = "black",
                                   size = 12, vjust = 0.6), 
        axis.text.y = element_text(color = 'black', size = 12),
        legend.title = element_text(size = 10))+
  labs(y = '', x = '')
heatmap
GD <- "Graphs"
tablename <- paste("Figs/Env_Contri",".pdf",sep = "")
ggsave(tablename, heatmap, width = 4, height = 5)

## overall contributions
RF_im <-  as.data.frame(read_excel("tem/tem.xlsx", sheet = 4, col_names = TRUE)) # one column, named 'Group'
RF_im$Con <-  RF_im$Con*100 
Exp_var <- ggplot(RF_im, aes(factor(M_D,levels = level_od_x), Con))+
  geom_bar(stat = "identity", fill = "steelblue")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5)) +
  labs(y = expression('Explained variation (%)'), x = expression(''))+
  scale_color_aaas()+
  scale_fill_aaas()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 0.9, size=10,color="black"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=10,color="black"))+
  theme(axis.title.x = element_text(size = rel(1)))+
  theme(axis.title.y = element_text(size = rel(1)))
Exp_var
tablename <- paste("Figs/explained_vari",".pdf",sep = "")
ggsave(tablename, Exp_var, width = 4, height = 2)


reCor <- Cortem
for (i in 1:4){
  for (e in 5:27){
    CoR <- cor.test(x = Ad[,i], y = Ad[,e])
    Cortem <- data.frame(CoR[c("estimate", "p.value")])
    Cortem$pro <- colnames(Ad)[i]
    Cortem$pro2 <- colnames(Ad)[e]
    reCor <- rbind(reCor, Cortem)
  }
}
write.csv(reCor[-1,], "tmp/alp.csv")

############## SIP data analysis ##################
library(reshape2)
library(vegan)
library(dplyr)
load(file = "data_raw/SIP.Rdata")

SIP_ASV <- SIP_ASV+1
df <- as.data.frame(t(decostand(t(SIP_ASV),"total")))  # decostand should site in rows and Spe/Envs in columns
df$SUM <- rowSums(df)
df2 <- arrange(df, desc(SUM))
df3 <- subset(df2, SUM>0.001)
colSums(df3)

dim(df3)

for (i in 1:537){
  tem1 <- log(df3[i,3]/df3[i,2])
  tem2 <- log(df3[i,5]/df3[i,4])
  tem3 <- log(df3[i,6]/df3[i,1])
  df3[i,8:10] <- c(tem1, tem2, tem3)
}
colnames(df3)[8:10] <- c("NQ", "AR", "GC")

Re <- merge(df3, SIP_tax, by = "row.names")

Re2 <- subset(Re, (NQ+AR+GC)>1.2)
write.csv(Re, "tem/SIP_RR1.csv")
write.csv(Re2, "tem/SIP_RR2.csv")

################################################# Enrichned Genus ##########
library(readxl)
SIP_Enr <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 7, col_names = TRUE))
rownames(SIP_Enr) <- SIP_Enr[,1]
SIP_Enr[,1] <- NULL
RE2 <- aggregate(.~ Genus, data=SIP_Enr[,c(16,1:6)], sum)


for (i in 1:23){
  tem1 <- log(RE2[i,4]/RE2[i,3])
  tem2 <- log(RE2[i,6]/RE2[i,5])
  tem3 <- log(RE2[i,7]/RE2[i,2])
  margin <- qt(0.025,df=2)*sd(c(tem1,tem2,tem3))/sqrt(3)
  RE2[i,8:11] <- c(tem1, tem2, tem3,margin)
}
colnames(RE2)[8:11] <- c("NQ", "AR", "GC","margin")

write.csv(RE2, "tem/SIP_RR3.csv")


##############################################################
library(ggplot2)
library(ggsci)
library(ggthemes)
SIP_Genus_Enr <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 8, col_names = TRUE))
SIP_Genus_Enr <- SIP_Genus_Enr[1:15,]
lo <- c("Mesorhizobium",	"Bradyrhizobium",	"Rhodovulum",	"Desulfovibrio",	"Skermanella",	"Geobacter",	"Rhodopseudomonas",	"Insolitispirillum",	"Methylobacterium",	"Roseospirillum",	"Novosphingobium",	"Pelosinus",	"Rubrivivax",	"Paenibacillus",	"Azohydromonas")
p <- ggplot(SIP_Genus_Enr, aes(Mean, factor(Genus, levels = lo))) +
  geom_point(aes(size=RA))+
  geom_errorbarh(aes(xmax = Mean +margin, xmin = Mean -margin),size= 0.5,height = 0.2) +
  scale_x_continuous(limits= c(-2, 6))+
  geom_vline(aes(xintercept = 0),color="gray",linetype="dashed", size = 0.6) +
  xlab('LnRR')+ 
  ylab(' ')+
  theme_few()+
  theme(axis.text.x = element_text(size = 24, color = "black"))+
  theme(axis.text.y = element_text(size = 24, color = "black"))+
  theme(title=element_text(size=24))+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))
tablename <- paste("Figs/SIP_Genus",".pdf",sep = "")
ggsave(tablename, p, width = 4, height = 3)


##################### SIP Genus composition #################
SIP_Genus_Enr <- as.data.frame(read_excel("data_raw/SIP_RE.xlsx", sheet = 9, col_names = TRUE))
ord <- c("Mesorhizobium",	"Bradyrhizobium",	"Rhodovulum",	"Desulfovibrio",	"Skermanella",	"Geobacter",	"Others")
p <- ggplot(SIP_Genus_Enr, aes(x=Tr, y=value, fill = factor(Genus,level = ord))) + 
  geom_col(aes(color=factor(Genus,level = ord)))+
  theme(strip.text.x = element_text(size = 32, colour = "black", angle = 0))+  # change the text size of facet
  #theme(strip.background = element_rect(fill="gray"))+ # color of facet 
  theme(axis.text.x=element_text(angle=60,hjust=0.5, vjust=0.5)) +
  labs(x='', y= 'Relative abundance (%)') +
  theme(axis.title.x = element_text(size = rel(2))) +
  theme(axis.title.y = element_text(size = rel(2))) +
  theme(axis.text.x = element_text(colour = 'black', size = 20, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(colour = 'black', size = 20, hjust = 0.5, vjust = 0.5))+
  scale_color_aaas()+
  scale_fill_aaas()
p
## save results
Gname <- paste("Figs/SIN_Genuscompsosition",".pdf",sep = "")
ggsave(Gname, p, width = 5, height = 6)

############## Phylogenetic tree for nifH genes ###########
library(seqinr)
ref <- read.fasta("data_raw/sequences.fasta")
list <- read.table("data_raw/seqid.txt", header = TRUE)
filter <- ref[names(ref) %in% list$ASV]
write.fasta(filter,names = names(filter), file.out ="Enirichednif.fasta")

#### sequence alignment ###########
library(msa)
Seqs <- readAAStringSet("data_raw/Enirichednif.fasta")
myClustalWAlignment <- msa(Seqs, "ClustalW")
myClustalWAlignment
hemoAln2 <- msaConvert(myClustalWAlignment, type="seqinr::alignment")
d <- dist.alignment(hemoAln2, "identity")
Tree <- nj(d)
plot(Tree, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")

write.tree(Tree, "data_raw/nifE.nwk")

######################## NifH Phylogenetic tree #####################
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(tidyverse)
library(readxl)
library(ggsci)

Tax <- read_excel("data_raw/SIP_Enr.xlsx", sheet = 1)
rownames(Tax) <- Tax$ASV
Anno <- read_excel("data_raw/SIP_Enr.xlsx", sheet = 2)
rownames(Anno) <- Anno$ASV
Tre <- read.tree("data_raw/nifE.nwk")



p <- ggtree(Tre, layout="circular", open.angle=10,branch.length="none")%<+%Tax + 
     geom_tippoint(aes(color=Phylum), size=1.5,position = "identity",
                show.legend=T)


p <- rotate_tree(p, -90)+  scale_color_npg()+
  scale_fill_npg()

p <- p +
  geom_fruit(data = Anno,
    geom=geom_col,
    mapping = aes(y=ASV, x=LnRR,
      group=Phylum2,
      fill=Phylum2,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 

p <- p+  scale_color_lancet()+
  scale_fill_lancet()

## save results
Gname <- paste("Figs/Enriched_Phylogenetic",".pdf",sep = "")
ggsave(Gname, p, width = 10, height = 10)

