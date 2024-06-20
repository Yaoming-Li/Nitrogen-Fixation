df<-read.csv(file.choose())

library(randomForest)
library(tidyverse)
library(rfUtilities)
library(rfPermute)
require(caret)

###raining set (70%) and the test set (30%)
set.seed(123)
train <- sample(nrow(df), nrow(df)*0.7)
mydata <- df[train, ]
df_test <- df[-train, ]
###construct Random forest model 
set.seed(123)
mydata.forest <- randomForest(ANF~., data =mydata ,importance = TRUE)
mydata.forest

plot(mydata.forest,main="Bagging OOB Errors")
###Importance of variable
importance(mydata.forest)
varImpPlot(mydata.forest,main="Variable Importance Plot")

set.seed(123)
treat_perm <- rf.significance(mydata.forest, mydata, nperm=99, ntree=500)
treat_perm

set.seed(123)
ANPP_rfP<- rfPermute(ANF~ ., data = mydata, ntree = 500,
                     na.action = na.omit, nrep = 100,num.cores = 1)
ANPP_dat <- importance(ANPP_rfP, sort.by = NULL, decreasing = TRUE)
ANPP_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(IncNodePurity.pval < 0.001,"***",
                         if_else(IncNodePurity.pval <0.01,"**",
                                 if_else(IncNodePurity.pval<0.05,"*","ns"))),
         IncNodePurity = as.numeric(IncNodePurity)) %>%
  arrange(IncNodePurity) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>%
  ggplot(aes(x = names, y = IncNodePurity))+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = IncNodePurity + 1,label = label))+
  labs(x = "", y = "IncNodePurity")+
  coord_flip()+
  theme_bw()

partialPlot(mydata.forest,mydata,x.var=MAT)
partialPlot(mydata.forest,mydata,x.var=MAP)
partialPlot(mydata.forest,mydata,x.var=TN)
partialPlot(mydata.forest,mydata,x.var=pH)

###Evaluate predictive performance using test sets
ANFR_predict <- predict(mydata.forest, df_test)
plot(df_test$ANF, ANFR_predict, main = 'Test',
     xlab = 'Actual', ylab = 'Predict')
abline(0, 1)

set.seed(123)
foldid <- sample(1:5,size=76,replace=TRUE) 

RMSE <- matrix(rep(0,60),ncol=5)   # initialize CV error for mtry 
MAE <- matrix(rep(0,60),ncol=5)
for (i in 1:12){
  for (j in 1:5){
    train_cv <- mydata[foldid!=j,]
    holdout <- mydata[foldid==j,]
    fit <- randomForest(ANF~.,data=train_cv,mtry=i)
    pred <- predict(fit,newdata=holdout)
    y.test <- holdout[,"ANF"]
    RMSE[i,j] <-sqrt(mean((pred-y.test)^2)) 
    MAE[i,j]<- mean(abs(pred-y.test))
  }
}
cv.error.RMSE <- apply(RMSE,1,mean)
cv.error.MAE <- apply(MAE,1,mean)
min(cv.error)
which.min(cv.error)

##R2
# define training control
control <- trainControl(method="repeatedcv", number=5, repeats=5, search="random")
set.seed(499)
# train the model
model <-  train(ANF~., data=mydata, method="rf", ntree=500, tuneGrid = data.frame(mtry = c(1,2,3,4,5,6,7,8,9,10,11,12)), trControl=control)
# summarize results
print(model)

###five repeats 5-fold cross validation
set.seed(123)
mydata.cv <- replicate(5, rfcv(mydata[-ncol(mydata)], mydata$ANF, cv.fold = 5, step = 1.1), simplify = FALSE)
mydata.cv

mydata.cv <- data.frame(sapply(mydata.cv, '[[', 'error.cv'))
mydata.cv$otus <- rownames(mydata.cv)
mydata.cv <- reshape2::melt(mydata.cv, id = 'otus')
mydata.cv$otus <- as.numeric(as.character(mydata.cv$otus))

mydata.cv.mean <- aggregate(mydata.cv$value, by = list(mydata.cv$otus), FUN = mean)
head(mydata.cv.mean, 10)

library(ggplot2)

ggplot(mydata.cv.mean, aes(Group.1, x)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 12, 2)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +  
  labs(title = '',x = 'Number of Variables', y = 'Cross-validation error')

###Retain variables as prompted
df.select<- df[c('MAT','MAP','TN','pH','ANF')]
set.seed(123)
train <- sample(nrow(df.select), nrow(df.select)*0.7)
mydata.select <- df.select[train, ]
df_test.select <- df.select[-train, ]

set.seed(123)
mydata.select.forest <- randomForest(ANF~., mtry = 2,data = mydata.select,importance = TRUE)
mydata.select.forest

pred <- predict(mydata.select.forest,newdata=df_test.select)
y.test <- df_test.select$"ANF"
RMSE <- sqrt(mean((pred-y.test)^2))
MAE<- mean(abs(pred-y.test))
#Model construction complete

#######Spatial prediction
library(raster) 
library(rasterVis) 
library(ggspatial)
library(ggplot2)   
library(sf)

Raster_Path <- dir("D:/download data/mod/fourth time/four predictors",full.names = TRUE)

MAP <- raster(Raster_Path[1])
MAT <- raster(Raster_Path[2])
PH <- raster(Raster_Path[3])
TN<- raster(Raster_Path[4])

###Harmonization of coordinates and resampling
MAP.gcs<- projectRaster(MAP, crs="+proj=longlat +datum=WGS84")
st_crs(MAP.gcs)
MAT.gcs<- projectRaster(MAT, crs="+proj=longlat +datum=WGS84")
PH.gcs<- projectRaster(PH, crs="+proj=longlat +datum=WGS84")
TN.gcs<- projectRaster(TN, crs="+proj=longlat +datum=WGS84")

MAP_r <- resample( MAP.gcs,PH.gcs,
                   method = "ngb")
MAT_r <- resample( MAT.gcs,PH.gcs,
                   method = "ngb")

tn=as.data.frame(TN.gcs,row.names = NULL,optional=T ,xy=T )
tn_data=na.omit(tn)

map=as.data.frame(MAP_r,row.names = NULL,optional=T ,xy=T )
map_data=na.omit(map)

mat=as.data.frame(MAT_r,row.names = NULL,optional=T ,xy=T )
mat_data=na.omit(mat)

ph=as.data.frame(PH.gcs,row.names = NULL,optional=T ,xy=T )
ph_data=na.omit(ph)

TWO=merge(ph_data,tn_data,by = c("x","y"))
THREE=merge(TWO,map_data,by = c("x","y"))
FOUR=merge(THREE,mat_data,by = c("x","y"))

library(tibble)
library(dplyr)
##prediction
ANFR_pred <- predict(mydata.select.forest, FOUR)
FOUR$ANFR <- ANFR_pred
###output result
write.table (FOUR, file ="D:/download data/mod/fourth time/ANFR.txt", sep =" ", row.names =F, col.names =TRUE, quote =TRUE)



