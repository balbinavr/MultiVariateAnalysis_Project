library("DMwR")
library("chemometrics")
library("VIM")
library("robustbase")
library("mice")
library("mvoutlier")
library("fBasics")
library("FactoMineR")
library("factoextra")
library("readr")
library("rgl")
library("pls")
library("rpart.plot")
library("randomForest")
library("stats")
library("ROCR")
source("mahalanobisDistance.R")
library("VGAM")
library(caret)

#setwd("~/Desktop/MVA/Practiques")

set.seed(7)

getAccuracy <- function(predictions, labels, threshold) {
     # ROC CURVE
     pred <- prediction(predictions, labels)
     perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
     plot(perf, col=rainbow(10))
     
     #ESTABLISH 0.6 as threshold
     predictions <- sapply(predictions, function (value) {
         if (value > threshold) return(TRUE)
         else return(FALSE)
     })
     
     count = 0
     i = 1
     for (label in labels) {
         if (label == predictions[i]) {
             count =  count + 1;
         }
         i = i + 1
     }
     
     acc <- 100.0*(count/length(labels))
     
     
     return(acc)
 }

# PREPROCESSING
# Remove barcode, ID of individual, Patient
# Interpret variables with levels / binary variables 
NKI_initial <- as.data.frame(
    read_csv("NKI_cleaned.csv", 
                        col_types = cols(amputation = col_logical(), 
                                         eventdeath = col_logical(), 
                                         angioinv = col_factor(levels = c("1","2", "3")), 
                                         chemo = col_logical(), 
                                         grade = col_factor(levels = c("1","2", "3")), 
                                         histtype = col_factor(levels = c("1", "2", "4", "5", "7")), 
                                         hormonal = col_logical(), 
                                         lymphinfil = col_factor(levels = c("1","2", "3")))))

# MISSING VALUES

sum(is.na(NKI_initial))

# FEATURE SELECTION
rowsToRemove = c("barcode", "ID", "Patient")
NKI <- NKI_initial[ , -which(names(NKI_initial) %in% rowsToRemove)]

corrplot(cor(NKI[,c(1, 3:4, 9)]))
columnsToDrop <- findCorrelation(cor(NKI[,c(1, 3:4, 9)]), cutoff = 0.9, names=TRUE)
NKI <- NKI[ , -which(names(NKI) %in% columnsToDrop)]

first_gene <- grep("esr1", colnames(NKI))

# OUTLIER DETECTION

# LOF
foundOutliersLOF.scores <- lofactor(NKI[, first_gene:nrow(NKI)], k=5)
plot(foundOutliersLOF.scores)
foundOutliersLOF = order(foundOutliersLOF.scores, decreasing=T)
rownames(NKI)[foundOutliersLOF]

# CANT BE DONE n < p
multivariateOutliers(NKI, 0.975, 0.1)

# FIRST APPROACH -> PRINCIPAL COMPONENTS
pcr_model <- pcr(eventdeath~., data = NKI[,c(2, first_gene:ncol(NKI))], scale=TRUE, validation = "CV") # validations = "LOO"
validationplot(pcr_model, val.type = "R2")
validationplot(pcr_model, val.type = "MSEP")

NKI_comp <- cbind(NKI[, 1:(first_gene-1)], pcr_model$scores[, 1:21])

# TRAIN / TEST SPLIT
test_threshold = 0.2
test_df = NKI_comp[1:(test_threshold*nrow(NKI_comp)),]
train_df = NKI_comp[(test_threshold*nrow(NKI)):nrow(NKI_comp),]

train_df_comp <- test_df
pred_df <- train_df

# New data
#train_df_comp <- cbind(train_df[, 1:(first_gene-1)], pcr_model$scores[(nrow(test_df)+1):nrow(NKI), 1:21])
#pred_df <- cbind(test_df[, 1:(first_gene-1)], pcr_model$scores[1:nrow(test_df), 1:21])

pred <- glm(eventdeath ~ ., data=as.data.frame(train_df_comp))
NKI_comp <- rbind(pred_df, train_df_comp)
    
predictedData <- predict(pred, newdata=pred_df, type="response")

#predictedDataExp <- 1/(1+exp(-predictedData))

acc <- getAccuracy(predictedData, pred_df[, 2], 0.65)
acc

# MCA 
famd.NKI <- FAMD(NKI_comp, ncp = ncol(NKI_comp), graph = TRUE)

eigen = as.data.frame(famd.NKI$eig)
plot(eigen$eigenvalue, type="b")
abline(h = mean(eigen$eigenvalue), lty = 2)

# DATA VISUALIZATION
pca.NKI <- PCA(NKI[, first_gene:ncol(NKI)], ind.sup = 1:nrow(test_df), ncp=nrow(NKI), graph=T)
eigen = as.data.frame(pca.NKI$eig)
plot(eigen$eigenvalue, type="b")
abline(h = mean(eigen$eigenvalue), lty = 2)


# DEC TREE
p2 = rpart(eventdeath~., data=train_df_comp, control=rpart.control(cp=0.001, xval=10))
printcp(p2)


p2$cptable = as.data.frame(p2$cptable) 
ind = which.min(p2$cptable$xerror) 
xerr <- p2$cptable$xerror[ind]
xstd <- p2$cptable$xstd[ind]
i= 1
while (p2$cptable$xerror[i] > xerr+xstd) {
    i = i+1
    alfa = p2$cptable$CP[i]
}
# AND PRUNE THE TREE ACCORDINGLY
p1 <- prune(p2,cp=alfa)

rpart.plot(p1)

p1$variable.importance

train_df_comp$eventdeath <- as.character(train_df_comp$eventdeath)
train_df_comp$eventdeath <- as.factor(train_df_comp$eventdeath)
pred_df[,2] <- as.character(pred_df[,2])
pred_df[,2] <- as.factor(pred_df[,2])

names(pred_df) <- make.names(names(train_df_comp))
names(train_df_comp) <- make.names(names(train_df_comp))
colnames(pred_df)
colnames(train_df_comp)

p1.rf = randomForest(formula = eventdeath ~ ., data = train_df_comp,
              importance = TRUE, xtest = pred_df[, c(1, 3:ncol(pred_df)) ], 
             ytest = pred_df[ , 2], ntree = 2000, nodesize = 10, maxnodes = first_gene)
p1.rf
# CLUSTERING 

d <- dist(rbind(famd.NKI$ind$coord), method = "euclidean")
hc <- hclust(d, method = "ward.D2")
plot(hc)
abline(h=27, col = 2)

barplot(hc$height, main="Heights")
abline(h=6, col = 2)

nc = 2
c1 <- cutree(hc,nc)
cdg<- aggregate(as.data.frame(famd.NKI$ind$coord),list(c1),mean)[,2:(ncol(famd.NKI$ind$coord)+1)]
Bss <- sum(rowSums(cdg^2)*as.numeric(table(c1)))
Tss <- sum(rowSums(famd.NKI$ind$coord^2))
Ib4 <- 100*Bss/Tss

kmeans = kmeans(as.data.frame(famd.NKI$ind$coord), centers=cdg)

plot3d(famd.NKI$ind$coord[,1], famd.NKI$ind$coord[,2], famd.NKI$ind$coord[,3], col = NKI$eventdeath+1)

for(cluster in 1:nc) {
    cluster_ind = famd.NKI$ind$coord[which(kmeans$cluster %in% cluster),]
    if (ncol(cluster_ind) > 1) {
        ellips <- ellipse3d(cov(cbind(cluster_ind[,1],cluster_ind[,2],cluster_ind[,3])),  centre = c(mean(cluster_ind[,1]), mean(cluster_ind[,2]), mean(cluster_ind[,3])), level = 0.95)
    } else {
        ellips <- ellipse3d(cov(cbind(cluster_ind[1],cluster_ind[2],cluster_ind[3])),  centre = c(mean(cluster_ind[1]), mean(cluster_ind[2]), mean(cluster_ind[3])), level = 0.95)
    }
   
    plot3d(ellips, col = "yellow", alpha = 0.5, add = TRUE, type = "shade")
}

catdesPlot = catdes(cbind(as.factor(kmeans$cluster),NKI_comp),1)
