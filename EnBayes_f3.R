#GENETIC ALGORITHM FOR MAXIMIZING THE GP ACCURACY ON THE BASIS OF OBJECTIVE FUNCTION f3(w)#######################
#setwd("G:/Drive (F)/Genomic_selection/Analysis/Final_data/dataset/Dataset_G3")
library(GA)
model <- c("BayesCpi", "BayesA", "BayesL", "BayesR", "BayesB", "BayesC", "BayesBpi", "BayesRR")
wn_pcc_mse_diff <- matrix(0, nrow=8, ncol=nexp)
rownames(wn_pcc_mse_diff) <- model
ccolnames(wn_pcc_mse_diff) <- paste("exp_", 1:nexp, sep="")

for (i in 1:nexp){

path <- paste(i,"_bayes_result.txt", sep="")
zz <- read.table(path, header=TRUE) #read the Bayeisian predicted values

Actual=zz[,1]
BCPi=zz[,2]
BA=zz[,3]
BL=zz[,4]
BR=zz[,5]
BB=zz[,6]
BC=zz[,7]
BBPi=zz[,8]
BRR=zz[,9]

f3<- function(x)
{
g <- cor(Actual,(x[1]*BCPi+x[2]*BA+x[3]*BL+x[4]*BR+x[5]*BB+x[6]*BC+x[7]*BBPi+x[8]*BRR))- mean((Actual-x[1]*BCPi-x[2]*BA-x[3]*BL-x[4]*BR-x[5]*BB-x[6]*BC-x[7]*BBPi-x[8]*BRR)^2)
  h<-x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]
  k<-ifelse(h==1,g,g)
}
GA3 <- ga(type = "real-valued",
         fitness =f3,
         lower = c(0,0,0,0,0,0,0,0), upper = c(1,1,1,1,1,1,1,1),
         popSize = 300, maxiter = 1000)

ww_pcc_mse_diff <- as.numeric(attr(GA3, "solution")[1,])
wn_pcc_mse_diff[,i] <- as.matrix(ww_pcc_mse_diff/sum(ww_pcc_mse_diff), nrow=8, ncol=1) #change wn[,i]

}

#write the weight matrix
write.table(wn_pcc_mse_diff, "wt_pcc_mse_diff.txt", quote=FALSE, sep="\t")
wt_pcc_mse_diff <- read.table("wt_pcc_mse_diff.txt", header=TRUE)

gp_pcc_mse_diff_pcc <- gp_pcc_mse_diff_mse <- matrix(0, nrow=9, ncol=nexp)

rownames(gp_pcc_mse_diff_pcc )<-rownames(gp_pcc_mse_diff_mse )<-c(model,"EnBayes")
colnames(gp_pcc_mse_diff_pcc )<-colnames(gp_pcc_mse_diff_mse )<-paste("exp_", 1:nexp, sep="")

for(j in 1:nexp){
path <- paste(j,"_bayes_result.txt", sep="")
zz <- read.table(path, header=TRUE)
mm <- as.matrix(zz)

#PCC
pcc_bayes <- cor(mm)[1,-1]
kk_gpc <- mm[,-1]%*%as.matrix(wt_pcc_mse_diff[,j], ncol=1)
pcc_ensem <- cor(mm[,1], kk_gpc[,1])
gp_pcc_mse_diff_pcc[,j] <- as.matrix(c(pcc_bayes, pcc_ensem), nrow=9, ncol=1)

#MSE

rr <- matrix(rep(mm[,1],8), ncol=8, nrow=nrow(mm))
mse_bayes <- apply((mm[,-1]-rr)^2, 2, sum)/nrow(mm)
kk_mse <- mm[,-1]%*%as.matrix(wt_pcc_mse_diff[,j], ncol=1)
mse_ensem <- sum((mm[,1]-kk_mse[,1])^2)/nrow(mm)
gp_pcc_mse_diff_mse[,j] <- as.matrix(c(mse_bayes, mse_ensem), nrow=9, ncol=1)

} 

m1 <- paste(round(apply(gp_pcc_mse_diff_pcc, 1, mean),4),"\u00B1", round(apply(gp_pcc_mse_diff_pcc, 1, sd)/sqrt(nexp),4))
m2 <- paste(round(apply(gp_pcc_mse_diff_mse, 1, mean),4),"\u00B1", round(apply(gp_pcc_mse_diff_mse, 1, sd)/sqrt(nexp),4))

res <- data.frame(Model=c(model,"EnBayes_f3(w)"),PCC=m1, MSE=m2) #Accuracy (in terms of PCC and MSE along with standard errors)for 8 Bayesian models and EnByaes models while objective function f3(w) is used
write.table(res,"GP_Accuracy_f3(w).txt", row.names=FALSE, quote=FALSE, sep="\t") 
#PCC_GP_Accuracy is the result file of the Bayesian models along with the ensemble model "EnBayes"

