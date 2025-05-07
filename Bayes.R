##################################################EXECUTION OF BAYESIAN MODEL FOR FINDING GENOMIC ESTIMATED BREEDING VALUES ####################################
library(hibayes) 
model <- c("BayesCpi", "BayesA", "BayesL", "BayesR", "BayesB", "BayesC", "BayesBpi", "BayesRR")
nrh <- floor(length(y)/5) 
y1 <- round(as.numeric(scale(y)),3) 
nexp <- 5 #number of times the number of experiments need to be repeated

for(j in 1:nexp){ #
set.seed(j)
sam <- sample(length(y)) # for randomly resuffling tht data
my_mod_g <- matrix(0, nrow=nrh, ncol=9) # the first column is the observed value 
y2 <- y1[sam]
y3 <- y2[1:nrh] #
y2[1:nrh]<-NA
x1 <- x[sam,] 
my_mod_g[,1] <- as.matrix(y3, ncol=1)
for(k in 1:8){
fit <- bayes(y=y2, M=x1, model = model[k], lambda=0.0001)
my_mod_g[,k+1] <- as.matrix(round(fit$g[1:nrh],3), ncol=1)

}
path1 <- paste(j,"_bayes_result.txt", sep="")
write.table(my_mod_g, path1, sep="\t", col.names=c("Observed",model), row.names=FALSE, quote=FALSE)
}
