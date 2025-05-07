####EXTRACT GENOTYPIC AND PHENOTYPIC DATASETS######################

#Load the genotypic and phenotypic dataset provided in the .RData file or you can use your
#own genotypic and phenotypic dataset but remember to impute all the missing values

load ("EnBayes.RData")

#set the working directory in which the result file will be saved
#setwd("path")# supply the path of the directory

#The .RData file contains 18 phenotypic traits as mentioned in the manuscript
# (M_SS, M_WW, WY1,WY2, WY3, WY4, GY1 , GY2 , GY3 , GY4 , RY11 , RY12 ,
# BHU_Zn , BHU_Fe , IIWBR_Zn , IIWBR_Fe , PAU_Zn and PAU_Fe) 






#and 5 genotypic datasets 
# (maize, rice, groundnut, wheat_yield and wheat_nutrient). The corresponding 
#phenotypic and genotypic are as follows:
# For phenotype M_SS and M_WW, the genotypic data is "maize"
# For phenotype WY1,WY2, WY3 and WY4 the genotypic data is "wheat_yield"
#For phenotype GY1 , GY2 , GY3 and GY4  the genotypic data is "groundnut"
#For phenotype RY11 and RY12 the genotypic data is "rice"
# For phenotype BHU_Zn , BHU_Fe , IIWBR_Zn , IIWBR_Fe , PAU_Zn and PAU_Fe the genotypic data is "wheat_nutrient" 

#Let's do the genomic prediction for M_SS trait
#read the phenotypic and genotypic data

y <- phenotypic_datasets$M_SS[,1] # phenotypic dataset M_SS 
x <- genotypic_datasets$maize #genotypic dataset of maize 
