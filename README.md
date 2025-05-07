# Overview
We propose an approach called EnBayes—an ensemble method based on Bayesian alphabet models—to improve genomic prediction accuracy. This method uses a constrained weight optimization strategy powered by a genetic algorithm to assign optimal weights to predicted values generated from eight Bayesian models: BayesCpi, BayesA, BayesL, BayesR, BayesB, BayesC, BayesBpi, and BayesRR. The optimal weights are determined using four different objective functions, denoted as f₁(w), f₂(w), f₃(w), and f₄(w).

- f₁(w): Maximizes the Pearson correlation coefficient (PCC) between the ensemble predictions and the observed trait values.
- f₂(w): Minimizes the mean squared error (MSE) between the ensemble predictions and the observed trait values.
- f₃(w): Maximizes the difference between PCC and MSE.
- f₄(w): Maximizes the ratio of PCC to MSE.

All four objective functions are subject to the constraint that the sum of all weights equals 1, and each weight must lie between 0 and 1. A step-by-step procedure for implementing the EnBayes method using R code is provided below.

# Requirements
R 4.0.2 version (The program runs only in the R version 4.0.2, no older or latest version will work)

## R packages
hibayes 1.0.0 (install only version 1.0.0, no older or latest version will work)
- Link to download hibayes 1.0.0 -> https://cran.r-project.org/src/contrib/Archive/hibayes/hibayes_1.0.0.tar.gz
- After downloading hibayes_1.0.0.tar.gz file, run following command in the R interface to install hibayes 1.0.0

      > install.packages("hibayes_1.0.0.tar.gz")
             
- To check if the hibayes package has been downloaded successfully run following command in the R interface

      > library(hibayes)
    
## Models
    viz.
    BayesCpi
    BayesA
    BayesL
    BayesR
    BayesB
    BayesC
    BayesBpi
    BayesRR

# Download files
- Bayes.R
- EnBayes_f1.R
- EnBayes_f2.R
- EnBayes_f3.R
- EnBayes_f4.R
- EnBayes.RData
- Pheno_Geno_Data.R

# File description
### Bayes.R
R script to run the eight Bayesian models i.e BayesCpi, BayesA, BayesL, BayesR, BayesB, BayesC, BayesBpi and BayesRR. Repeated cross-validation approach is adopted, where 20% of the dataset was used for testing and rest 80% for training in each repeatation (nexp)

### EnBayes_f1.R
R script contains genetic algorithm for maximizing the Genomic Prediction (GP) accuracy on the basis of objective function f1(w) by maximixing the pearson correlation coefficient (PCC) between the predicted values of ensemble model and their respective observed values

### EnBayes_f2.R
R script contains genetic algorithm for maximizing the GP accuracy on the basis of objective function f2(w) by minimizing the mean square error (MSE) between the predicted values of ensemble model and their respective observed values

### EnBayes_f3.R
R script contains genetic algorithm for maximizing the GP accuracy on the basis of objective function f3(w) by maximizing the difference of PCC and MSE i.e. PCC-MSE

### EnBayes_f4.R
R script contains genetic algorithm for maximizing the GP accuracy on the basis of objective function f4(w) by maximizing the ratio of PCC to MSE i.e. PCC/MSE

### EnBayes.RData
contains phenotypic and genotypic dataset for maize, rice, groundnut and wheat

### Pheno_Geno_Data.R
R script to load the genotypic and phenotypic dataset provided in the EnBayes.RData

# Usage 
- Create a working directory.
- Place all the downloaded files in the working directory.
- First, execute the Pheno_Geno_Data.R script from the R interface to load the dataset provided in EnBayes.RData. We have provided a total of 5 genotypic datasets and there corresponding 18 phenotypic traits. User can choose any dataset by using the symbols listed in the table below. Alternatively, user can provide their own phenotypic and genotypic files—just ensure that all missing values are imputed. Place the phenotypic file against y and genotypic file against x in the Pheno_Geno_Data.R file.

    - y (phenotypic file) represents the numeric vector of trait values for n lines/individuals
    - x (genotypic file) represents the data frame of n rows and m columns, where n is the number of lines/individuals and m is the number of markers.

      | |Genotypic Dataset|Phenotypic Dataset|
      |---|---|---|
      |1|maize|M_SS, M_WW|
      |2|rice|RY11, RY12|
      |3|groundnut|GY1, GY2, GY3, GY4|
      |4|wheat_yield|WY1, WY2, WY3, WY4|
      |4|wheat_nutrient|BHU_Zn, BHU_Fe, IIWBR_Zn, IIWBR_Fe, PAU_Zn, PAU_Fe|

          [6] load ("EnBayes.RData")
          [32] y <- phenotypic_datasets$M_SS[,1] #phenotypic dataset M_SS 
          [33] x <- genotypic_datasets$maize #genotypic dataset of maize
  
- After running the Pheno_Geno_Data.R script, execute the Bayes.R script to run the eight Bayesian models from the same R interface. The user can modify how many times the experiments are repeated by replacing the value against nexp, the default value is 5.

      [6] nexp <- 5 #number of times the number of experiments need to be repeated

- The Bayes.R script will generate an output file which will be used in the R script for genetic algorithms i.e EnBayes_f1.R, EnBayes_f2.R, EnBayes_f3.R and EnBayes_f4.R.
- Now execute EnBayes_f1.R, EnBayes_f2.R, EnBayes_f3.R and EnBayes_f4.R from the same R interface.
   
NOTE: All R codes should be executed within the same R interface to avoid potential errors. Alternatively, if users prefer to run the scripts via the command line (cmd), they can combine all R files in a format—Pheno_Geno_Data.R → Bayes.R → EnBayes_f1.R → EnBayes_f2.R → EnBayes_f3.R → EnBayes_f4.R into a single script and execute it using a single command in the cmd.

    path/to/R-4.0.2/Rscript <script.R>
        
Example 
      
    C:/PROGRA~1/R/R-4.0.2/bin/Rscript combined_code.R

# Output description
Once the script is executed, an individual table will be generated for each of the four objective functions i.e f1(w), f2(w), f3(w) and f4(w) as shown below. Each table’s first column lists the eight Bayesian models, with the ensemble model in the last row. The second column displays the PCC values along with their standard errors, while the third column presents the MSE values and their respective standard errors.

|Model|PCC|MSE|
|---|---|---| 
|BayesCpi|0.4545 ± 0.0306|0.9325 ± 0.1056|
|BayesA|0.4565 ± 0.0197|0.8777 ± 0.089|
|BayesL|0.4591 ± 0.0261|0.8814 ± 0.0913|
|BayesR|0.4596 ± 0.0242|0.8803 ± 0.0867|
|BayesB|0.4074 ± 0.0194|0.9032 ± 0.0893|
|BayesC|0.4134 ± 0.0231|0.9013 ± 0.0939|
|BayesBpi|0.4224 ± 0.0373|0.9256 ± 0.1143|
|BayesRR|0.466 ± 0.0237|0.8747 ± 0.077|
|EnBayes_f1(w)|0.4796 ± 0.0217|0.8565 ± 0.0845|
