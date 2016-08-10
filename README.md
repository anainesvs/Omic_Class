## Integrating multiple Omics for Prediction of BC Outcome

The following scripts illustrate how to fit shrinkage and variable selection models to integrate omics. 
A full article using omic extended models is  presented in *Vazquez et al., Genetics, 2016*, the scripts used in the article are also provided at: [https://github.com/anainesvs/VAZQUEZ_etal_GENETICS_2016](https://github.com/anainesvs/VAZQUEZ_etal_GENETICS_2016).

**Contact**: avazquez@msu.edu

#### (1) Installing BGLR
The code below illustrates how to install and load the necessary package from CRAN using `install.packages()`.
```R
   install.packages(pkg='BGLR')    # install BGLR
   library(BGLR); 
```  

#### (2) Loading data
**Data**: The data can be obtained from NIH/NCI at: [https://gdc-portal.nci.nih.gov]. Download GE and Methylation chromosome 21, and clinical data. 
The code assumes that the user starts with data organized in a file called: `OMIC_DATA_class.rda` 
the objects that contain the phenotypic information, clinical covariates, and omic data. 
The code assumes that the file `OMIC_DATA_class.rda` contain the following objects:
   * `XF`: an incidence matrix for clinical covariates.
   * `Xge`: an incidence matrix for gene expression. 
   * `Xmt`: an incidence matrix for methylation values at various sites of chromosome 21 (only, since the full array has 450,000 CpG sites).
   * `y`: a matrix with an id, pathologic N, time to last follow up and alive status at the last follow up (0: alive, 1:death).
   * `XR`: a numeric vector with batches. 
The code below assumes that all the predictors were edited by removing outliers 
and predictors that did not vary in the sample, transformed if needed, and 
missing values were imputed.

#### (3) Clinical Covariates.
The following covariates are available: 
Race/ethnic group:  latino, african 
Age at diagnosis: age 
Cancer stage: PatStage 
Cancer histological type, specifically if the cancer is located in the lobular tissue: is_lobCarc 
"Cancer staging". National Cancer Institute. Retrieved 4 January 2013.

Response variable: y. In this class we will be using Pathological Nodal status: patN 
Pathologic N scores the degree of spread to regional lymph nodes at diagnosis. 
    N0: tumor cells absent from regional lymph nodes
    N1: regional lymph node metastasis present; (at some sites: tumor spread to closest or small number of regional lymph nodes)
    N2: tumor spread to an extent between N1 and N3 (N2 is not used at all sites)
    N3: tumor spread to more distant or numerous regional lymph nodes (N3 is not used at all sites)

Breast cancer subtypes can be defined as luminal types, Her2+ and Triple Negatives based on histochemstry tests. 
The last three columns of XF contain indicator variables showing the cancer subtype of each tumor.  

#### (4) Computing similarity matrices
 In this lab we will incorporate omics by incorporating correlated random effects. 
 For this, we will compute similarity between tumors in the methylation omic set and similarites in GE patterns with matrices of the form G= X * t(X) 
  computed from omics. The following code illustrates (3.a) editions necessary to perform before computing the similarity matrices, (3.b) code to compute this matrix.
 A similar code could be use to compute a G-matrix for gene expression, microRNA or other omics (see (6)).
 
 (4.a)  Exploring the methylation data
  1. Look at the data
  2. Study if there are missing values
  3. See if there is methylation sites that are constants
  4. Compute number of missing CpG sites per subject, remove sub >20% missing values
  5. Imputation of missing values
  6. Scale, centering, of distance between subjects (G)

```R 
#1
plot(Xmt[,1], type='l', col='gray', ylim=c(0,1))
for(i in 2:20){ lines(Xmt[,i], type='l', col=sample(1:20)) ; Sys.sleep(1)}

#2
namth<-sdmth<-numeric(ncol(Xmt))
for(i in 1:ncol(Xmt)){namth[i]<-sum(is.na(Xmt[,i]))}
hist(namth)

#3 Rm CpG sites that are constant
for(i in 1:ncol(Xmt)){sdmth[i]<-sd(Xmt[,i],na.rm=TRUE)} 
hist(sdmth)
Xmt<-Xmt[,sdmth>0.01]

#4 identify subjects with problematic assays
namth<-numeric(nrow(Xmt))
for(i in 1:nrow(Xmt)){namth[i]<-sum(is.na(Xmt[,i]))}
any(  (namth/ncol(Xmt)) >0.2 )

#5 naive imputation imputing mean of the CpG sites
for(i in 1:ncol(Xmt)){Xmt[,i]<-ifelse(is.na(Xmt[,i]), mean(Xmt[,i],na.rm=TRUE), Xmt[,i])}
any(is.na(Xmt))
```

(4.b)  Look and edit if needed Gene Expression data. The following code shows how to compute a similarity matrix for gene expression, please also do one for methylation called Gmt.

 ```R 
  #Computing a similarity matrix for gene-expression data
   Xge.s<- scale(Xge, scale=TRUE, center=TRUE) #centering and scaling
   Gge<-tcrossprod(Xge.s)                      #computing crossproductcts
   Gge<-Gge/mean(diag(Gge))                    #scales to an average diagonal value of 1.
   plot(Gge[,1])                               #observe Gge
   Gge[,100]
```
**NOTE**: for larger data sets it may be more convinient to use the `getG()` function of the [BGData](https://github.com/quantgen/BGData) R-package. This function allows computing G without loading all the data in RAM and offers methods for multi-core computing. 

#### (5)  Looking at the omic data to have insight of factors affecting its variability.
Observe PC (derived from methylation) and cancer subtypes. Explore which loading in PC 1 and 2 (PCs from both, methylations sets and gene expression). 

 ```R 
evGmt<- eigen(Gmt)
plot(evGmt$vector[,1:2], xlab='PC 1', ylab='PC 2', pch=20, col='gray80', cex=0.7, main='PC 1&2, CpG sites Ch 21, and Luminal Subtypes')
points(pcGmt$vectors[luminal,1:2], col='red', pch=8)
```

Omics can have important batch effects. Assume XR contains the batches where samples were analyzed for Methylation. The number of samples per bach can be seen using `table(XR)`
Then we look at the distribution of batches on the loading of the PC, and a regression of the batches on the PC 1 to 10. After this analysis we could consider the batch as a random effect in the model.

 ```R 
for(i in levels(XR)){
  plot(pcGmt$vectors[,1:2], xlab='pc1', ylab='pc2', pch=20, col='gray80', cex=0.7)
  points(pcGmt$vectors[XR==i,1:2], col='red', pch=8)
  Sys.sleep(1)
}

for(i in 1:10){print(paste('PC ',i)); print(summary(lm(pcGmt$vectors[,i]~XR)))}
```

#### (6)  Fitting a linear regression for (the "fixed effects" of) Clinical Coavariates using BGLR (COV)
The following code illustrates how to use BGLR to fit a fixed effects model. The matrix XF is an incidence matrix for clinical covariates. There is no column for intercept in XF because BGLR adds the intercept automatically. The response variable will be the number of nodes affected and will be assumed continuous. The argument `response_type` is used to indicate to BGLR that the response is gaussian. Predictors are given to BGLR in the form a two-level list. The argument `save_at` can be used to provide a path and a pre-fix to be added to the files saved by BGLR. For further details see [Pérez-Rodriguez and de los Campos, Genetics, 2014](http://www.genetics.org/content/genetics/198/2/483.full.pdf). The code also shows how to retrieve estimates of effects. In the examples below we fit the model using the default number of iterations (1,500) and burn-in (500). In practice longer chains are needed, the user can increase the numbrer of iterations or the burn-in using the arguments `nIter` and `burnIn` of `BGLR`.

```R
### Inputs
# centering and scaling the incidence matrix for fixed effects.
 XFc<- scale(XF, scale=FALSE, center=TRUE) 
 ETA.COV<-list( COV=list(X=XFc, model='FIXED') )
# defining the response variable
 yNM<-y[,2]
# Fitting the model
 fm=BGLR(y=yNM, ETA=ETA.COV, saveAt='cov_', response_type='gaussian')
# Retrieving estimates
 fm$ETA$COV$b      # posterior means of fixed effects
 fm$ETA$COV$SD.b   # posteriro SD of fixed effects
# Comparing posterior means of the fixed effects with OLS coefficients estimated with lm() function.
plot(x=coefficients(lm(yNM~XF))[-1], y=fm$ETA$COV$b , ylab='BGLR posterior means', xlab='lm solutions')
summary(lm(yNM~XF[,1:5]+triplenegative+luminal))
```

#### (6)  Fitting a binary/ordinal regression to the model above described. 
The following code illustrates how to use BGLR to fit the same fixed effects model but considering the response as a binary response, now yNM is 0 if there is not nodal metastasis and 1 if there is. 
```R
### Inputs
# generating the binary response variable.
 table(yNM, yNM01<- ifelse(yNM==0,0,1))
# Fitting a binary model
 fm01=BGLR(y=yNM01, ETA=ETA.COV, saveAt='cov_', response_type='ordinal')
# Retrieving estimates
 fm$ETA$COV$b      # posterior means of fixed effects
 fm$ETA$COV$SD.b   # posteriro SD of fixed effects
 head(fm$probs)    # estimated probabilities for the 0/1 outcomes.
 
# Fitting an ordinal model for a response 0,1,2,3
# fm=BGLR(y=yNM, ETA=ETA.COV, saveAt='cov_', response_type='ordinal')
```

#### (7)  Fitting a linear regression model for fixed effects and whole genome methylation (Metyl) using BGLR (COV+Metyl)
The following code illustrates how to use BGLR to fit a mixed effects model that accomodates both clinical covariates and whole-genome-methylation. 
```R
# Setting the linear predictor
  ETA.COV.Methyl<-list( COV=list(X=XFc, model='FIXED'), Methyl=list(K=Gmt, model='RKHS'))
# Fitting the model
  fm.COV.Methyl<- BGLR(y=yNM, ETA=ETA.COV.Methyl, response_type='gaussian',saveAt='cov_mt_')
#  Retrieving predictors
  fm.COV.Methyl$mu            # intercept
  fm.COV.Methyl$ETA$COV$b     # effects of covariates
  fm.COV.Methyl$ETA$Methyl$varU   # variance associated to methylation SD.varU gives posterior SD
  fm.COV.Methyl$ETA$Methyl$u      # random effects associated to gene expression
  plot(scan('cov_mt_ETA_Methyl_varU.dat'),type='o',col=4) # trace plot of variance of methylation in chromosome 21.
```

#### (8)  Covariates+Gene Expression information.
Please, fit an ordinal model to the yNM variable (0,1,2,3) for fixed effects covariates and gene expression (COV+GE)
**NOTE**: to fit a similar model for COV+GE one just needs to change the inputs in the defintiion of the linear predictor by providing Gge instead of Gmt.

### (9)
The model above is not accounting for batch, thus we can either pre-correct GE and methylation by batch effects or we can incorporate batches to the model. Next there is an example in how to incorporate batches of methylation in the model:
```R
# Setting the linear predictor
  ETA.COV.Methyl<-list( COV=  list(X=XFc, model='FIXED'), 
                        Batch=list(~factor(MB), model='BRR'),
                        Methyl=list(K=Gmt, model='RKHS')
                       )
# Fitting the model
  fm.COV.Methyl<- BGLR(y=yNM, ETA=ETA.COV.Methyl, response_type='gaussian',saveAt='cov_mt_')
#  Retrieving predictors
  fm.COV.Methyl$mu            # intercept
  fm.COV.Methyl$ETA$COV$b     # effects of covariates
  fm.COV.Methyl$ETA$Methyl$varU   # variance associated to Methyl SD.varU gives posterior SD
  fm.COV.Methyl$ETA$Methyl$u      # random effects associated to methylation
  plot(scan('cov_mt_ETA_Methyl_varU.dat'),type='o',col=4) # trace plot of variance of the methylation sites at chromosome 21.
  fm.COV.Methyl$ETA$Batch$varB    #v ariance of the batches and batch effects     
  fm.COV.Methyl$ETA$Batch$b
```

#### (10)  The following code demostrates how to extend the model to fixed effects covariates and 2 omics (COV+GE+METH)
The model `COV+GE` was extended to incorporate methylation data.
```R
#Computing a similarity matrix for methylation data
ETA.COV.GE.MT<-list( COV=  list(X=XF, model='FIXED'),
                     GE=   list(K=Gge, model='RKHS'),
                     METH= list(K=Gmt, model='RKHS'))
# Fitting models 
fm.COV.GE.MT<- BGLR(y=y, ETA=ETA.COV.GE.MT, 
                 response_type='ordinal',saveAt='cov_ge_mt_')
```

#### Extensions
Code to model omic by omic interactions, and validation to evaluate prediction accuracy are provided for a similar example at 
[https://github.com/anainesvs/VAZQUEZ_etal_GENETICS_2016](https://github.com/anainesvs/VAZQUEZ_etal_GENETICS_2016).
Omic by a systematic effect interaction (e.g. omic by treatment) can also be accomodated. See: Gonzalez et al., 2016 submitted. 

