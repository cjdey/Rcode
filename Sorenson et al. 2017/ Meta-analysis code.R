library(ape)
library(MCMCglmm)


###import the data set##
data<-ES.All

###making binomial compatible with tree###
data$Latin<-as.factor(gsub(" ","_", data$Latin) )
data$Study<-as.factor(paste(data$Author, data$Year, sep = " "))

##Import the complete tree
tree<-read.nexus(file.choose())

###Priors
prior2 <- list(G = list(G1 = list(V = 1, nu = 0.02), G2 = list(V = 1, nu = 0.02)), R = list(V = 1, nu = 0.02))

##################################
#INTERCEPT MODELS
##############################

## phylogenetic
tree.inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
M1 <- MCMCglmm(Zr ~ 1, random = ~Latin + Study, 
               family = "gaussian", #mev = 1/(data$N - 3), 
               ginverse = list(Latin = tree.inv.phylo$Ainv), prior = prior2, 
               data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  
summary(M1)
#plot(M1)

predict(M1, marginal = NULL)

## non phylogenetic
M2 <- MCMCglmm(Zr ~ 1 , random = ~Latin + Study, 
               family = "gaussian", mev = 1/(data$N - 3), prior = prior2, 
              data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  ####5e+06
summary(M2) ## non-phylogenetic
#plot(M2)


###################
# Models including Category
#each category compared to 0
############

M3 <- MCMCglmm(Zr ~ Category - 1, random = ~Latin + Study, 
               family = "gaussian", mev = 1/(data$N - 3), 
               ginverse = list(Latin = tree.inv.phylo$Ainv), prior = prior2, 
               data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  
summary(M3)
#plot(M3)

## non phylogenetic
M4 <- MCMCglmm(Zr ~ Category - 1 , random = ~Latin + Study, 
               family = "gaussian", mev = 1/(data$N - 3), prior = prior2, 
               data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  ####5e+06
summary(M4) ## non-phylogenetic
#plot(M4)


###################
# Models including Category and Method
############

M5 <- MCMCglmm(Zr ~  Category + Method -1 , random = ~Latin + Study, 
               family = "gaussian", mev = 1/(data$N - 3), 
               ginverse = list(Latin = tree.inv.phylo$Ainv), prior = prior2, 
               data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  
summary(M5)
#plot(M5)

## non phylogenetic
M6 <- MCMCglmm(Zr ~ Category + Method - 1  , random = ~Latin + Study, 
               family = "gaussian", mev = 1/(data$N - 3), prior = prior2, 
               data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  ####5e+06
summary(M6) ## non-phylogenetic
#plot(M6)


DIC<-data.frame(Model = c("M1","M2","M3","M4","M5","M6"), DIC =c(M1$DIC,M2$DIC,M3$DIC,M4$DIC,M5$DIC,M6$DIC ))
##M4 has lowest DIC 


###################
# Models including Category and Method and Stage
############

M7 <- MCMCglmm(Zr ~  Category + Method + Stage -1 , random = ~Latin + Study, 
               family = "gaussian", mev = 1/(data$N - 3), 
               ginverse = list(Latin = tree.inv.phylo$Ainv), prior = prior2, 
               data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  
summary(M7)
#plot(M7)

## non phylogenetic
M8 <- MCMCglmm(Zr ~ Category + Method + Stage - 1  , random = ~Latin + Study, 
               family = "gaussian", mev = 1/(data$N - 3), prior = prior2, 
               data = data, nitt = 5e+05, thin = 500, burnin = 1e+05, pr=TRUE)  ####5e+06
summary(M8) ## non-phylogenetic
#plot(M8)


DICs<-data.frame(Model = c("M1","M2","M3","M4","M5","M6","M7","M8"), DIC =c(M1$DIC,M2$DIC,M3$DIC,M4$DIC,M5$DIC,M6$DIC,M7$DIC,M8$DIC ))
DICs<-DICs[order(DICs$DIC),]
DICs
##M4 has lowest DIC 




























## Egger's Regression
##use 'best' model
M4$Random$formula<-update(M4$Random$formula, ~.+leg(mev, -1, FALSE):units)

Pred<-predict(M4, marginal=~leg(mev, -1, FALSE):units)
Prec<-1/sqrt(data$ZrV) # Precision = 1/sqrt(V)
es<-data$Zr-Pred
data2<-data.frame(ES=es,Prec=Prec)

library(metafor)
res2<-rma(yi= es, sei= 1/Prec, method="REML", data=data2)

Egger<-regtest(res2,model="lm")
Egger


##Funnel plots
##on Raw data
plot(1/sqrt(data$ZrV)~data$Zr, ylab="Precision", xlab="Fisher's Z")
abline(v = 0)
abline(v = -0.15, lty=2)

#on residuals
plot(data2$Prec~data2$ES, ylab="Precision", xlab="Residual Effect Size (Zr)")
abline(v = 0)











