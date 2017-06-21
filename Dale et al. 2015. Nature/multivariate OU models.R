library(mvMORPH)
library(ape)

data<-read.csv("plumage_scores.csv")  ##full data set
tree<-read.tree("5831_species.tre")  #full tree
rownames(data)<-data$TipLabel
data<-data[,c(4:5)]   # just the color data
data <- data[match(tree$tip.label,rownames(data)),] 
ifelse(rownames(data)==tree$tip.label, "OK", "ERROR")   ##double check the tree and the data set are in the same order

####the mvOU models###
####these take a very long time to run (days or weeks) with the dataset from Dale et al. 2015 Nature
#####

M1<-mvOU(tree, data, model="OU1", param=list(sigma="constraint", method = "sparse", decomp="diagonalPositive" ))  
###independent evolution model
###sigma = constraint actually constrains the model to have a diagonal sigma matrix, which means the evolution of the two traits doesn't interact
###the decomp arguement specifies the decomposition of the alpha matrix, diagonalPositive restricts it to be diagonal (no interaction) and positive

M2<-mvOU(tree, data, model="OU1", method="sparse", param=list(decomp="diagonalPositive"))   
#here we don't restrict sigma to be diagonal, allowing for correlation in drift

M3<-mvOU(tree, data, model="OU1", method="sparse", param=list(decomp="symmetricPositive")) 
#here alpha is also not restricted to be diagonal, allowing correlated evolution (i.e. constrained evolution) in additional to correlated drift

###These are the brownian motion models###
###they dont have alpha matrices, only sigma matrices
###they can be run very quickly
####

M4<-mvBM(tree, data, model="BM1", method="pic", param=list(constraint="diagonal"))  
###each trait evolves independently under BM

M5<-mvBM(tree, data, model="BM1", method="pic")    
##allows correlation in the drift parameter
