require(ape)
require(igraph)

data<-read.csv(file.choose())  ##discrete codes csv

data<-na.omit(data)
data$Species<-as.character(paste(data$Species))
data[which(data[1]=="Neolamprologus_moori"),1]<-"Neolamprologus_moorii"

tree<-read.nexus(file.choose()) ##consensus tree

# Prune the tree 
allNames=tree$tip.label
keeps=data$Species
drops=allNames[!(allNames %in% keeps)]
tree=drop.tip(tree,drops)

## Check species are in both tree and list
list1<-as.character(data$Species)  # Extracts a vector of species names from the data
list2<-tree$tip.label  # Extracts a vector of species names from the tree
list1[which((list1 %in% list2) == FALSE)]  # Species in list1 missing in list2
list2[which((list2 %in% list1) == FALSE)]  # Species in list2 missing in list1

## Sort data in same order as tree
rownames(data)<-data$Species
data<-data[match(tree$tip.label,rownames(data)),] 
names(data)
colnames(data)<-c("Species", "C", "M", "G", "P", "D")


## Create a function for calculating the number of conditional independencies given the number of vertices and links 
# in each conceptual model
condNum <- function(V, A) {   ##V is number of vertices, A is number of causal links
  (factorial(V)/(2 * factorial(V - 2))) - A
}

## Create a function for calculating CIC values
CICc <- function(C, q, n) {
  C + 2 * q * (n/(n - 1 - q))
}

#######ENTER THE SET OF VERTICES YOU WILL USE (SAME FOR ALL MODELS)
variables <- data.frame(name=c("C", "M", "G", "P", "D"))
v = nrow(variables) ##number of vertices

#### Here's a function to determine the statistical models you need to test the conditional independencies in each model
#### It also runs the models (binaryPGLMM)


TestCondInd<-function(vertices, graph) {
  
  g=graph.data.frame(graph, directed=TRUE, vertices=vertices)
  plot(g)  #double check this is what you want
  
  x = shortest.paths(g)
  x[upper.tri(x, diag=TRUE)]<-NA
  y = which(x > 1, arr.ind=TRUE, useNames=TRUE)  #if the shortest path is > 1, then they are conditionally independent
  y<-cbind(rownames(y), colnames(x)[y[,2]])   #gives you the set of conditional independencies
  
  ifelse( condNum(v,length(E(g))) == nrow(y) , "OK",  stop("Something went wrong") )  ##double check that you have the correct number of conditional independencies
  
  models<-c()  ##a list of the statistical models used to test the conditional independencies
  for (i in 1:nrow(y)){
    civ1<-which( V(g)$name  == y[i,1] )
    civ2<-which( V(g)$name  == y[i,2] )
    
    
    a<-  V(g)$name [neighbors(g, civ1, mode="in") ]  ##causal parents of variable 1
    b<-  V(g)$name [neighbors(g, civ2, mode="in") ]  ## causal parents of variable 2
    
    #to determine which of the conditionally independent variables should be the 'response' variable in the model
    #you need to follow the direction of causality in the model
    #in some cases it is arbitrary (i.e. if both variables are at the end of causal chains), so just pick the first variable
    
    ifelse( shortest.paths(g, v = civ1, to = civ2 , mode = "all") != Inf,           ## if they are connected in a causal chain
            ifelse( shortest.paths(g, v = civ1, to = civ2 , mode = "out") < shortest.paths(g, v = civ1, to = civ2 , mode = "in"),
                    models [i] <- paste( y[i,2],"~", y[i,1] ), models [i] <- paste( y[i,1],"~", y[i,2] )), #the response var will be the most downstream 
            ifelse (length(b) > length(a), models [i] <- paste( y[i,2],"~", y[i,1] ), models [i] <- paste( y[i,1],"~", y[i,2] ))  ##otherwise, the response var will be the one with the most causal parents
    )
    
    #add the causal parents of both variables as independent variables
    if (length(a) > 0) for (j in 1:length(a)){ models[i] <- paste(models[i], "+", a[j]) } 
    if (length(b) > 0) for (k in 1:length(b)){ models[i] <- paste(models[i], "+", b[k]) }   ##in some cases a variable will be added twice on the right side of the formula. This doesn't effect model output so I just left it
  }
  
  pvals<-c()
  
  ##each of the models is run as a binaryPGLMM (ape) model to account for phylogenetic relationships
  for(w in 1:length(models)){
    nterm<-length(lm(models[w], data=data)$coefficients)  ##easiest way to extract how many terms are in the model
    
    z<-binaryPGLMM( formula(models[w]), B.init=matrix(c(rep(0,nterm)),nterm,1),   data=data, phy=tree ) 
    pvals[w]<-z$B.pvalue[2]    ##extract the p-value for the variable that is conditionally independent
}
  x<-c(pvals) 
  C <- -2 * (sum(log(x)))  # C - statistic
  CICc <- CICc (C, (v+length(E(g))), nrow(data))   # CIC value for model 
  list(models=models, C=C, CICc=CICc)
}


##################################
## Model A 
###Create the causal links in the model. 
modelA<-data.frame(rbind(
  c("M", "C"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MA<-TestCondInd(variables, modelA)   


## Model B 
###Create the causal links in the model. 
modelB<-data.frame(rbind(
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MB<-TestCondInd(variables, modelB)   

## Model C 
###Create the causal links in the model. 
modelC<-data.frame(rbind(
  c("M", "P"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MC<-TestCondInd(variables, modelC)   

## Model D 
###Create the causal links in the model. 
modelD<-data.frame(rbind(
  c("P", "M"),
  c("P", "G"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MD<-TestCondInd(variables, modelD)   

## Model E 
###Create the causal links in the model. 
modelE<-data.frame(rbind(
  c("M", "P"),
  c("P", "G"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

ME<-TestCondInd(variables, modelE)   

## Model F 
###Create the causal links in the model. 
modelF<-data.frame(rbind(
  c("M", "P"),
  c("G", "P"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MF<-TestCondInd(variables, modelF)   

## Model G 
###Create the causal links in the model. 
modelG<-data.frame(rbind(
  c("P", "M"),
  c("G", "P"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MG<-TestCondInd(variables, modelG)

## Model H 
###Create the causal links in the model. 
modelH<-data.frame(rbind(
  c("P", "M"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MH<-TestCondInd(variables, modelH)

## Model I 
###Create the causal links in the model. 
modelI<-data.frame(rbind(
  c("P", "G"),
  c("P", "C"),
  c("G", "C"),
  c("D", "C")
))

MI<-TestCondInd(variables, modelI)

## Model J 
###Create the causal links in the model. 
modelJ<-data.frame(rbind(
  c("P", "M"),
  c("P", "C"),
  c("G", "C"),
  c("D", "G")
))

MJ<-TestCondInd(variables, modelJ)

## Model K 
###Create the causal links in the model. 
modelK<-data.frame(rbind(
  c("M", "P"),
  c("P", "C"),
  c("G", "C"),
  c("D", "G")
))

MK<-TestCondInd(variables, modelK)


## Model L 
###Create the causal links in the model. 
modelL<-data.frame(rbind(
  c("M", "C"),
  c("M", "P"),
  c("P", "C"),
  c("G", "P"),
  c("G", "C"),
  c("D", "C")
))

ML<-TestCondInd(variables, modelL)


models<-list(MA,MB,MC,MD,ME,MF,MG,MH,MI,MJ,MK,ML)

cvals<-c() 
for(i in 1:12){cvals[i]<-models[[i]]$C}

CICvals<-c() 
for(i in 1:12){CICvals[i]<-models[[i]]$CICc}
results.dat<-data.frame(c(LETTERS[1:12]),cvals, CICvals )
colnames(results.dat)<-c("Model", "C", "CICc")
results.dat <- results.dat[order(results.dat$CICc), ]


deltaCIC<-c()
deltaCIC[1]<-0
for (i in 2:12){deltaCIC[i]<-results.dat$CICc[i]-31.79669}
results.dat$deltaCICc<-round(deltaCIC,2)
results.dat$likelihood<-round(exp(-1/2*results.dat$deltaCICc) ,3)
results.dat$weight<-round(results.dat$likelihood/sum(results.dat$likelihood),2)

results.dat

0.51/0.17


###Path coefficients for best model

ms6.1<-binaryPGLMM(C~P+G+D, B.init=matrix(c(rep(0,4)),4,1), phy=tree, data=data)
ms6.2<-binaryPGLMM(P~G+M, B.init=matrix(c(rep(0,3)),3,1), phy=tree, data=data)

