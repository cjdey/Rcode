###This is a new (& improved?) function to overwrite the Discrete() function in the
##BayesTrait wrapper package created by Randi Griffin
##This package wraps the free software BayesTrait for Bayesian Phylogenetic Comparative Analysis


##This improved function formats the output in a way that is convenient for using the data in 
## TRACER (another free software useful for interpreting the results from BayesTrait).
##It also creates the 'Schedule' file in the BayesTrait output, which is not created by default in BayesTraitV2
## but is useful for examining acceptance rates.


##To Install the package
#library(devtools)
#install_github("rgriff23/btw")



Discrete<-function (tree, data, mode = "ML", dependent = FALSE, res = NULL, 
          resall = NULL, mrca = NULL, fo = NULL, mlt = 10, it = 1e+05, 
          bi = 5000, sa = 100, pr = NULL, pa = NULL, hp = NULL, hpall = NULL, 
          rj = NULL, rjhp = NULL, RD = NULL, silent = TRUE) 
{
  if (class(tree) == "phylo") {
    tree$node.label = NULL
    treelabs = tree$tip.label
  }
  else if (class(tree) == "multiPhylo") {
    treelabs = attributes(tree)$TipLabel
  }
  else {
    stop("Tree must be of class phylo or multiPhylo")
  }
  if (!(class(data[, 1]) %in% c("character", "factor"))) {
    stop("First column of data should contain species names.")
  }
  if (length(setdiff(treelabs, data[, 1])) > 0) {
    stop(paste("No match found in the data:", paste(setdiff(tree$tip.label, 
                                                            data[, 1]), collapse = ", ")))
  }
  if (length(setdiff(data[, 1], treelabs)) > 0) {
    stop(paste("No match found in the phylogeny:", paste(setdiff(data[, 
                                                                      1], tree$tip.label), collapse = ", ")))
  }
  if (length(setdiff(treelabs, data[, 1])) > 0 | length(setdiff(data[, 
                                                                     1], treelabs)) > 0) {
    stop("Species in your phylogeny and data must match up exactly.")
  }
  if (ncol(data) > 3) {
    stop("Too many columns in data: BayesTraits can only analyze one or two discrete traits.")
  }
  if (!exists(".BayesTraitsPath") | !file.exists(.BayesTraitsPath)) {
    stop("Must define '.BayesTraitsPath' to be the path to BayesTraitsV2 on your computer. For example: .BayesTraitsPath <- User/Desktop/BayesTraitsV2")
  }
  if (mode == "Bayesian") {
    mode = 2
  }
  else {
    mode = 1
  }
  if (dependent == FALSE) {
   model = 2
  }
  else (model = 3)
  input = c(model, mode)
  
  
  if (!is.null(res)) {
    for (i in 1:length(res)) {
      input = c(input, paste("Restrict", res[i]))
    }
  }
  if (!is.null(resall)) {
    input = c(input, paste("resall", resall))
  }
  if (!is.null(mrca)) {
    for (i in 1:length(mrca)) {
      input = c(input, paste("mrca", paste("mrcaNode", 
                                           i, sep = ""), mrca[i]))
    }
  }
  if (!is.null(fo)) {
    for (i in 1:length(fo)) {
      input = c(input, paste("Fossil", paste("fossilNode", 
                                             i, sep = ""), fo[i]))
    }
  }
  if (mode == 1) {
    input = c(input, paste("mlt", as.numeric(mlt)))
  }
  if (mode == 2) {
    input = c(input, paste("it", format(it, scientific = F)))
    input = c(input, paste("bi", format(bi, scientific = F)))
    input = c(input, paste("sa", format(sa, scientific = F)))
    if (!is.null(pr)) {
      for (i in 1:length(pr)) {
        input = c(input, paste("prior", pr[i]))
      }
    }
    if (!is.null(pa)) {
      input = c(input, paste("pa", pa))
    }
    if (!is.null(rj)) {
      input = c(input, paste("rj", rj))
    }
    if (!is.null(hp)) {
      for (i in 1:length(hp)) {
        input = c(input, paste("hp", hp[i]))
      }
    }
    if (!is.null(hpall)) {
      input = c(input, paste("Hpall", hpall))
    }
    if (!is.null(rjhp)) {
      input = c(input, paste("rjhp", rjhp))
    }
    if (!is.null(RD)) {
      input = c(input, paste("RD", RD))
    }
    
  }
  
  filename = paste("./", var1, var2, "_", ifelse (dependent == FALSE, "ind", "dep"), ".log.txt", sep="")
  input = c(input, paste("lf", filename))
  input = c(input, paste("Schedule")) 
  input = c(input, "run")
  write(input, file = "./inputfile.txt")
  ape::write.nexus(tree, file = "./tree.nex", translate = T)
  write.table(data, file = "./data.txt", quote = F, col.names = F, 
              row.names = F)
  system(paste(.BayesTraitsPath, "./tree.nex", "./data.txt", 
               "< ./inputfile.txt"), ignore.stdout = silent)
  Skip = grep("Tree No", scan(file = filename, what = "c", 
                              quiet = T, sep = "\n", blank.lines.skip = FALSE)) - 1
  Results = read.table(filename, skip = Skip, sep = "\t", 
                       quote = "\"", header = TRUE)
  Results = Results[, -ncol(Results)]
  system(paste("rm", filename))
  system(paste("rm ./inputfile.txt"))
  system(paste("rm", "./tree.nex"))
  system(paste("rm", "./data.txt"))
  return(Results)
}
