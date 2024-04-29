library(treeWAS)
library(phangorn)
##############################################
#test data provided by treeWAS
#used to check what kind of data structure etc.
data(snps)
str(snps)
#data(tree)
#str(tree)
#test <- snps
#test <- t(na.omit(t(test)))

#test[test == 1] <- "x"
#test[test == 0] <- ""
#test <- test[,1:30]
#test <- test[,13:30]
#own_data <- test
#l.edge <- length(tree$edge.length)
#print(tree$edge.length)

#treeWAS(snps)
#t = tree.reconstruct(snps, method="NJ")
#get.fitch.n.mts(snps, t)
##############################################

mypath = "/home/chris/Workspace/coin_treewas_mp/"
own_data <- read.csv(paste0(mypath,"test_files/mini_example.csv"), row.names = 1, stringsAsFactors = FALSE)
own_data <- read.csv("/home/chris/Workspace/coin_treewas_mp/coinfinder_data/gene_presence_absence.csv",row.names = 1)
own_data <- read.csv(paste0(mypath,"test_files/new_test_2.csv"), row.names = 1, stringsAsFactors = FALSE)
own_data <- read.csv(paste0(mypath,"test_files/mini_example_new.csv"), row.names = 1, stringsAsFactors = FALSE)


own_data <- Filter(function(x)!all(is.na(x)), own_data) #remove all NaN columns
own_data[own_data == "x"] <- 1
own_data[own_data == ""] <- 0
own_data <- as.matrix(sapply(own_data, as.numeric))

dna <- t(own_data) #transposing matrix to get data ready (same format) for treewas
print(dna)
#dna <- dna[,-c(1,2)]
#print(dna)
#dna <- test

######################################################################matrix where 1 and 0 changed to A and C
if(is.matrix(dna)){
  # dna <- as.DNAbin(dna)
  sp <- matrix(as.character(dna), nrow=nrow(dna), ncol=ncol(dna))
  rownames(sp) <- rownames(dna)
  colnames(sp) <- colnames(dna)
  
  ## Check/convert levels:
  levs <- unique(as.vector(unlist(sp)))
  nts <- c("a", "c", "g", "t")
  if(length(levs[!is.na(levs)]) > 4){
    stop("There must be no more than 4 unique values in dna, excluding NAs.")
  }
  if(!all(levs %in% nts)){
    for(i in 1:length(levs)){
      sp <- replace(sp, which(sp == levs[i]), nts[i])
    } # end for loop
  } # end levs conversion
  dna <- as.DNAbin(sp)
  rownames(dna) <- rownames(sp)
  colnames(dna) <- colnames(sp)
}
###########################################################distance matrix
#dna <- as.data.frame(dna)
D <- dist.dna(dna, model= "JC69")
print(D)
###########################################################

#neighbor joining and print tree
tree <- nj(D)
print(tree[["edge.length"]]) #negative distances, but apparently this
#happens in R??? Solution, set them to 0???
warning("negative branch lengths from neighbor joining have been set to zero")
tree$edge.length[tree$edge.length < 0] = 0
print(tree[["edge.length"]])

#midpoint = locating midpoint of longest path between any two tips and putting
#root in that location.
#ladderize = reorganizes internal structure of tree to get 
#ladderized effect when plotted
tree <- midpoint(ladderize(tree))
## Always work with tree in pruningwise order:
#tree <- reorder.phylo(tree, order="pruningwise")
#tree <- reorder.phylo(tree, order = "postorder")
print(tree[["edge.length"]])
## Trees must be rooted:
if(!is.rooted(tree)) tree <- midpoint(tree)
plot(tree)
plot(tree, edge.width=2, cex=0.5)
title("Neighbour-joining tree")
axisPhylo()

write.tree(tree)

print(tree[["edge.length"]])
## Get Anc-Des EDGE ORDER ##
############################
## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
## Get for loop index of rows in tree$edge[,1], in pairs (if binary tree), from lowest to highest:
x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
####################################################################
#########################################################
## get parsimomy cost for each SNP locus using fitch:
n.subs <- get.fitch.n.mts(x=dna, tree=tree)
n.subs <- table(n.subs) #same fitch scores
## handle n.subs "levels" with 0 SNP loci at those levels:
noms <- as.numeric(names(n.subs))
temp <- rep(0, max(noms))
for(i in 1:max(noms)){
  if(i %in% noms) temp[i] <- n.subs[which(noms==i)]
}
n.subs <- temp
names(n.subs) <- 1:length(n.subs)

#t <- tree.reconstruct(dna, method="NJ")
#fitch <- get.fitch.n.mts(dna, tree)  

## check:
#barplot(n.subs, col=transp("blue", 0.5), names=c(1:length(n.subs)))
#title("Homoplasy distribution \n(treeWAS Fitch)")
##################################################################

#################Distribution
dist <- n.subs     ########<- this one is important, rest not so much
## check for names first!
#if(!is.null(names(dist))){
  ## only modify if names are numeric
#  if(all.is.numeric(names(dist))){ ###############worked but now doesnt?????
#    noms <- as.numeric(names(dist))
#    aligned <- sapply(c(1:length(dist)), function(e) noms[e] == e)
    ## if any names do not correspond to their index, add zeros where missing:
#    if(any(aligned == FALSE)){
#      dist.new <- rep(0, max(noms))
#      dist.new[noms] <- dist
#      names(dist.new) <- c(1:length(dist.new))
#      dist <- dist.new
#    }
#  }
#} # end check for missing places
###############################################################################Adjust mutation to simulation size 10000

## get dist.prop, a distribution containing the counts
## of the number of SNPs to be simulated that will have
## i many substitutions
gen.size <- 10000
dist.sum <- sum(dist)
dist.prop <- round((dist/dist.sum)*gen.size)
## check that these counts sum to gen.size,
## else add the remainder to the largest n.subs count
## (OR should we just add these to the n.subs=1 set ??? ###
## likely to be the same thing, but could not be...)
if(sum(dist.prop) != gen.size){
  m <- which.max(dist.prop)
  #m <- 1
  if(sum(dist.prop) < gen.size){
    dist.prop[m] <- dist.prop[m] + (gen.size - sum(dist.prop))
  }
  if(sum(dist.prop) > gen.size){
    dist.prop[m] <- dist.prop[m] - (sum(dist.prop) - gen.size)
  }
}
## get rid of useless trailing 0s
while(dist.prop[length(dist.prop)] == 0){
  dist.prop <- dist.prop[c(1:(length(dist.prop)-1))]
}
##############################################################determine for each site how many mutations will happen
seed <- 1
## make n.mts, a vector of length ncol(snps)
n.mts <- rep(1111, gen.size)
loci.available <- c(1:gen.size)
## assign dist.prop[i] elements of n.mts
## to be the same as the n.subs
## indicated by i, the given element of dist.prop
for(j in 1:length(dist.prop)){
  ## provided there are not 0 sites to have this number of substitutions...
  if(dist.prop[j] > 0){
    if(length(loci.available) > 1){
      if(!is.null(seed)) set.seed(seed)
      ## assign dist.prop[i] elements of n.mts to be i
      loci.selected <- sample(loci.available, dist.prop[j], replace = FALSE)
      loci.available <- loci.available[-which(loci.available %in% loci.selected)]
    }else{
      ## if there is only 1 (the last) loci available,
      ## we select this one:
      loci.selected <- loci.available
    }
    n.mts[loci.selected] <- j
  }
}
####################################################################assign mutations to branches
## for each site, draw the branches to which
## you will assign the mts for this site
## (~ branch length):

l.edge <- length(tree$edge.length)
## Get vector of FALSEs of length tree$edge.length:
null.vect <- rep(FALSE, l.edge)

if(!is.null(seed)) set.seed(seed)
snps.loci <- sapply(c(1:length(n.mts)),
                    function(e)
                      replace(null.vect,
                              sample(c(1:l.edge),
                                     n.mts[e],
                                     replace=FALSE,
                                     prob=tree$edge.length), TRUE))


## rearrange snps.loci s.t it becomes a
## list of length tree$edge.length,
## each element of which contains the
## locations of the mutations that will
## occur on that branch
snps.loci <- sapply(c(1:nrow(snps.loci)),
                    function(e) which(snps.loci[e,] == TRUE))

## get the node names for all individuals (terminal and internal)
all.inds <- sort(unique(as.vector(unlist(tree$edge))))


# we will store the output in a list called snps:
snps <- list()

## Simulate genotype for root individual: ##

## For n.subs = n or = dist approaches:
gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
## we start w all inds having same genotype as root:
snps[all.inds] <- rep(list(gen.root), length(all.inds))

## store replacement nts in list new.nts:
new.nts <- list()

## distinguish btw list of loci and unique list
snps.loci.ori <- snps.loci

## will need to treat repeat loci differently...
snps.loci.unique <- lapply(snps.loci, unique)

#############################
## For Loop to get new nts ##
#############################
for(i in x){
  ## for all snps other than root, we mutate the
  ## genome of the node preceding it, according to snps.loci.
  ## Draw new nts for each locus selected for mutation:
  if(!.is.integer0(snps.loci.unique[[i]])){
    new.nts[[i]] <- !snps[[tree$edge[i,1]]][snps.loci.unique[[i]]]
    
    
    ## if any loci are selected for multiple mutations
    ## within their given branch length:
    if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
      ## identify which loci are repeaters
      repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
      ## how many times they repeat
      n.reps <- repeats - 1
      ## the positions of these loci in the vector of snps loci
      toRepeat <- which(snps.loci.unique[[i]] %in% names(repeats))
      ## run chain of re-sampling to end in our new nt for repeater loci:
      foo <- list()
      for(j in 1:length(toRepeat)){
        foo[[j]] <- new.nts[[i]][toRepeat[j]]
        for(k in 1:n.reps[j]){
          if(k==1){
            foo[[j]][k] <- !foo[[j]][1]
            
          }else{
            foo[[j]][k] <- !foo[[j]][k-1]
          }
        }
        ## retain only the last nt selected
        out <- sapply(c(1:length(foo)),
                      function(e) foo[[e]][length(foo[[e]])])
      }
      ## for the loci with repeated mts, replace these positions
      ## in new.nts with the corresponding elements of out, above.
      new.nts[[i]][toRepeat] <- out
    } # end of if statement for repeaters
    
    ## update ancestral genotype with new.nts:
    temp <- snps[[tree$edge[i,1]]]
    temp[snps.loci.unique[[i]]] <- new.nts[[i]]
    snps[[tree$edge[i,2]]] <- temp
    
  }else{
    ## if no mts occur on branch, set genotype of
    ## downstream individual to be equal to ancestor's
    snps[[tree$edge[i,2]]] <- snps[[tree$edge[i,1]]]
  }
} # end of for loop selecting new nts at mutator loci
###############################################determining
#how often a specific branch mutates
#collect <- c()
#i <- 1
#while (i<10001){
#  temporary <- which(snps.loci[,i])
#  collect <- c(collect, temporary)
#  i <- i+1
#}
#print(collect)

#branch_number_selected <- table(collect)
#barplot(branch_number_selected, col=transp("blue", 0.5), names=rownames(branch_number_selected))
#title("Branches selected to mutate")

#print(branch_number_selected)
#print(sum(branch_number_selected))




#test
#myTree <- ape::read.tree(text='((PagglomeransE325:0.000000, (Epyrifoliae:0.000000, (Etasmaniensis:0.000000, PvagansC91:-0.000000):0.000000):0.232616):0.006648, ((PananatisLMG:0.000000, PantoeaspaBvalens:0.232616):0.013297, (EbillingiaeEb661:0.250345, PantoeaspAt9b:0.000000):0.013297):0.006648, EamylovoraCFBP1430:0.000000);')
#plot(myTree)

#myTree <- ape::read.tree(text='((s3:0.527512, s1:0.000000):0.222335, (s4:0.000000, s2:0.823959):0.222335, s5:0.601624);')
#plot(myTree)


myTree <- ape::read.tree(text='(s4:0.2650185055746478,((s3:34.538776394910684,s2:34.538776394910684):34.538776394910684,s5:34.538776394910684):34.538776394910684,s1:0.415049208919844):0.0;')
myTree <- midpoint(ladderize(myTree))
## Always work with tree in pruningwise order:
#tree <- reorder.phylo(tree, order="pruningwise")
myTree <- reorder.phylo(tree, order = "postorder")
## Trees must be rooted:
if(!is.rooted(myTree)) myTree <- midpoint(myTree)

plot(myTree, edge.width=2, cex=0.5)





myTree <- ape::read.tree(text='((EbillingiaeEb661:0.568546,PantoeaspAt9b:0.0):0.03192700000000004,((PananatisLMG:0.0,PantoeaspaBvalens:0.44084):0.09578,((PagglomeransE325:0.0,(Epyrifoliae:0.0,(Etasmaniensis:0.0,PvagansC91:-0.0):0.0):0.44084):0.04789,EamylovoraCFBP1430:0.0):0.04789):0.06385299999999997);')
plot(myTree, edge.width=2, cex=0.5)

myTree <- ape::read.tree(text='((EbillingiaeEb661:0.568546,PantoeaspAt9b:0.0)1:0.03192700000000004,((PananatisLMG:0.0,PantoeaspaBvalens:0.44084)1:0.09578,((PagglomeransE325:0.0,(Epyrifoliae:0.0,(Etasmaniensis:0.0,PvagansC91:-0.0)1:0.0)1:0.44084)1:0.04789,EamylovoraCFBP1430:0.0)1:0.04789):0.06385299999999997);')
plot(myTree, edge.width=2, cex=0.5)