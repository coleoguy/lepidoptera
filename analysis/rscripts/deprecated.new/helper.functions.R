# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# helper functions

# whenwver there are multiple chromosome values or a range of chromosome values
# are given this function will sample a single chromosome number from the
# given range. Infput for this function is a data table of chromosome numbers

chromSampler <- function(dat){
  # this function will sample a single chromosome number from a given
  # range of chromosome numbers.
  for(i in 1:nrow(dat)){
    if(dat$notesMale2N[i] != ""){
      foo.a <- as.list(unlist(strsplit(dat$notesMale2N[i], split = ",")))
      foo.x <- vector(mode = "list", length = length(foo.a))
      for(j in 1:length(foo.a)){
        if((0 < gregexpr(pattern = "-", foo.a[[j]])[[1]])){
          x <- unlist(strsplit(x = foo.a[[j]], split = "-"))
          foo.x[[j]] <- seq(from = as.numeric(x[1]),
                            to = as.numeric(x[2]),
                            by = 2)
        }else{
          foo.x[[j]] <- as.numeric(foo.a[[j]])
        }
        dat$male2N[i] <-  sample(unlist(foo.x), size = 1)
      }
    }
  }
  return(dat)
}

sp.matches <- function(dat, trees, use.sub.species = F){
  dat.new <- dat
  if(use.sub.species == F){
    for(i in 1:nrow(dat)){
      dat.new$binomial[i] <- paste(unlist(strsplit(dat.new$binomial[i], split = " "))[1],
                                   unlist(strsplit(dat.new$binomial[i], split = " "))[2], 
                                   sep = "_")
    }
  }
  if(class(trees) == "multiPhylo"){
    tree <- trees[[1]]
  }
  
  if(class(trees) == "Phylo"){
    tree <- trees
  }
  # get the species level matches
  dat.new.s.m <- dat.new[dat.new$binomial %in%  tree$tip.label,]
  
  # from these species level matches if there are multiple hits for a single
  # species get only 1 hit
  unique.sp.matches <- unique(dat.new.s.m$binomial)
  # make an empty data frame to store values unique species matches
  dat.new.s.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.s.m$binomial)),
                                             ncol = ncol(dat.new)))
  # give column names
  colnames(dat.new.s.m.unique) <- colnames(dat.new)
  # fill in the data frame
  for(i in 1:length(unique.sp.matches)){
    dat.new.temp <- dat.new.s.m[dat.new.s.m$binomial %in% unique.sp.matches[i],]
    dat.new.s.m.unique[i,] <- dat.new.temp[sample(1:nrow(dat.new.temp), 1),]
  }
  # subset the original dataset by the species which had no hit to the tree 
  dat.new.unhit <- dat.new[!(dat.new$binomial %in% tree$tip.label),]
  # get the genera names of the species had matches to the tree dataset
  hit.genera <- c(unique(dat.new.s.m$Genus))
  # get the genera names that have no species level matches to the tree dataset
  unhit.gen <- unique(dat.new$Genus[!(dat.new$Genus %in% hit.genera)])
  # from the dataset that have species which have no species level matches remove
  # species that share the genus with the species that have a match with the tree
  dat.new.unhit <- dat.new.unhit[dat.new.unhit$Genus %in% unhit.gen,]
  # get the names of the species from the tree that have no species level match
  tree.unhit <- tree$tip.label[!(tree$tip.label %in% unique.sp.matches)]
  # get the genus names of the species hat have no species level match with the
  # trait dataset
  tree.unhit.gen <- c()
  for(i in 1:length(tree.unhit)){
    tree.unhit.gen[i] <- unlist(strsplit(tree.unhit[i], split = "_"))[[1]]
  }
  # get the index of the genera which have species level matches
  tree.hit.gen.num <- which(tree.unhit.gen %in% hit.genera)
  # get the names of the genera that appear on the phylogeny that have no species level match
  tree.true.unhit.gen <- tree.unhit.gen[!(tree.unhit.gen %in% hit.genera)]
  # get genera level matches
  dat.new.g.m <- dat.new.unhit[dat.new.unhit$Genus %in%  tree.true.unhit.gen,]
  # get the names of unique genera level matches
  unique.gn.matches <- unique(dat.new.g.m$Genus)
  # from these genera level matches if there are multiple hits for a single
  # genera get only 1 hit
  dat.new.g.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.g.m$Genus)),
                                             ncol = ncol(dat.new)))
  # give column names
  colnames(dat.new.g.m.unique) <- colnames(dat.new)
  # fill in the data frame
  for(i in 1:length(unique.gn.matches)){
    dat.new.temp <- dat.new.g.m[dat.new.g.m$Genus %in% unique.gn.matches[i],]
    dat.new.g.m.unique[i,] <- dat.new.temp[sample(1:nrow(dat.new.temp), 1),]
    dat.new.g.m.unique$binomial[i] <- paste(dat.new.g.m.unique$Genus[i], "_sp.", sep = "")
  }
  # get the index of genera level matches
  tree..unhit.hit.gen.num <- which(tree.unhit.gen %in% dat.new.g.m.unique$Genus)
  # get the names of the species as they appear on the phylogeny which have a
  # genera level match
  tree..unhit.hit.sp.names <- tree.unhit[tree..unhit.hit.gen.num]
  # this will be used to rename the names in the tree
  # make a data frame to store the species names of those species that have a
  # genera level match
  unhit.temp.dat.new <- as.data.frame(matrix(data = NA,
                                             nrow = length(tree..unhit.hit.sp.names),
                                             ncol = 2))
  # give column names
  colnames(unhit.temp.dat.new) <- c("gen", "sp")
  # fill in the species names 
  unhit.temp.dat.new$sp <- tree..unhit.hit.sp.names
  # give them the proper naming for genera level match
  for(i in 1:length(tree..unhit.hit.sp.names)){
    unhit.temp.dat.new$gen[i] <- unlist(strsplit(tree..unhit.hit.sp.names[i], split = "_"))[[1]]
  }
  # some genera are duplicated. get the unique genera name and randomly sample
  # a single species to represent these genera
  tree.gen.matches <- as.data.frame(matrix(data = NA,
                                           nrow = length(unique(unhit.temp.dat.new$gen)),
                                           ncol = 2))
  # give column names
  colnames(tree.gen.matches) <- c("gen", "sp")
  # fill in the data table  
  tree.gen.matches$gen <- unique(unhit.temp.dat.new$gen)
  for(i in 1:nrow(tree.gen.matches)){
    tree.gen.matches$sp[i] <- sample(unhit.temp.dat.new$sp[unhit.temp.dat.new$gen == tree.gen.matches$gen[i]],1)
  }
  # rename the genera
  tree.gen.matches$gen <- paste(tree.gen.matches$gen, "_sp.", sep = "")
  # for tree name correction
  tree.names <- as.data.frame(matrix(data = NA,
                                     nrow = (nrow(dat.new.s.m.unique) + nrow(dat.new.g.m.unique)),
                                     ncol = 2))
  colnames(tree.names) <- c("species_name", "name_on_tree")
  # fill in
  tree.names$species_name[1:nrow(dat.new.s.m.unique)] <- tree.names$name_on_tree[1:nrow(dat.new.s.m.unique)] <- dat.new.s.m.unique$binomial
  tree.names$species_name[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$gen
  tree.names$name_on_tree[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$sp
  # get the finalized dataset
  finaldat.new <- rbind(dat.new.s.m.unique, dat.new.g.m.unique)
  # get the results
  results <- list(finaldat.new, tree.names)
  names(results) <- c("chroms", "name_corrections")
  return(results)
}


buildQ <- function (pars, limits){ 
  print("building q-matrix")
  if (length(pars) != 10) 
    stop("pars should have length of 10")
  if (limits[2] < 100) 
    pad <- 2
  if (limits[2] >= 100) 
    pad <- 3
  if (limits[2] < 10) 
    pad <- 1
  parMat <- matrix(0, 2 * length(limits[1]:limits[2]), 
                   2 * length(limits[1]:limits[2]))
  colnames(parMat) <- sprintf(paste("%0", pad, "d", 
                                    sep = ""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)
  split <- ncol(parMat)/2
  chroms <- limits[1]:limits[2]
  for (i in 1:(split - 1)) {
    parMat[i, (i + 1)] <- pars[1]
    parMat[(i + 1), i] <- pars[3]
    if ((chroms[i] * 2) <= max(chroms)) 
      parMat[i, which(chroms == (chroms[i] * 2))] <- parMat[i, 
                                                            which(chroms == (chroms[i] * 2))] + pars[7]
    if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
      x <- chroms[i] * 1.5
      if (x%%1 == 0) 
        parMat[i, which(chroms == x)] <- parMat[i, 
                                                which(chroms == x)] + pars[5]
      if (x%%1 != 0) 
        parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- parMat[i, 
                                                                        which(chroms %in% c(floor(x), ceiling(x)))] + 
          (pars[5]/2)
    }
    parMat[i, (i + split)] <- pars[9]
    if (i == (split - 1)) 
      parMat[(i + 1), (i + 1 + split)] <- pars[9]
  }
  for (i in (split + 1):(nrow(parMat) - 1)) {
    parMat[i, (i + 1)] <- pars[2]
    parMat[(i + 1), i] <- pars[4]
    if ((chroms[i - split] * 2) <= max(chroms)) 
      parMat[i, (which(chroms[i - split] * 2 == chroms) + 
                   split)] <- parMat[i, (which(chroms[i - split] * 
                                                 2 == chroms) + split)] + pars[8]
    if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
      x <- chroms[i - split] * 1.5
      if (x%%1 == 0) 
        parMat[i, (which(chroms == x) + split)] <- parMat[i, 
                                                          (which(chroms == x) + split)] + pars[6]
      if (x%%1 != 0) 
        parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + 
                     split)] <- parMat[i, (which(chroms %in% c(floor(x), 
                                                               ceiling(x))) + split)] + (pars[6]/2)
    }
    parMat[i, (i - split)] <- pars[10]
    if (i == (nrow(parMat) - 1)) 
      parMat[(i + 1), (i + 1 - split)] <- pars[10]
  }
  return(parMat)
}