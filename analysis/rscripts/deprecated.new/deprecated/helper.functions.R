# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# helper functions

# whenwver there are multiple chromosome values or a range of chromosome values
# are given this function will sample a single chromosome number from the
# given range. Infput for this function is a data table of chromosome numbers

chromSampler <- function(dat){
  # read in data
  for(i in 1:nrow(dat)){
    # our data contains three types of chromosome ranges which are shown as
    # follows.
    
    # xx - yy (only extreme values are given seperated by a dash)
    # xx , yy (all possible values are given seperated by a comma)
    # xx - yy , zz, aa - bb (both extreme values and possible values are given seperated by both commas and dashes)
    
    # fist discard all the rows that contain information in multimple chromosome
    # numbers
    if(dat$notesMale2N[i] != ""){
      # this linle will select the rows that have a dash
      if(!(-1 %in% gregexpr(pattern = "-", dat$notesMale2N[i])[[1]])){
        # this will select rows that have a single dash
        if(length(gregexpr(pattern = "-", dat$notesMale2N[i])[[1]]) == 1){
          x <- unlist(strsplit(x = dat$notesMale2N[i], split = "-"))
          # here we will get all possible chromosome values between the two 
          # extremes and sample a single value. 
          # each possible chromosome value is given equal probability
          pos.chrom.vals <- seq(from = as.numeric(x[1]),
                                to = as.numeric(x[2]),
                                by = 2)
          dat$male2N[i] <-as.numeric(sample(pos.chrom.vals, size = 1))
        }
      }
      # this line will select the rows that contains a comma
      if(!(-1 %in% gregexpr(pattern = ",", dat$notesMale2N[i])[[1]])){
        # fist we seperate the lines that have a single comma and randomely 
        # sample a single value from the two possibilities
        if(length(gregexpr(pattern = ",", dat$notesMale2N[i])[[1]]) == 1){
          pos.chrom.vals <- unlist(strsplit(x = dat$notesMale2N[i], split = ","))
          dat$male2N[i] <- as.numeric(sample(pos.chrom.vals, size = 1))
        }
        # then we select rows that have multiple commas. Now in our dataset
        # rows that include multiple commas contains dashes as well.
        # fist we split all terms separated by a comma
        if(length(gregexpr(pattern = ",", dat$notesMale2N[i])[[1]]) >= 2){
          hit <- unlist(strsplit(x = dat$notesMale2N[i], split = ","))
          pos.chrom.vals <- c()
          # now we split terms that are seperated by a dash
          # if there are no dashes present then that is a single value which
          # was seperated in the previous step and we will store that as it is.
          # if a dash is present this suggests a range and we will get all possibilities
          # and finally sample a single value
          for(j in 1:length(hit)){
            x <- unlist(strsplit(x = hit[j], split = "-"))
            if(length(x) == 2){
              pos.chrom.vals <- c(pos.chrom.vals, seq(from = as.numeric(x[1]),
                                                      to = as.numeric(x[2]),
                                                      by = 2))
            }else{
              pos.chrom.vals <- c(pos.chrom.vals, x)
            }
            dat$male2N[i] <- as.numeric(sample(pos.chrom.vals, size = 1))
          }
        }
      }
    }
  }
  return(dat)
}

# whenwver there are multiple chromosome values or a range of chromosome values
# are given this function will sample a single chromosome number from the
# given range accounting for unknown values as well. 
# this finction perfoms as same as the above function. The code is similar.

chromSamplerMuSSE <- function(dat){
  for(i in 1:nrow(dat)){
    if(dat$notesMale2N[i] != ""){
      if(!(-1 %in% gregexpr(pattern = "-", dat$notesMale2N[i])[[1]])){
        if(length(gregexpr(pattern = "-", dat$notesMale2N[i])[[1]]) == 1){
          x <- unlist(strsplit(x = dat$notesMale2N[i], split = "-"))
          if((as.numeric(x[2]) - as.numeric(x[1])) == 2){
            pos.chrom.vals <- c(as.numeric(x[1]),
                                as.numeric(x[2]))
          }
          if((as.numeric(x[2]) - as.numeric(x[1])) > 2){
            pos.chrom.vals <- c(as.numeric(x[1]),
                                as.numeric(x[2]),
                                NA)
          }
          dat$male2N[i] <-as.numeric(sample(pos.chrom.vals, size = 1))
        }
      }
      if(!(-1 %in% gregexpr(pattern = ",", dat$notesMale2N[i])[[1]])){
        if(length(gregexpr(pattern = ",", dat$notesMale2N[i])[[1]]) == 1){
          pos.chrom.vals <- unlist(strsplit(x = dat$notesMale2N[i], split = ","))
          dat$male2N[i] <- as.numeric(sample(pos.chrom.vals, size = 1))
        }
        if(length(gregexpr(pattern = ",", dat$notesMale2N[i])[[1]]) >= 2){
          hit <- unlist(strsplit(x = dat$notesMale2N[i], split = ","))
          pos.chrom.vals <- c()
          for(j in 1:length(hit)){
            x <- unlist(strsplit(x = hit[j], split = "-"))
            if(length(x) == 2){
              pos.chrom.vals <- c(pos.chrom.vals,
                                  as.numeric(x[1]),
                                  as.numeric(x[2]),
                                  NA)
            }else{
              pos.chrom.vals <- c(pos.chrom.vals, x)
            }
            dat$male2N[i] <- as.numeric(sample(pos.chrom.vals, size = 1))
          }
        }
      }
    }
  }
  return(dat)
}





# this function will make chimaric taxa. Here this will look at species 
# level matches and then look at genera level matches. inputs are a data
# table with chromosome numbers and a phylogenetic tree

TraitOverlap <- function(dat = "../data/chroms/lep.chroms.csv", 
                         trees = "../data/trees/papilionoidea-10.trees",
                         hosts = "../data/hosts/papilionoidea-bisse-data.txt"){
  # first get species level matches
  
  # get packages
  # library(chromePlus)
  library(diversitree)
  
  # # read in chromosome data
  # dat <- read.csv(dat, as.is = T)
  # 
  # # read in trait data
  # hosts <- read.csv(hosts, as.is = T)
  # 
  # # read trees
  # trees <- read.tree(trees)
  
  # make a new column to store genera and family names
  dat$genera <- dat$family <- dat$hosts <- NA
  
  for(i in 1:nrow(dat)){
    dat$family[i] <- unlist(strsplit(x = dat$species[i], split = "_"))[1]
    dat$genera[i] <- unlist(strsplit(x = dat$species[i], split = "_"))[2]
  }
  
  dat <- dat[,c(11,12,1:10)]
  
  # make a new data table. Here we only want the species name, 
  # hapliod chromosome number the hosts and the tree tip name.
  # in hosts dataset 0 means generalists and 1 means specialists
  
  new.dat <- as.data.frame(matrix(data = NA,
                                  nrow = nrow(dat),
                                  ncol = 4))
  counter <- 1
  hit.genera <- c()
  hit.family <- c()
  colnames(new.dat) <- c("species", "haploid", "hosts","TreeTip")
  
  for(i in 1:length(trees[[1]]$tip.label)){
    current <- trees[[1]]$tip.label[i]
    if(current %in% dat$species){
      hit <- which(dat$species == current)
      hit.trait <- hosts$hosts[hosts$sp == current]
      if(length(hit) > 1) hit <- sample(hit, 1)
      if(length(hit.trait) > 1) hit.trait <- sample(hit.trait, 1)
      new.dat[counter, ] <- c(dat[hit, c(3,10)], hit.trait)
      hit.genera <- c(hit.genera, dat[hit,2])
      hit.family <- c(hit.family, dat[hit, 1])
      counter <- counter + 1
    }
  }
  
  # find genera found
  hit.genera <- unique(hit.genera)
  # find genera not found
  unhit.genera <- unique(dat$genera[!dat$genera %in% hit.genera])
  # split names at underscores
  # making a table with genus and species epethet split out to columns
  tree.taxa <- as.data.frame(matrix(data = NA,
                                    nrow = length(trees[[1]]$tip.label),
                                    ncol = 3))
  
  colnames(tree.taxa) <- c("family", "genera", "species")
  
  for(i in 1:length(trees[[1]]$tip.label)){
    tree.taxa$family[i] <- unlist(strsplit(x = trees[[1]]$tip.label[i], split = "_"))[1]
    tree.taxa$genera[i] <- unlist(strsplit(x = trees[[1]]$tip.label[i], split = "_"))[2]
    tree.taxa$species[i] <- trees[[1]]$tip.label[i]
  }
  
  for(i in 1:length(unhit.genera)){
    if(unhit.genera[i] %in% tree.taxa$genera){
      hit <- which(tree.taxa$genera == unhit.genera[i])
      if(length(hit) > 1) hit <- sample(hit, 1)
      hit.chroms <- which(dat$genera == tree.taxa$genera[hit])
      if(length(hit.chroms ) > 1) hit.chroms  <- sample(hit.chroms,  1)
      new.dat[counter, ] <- c(dat[hit.chroms , c(3,10)],
                              sample(hosts$hosts[hosts$sp %in% tree.taxa$species[which(tree.taxa$genera == unhit.genera[i])]],1),
                              sample(tree.taxa$species[which(tree.taxa$genera == unhit.genera[i])],1))
      new.dat$species[counter] <- paste(strsplit(new.dat$species[counter], "_")[[1]][2], "_sp",sep="")
      counter <- counter + 1
      hit.family <- c(hit.family, dat$family[hit.chroms])
    }
    # remove emplty rows from new.dat
    new.dat <- new.dat[!is.na(new.dat$species),]
  }
  return(new.dat)
}

# this function will prune the tree so that number of taxa in the tree matches
# with the number of species in our dataset. Also this will rename those tips
# that represent chimaric taxa accordingly
changeTipnames <- function(tree, dat){
  hit <- grep(pattern = "_sp", dat$species)
  new.tree <- keep.tip(phy = tree ,tip = dat$TreeTip)
  new.tree$tip.label[new.tree$tip.label %in% dat$TreeTip[hit]] <- dat$species[hit]
  return(new.tree)
}
