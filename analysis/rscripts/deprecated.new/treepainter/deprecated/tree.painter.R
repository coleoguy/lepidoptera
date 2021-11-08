tree.paintR.ver.3 <- function(tree = NULL,
                              tip_states = NULL,
                              qmat = NULL,
                              iter = NULL,
                              rate.class = NULL,
                              rate = NULL,
                              sample_mode = NULL,
                              return_br_table = T){
  # store trees in a new object
  tree1 <- tree
  # store traits in a new object
  trait <- tip_states
  # get the edge table
  edge.tab <- tree1$edge
  # get the likelihood of the original tree
  tree1.lik <- asr_mk_model(tree = tree1,
                            tip_states = trait,
                            transition_matrix = qmat,
                            Nstates = ncol(qmat))$loglikelihood
  # make a table to hold branch rates
  br.rate.class <- matrix(data = NA,
                          nrow = iter,
                          ncol = Nedge(tree))
  # fill the table
  for(j in 1:nrow(br.rate.class)){
    # get the second tree which we will change the branch lengths
    tree2 <- tree
    tree2.lik <- tree1.lik
    # tracker 1
    print(paste("working on iteration: ", j, sep = ""))
    # keep track of the sampled edges
    edges.sampled <- c()
    x <- T
    # randomly sample and edge and see if by changing the edge length the 
    # likelihood improves
    for(i in 1){
      while (x) {
        # randomly sample a single edge at a given iteration
        sampled.edge <- sample(1:Nedge(tree2), 1)
        # place holder for sampled edges
        edges.sampled <- c(sampled.edge, edges.sampled)
        # get the likelihood of the current tree
        current.tree2.lik <- tree2.lik
        # print to keep track
        # tracker 2
        print("A branch have been sampled successfully")
        print(paste("current branch:", sampled.edge))
        print(paste("remaining number of iteratons:", Nedge(tree2) - length(edges.sampled)))
        print(paste("Likelihood of the original tree:", tree1.lik))
        print(paste("Likelihood of the modified tree:", current.tree2.lik))
        # store the tree2 in two objects. One to test higher rate class and 
        # get the current rate class of the sampled edge
        pos.rates <- get.rates(tree = tree2,
                               edge = sampled.edge,
                               rate = rate)
        # get likelihoods
        liks <- vector(mode = "numeric", length = length(pos.rates))
        # store trees in a new vector
        trees <- vector(mode = "list", length = length(pos.rates))
        for(k in 1:length(pos.rates)){
          trees[[k]] <- tree2
        }
        class(trees) <- "multiPhylo"
        # calculate the likelihod
        for (k in 1:length(pos.rates)) {
          if(pos.rates[k] == 3){
            liks[k] <- current.tree2.lik
          }
          if(pos.rates[k] != 3){
            trees[[k]]$edge.length[sampled.edge] <- trees[[k]]$edge.length[sampled.edge]*rate[which(rate.class == pos.rates[k])]
            trees[[k]]$rates[sampled.edge] <- pos.rates[k]
            liks[k] <- asr_mk_model(tree = trees[[k]],
                                    tip_states = trait,
                                    transition_matrix = qmat,
                                    Nstates = ncol(qmat))$loglikelihood
            print(paste("Likelihood of the modified tree:", liks[k]))
            print(pos.rates[k])
          }
        }
        # get the best likelihood and make the tree with the highest likelihood
        # as the current tree.
        max.lik <- which(liks == max(liks))
        if(length(max.lik) == 1){
          tree2 <- trees[[max.lik]]
          tree2.lik <- liks[max.lik]
        }else{
          tree2 <- trees[[sample(max.lik),1]]
          tree2.lik <- liks[sample(max.lik,1)]
        }
        # stop once all edges are sampled
        if(length(edges.sampled) == Nedge(tree2))
          x <- F
      }
    }
    br.rate.class[j,] <- tree2$rates
  }
  if(sample_mode == "mean"){
    tree2$rates <- colMeans(br.rate.class)
  }
  if(sample_mode == "median"){
    for(i in 1:ncol(br.rate.class)){
      tree2$rates[i] <- median(br.rate.class[,i])
    }
  }
  if(return_br_table == T){
    results <- vector(mode = "list", length = 2)
    results[[1]] <- tree2
    results[[2]] <- br.rate.class
    
    names(results) <- c("tree", "branch.rate.class")
  }else{
    results <- tree2  
  }
  return(results)
}


### older functions ###

tree.paintR.ver.2 <- function(tree = NULL,
                              tip_states = NULL,
                              qmat = NULL,
                              iter = NULL,
                              rate.class = NULL,
                              rate = NULL,
                              sample_mode = NULL,
                              return_br_table = T){
  # store trees in a new object
  tree1 <- tree
  # store traits in a new object
  trait <- tip_states
  # get the edge table
  edge.tab <- tree1$edge
  # get the likelihood of the original tree
  tree1.lik <- asr_mk_model(tree = tree1,
                            tip_states = trait,
                            transition_matrix = qmat,
                            Nstates = ncol(qmat))$loglikelihood
  # make a table to hold branch rates
  br.rate.class <- matrix(data = NA,
                          nrow = iter,
                          ncol = Nedge(tree))
  # fill the table
  for(j in 1:nrow(br.rate.class)){
    # get the second tree which we will change the branch lengths
    tree2 <- tree
    tree2.lik <- tree1.lik
    # tracker 1
    print(paste("working on iteration: ", j, sep = ""))
    # keep track of the sampled edges
    edges.sampled <- c()
    x <- T
    # randomly sample and edge and see if by changing the edge length the 
    # likelihood improves
    for(i in 1){
      while (x) {
        # randomly sample a single edge at a given iteration
        sampled.edge <- sample(1:Nedge(tree2), 1)
        if(sampled.edge %in% edges.sampled){
          if(length(edges.sampled) < Nedge(tree1))
            while(sampled.edge %in% edges.sampled){
              sampled.edge <- sample(1:Nedge(tree2), 1)
            }
        }
        # place holder for sampled edges
        edges.sampled <- c(sampled.edge, edges.sampled)
        # get the likelihood of the current tree
        current.tree2.lik <- tree2.lik
        # print to keep track
        # tracker 2
        print("A branch have been sampled successfully")
        print(paste("current branch:", sampled.edge))
        print(paste("remaining number of edges:", Nedge(tree2) - length(edges.sampled)))
        print(paste("Likelihood of the original tree:", tree1.lik))
        print(paste("Likelihood of the modified tree:", current.tree2.lik))
        # store the tree2 in two objects. One to test higher rate class and 
        # the other to test lower rate class
        tree2.high <- tree2
        tree2.low <- tree2
        # get the current rate class of the sampled edge
        # current rate class of the root will be set to the middle class 
        # get the basal branches
        hit.base <- which(edge.tab[,1] == (Nnode(tree1)+2))
        if(sampled.edge == hit.base[1] | sampled.edge == hit.base[2]){
          current.rate.class <- tree2$rates[sampled.edge]
        }else{
          parent.of.sampled.edge <- edge.tab[sampled.edge,1]
          ancestral.edge <- which(edge.tab[,2] == parent.of.sampled.edge)
          current.rate.class <- tree2$rates[ancestral.edge]
          # get the position of the current rate class of the ancestral edge of
          # the sampled edge
          pos <- which(rate.class == current.rate.class)
          # calculate the likelihood one rate class above the ancestral rate class
          if((pos + 1) <= max(rate.class)){
            tree2.high$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos + 1]
            tree2.high.lik <- asr_mk_model(tree = tree2.high,
                                           tip_states = trait,
                                           transition_matrix = qmat,
                                           Nstates = ncol(qmat))$loglikelihood
            tree2.high$rates[sampled.edge] <- pos+1
            print(paste("Likelihood of the modified tree high rate class:", tree2.high.lik))
          }
          # if position of the rate class one step above the ancestral edge rate
          # class is greater than the max rate class the do not go beyond the 
          # ancestral rate class.
          if((pos + 1) > max(rate.class)){
            tree2.high$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos]
            tree2.high.lik <- asr_mk_model(tree = tree2.high,
                                           tip_states = trait,
                                           transition_matrix = qmat, 
                                           Nstates = ncol(qmat))$loglikelihood
            tree2.high$rates[sampled.edge] <- pos
            print(paste("Likelihood of the modified tree high rate class - max rate class reached:", tree2.high.lik))
          }
          # calculate the likelihood one rate class below the ancestral rate class
          if((pos - 1) >= min(rate.class)){
            tree2.low$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos - 1]
            tree2.low.lik <- asr_mk_model(tree = tree2.low,
                                          tip_states = trait,
                                          transition_matrix = qmat, 
                                          Nstates = ncol(qmat))$loglikelihood
            tree2.low$rates[sampled.edge] <- pos-1
            print(paste("Likelihood of the modified tree low rate class:", tree2.low.lik))
          }
          # if position of the rate class one step below the ancestral edge rate
          # class is smaller than the min rate class the do not go beyond the 
          # ancestral rate class.
          if((pos - 1) < min(rate.class)){
            tree2.low$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos]
            tree2.low.lik <- asr_mk_model(tree = tree2.low,
                                          tip_states = trait,
                                          transition_matrix = qmat, 
                                          Nstates = ncol(qmat))$loglikelihood
            tree2.low$rates[sampled.edge] <- pos
            print(paste("Likelihood of the modified tree low rate class - min rate class reached:", tree2.low.lik))
          }
        }
        # store the likelihood values
        liks <- vector(mode = "numeric", length = 3)
        liks <- c(current.tree2.lik, tree2.high.lik, tree2.low.lik)
        # get the best likelihood and make the tree with the highest likelihood
        # as the current tree.
        max.lik <- which(liks == max(liks))
        if(length(max.lik) == 1){
          if(length(max.lik)>1){
            print("yes")
          }
          if(max.lik == 1){
            tree2 <- tree2
            tree2.lik <- tree2.lik
          }
          if(max.lik == 2){
            tree2 <- tree2.high
            tree2.lik <- tree2.high.lik
          }
          if(max.lik == 3){
            tree2 <- tree2.low
            tree2.lik <- tree2.low.lik
          }
        }else{
          tree2 <- sample(c(tree2, tree2.high, tree2.low)[max.lik], 1)
          tree2 <- tree2[[1]]
          tree2.lik <- sample(liks, 1)
        }
        # stop once all edges are sampled
        if(length(edges.sampled) == Nedge(tree2))
          x <- F
      }
    }
    br.rate.class[j,] <- tree2$rates
  }
  if(sample_mode == "mean"){
    tree2$rates <- colMeans(br.rate.class)
  }
  if(sample_mode == "median"){
    for(i in 1:ncol(br.rate.class)){
      tree2$rates[i] <- median(br.rate.class[,i])
    }
  }
  if(return_br_table == T){
    results <- vector(mode = "list", length = 2)
    results[[1]] <- tree2
    results[[2]] <- br.rate.class
    
    names(results) <- c("tree", "branch.rate.class")
  }else{
    results <- tree2  
  }
  return(results)
}

#### old version ####

tree.paintR <- function(tree = NULL,
                        tip_states = NULL,
                        qmat = NULL,
                        iter = NULL,
                        rate.class = NULL,
                        rate = NULL,
                        sample_mode = NULL,
                        return_br_table = T){
  # store trees in a new object
  tree1 <- tree
  trait <- tip_states
  # get the likelihood of the original tree
  tree1.lik <- asr_mk_model(tree = tree1,
                            tip_states = trait,
                            transition_matrix = qmat,
                            Nstates = ncol(qmat))$loglikelihood
  # make a table to hold branch rates
  br.rate.class <- matrix(data = NA,
                          nrow = iter,
                          ncol = Nedge(tree))
  for(j in 1:nrow(br.rate.class)){
    tree2 <- tree
    print(paste("working on iteration: ", j, sep = ""))
    # set the initial likelihood of the tree 2 to same as that of tree 1
    tree2.lik <- tree1.lik
    # keep track of the sampled edges
    edges.sampled <- c()
    x <- T
    # randomly sample and edge and see if by changing the edge length the 
    # likelihood improves
    for(i in 1){
      while (x) {
        # randomly sample a single edge at a given iteration
        sampled.edge <- sample(1:Nedge(tree2), 1)
        if(sampled.edge %in% edges.sampled){
          if(length(edges.sampled) < Nedge(tree1))
            while(sampled.edge %in% edges.sampled){
              sampled.edge <- sample(1:Nedge(tree2), 1)
            }
        }
        # place holder for sampled edges
        edges.sampled <- c(sampled.edge, edges.sampled)
        # print to keep track
        print("A branch have been sampled successfully")
        print(paste("current branch:", sampled.edge))
        print(paste("remaining number of edges:", Nedge(tree2) - length(edges.sampled)))
        print(paste("Likelihood of tree 1:", tree1.lik))
        # store the tree2 in two objects. One to test higher rate class and 
        # the other to test lower rate class
        tree2.high <- tree2
        tree2.low <- tree2
        # get the current rate class of the sampled edge
        # current rate class of the root will be set to the middle class 
        if((sampled.edge-1) == 0 | (sampled.edge-1) == (which(tree2$edge[,1] %in% (Ntip(tree2) + 1))[2])-1){
          current.rate.class <- tree2$rates[sampled.edge]
        }else{
          current.rate.class <- tree2$rates[sampled.edge-1]
          # get the possition of the current rate class
          pos <- which(rate.class == current.rate.class)
          if((pos + 1) < max(rate.class)){
            tree2.high$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos + 1]
            tree2.high.lik <- asr_mk_model(tree = tree2.high,
                                           tip_states = trait,
                                           transition_matrix = qmat,
                                           Nstates = ncol(qmat))$loglikelihood
            tree2.high$rates[sampled.edge] <- pos+1
            print(paste("Likelihood of tree 2 high rate class:", tree2.high.lik))
          }
          if((pos + 1) > max(rate.class)){
            tree2.high$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos]
            tree2.high.lik <- asr_mk_model(tree = tree2.high,
                                           tip_states = trait,
                                           transition_matrix = qmat, 
                                           Nstates = ncol(qmat))$loglikelihood
            tree2.high$rates[sampled.edge] <- pos
            print(paste("Likelihood of tree 2 high rate class:", tree2.high.lik))
          }
          if((pos - 1) < min(rate.class)){
            tree2.low$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos - 1]
            tree2.low.lik <- asr_mk_model(tree = tree2.low,
                                          tip_states = trait,
                                          transition_matrix = qmat, 
                                          Nstates = ncol(qmat))$loglikelihood
            tree2.low$rates[sampled.edge] <- pos-1
            print(paste("Likelihood of tree 2 low rate class:", tree2.low.lik))
          }
          if((pos - 1) > min(rate.class)){
            tree2.low$edge.length[sampled.edge] <- tree2.high$edge.length[sampled.edge]*rate.class[pos]
            tree2.low.lik <- asr_mk_model(tree = tree2.low,
                                          tip_states = trait,
                                          transition_matrix = qmat, 
                                          Nstates = ncol(qmat))$loglikelihood
            tree2.low$rates[sampled.edge] <- pos
            print(paste("Likelihood of tree 2 low rate class:", tree2.low.lik))
          }
        }
        if(tree2.high.lik > tree1.lik | tree2.low.lik > tree1.lik ){
          if(tree2.high.lik > tree1.lik & tree2.low.lik > tree1.lik){
            if(tree2.high.lik > tree2.low.lik){
              tree2 <- tree2.high
            }
            if(tree2.low.lik > tree2.high.lik){
              tree2 <- tree2.low
            }
          }
        }
        if(length(edges.sampled) == Nedge(tree2))
          x <- F
      }
    }
    br.rate.class[j,] <- tree2$rates
  }
  if(sample_mode == "mean"){
    tree2$rates <- colMeans(br.rate.class)
  }
  if(sample_mode == "median"){
    for(i in 1:ncol(br.rate.class)){
      tree2$rates[i] <- median(br.rate.class[,i])
    }
  }
  if(return_br_table == T){
    results <- vector(mode = "list", length = 2)
    results[[1]] <- tree2
    results[[2]] <- br.rate.class
    
    names(results) <- c("tree", "branch.rate.class")
  }else{
    results <- tree2  
  }
  return(results)
}

