get.rates <- function(tree = NULL, edge = NULL, rate = NULL){
  # store edge table
  edge.tab <- tree$edge
  # get the start node of the chosen edge
  start.node <- edge.tab[edge,1]
  # get the end node of the chosen edge
  end.node <- edge.tab[edge,2]
  # get the parent edge of the chosen edge
  parent.edge <- which(edge.tab[,2] == start.node)
  # if there is no parent edge (eg: basal banches) then set the parent edge
  # which is the root of the tree as the class that does not multiply the 
  # branch length
  # get all daughter edges of the chosen edge
  daugter.edges <- which(edge.tab[,1] == end.node)
  # get the parent edge rate clasee
  # if there is no parent edge (eg: basal banches) then set the parent edge
  # which is the root of the tree as the class that does not multiply the 
  # branch length
  if(length(parent.edge) == 0){
    parent.edge.rate.class <-  which(rate == 1)
  }else{
    parent.edge.rate.class <- tree$rates[parent.edge]  
  }
  # get the rate classes of all daughter edges
  daughter.edges.rate.class <- tree$rates[daugter.edges]
  # possible rate classes based on the parent edge
  pos.rc.1 <- c(parent.edge.rate.class-1,
                parent.edge.rate.class,
                parent.edge.rate.class+1)
  # possible rate classes based on the daughter edgees
  if(length(daugter.edges) == 0){
    pos.rc.2 <- NULL
  }else{
    pos.rc.2 <- c(daughter.edges.rate.class-1,
                  daughter.edges.rate.class,
                  daughter.edges.rate.class+1)
  }
  # get the intersect of parent and daughter edges
  # if a terminal branch is sampled then get rates from the parent brach
  if(is.null(pos.rc.2)){
    res <- pos.rc.1
  }else{
    res <- intersect(pos.rc.1, pos.rc.2)  
  }
  return(res)
}


