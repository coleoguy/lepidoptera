# Terrence Sylvester
# 28th October 2020
# pradakshanas@gmail.com
# process and visualise the mcmc output

# load results
load("../results/w.poly.all.matches.musse.MCMC.RData")
results <- x
# keep results and remove other data from the environment
rm(list = ls()[-which(ls() %in% "results")])

i <- 1
results[[i]]$p <- results[[i]]$p * (1/results[[i]]$i[1])

BurninCalc <- function(dat = NULL,
                       window.size = NULL,
                       nGen = NULL,
                       likelihood.Column = NULL,
                       get.post.burnin = NULL){
  
  results <- dat
  
  window.size <- 5
  nGen <- 100
  
  goodRuns <- badRuns <- as.data.frame(matrix(data = NA,
                                              nrow = length(results),
                                              ncol = 2))
  
  colnames(goodRuns) <- colnames(badRuns) <- c("iteration", "generation")
  EndofBurn <- c()
  
  tTestVar <- tTestSD <- tTestPecent <- as.data.frame(matrix(data = NA,
                                                             nrow = length(results) * nGen,
                                                             ncol = 3))
  
  colnames(tTestVar) <- colnames(tTestSD) <- colnames(tTestPecent) <- c("tree", "gen", "p-value")
  
  counter <- 1
  for(i in 1:length(results)){
    varVector <- vector(mode = "numeric", length = nrow(results[[i]]))
    sdVector <- vector(mode = "numeric", length = nrow(results[[i]]))
    # meanVector <- vector(mode = "numeric", length = nrow(results[[i]]))
    pecentChange <- vector(mode = "numeric", length = nrow(results[[i]]))
    degreeOfChange <- vector(mode = "numeric", length = nrow(results[[i]]))
    
    atanVectorLikelihood <- vector(mode = "numeric", length = nrow(results[[i]]))
    atanVectorVar <- vector(mode = "numeric", length = nrow(results[[i]]))
    atanVectorSD <- vector(mode = "numeric", length = nrow(results[[i]]))
    atanVectorMean <- vector(mode = "numeric", length = nrow(results[[i]]))
    
    
    for(j in 2:(nrow(results[[i]]))){
      pecentChange[j] <- ((results[[i]]$p[j] - results[[i]]$p[j-1]) / results[[i]]$p[j-1]) * 100
      degreeOfChange[j] <- (atan(results[[i]]$p[j] - results[[i]]$p[j - 1]) * (180/pi)) 
      
      if(nrow(results[[i]]) - j > window.size){
        varVector[j] <- var(results[[i]]$p[j:(j+window.size)]) 
        sdVector[j] <- sd(results[[i]]$p[j:(j+window.size)]) 
        meanVector[j] <- mean(results[[i]]$p[j:(j+window.size)])
        
      }else{
        varVector[j] <- var(results[[i]]$p[j:nrow(results[[i]])])
        sdVector[j] <- sd(results[[i]]$p[j:nrow(results[[i]])])
        meanVector[j] <- mean(results[[i]]$p[j:nrow(results[[i]])])
      }
    }
    
    
    for(j in 1:nGen){
      
      if((j+(window.size *2)) < nGen){
      
      tTestVar$tree[counter] <- i
      tTestVar$gen[counter] <- j
      tTestVar$`p-value`[counter] <- t.test(varVector[j:(j+window.size)], varVector[(j+window.size):(j+(window.size *2))])$p.value
      
      tTestSD$tree[counter] <- i
      tTestSD$gen[counter] <- j
      tTestSD$`p-value`[counter] <- t.test(sdVector[j:(j+window.size)], sdVector[(j+window.size):(j+(window.size *2))])$p.value
      
      tTestPecent$tree[counter] <- i
      tTestPecent$gen[counter] <- j
      tTestPecent$`p-value`[counter] <- t.test(pecentChange[j:(j+window.size)], pecentChange[(j+window.size):(j+(window.size *2))])$p.value
      }else{
        tTestVar$tree[counter] <- i
        tTestVar$gen[counter] <- j
        tTestVar$`p-value`[counter] <- t.test(varVector[j:(j+window.size)], varVector[(j+window.size):nGen])$p.value
        
        tTestSD$tree[counter] <- i
        tTestSD$gen[counter] <- j
        tTestSD$`p-value`[counter] <- t.test(sdVector[j:(j+window.size)], sdVector[(j+window.size):nGen])$p.value
        
        tTestPecent$tree[counter] <- i
        tTestPecent$gen[counter] <- j
        tTestPecent$`p-value`[counter] <- t.test(pecentChange[j:(j+window.size)], pecentChange[(j+window.size):nGen])$p.value
      }
      
      counter <- counter + 1
    }
    
    
    EndofBurn <- min(which(degreeOfChange[-1] < 0)) + 1
    if(EndofBurn * 2 >= length(results[[i]]$p * 0.5)){
      # runs that took long time to burn
      badRuns$iteration[i] <- i
      badRuns$generation <- EndofBurn
    }else{
      # runs that burn well
      goodRuns$iteration[i] <- i
      goodRuns$generation[i] <- EndofBurn
    }
  }
  # get the highest generation when the burn period had stopped
  maxBurn <- max(goodRuns$generation)
  # get suggested burnin based on all runs
  suggestedBurnin <-  (maxBurn / nGen) * 2
  # print this suggestion
  # print(paste("Suggested burn-in ", suggestedBurnin, "%", sep = ""))
  # get the output
  return(suggestedBurnin)
}

BurninCalc(dat = results)


getPostBurnin <- function(dat = NULL,
                          burn.in = NULL){
  
}




# # remove the first entry of XXXXXXXX because XXXXXXXX
# pecentChange <- pecentChange[-1]
# 
# for(j in 1:99){
#   atanVectorLikelihood[j] <-  abs(atan(results[[i]]$p[i + 1] - results[[i]]$p[j]) * (180/pi)) 
#   atanVectorVar[j] <- abs(atan(varVector[j + 1] - varVector[j]) * (180/pi))
#   atanVectorSD[j] <- abs(atan(sdVector[j + 1] - sdVector[j]) * (180/pi))
#   atanVectorMean[j] <- abs(atan(meanVector[j + 1] - meanVector[j]) * (180/pi))
#   
# }
# 
# plot(pecentChange, type = "l")
