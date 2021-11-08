# Terrence Sylvester
# 26 October 2020
# pradakshanas@gmail.com

# load libraries
library(doMC)

# make the function to calculate the genome diversity
processFile = function(filepath) {
  # define parameters
  # total length of the genome excluding N's
  tot.length <- 0
  # total ambiguous bases
  amb <- 0
  # connection to open
  # lets make it so that each file takes 3 min to run
  # start.time <- unclass(as.POSIXlt(Sys.time()))$min
  con = file(filepath, "r")
  while ( TRUE ) {
    # read in 1 line at a time
    y = readLines(con, n = 1)
    # if the read connection is length 0 then end this loop
    if ( length(y) == 0 ) {
      break
    }
    # if (unclass(as.POSIXlt(Sys.time()))$min - start.time == 1) {
    #   break
    # }
    # here we ignore the header of fasta sequence
    if((gregexpr(pattern = ">", text =  toupper(y))[[1]]) == -1){
      # this line will get the number of bases per line ignoring N's
      if(sample(gregexpr(pattern = "N", text =  toupper(y))[[1]],1) == -1){
        nonNlength = nchar(y)
        tot.length <- tot.length + nonNlength
      }else{
        nonNlength <- (nchar(y) - length(gregexpr(pattern = "N", text =  toupper(y))[[1]]))
        tot.length <- tot.length + nonNlength
      }
      # length of ambiguity letters
      #Y
      if(sample(gregexpr(pattern = "Y", text =  toupper(y))[[1]],1) != -1){
        length.y <- length(gregexpr(pattern = "Y", text =  toupper(y))[[1]])
      }else{
        length.y <- 0
      }
      #R
      if(sample(gregexpr(pattern = "R", text =  toupper(y))[[1]],1) != -1){
        length.r <- length(gregexpr(pattern = "R", text =  toupper(y))[[1]])
      }else{
        length.r <- 0
      }
      #W
      if(sample(gregexpr(pattern = "W", text =  toupper(y))[[1]],1) != -1){
        length.w <- length(gregexpr(pattern = "W", text =  toupper(y))[[1]])
      }else{
        length.w <- 0
      }
      #S
      if(sample(gregexpr(pattern = "S", text =  toupper(y))[[1]],1) != -1){
        length.s <- length(gregexpr(pattern = "S", text =  toupper(y))[[1]])
      }else{
        length.s <- 0
      }
      #K
      if(sample(gregexpr(pattern = "K", text =  toupper(y))[[1]],1) != -1){
        length.k <- length(gregexpr(pattern = "K", text =  toupper(y))[[1]])
      }else{
        length.k <- 0
      }
      #M
      if(sample(gregexpr(pattern = "M", text =  toupper(y))[[1]],1) != -1){
        length.m <- length(gregexpr(pattern = "M", text =  toupper(y))[[1]])
      }else{
        length.m <- 0
      }
      # Here we calculate the sum of ambiguouse bases
      A <-  sum(length.y, length.r, length.w, length.s, length.k, length.m)
      amb <- amb + A
    }
  }
  # we calculate the divergence of the given genome
  distance <- amb / tot.length
  return(distance)
  close(con)
}

files.specialists <- dir("../Genomes/FASTA-files/specialists/")
files.generalists <- dir("../Genomes/FASTA-files/generalists/")

dat.specialists <- matrix(data =NA, nrow = length(files.specialists), ncol = 2)
dat.generalists <- matrix(data =NA, nrow = length(files.generalists), ncol = 2)

colnames(dat.generalists) <- colnames(dat.specialists) <- c("Name", "genomeDistance")

# process specialist genomes
registerDoMC(3)
specialists <- foreach(i  = 1:24)%dopar%{
  spec <- (processFile(paste("../Genomes/FASTA-files/specialists/", files.specialists[i], sep = "")))
  names(spec) <- files.specialists[i]
  spec
  
  # print(files.specialists[i])
  # x <- Sys.time()
  # dat.specialists[i,1] <- files.specialists[i]
  # dat.specialists[i,2] <- (processFile(paste("../Genomes/FASTA-files/specialists/", files.specialists[i], sep = "")))
  # print(paste("time taken to process:"))
  # print(Sys.time() - x)
}

for(i in 1:length(files.specialists)){
  dat.specialists[i,1] <- names(specialists[[i]])
  dat.specialists[i,2] <- specialists[[i]]
}

# process generalist genomes
registerDoMC(3)
generalists <- foreach(i = 1:14)%dopar%{
  gen <- (processFile(paste("../Genomes/FASTA-files/generalists/", files.generalists[i], sep = "")))
  names(gen) <- files.generalists[i]
  gen
  
  # x <- Sys.time()
  # dat.specialists[i,1] <- files.generalists[i]
  # dat.specialists[i,2] <- (processFile(paste("../Genomes/FASTA-files/generalists/", files.generalists[i], sep = "")))
  # print(paste("time taken to process:"))
  # print(Sys.time() - x)
}

for(i in 1:length(files.generalists)){
  dat.generalists[i,1] <- names(generalists[[i]])
  dat.generalists[i,2] <- generalists[[i]]
}

# write the csv files
write.csv(x = dat.generalists,file = "../tables/diversity.generalists.csv",row.names = F)
write.csv(x = dat.specialists,file = "../tables/diversity.specialists.csv",row.names = F)

save.image("../genome.diversity.RData")

