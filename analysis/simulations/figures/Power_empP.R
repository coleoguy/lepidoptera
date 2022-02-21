# Terrence Sylvester
# get file names and make data table to store results
files <- dir("../results/empP-MCMC/tree.67/")
dat <- as.data.frame(matrix(data = NA, nrow = length(files),ncol = 4))
colnames(dat) <- c("run", "fission", "fusion","aneuploidy")
# get results
for(i in 1:length(files)){
  # load results
  load(paste("../results/empP-MCMC/tree.67/", files[i], sep = ""))
  dat$run[i] <- files[i]
  dat$fission[i] <- empP$empiricalP$EmpPvalueFission
  dat$fusion[i] <- empP$empiricalP$EmpPvalueFusion
  dat$aneuploidy[i] <- empP$empiricalP$EmpPvalueAneuploidy
  # remove unwanted results
  rm(empP)
}
# make new columns to indicate condition rr and tips
dat$tips <- rep(rep(c(100,200,50), each = 4), 3)
dat$rr <- rep(c(1,10,2,5),9)
dat$cond <- rep(c(1,2,3), each = 12)
# sort by tips
dat <- dat[order(dat$tips),]
# plot
par(mfrow = c(3,3))
for(i in 1:3){
  # fission
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "P value",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$fission[dat$cond == i & dat$rr == 1],type = "o",col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$fission[dat$cond == i & dat$rr == 2],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$fission[dat$cond == i & dat$rr == 5],type = "o",col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$fission[dat$cond == i & dat$rr == 10],type = "o",col = "#5e3c99", pch = 16)
  mtext(text = paste("Condition", i),side = 3,adj = 0,
        cex = 0.7)
  mtext(text = "Fission",side = 3,adj = 1,
        cex = 0.7)
  abline(h = 0.05, col = "red", lty = 2)
  # fusion
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "P value",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$fusion[dat$cond == i & dat$rr == 1],type = "o", col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$fusion[dat$cond == i & dat$rr == 2],type = "o", col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$fusion[dat$cond == i & dat$rr == 5],type = "o", col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$fusion[dat$cond == i & dat$rr == 10],type = "o", col = "#5e3c99", pch = 16)
  mtext(text = "Fusion",side = 3,adj = 1,
        cex = 0.7)
  abline(h = 0.05, col = "red", lty = 2)
  # Aneuploidy
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "P value",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$aneuploidy[dat$cond == i & dat$rr == 1],type = "o",col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$aneuploidy[dat$cond == i & dat$rr == 2],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$aneuploidy[dat$cond == i & dat$rr == 5],type = "o", col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$aneuploidy[dat$cond == i & dat$rr == 10],type = "o",col = "#5e3c99", pch = 16)
  mtext(text = "Aneuploidy",side = 3,adj = 1,
        cex = 0.7)
  abline(h = 0.05, col = "red", lty = 2)
  # # legend
  # if(i == 3){
  # # add legend
  # legend("topright",inset=.02,
  #        legend=c("1",
  #                 "2",
  #                 "5",
  #                 "10"),
  #        col=c("gray","#e66101" ,"#fdb863" ,"#5e3c99"),
  #        lty=1,
  #        cex=1,
  #        lwd = 2,title = "Rates ratio")
  # }
}

