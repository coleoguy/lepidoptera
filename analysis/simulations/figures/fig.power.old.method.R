# Terrence Sylvester
# pradakshanas@gmail.com
# 3rd February 2022

# read data
dat <- read.csv("../results/power.old.method.csv", as.is = T)

# plot
par(mfrow = c(3,3))
for(i in 1:3){
  # Fission
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$sigFission[dat$cond == i & dat$rr == 1],type = "o",col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$sigFission[dat$cond == i & dat$rr == 2],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$sigFission[dat$cond == i & dat$rr == 5],type = "o",col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$sigFission[dat$cond == i & dat$rr == 10],type = "o",col = "#5e3c99", pch = 16)
  mtext(text = paste("Condition", i),side = 3,adj = 0,
        cex = 0.7)
  mtext(text = "Fission",side = 3,adj = 1,
        cex = 0.7)
  # if(i == 1){
  # # # add legend
  # legend("topleft",inset=.02,
  #        legend=c("1",
  #                 "2",
  #                 "5",
  #                 "10"),
  #        col=c("gray","#e66101" ,"#fdb863" ,"#5e3c99"),
  #        lty=1,
  #        cex=1,
  #        lwd = 2,title = "Rates ratio")
  # }
  # Fusuin
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$sigFusion[dat$cond == i & dat$rr == 1],type = "o", col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$sigFusion[dat$cond == i & dat$rr == 2],type = "o", col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$sigFusion[dat$cond == i & dat$rr == 5],type = "o", col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$sigFusion[dat$cond == i & dat$rr == 10],type = "o", col = "#5e3c99", pch = 16)
  mtext(text = "Fusion",side = 3,adj = 1,
        cex = 0.7)
  # Aneuploidy
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 1],type = "o",col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 2],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 5],type = "o", col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 10],type = "o",col = "#5e3c99", pch = 16)
  mtext(text = "Aneuploidy",side = 3,adj = 1,
        cex = 0.7)
  
}


