# Terrence Sylvester
# pradakshanas@gmail.com
# 3rd February 2022

# read data
dat <- read.csv("../results/power.old.method.csv", as.is = T)

# plot
par(mfrow = c(3,3))
for(i in c(2,5,10)){
  # Fission
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  mtext(text = paste("Rates ratio", i),side = 3,adj = 0,
        cex = 0.7)
  mtext(text = "Fission",side = 3,adj = 1,
        cex = 0.7)
  # add power lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == i], y = dat$sigFission[dat$cond == 1 & dat$rr == i],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == i], y = dat$sigFission[dat$cond == 2 & dat$rr == i],type = "o",col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == i], y = dat$sigFission[dat$cond == 3 & dat$rr == i],type = "o",col = "#5e3c99", pch = 16)
  # add fasle positive lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == 1], y = dat$sigFission[dat$cond == 1 & dat$rr == 1],type = "o",col = "#e66101", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == 1], y = dat$sigFission[dat$cond == 2 & dat$rr == 1],type = "o",col = "#fdb863", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == 1], y = dat$sigFission[dat$cond == 3 & dat$rr == 1],type = "o",col = "#5e3c99", pch = 16, lty = 2)
  # Fusion
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  mtext(text = "Fusion",side = 3,adj = 1,
        cex = 0.7)
  # add power lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == i], y = dat$sigFusion[dat$cond == 1 & dat$rr == i],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == i], y = dat$sigFusion[dat$cond == 2 & dat$rr == i],type = "o",col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == i], y = dat$sigFusion[dat$cond == 3 & dat$rr == i],type = "o",col = "#5e3c99", pch = 16)
  # add false positive lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == 1], y = dat$sigFusion[dat$cond == 1 & dat$rr == 1],type = "o",col = "#e66101", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == 1], y = dat$sigFusion[dat$cond == 2 & dat$rr == 1],type = "o",col = "#fdb863", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == 1], y = dat$sigFusion[dat$cond == 3 & dat$rr == 1],type = "o",col = "#5e3c99", pch = 16, lty = 2)
  # Aneuploidy
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  mtext(text = "Aneuploidy",side = 3,adj = 1,
        cex = 0.7)
  # add power lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == i], y = dat$sigAneuploidy[dat$cond == 1 & dat$rr == i],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == i], y = dat$sigAneuploidy[dat$cond == 2 & dat$rr == i],type = "o",col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == i], y = dat$sigAneuploidy[dat$cond == 3 & dat$rr == i],type = "o",col = "#5e3c99", pch = 16)
  # add false positive lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == 1 & dat$rr == 1],type = "o",col = "#e66101", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == 2 & dat$rr == 1],type = "o",col = "#fdb863", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == 3 & dat$rr == 1],type = "o",col = "#5e3c99", pch = 16, lty = 2)
}


