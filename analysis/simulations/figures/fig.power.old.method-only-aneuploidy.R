# Terrence Sylvester
# pradakshanas@gmail.com
# 3rd February 2022

# read data
dat <- read.csv("../results/power.old.method.csv", as.is = T)
# plot
par(mfrow = c(1,3))
for(i in c(2,5,10)){
  # Aneuploidy
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  mtext(text = paste("Rates ratio", i),side = 3,adj = 0,
        cex = 0.7)
  # add power lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == i], y = dat$sigAneuploidy[dat$cond == 1 & dat$rr == i],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == i], y = dat$sigAneuploidy[dat$cond == 2 & dat$rr == i],type = "o",col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == i], y = dat$sigAneuploidy[dat$cond == 3 & dat$rr == i],type = "o",col = "#5e3c99", pch = 16)
  # add false positive lines
  lines(x = dat$tips[dat$cond == 1 & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == 1 & dat$rr == 1],type = "o",col = "#e66101", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 2 & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == 2 & dat$rr == 1],type = "o",col = "#fdb863", pch = 16, lty = 2)
  lines(x = dat$tips[dat$cond == 3 & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == 3 & dat$rr == 1],type = "o",col = "#5e3c99", pch = 16, lty = 2)
  
  if(i == 2){
    # add legend
    legend("topleft",inset=.02,
           legend=c("condition 1",
                    "condition 2",
                    "condition 3",
                    "rates ratio 1"),
           col=c("#e66101" ,"#fdb863" ,"#5e3c99", "black"),
           lty=c(1,1,1,2),
           cex=1,
           lwd = 2,title = "Legend")
  }
}


