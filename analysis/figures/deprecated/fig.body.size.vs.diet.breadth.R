# Terrence Sylvester
# pradakshanas@gmail.com
# 24 November 2020

# read in data
bsz <- read.csv("../data/body.size/body.size.csv", as.is = T)
bsz$speciesNames <- paste(bsz$Genus, bsz$Species, sep = "_")
# remove empty rows
bsz <- bsz[!is.na(bsz$wingspan.mm.),]
bsz <- bsz[!is.na(bsz$hostType),]
# remove duplicated species
bsz <- bsz[!duplicated(bsz$speciesNames),]
# do a t.test
t.test <- t.test(bsz$wingspan.mm. ~ bsz$hostType)
# plot data
# make an empty plot
plot(x = NULL,
     y = NULL,
     xlim = c(0.5,2.5),
     ylim = c(0,140),
     axes = F,
     xlab = "Larval host breadth",
     ylab = "Wingspan (mm)")
#plot the points
points(x = jitter(bsz$hostType + 1, factor = 2),
       y = bsz$wingspan.mm.,
       pch = 16,
       col = rgb(0,0,0,.2))
# plot a boxplot
boxplot(bsz$wingspan.mm. ~ bsz$hostType,
        outline = F,
        ylim = c(0,140),
        add = T, col = rgb(1,1,1,.7),
        axes = F)
# plot axiz
axis(side = 1, at = c(1,2), labels = c("Generalist", "Specialist"))
axis(side = 2)
box(bty = "l")
# plot the p-value
text(x = 1.7,
     y = 0,
     labels = expression(paste("p-value ", 2 %*% 10^-5)),
     pos = 4)
