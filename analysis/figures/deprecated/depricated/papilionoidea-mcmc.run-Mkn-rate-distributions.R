# load data
dat <- read.csv("../results/papilionoidea-mcmc.run-Mkn-post-burn-in.csv",
                as.is = T)

# specialists
plot(x = NULL,
     y = NULL,
     ylim  = c(0,70),
     xlim = c(-0,.22),
     xlab = "Rate",
     ylab = "probability density",
     main = "Distribution of rates of fusion and fission in specialist species")

polygon(density(dat$asc2),
        col = rgb(1,0,0,.2))

polygon(density(dat$desc2),
        col = rgb(0,0,1,.2))

polygon(density(dat$pol2),
        col = rgb(0,1,0,.2))

text(x = rep(0.05,3),
     y = c(60,55,50),
     labels = c("Rates of fission", 
                "Rates of fusion",
                "Rates of polyploidy"),
     col = c(rgb(1,0,0,.7),
             rgb(0,0,1,.7),
             rgb(0,1,0,.7)),
     pos = 4)

# generalists
plot(x = NULL,
     y = NULL,
     ylim  = c(0,50),
     xlim = c(-.01,.22),
     xlab = "Rate",
     ylab = "probability density",
     main = "Distribution of rates of fusion and fission in generalist species")

polygon(density(dat$asc1),
        col = rgb(1,0,0,.2))

polygon(density(dat$desc1),
        col = rgb(0,0,1,.2))

polygon(density(dat$pol1),
        col = rgb(0,1,0,.2))

text(x = rep(0.15,3),
     y = c(50,45,40),
     labels = c("Rates of fission", 
                "Rates of fusion",
                "Rates of polyploidy"),
     col = c(rgb(1,0,0,.7),
             rgb(0,0,1,.7),
             rgb(0,1,0,.7)),
     pos = 4)

# transision rates
plot(x = NULL,
     y = NULL,
     ylim  = c(0,20),
     xlim = c(0,.22),
     xlab = "Rate",
     ylab = "probability density",
     main = "Distribution of transision rates between generalis and specialist states")

polygon(density(dat$tran12),
        col = rgb(1,0,0,.2))

polygon(density(dat$tran21),
        col = rgb(0,0,1,.2))

text(x = rep(0.12,3),
     y = c(20,19),
     labels = c("transision from generalist to specialist", 
                "transision from specialist to generalist"),
     col = c(rgb(1,0,0,.7),
             rgb(0,0,1,.7)),
     pos = 4)
