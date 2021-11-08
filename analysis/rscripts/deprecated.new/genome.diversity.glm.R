# Terrence Sylvester
# 12 November 2020
# pradakshanas@gmail.com

# load in data
dat <- read.csv("../data/genetic.diversity/genetic.diversity.csv", as.is = T)
# get mean karyotype when there are multiple records
for(i in 1:nrow(dat)){
  dat$Karyotype[i] <- mean(as.numeric(unlist(strsplit(dat$Karyotype[i], split = "-"))))
}
dat$Karyotype <- as.numeric(dat$Karyotype)
dat$Karyotype[is.nan(dat$Karyotype)] <- NA
# make colours
col <- c("red", "blue")[as.factor(dat$LHP.breadth)]
# get the census pop size
pop.size <- dat$Abundance * dat$Distribution..km2.
# get a second dataset where the two outliers for karyotype are removed
dat.2 <- dat[!(dat$Karyotype %in% c(13,15)),]
col.2 <- c("red", "blue")[as.factor(dat.2$LHP.breadth)]
# make glm models
# single parameter models
m.1.01 <- glm(log(dat$X._4D) ~ dat$LHP.breadth)
m.1.02 <- glm(log(dat$X._4D) ~ dat$Forewing.length..mm.)
m.1.03 <- glm(log(dat$X._4D) ~ dat$Egg.volume..mm.3.)
m.1.04 <- glm(log(dat$X._4D) ~ dat$Abundance)
m.1.05 <- glm(log(dat$X._4D) ~ dat$Karyotype)
m.1.06 <- glm(log(dat.2$X._4D) ~ dat.2$Karyotype)
m.1.07 <- glm(log(dat$X._4D) ~ dat$Genome.size)
m.1.08 <- glm(log(dat$X._4D) ~ pop.size)
m.1.09 <- glm(log(dat$X._4D) ~ dat$Distribution..km2.)
# get summary for each model
m.1.summaries <- list(sum.m.1.01 <- summary(m.1.01),
                      sum.m.1.02 <- summary(m.1.02),
                      sum.m.1.03 <- summary(m.1.03),
                      sum.m.1.04 <- summary(m.1.04),
                      sum.m.1.05 <- summary(m.1.05),
                      sum.m.1.06 <- summary(m.1.06),
                      sum.m.1.07 <- summary(m.1.07),
                      sum.m.1.08 <- summary(m.1.08),
                      sum.m.1.09 <- summary(m.1.09))

# two parameter models where the second parameter is the larval host plant breadth
m.2.1.01 <- glm(log(dat$X._4D) ~ dat$Forewing.length..mm. + dat$LHP.breadth)
m.2.1.02 <- glm(log(dat$X._4D) ~ dat$Egg.volume..mm.3.+ dat$LHP.breadth)
m.2.1.03 <- glm(log(dat$X._4D) ~ dat$Abundance + dat$LHP.breadth)
m.2.1.04 <- glm(log(dat$X._4D) ~ dat$Karyotype + dat$LHP.breadth)
m.2.1.05 <- glm(log(dat.2$X._4D) ~ dat.2$Karyotype + dat.2$LHP.breadth)
m.2.1.06 <- glm(log(dat$X._4D) ~ dat$Genome.size + dat$LHP.breadth)
m.2.1.07 <- glm(log(dat$X._4D) ~ pop.size + dat$LHP.breadth)
m.2.1.08 <- glm(log(dat$X._4D) ~ dat$Distribution..km2. + dat$LHP.breadth)
# get summary for each model
m.2.1.summaries <- list(sum.m.2.1.01 <- summary(m.2.1.01),
                        sum.m.2.1.02 <- summary(m.2.1.02),
                        sum.m.2.1.03 <- summary(m.2.1.03),
                        sum.m.2.1.04 <- summary(m.2.1.04),
                        sum.m.2.1.05 <- summary(m.2.1.05),
                        sum.m.2.1.06 <- summary(m.2.1.06),
                        sum.m.2.1.07 <- summary(m.2.1.07),
                        sum.m.2.1.08 <- summary(m.2.1.08))

# two parameter models which look only at the individual impact interaction of each parameter
m.2.2.01 <- glm(log(dat$X._4D) ~ dat$Forewing.length..mm. : dat$LHP.breadth)
m.2.2.02 <- glm(log(dat$X._4D) ~ dat$Egg.volume..mm.3.: dat$LHP.breadth)
m.2.2.03 <- glm(log(dat$X._4D) ~ dat$Abundance : dat$LHP.breadth)
m.2.2.04 <- glm(log(dat$X._4D) ~ dat$Karyotype : dat$LHP.breadth)
m.2.2.05 <- glm(log(dat.2$X._4D) ~ dat.2$Karyotype : dat.2$LHP.breadth)
m.2.2.06 <- glm(log(dat$X._4D) ~ dat$Genome.size : dat$LHP.breadth)
m.2.2.07 <- glm(log(dat$X._4D) ~ pop.size : dat$LHP.breadth)
m.2.2.08 <- glm(log(dat$X._4D) ~ dat$Distribution..km2. : dat$LHP.breadth)
# get summary for each model
m.2.2.summaries <- list(sum.m.2.2.01 <- summary(m.2.2.01),
                        sum.m.2.2.02 <- summary(m.2.2.02),
                        sum.m.2.2.03 <- summary(m.2.2.03),
                        sum.m.2.2.04 <- summary(m.2.2.04),
                        sum.m.2.2.05 <- summary(m.2.2.05),
                        sum.m.2.2.06 <- summary(m.2.2.06),
                        sum.m.2.2.07 <- summary(m.2.2.07),
                        sum.m.2.2.08 <- summary(m.2.2.08))

# two parameter models which look at the individual impact as well as
# the interaction of the parameters
m.2.3.01 <- glm(log(dat$X._4D) ~ dat$Forewing.length..mm. * dat$LHP.breadth)
m.2.3.02 <- glm(log(dat$X._4D) ~ dat$Egg.volume..mm.3.* dat$LHP.breadth)
m.2.3.03 <- glm(log(dat$X._4D) ~ dat$Abundance * dat$LHP.breadth)
m.2.3.04 <- glm(log(dat$X._4D) ~ dat$Karyotype * dat$LHP.breadth)
m.2.3.05 <- glm(log(dat.2$X._4D) ~ dat.2$Karyotype * dat.2$LHP.breadth)
m.2.3.06 <- glm(log(dat$X._4D) ~ dat$Genome.size * dat$LHP.breadth)
m.2.3.07 <- glm(log(dat$X._4D) ~ pop.size * dat$LHP.breadth)
m.2.3.08 <- glm(log(dat$X._4D) ~ dat$Distribution..km2. * dat$LHP.breadth)
# get summary for each model
m.2.3.summaries <- list(sum.m.2.3.01 <- summary(m.2.3.01),
                        sum.m.2.3.02 <- summary(m.2.3.02),
                        sum.m.2.3.03 <- summary(m.2.3.03),
                        sum.m.2.3.04 <- summary(m.2.3.04),
                        sum.m.2.3.05 <- summary(m.2.3.05),
                        sum.m.2.3.06 <- summary(m.2.3.06),
                        sum.m.2.3.07 <- summary(m.2.3.07),
                        sum.m.2.3.08 <- summary(m.2.3.08))


# get the slop and the p-value on a data table
m1.tab <- c()
r.name <- c()
for(i in 1:length(m.1.summaries)){
  m1.tab <- rbind(m1.tab, m.1.summaries[[i]]$coefficients[-1,])
  r.name[i] <- rownames(m.1.summaries[[i]]$coefficients)[-1]
}
rownames(m1.tab) <- r.name

# get the slop and the p-value on a data table
m2.1tab <- c()
for(i in 1:length(m.2.1.summaries)){
  m2.1tab <- rbind(m2.1tab, m.2.1.summaries[[i]]$coefficients[-1,])
}

# get the slop and the p-value on a data table
m2.2tab <- c()
for(i in 1:length(m.2.2.summaries)){
  m2.2tab <- rbind(m2.2tab, m.2.2.summaries[[i]]$coefficients[-1,])
}

# get the slop and the p-value on a data table
m2.3tab <- c()
for(i in 1:length(m.2.3.summaries)){
  m2.3tab <- rbind(m2.3tab, m.2.3.summaries[[i]]$coefficients[-1,])
}

write.csv(m1.tab, "../tables/genome.diversity.m1.tab.csv")
write.csv(m2.1tab, "../tables/genome.diversity.m2.1.tab.csv")
write.csv(m2.2tab, "../tables/genome.diversity.m2.2.tab.csv")
write.csv(m2.3tab, "../tables/genome.diversity.m2.3.tab.csv")

# t.tests
t.tests <- list(t.test01 <- t.test(dat$LHP.species ~ dat$LHP.breadth),
                t.test02 <- t.test(dat$Forewing.length..mm. ~ dat$LHP.breadth),
                t.test03 <- t.test(dat$Egg.volume..mm.3. ~ dat$LHP.breadth),
                t.test04 <- t.test(dat$Distribution..km2. ~ dat$LHP.breadth),
                t.test05 <- t.test(dat$Abundance ~ dat$LHP.breadth),
                t.test06 <- t.test(dat$Karyotype ~ dat$LHP.breadth),
                t.test07 <- t.test(dat.2$Karyotype ~ dat.2$LHP.breadth),
                t.test08 <- t.test(dat$Genome.size ~ dat$LHP.breadth),
                t.test09 <- t.test(dat$X._4D ~ dat$LHP.breadth),
                t.test10 <- t.test(dat$X._0D ~ dat$LHP.breadth))

t.test.results <- as.data.frame(matrix(data = NA,
                                       nrow = 10,
                                       ncol = 4))
colnames(t.test.results) <- c("parameter",
                              "Mean in generalists", 
                              "Mean in specialists",
                              "p-value")

t.test.results[,1] <- c("Number of Host species",
                        "Body size (mm)",
                        "Egg volume (mm3)",
                        "Range (km2)",
                        "Density",
                        "Karyotype",
                        "Karyotype - outliers removed",
                        "Genome size",
                        "Neutral genome diversity",
                        "Nonsynonimous diversity")
for( i in 1:length(t.tests)){
  t.test.results[i,2] <-  round(t.tests[[i]]$estimate[2], 3)
  t.test.results[i,3] <-  round(t.tests[[i]]$estimate[1], 3)
  t.test.results[i,4] <-  round(t.tests[[i]]$p.value, 3)
}
write.csv(t.test.results, "../tables/t.test.results.csv",row.names = F)
