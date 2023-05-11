library(phytools)
library(l1ou)
library(PhylogeneticEM)
library(dplyr)
library(ratematrix)
source("~/Google.Drive/R.Analyses/Convenient Scripts/make.OUwie.input_Script.R")
library(factoextra)
library(cluster)
library(NbClust)
library(Rmisc)
library(corrplot)
library(patchwork)


tree <- read.nexus("~/Google.Drive/R.Analyses/Microhylids/T428_50loci.tre")
morph.data <- read.csv("~/Google.Drive/R.Analyses/Microhylids/T428_RAWdata.csv", header=T)
    morph.data <- morph.data[,c("Name_in_Tree", "SVL", "HLL", "HW", "THIRD")]
        morph.data <- morph.data[complete.cases(morph.data),]
            #morph.data[2:5] <- log(morph.data[2:5])

# create a tibble to get the species means
sp.means <- morph.data %>%
  group_by(Name_in_Tree) %>%
  summarise_at(vars(SVL:THIRD), mean)
sp.means <- sp.means[-c(1:4,27,29,34),] # drop some samples, look View(sp.means)
# make sure to drop riparius, kaindensis, nubicola_2!
write.csv(sp.means, row.names=FALSE,
          file="/Users/Ian/Google.Drive/R.Analyses/Microhylids/T428_spMEANS.csv")
            
tree <- drop.tip(tree, tip=setdiff(tree$tip.label, sp.means$Name_in_Tree))
 # write.nexus(tree, file="/Users/Ian/Google.Drive/R.Analyses/Microhylids/T428_50loci_TRIMMED.tre")

data <- as.data.frame(sp.means[2:5]); rownames(data) <- sp.means$Name_in_Tree
  tree.pca <- phyl.pca(tree, data, method="lambda")
  plot(tree.pca)
  biplot(tree.pca)
  data <- cbind(data, tree.pca$S)
    write.csv(data, file="/Users/Ian/Google.Drive/R.Analyses/Microhylids/T428_PCAscores.csv")
  #normal.pca <- prcomp(data[1:5])
    #normal.pcs <- normal.pca$x

library(factoextra)
    res.pca <- prcomp(data, scale = TRUE)
    testo <- fviz_eig(res.pca)
    fviz_pca_ind(res.pca,
                 col.ind = "cos2", # Color by the quality of representation
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE)     # Avoid text overlapping
    fviz_pca_var(res.pca,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE)     # Avoid text overlapping
    get_eigenvalue(res.pca)
    
    
pcs <- data[6:8] 
  l1ou.data <- adjust_data(test.tree, test.data[,1])
shift.fit <- estimate_shift_configuration(l1ou.data$tree, l1ou.data$Y,
                                          nCores=8, criterion="AICc", quietly=F)
plot(shift.fit) # if you get 'figure margins' error, do: par(mar=c(1,1,1,1))

conv.fit <- estimate_convergent_regimes(shift.fit)
plot(conv.fit)

shift.config <- c(27, 23, 10)
ou.fit <- fit_OU(l1ou.data$tree, l1ou.data$Y, shift.config, criterion="pBIC")
shift.config <- c(28, 29, 25, 24, 10)
ou.fit <- fit_OU(l1ou.data$tree, l1ou.data$Y, shift.config, criterion="AICc")
shift.config <- c(10)
ou.fit <- fit_OU(l1ou.data$tree, l1ou.data$Y, shift.config, criterion="pBIC")


  
traits <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Microhylids/T428_spMEANS.csv")
    traits <- traits[-16,]
  
trait.in <- make.OUwie.input(traits, traits$Habitat, c(1,2), traits$Name_in_Tree,
                 tree, traits$SVL) 
trait.in <- make.OUwie.input(data, data$Habitat, c(1,2), rownames(data),
                             tree, data$PC1)

bm.fit <-    OUwie(trait.in[[4]], trait.in[[3]], model="BM1",  root.station=T, simmap.tree=T)
ou.fit <-    OUwie(trait.in[[4]], trait.in[[3]], model="OU1",  root.station=T, simmap.tree=T)
oum.fit <-   OUwie(trait.in[[4]], trait.in[[3]], model="OUM",  root.station=T, simmap.tree=T)
ouma.fit <-  OUwie(trait.in[[4]], trait.in[[3]], model="OUMA", root.station=T, simmap.tree=T)
oumv.fit <-  OUwie(trait.in[[4]], trait.in[[3]], model="OUMV", root.station=T, simmap.tree=T)
oumva.fit <- OUwie(trait.in[[4]], trait.in[[3]], model="OUMVA", root.station=T, simmap.tree=T)

# CALCULATE DISCRIMINABILITY RATIO, OUTLINE DISCRETE HYPOTHESES TO TEST

#test.data <- traits[c(2,3)]; rownames(test.data) <- traits$Name_in_Tree
test.data <- data[5];
test.data <- t(test.data) # transpose the data (taxa as columns, traits as rows)
test.data <- test.data[ , match(tree$tip.label, colnames(test.data))] # make the order of the data match the tree tip labels
frogs <- NULL; frogs$phy <- tree; frogs$dat <- test.data

#   res <- PhyloEM(phylo=tree, 
#                 Y_data=test.data,      # read in the trait data         
#                 process="scOU",         # identify the process to analyse    
#                 random.root=T,         #    
#                 stationary.root=T,
#                 nbr_alpha=2,
#                 K_max=20,               # set a maximum limit on the number of shifts to search
#                 check.tips.names=F,     # check to make sure names/trait data match          
#                 parallel_alpha=T,       # we want to parallelize the analysis        
#                 Ncores=8)               # with how many cores?
#   independent=F)          # if using multiple traits, are they independent?   
#   plot(res, show.tip.label=T, label_cex=0.1)

res <- PhyloEM(Y_data = frogs$dat, # read in the trait data
               phylo = frogs$phy,
               process = "BM",
               random.root = TRUE,
               stationary.root = TRUE,
               K_max = 10,
               #nbr_alpha = 20,
               parallel_alpha = TRUE,
               Ncores = 8)
plot(res, show.tip.label = T, label_cex=0.4)
plot(shifts_to_simmap(tree, c(27, 23, 10)))


data(monkeys)
res <- PhyloEM(Y_data = monkeys$dat,        ## data
               phylo = monkeys$phy,         ## phylogeny
               process = "scOU",            ## scalar OU
               random.root = TRUE,          ## root is stationary
               stationary.root = TRUE,
               K_max = 10,                  ## maximal number of shifts
               nbr_alpha = 4,               ## number of alpha values
               parallel_alpha = TRUE,       ## parallelize on alpha values
               Ncores = 8)
plot(res)




# LOG-SHAPE RATIOS

# this is basically calculating size as the geometric mean of all the morphometric measurements of each specimen. 
# Thsi is a great measure of size because it incorporates all the variables. 
# Using only SVL might give you a misleading measure of size (a tree snake being larger than a monitor)
# Then each measurement is divided by this "size" to give you a shape ratio
# The you log this no normalize it

# grabbed the data from down below (phylomorphospace plot)
data <- exp(pdata[,2:5]) #! What is this? Why is there "EXP"?????

#computing the geometric mean for obtaining size
#size <- apply(data[1:4], 1, prod)^(1/ncol(data[1:4]))
size <- apply(data, 1, prod)^(1/ncol(data))

#computing the log shape ratios
#LSR <- log(data[2:4]/size) # use if data isn't already log transformed
#LSR <- data[1:4]/size
LSR <- data/size
LSR$size <- size

# See how the autocorrelation between your variables due to size changes before and after correction!
#pairs(data[1:4])
pairs(data)
pairs(LSR)

M <- cor(LSR[,1:5])
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"), method="ellipse")
corrplot.mixed(M, upper="ellipse")

Q <- cor(data)
corrplot(Q, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
corrplot.mixed(Q, order="hclust", upper="ellipse")

# color the species
LSR$color <- "#E2E2E2"
LSR[rownames(LSR)=="Cophixalus_pakayakulangun","color"] <- "#D7191C"
LSR[rownames(LSR)=="Cophixalus_kulakula_1", "color"] <- "#F99898"
LSR[rownames(LSR)=="Cophixalus_zweifeli_1", "color"] <- "#2C7BB6"
LSR[rownames(LSR)=="Cophixalus_petrophilus_1", "color"] <- "#ABD9E9"
LSR[rownames(LSR)=="Cophixalus_saxatilis", "color"] <- "#F9F278"

# plot variables
size_SVL <-  ggplot(LSR) + geom_point(aes(x=size, y=SVL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
size_HLL <-  ggplot(LSR) + geom_point(aes(x=size, y=HLL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
size_HW <-   ggplot(LSR) + geom_point(aes(x=size, y=HW),    color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
size_THR <-  ggplot(LSR) + geom_point(aes(x=size, y=THIRD), color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
size_size <- ggplot(LSR) + geom_point(aes(x=size, y=size),  color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
# plot all size figures together
sizes <- size_SVL / size_HLL / size_HW / size_THR / size_size


SVL_SVL <-  ggplot(LSR) + geom_point(aes(x=SVL, y=SVL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
SVL_HLL <-  ggplot(LSR) + geom_point(aes(x=SVL, y=HLL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
SVL_HW <-   ggplot(LSR) + geom_point(aes(x=SVL, y=HW),    color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
SVL_THR <-  ggplot(LSR) + geom_point(aes(x=SVL, y=THIRD), color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
SVL_size <- ggplot(LSR) + geom_point(aes(x=SVL, y=size),  color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
# plot all SVL figures together
SVL <- SVL_SVL / SVL_HLL / SVL_HW / SVL_THR / SVL_size

HLL_SVL <-  ggplot(LSR) + geom_point(aes(x=HLL, y=SVL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HLL_HLL <-  ggplot(LSR) + geom_point(aes(x=HLL, y=HLL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HLL_HW <-   ggplot(LSR) + geom_point(aes(x=HLL, y=HW),    color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HLL_THR <-  ggplot(LSR) + geom_point(aes(x=HLL, y=THIRD), color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HLL_size <- ggplot(LSR) + geom_point(aes(x=HLL, y=size),  color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
#plot all HLL figures together
HLL <- HLL_SVL / HLL_HLL /  HLL_HW / HLL_THR / HLL_size

HW_SVL <-  ggplot(LSR) + geom_point(aes(x=HW, y=SVL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HW_HLL <-  ggplot(LSR) + geom_point(aes(x=HW, y=HLL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HW_HW <-   ggplot(LSR) + geom_point(aes(x=HW, y=HW),    color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HW_THR <-  ggplot(LSR) + geom_point(aes(x=HW, y=THIRD), color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
HW_size <- ggplot(LSR) + geom_point(aes(x=HW, y=size),  color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
# plot all HW figures together
HW <- HW_SVL / HW_HLL /  HW_HW / HW_THR / HW_size

THR_SVL <-  ggplot(LSR) + geom_point(aes(x=THIRD, y=SVL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
THR_HLL <-  ggplot(LSR) + geom_point(aes(x=THIRD, y=HLL),   color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
THR_HW <-   ggplot(LSR) + geom_point(aes(x=THIRD, y=HW),    color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
THR_THR <-  ggplot(LSR) + geom_point(aes(x=THIRD, y=THIRD), color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
THR_size <- ggplot(LSR) + geom_point(aes(x=THIRD, y=size),  color="black", fill=LSR$color, size=3, shape=21) + theme_bw()
# plot all THIRD figures together
THR <- THR_SVL / THR_HLL /  THR_HW / THR_THR / THR_size

# combine all plots into a pairwise comparison figure
SVL | HLL | HW | THR | size


# Now do the same for the non size-corrected data
# color the species
data$color <- "#E2E2E2"
data[rownames(data)=="Cophixalus_pakayakulangun","color"] <- "#D7191C"
data[rownames(data)=="Cophixalus_kulakula_1", "color"] <- "#F99898"
data[rownames(data)=="Cophixalus_zweifeli_1", "color"] <- "#2C7BB6"
data[rownames(data)=="Cophixalus_petrophilus_1", "color"] <- "#ABD9E9"
data[rownames(data)=="Cophixalus_saxatilis", "color"] <- "#F9F278"

# plot variables
raw_SVL_SVL <-  ggplot(data) + geom_point(aes(x=SVL, y=SVL),   color="black", fill=data$color, size=3, shape=21) + theme_bw()
raw_SVL_HLL <-  ggplot(data) + geom_point(aes(x=SVL, y=HLL),   color="black", fill=data$color, size=3, shape=21) + theme_bw()
raw_SVL_HW <-   ggplot(data) + geom_point(aes(x=SVL, y=HW),    color="black", fill=data$color, size=3, shape=21) + theme_bw()
raw_SVL_THR <-  ggplot(data) + geom_point(aes(x=SVL, y=THIRD), color="black", fill=data$color, size=3, shape=21) + theme_bw()
raw_SVL_size <- ggplot(data) + geom_point(aes(x=SVL, y=size),  color="black", fill=data$color, size=3, shape=21) + theme_bw()
# plot all SVL figures together
raw_SVL <- raw_SVL_SVL / raw_SVL_HLL / raw_SVL_HW / raw_SVL_THR / raw_SVL_size


# Test for isometry:

#get geomorph
#devtools::install_github("geomorphR/geomorph", ref="Stable") # for the stable version
library(geomorph)


# allo <- procD.lm(LSR[1:4] ~ LSR$size, iter=999, RRPP=T, plot=T)
# plotAllometry(allo, LSR$size, logsz=F, method="PredLine")

allo.SVL <- procD.pgls(SVL ~ size, data=LSR, phy=tree, iter=999, print.progress = FALSE, plot=T); summary(allo.SVL)
    plotAllometry(allo.SVL, LSR$size, logsz=F, method="PredLine")
allo.HLL <- procD.pgls(HLL ~ size, data=LSR, phy=tree, iter=999, print.progress = FALSE, plot=T); summary(allo.HLL)
    plotAllometry(allo.HLL, LSR$size, logsz=F, method="PredLine")
allo.HW <-  procD.pgls(HW ~ size, data=LSR, phy=tree, iter=999, print.progress = FALSE, plot=T); summary(allo.HW)
    plotAllometry(allo.HW, LSR$size, logsz=F, method="PredLine")
allo.THR <- procD.pgls(THIRD ~ size, data=LSR, phy=tree, iter=999, print.progress = FALSE, plot=T); summary(allo.THR)
    plotAllometry(allo.THR, LSR$size, logsz=F, method="PredLine")

allo.all <- procD.pgls(LSR[1:4] ~ size, data=LSR, phy=tree, iter=999, print.progress = FALSE, plot=T); summary(allo.all)
    plotAllometry(allo.all, LSR$size, logsz=F, method="PredLine")

allo.all.raw <- procD.pgls(raw.data[1:4] ~ LSR$size, phy=tree, iter=999); summary(allo.all.raw)
    plotAllometry(allo.all.raw, LSR$size, logsz=F, method="PredLine")

allo.test <- procD.pgls(test.LSR~test.size, test.tree)
#procD.allometry(LSR~size)

small.tree <- extract.clade(tree, 31); plot(small.tree)
small.LSR <- LSR[which(rownames(LSR)%in%small.tree$tip.label),]

allo.small <- procD.pgls(small.LSR[1:4] ~ size, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(allo.small)
    plotAllometry(allo.small, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
small.SVL <- procD.pgls(SVL ~ size, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(small.SVL)
    plotAllometry(small.SVL, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
small.HLL <- procD.pgls(HLL ~ size, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(small.HLL)
    plotAllometry(small.HLL, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
small.HW <-  procD.pgls(HW ~ size, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(small.HW)
    plotAllometry(small.HW, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
small.THR <- procD.pgls(THIRD ~ size, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(small.THR)
    plotAllometry(small.THR, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))    
    
small.LSR$ecology <- "terrestrial"
boulders <- c("Cophixalus_pakayakulangun",
              "Cophixalus_kulakula_1",
              "Cophixalus_zweifeli_1",
              "Cophixalus_petrophilus_1",
              "Cophixalus_saxatilis")
small.LSR[which(rownames(small.LSR)%in%boulders),"ecology"] <- "boulder"

allo.eco <- procD.pgls(small.LSR[1:4] ~ size * ecology, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(allo.eco)
    plotAllometry(allo.eco, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
eco.SVL <- procD.pgls(SVL ~ size * ecology, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(eco.SVL)
    plotAllometry(eco.SVL, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
eco.HLL <- procD.pgls(HLL ~ size * ecology, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(eco.HLL)
    plotAllometry(eco.HLL, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
eco.HW <- procD.pgls(HW ~ size * ecology, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(eco.HW)
    plotAllometry(eco.HW, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
eco.THR <- procD.pgls(THIRD ~ size * ecology, data=small.LSR, phy=small.tree, iter=999, print.progress = FALSE, plot=T); summary(eco.THR)
    plotAllometry(eco.THR, small.LSR$size, logsz=F, method="PredLine", ylim=c(-1,1))
    
# compare the simple models to complex ones
anova(allo.small, allo.eco)
anova(small.SVL, eco.SVL)
anova(small.HLL, eco.HLL)
anova(small.HW, eco.HW)
anova(small.THR, eco.THR)

    
# or course, a significant relationship reveals an allometric trend in your data 
# (shape changes with size), whereas no effect could be translated into isometry


# Use fitted values from the model to make prediction lines
plot(allo.pgls, type = "regression", 
     predictor = LSR$size, reg.type = "RegScore", 
     pch = 19, col = "green")



regime <- rep("terrestrial", 28)
regime[c(28,25,24,23,22,15)] <- "boulder" # terrestrial=0, boulder=1
names(regime) <- rownames(data)
tree.simm <- make.simmap(tree, regime)
frogs <- NULL; frogs$data <- data; frogs$phy <- tree.simm

estimateTimeMCMC(data=data[,c(1:2)], phy=tree.simm, gen=500000)
run1 <- ratematrixMCMC(data=data[,c(1:2)], phy=tree.simm, prior="empirical_mean",
                       start="mle", gen=500, dir="/Users/Ian/Google.Drive/R.Analyses/Microhylids/ratematrix")
run1 <- ratematrixMCMC(data=frogs$data[1:2], phy=frogs$phy, gen=500, start="mle", prior="uniform",
                       dir="/Users/Ian/Google.Drive/R.Analyses/Microhylids/ratematrix")
#run1 <- continueMCMC(run1, add.gen = 100000)
logAnalyzer(run1, burn=0.25, thin=1)
chain1 <- readMCMC(run1, burn=0.25)
checkConvergence(chain1)
plotRatematrix(chain=chain1)
testRatematrix(chain1, par="rates")
testRatematrix(chain1, par="correlation")

run2 <- ratematrixMCMC(data=data[,c(1:4)], phy=tree.simm, prior="empirical_mean", 
                       start="mle", gen=100000)
logAnalyzer(run2, burn=0.25, thin=1)
chain2 <- readMCMC(run2, burn=0.25)
plotRatematrix(chain=chain2)
testRatematrix(chain2, par="correlation")
testRatematrix(chain2, par="rates")

checkConvergence(chain2)


data("centrarchidae")
handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000,
                         dir="/Users/Ian/Google.Drive/R.Analyses/Microhylids/ratematrix")
handle <- continueMCMC(handle, add.gen = 10000)
logAnalyzer(handle, burn=0.25, thin=1)
posterior <- readMCMC(handle, burn=0.25, thin=1)
plotRatematrix(posterior)
#plotRootValue(posterior)
#plotPrior(handle)
#plotPrior(handle, root=TRUE)
testRatematrix(posterior, par="correlation")




chain1 <- ratematrixMCMC(data=frogs$data, phy=frogs$phy, prior=prior,
                         gen=100000, w_mu=w_mu, w_sd=w_sd,
                         dir="/Users/Ian/Google.Drive/R.Analyses/Microhylids/ratematrix")
logAnalyzer(chain1, burn=0.25, thin=1)
run1 <- readMCMC(chain1, burn=0.25)
checkConvergence(run1)
plotRatematrix(chain=chain1)
testRatematrix(chain1, par="rates")
testRatematrix(chain1, par="correlation")

data.range <- t(apply(frogs$data, 2, range))
w_mu <- (data.range[,2] - data.range[,1]) / 10
par.sd <- cbind(c(0,0), sqrt(c(10,10)))
w_sd <- matrix(0.2, ncol=2, nrow=4)
prior <- makePrior(r=4, p=2, den.mu="unif", par.mu=data.range,
                   den.sd="unif", par.sd=par.sd)

# run multiple mcmc chains by using lapply
res.list <- lapply(1:4, function(x) ratematrixMCMC(data=frogs$data, phy=frogs$phy, prior=prior,
                         gen=1000000, w_mu=w_mu, w_sd=w_sd,
                         dir="/Users/Ian/Desktop/output"))
# if you need to add extra generations to get the ESS up
for (i in 1:length(res.list)) {
  res.list[[i]] <- continueMCMC(res.list[[i]], add.gen=1000000,
                                dir="/Users/Ian/Desktop/output")
}

# read all the results into a list
posterior.list <- lapply(res.list, readMCMC)

#check for convergence within and among the runs
checkConvergence(posterior.list)

# merge the posteriors in the list
merged.posterior <- mergePosterior(posterior.list)

# plot the results
plotRatematrix(merged.posterior)
plotRootValue(merged.posterior)

# determine significance of correlations and rates:
testRatematrix(merged.posterior, par="all", plot=T) # or change 'par=c("rates", "correlation")'




logAnalyzer(chain2, burn=0.25, thin=1)
run2 <- readMCMC(chain2, burn=0.25, thin=1)
checkConvergence(run1, run2)

plotRatematrix(chain=run1)

ratematrixMCMC(data=frogs$data[1:2], phy=frogs$phy, gen=1000, 
               dir="/Users/Ian/Google.Drive/R.Analyses/Microhylids/ratematrix")

plotTree.barplot(tree,test4,list(fsize=0.4),
                 list(col="blue",space=1,log="x"))



tree.pca <- phyl.pca(tree, log(data[1:4]))
plot(tree.pca)
biplot(tree.pca)

tree.lsr.pca <- phyl.pca(tree, LSR[1:4])
plot(tree.lsr.pca)
biplot(tree.lsr.pca)
LSR <- cbind(LSR, tree.lsr.pca$S)

data <- as.data.frame(testdata[2:5]); rownames(data) <- testdata$Name_in_Tree
tree.pca <- phyl.pca(tree, data)
plot(tree.pca)
biplot(tree.pca)
data <- cbind(data, tree.pca$S)
write.csv(data, file="/Users/Ian/Google.Drive/R.Analyses/Microhylids/T428_PCAscores_reduced.csv")
#normal.pca <- prcomp(data[1:5])
#normal.pcs <- normal.pca$x


# Let's make a phylomorphospace plot to look at distribution of morphotypes
###########################################################################

# start by reading in the tree of interest 
all.tree <- read.nexus("~/Google.Drive/R.Analyses/Microhylids/T428_50loci_TRIMMED.tre")

# now read in the data you'd like to use
#all.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Microhylids/T428_PCAscores_reduced.csv", header=T) # excludes riparius, kaindiensis and nubicola
all.data <- read.csv("~/Google.Drive/R.Analyses/Microhylids/T428_PCAscores.csv", header=T)
  all.data <- all.data[which(!all.data$Name_in_Tree %in% c("Cophixalus_misimae",
                                                           "Cophixalus_nubicola_2",
                                                           "Cophixalus_kaindiensis",
                                                           "Cophixalus_riparius")),]


# You might notice the tree is for all the microhylid frogs
# We'll trim the data and tree down to match one another.
keepers <- intersect(all.tree$tip.label, all.data$Name_in_Tree)

# start by using 'drop.tip' and 'setdiff' to trim the tree
tree <- drop.tip(all.tree, setdiff(all.tree$tip.label, keepers))
# and then use 'filter' from 'dplyr' to reduce the data frame, 
# '%in%' searches for items of one matrix in another
pdata <- filter(all.data, Name_in_Tree %in% keepers)

# log transform the data to remove effects of massive differences in size, and normality
#pdata[,5:7] <- log(pdata[,5:7]) # pick the traits of interest, here you'd use your first two PCs
rownames(pdata) <- pdata$Name_in_Tree # provide names for the data

morpho <- pdata[,c("PC1", "PC2")]
ecology <- pdata$Ecology; names(ecology) <- pdata$Name_in_Tree

mycol <- character(length(ecology))
mycol[ecology=="saxicolous"] <- mpalette[3]
mycol[ecology=="terrestrial"] <- mpalette[1]
#mycol[ecology=="arboreal"] <- mpalette[2]



pdata$color <- mycol
pdata <- pdata[match(tree$tip.label, pdata$Name_in_Tree),]
group.colors <- pdata$color
names(group.colors) <- 1:Ntip(tree)
nodecols <- rep("black", tree$Nnode)
names(nodecols) <- (Ntip(tree)+1) : (Ntip(tree)+tree$Nnode)
colorz <- c(group.colors, nodecols)

phylomorphospace(tree, morpho, 
                 label = "horizontal",
                 node.size = c(1, 3),
                 control=list(col.node=colorz),
                 xlab = "PC1", ylab = "PC2")

phylomorphospace(tree, pdata[,c("PC1","PC2")], 
                 label = "horizontal",
                 node.size = c(1, 3),
                 control=list(col.node=lsr.colors),
                 xlab = "PC1", ylab = "PC2")

lsr.colors <- c(LSR$color, rep("black", Ntip(tree)-1)); names(lsr.colors) <- 1:length(lsr.colors)
phylomorphospace(tree, LSR[,c("PC1","PC2")],
                 label = "horizontal",
                 node.size = c(1, 3),
                 control=list(col.node=c(lsr.colors, rep("black", 23))),
                 xlab = "LSR_PC1", ylab="LSR_PC2")


# Let's determine the optimum number of phenotypic clusters
###########################################################################

cluster.data <- pdata[,2:5] # try using the PCA data or raw data
scaled <- scale(cluster.data)

fviz_nbclust(scaled, kmeans, method="wss")
fviz_nbclust(scaled, kmeans, method="silhouette")
fviz_nbclust(scaled, pam, method="wss")
fviz_nbclust(scaled, pam, method="silhouette")

nb <- NbClust(scaled, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "complete", index ="all")

pam.res <- pam(cluster.data, 2)
output <- fviz_cluster(pam.res, stand = F, geom = "point",
                       ellipse.type = "norm", show.clust.cent=T)
output+theme_classic()
cluster.data[,"clustering2"] <- pam.res$clustering
write.csv(cluster.frame, file="/Users/Ian/Desktop/Limbless.cluster.csv")



km.res <- kmeans(ln.test, 5, nstart = 25)
fviz_cluster(km.res, data=ln.test, geom = "point",
             stand = FALSE, frame.type = "norm")


km.res$cluster
ln.test[,"cluster"] <- km.res$cluster
ln.test[,"limbless?"] <- func.limbless
ln.test[,"taxon"] <- squam.species



#
## Another way to visualize correlation among traits
# info here: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html











p <- ggplot(diamonds, aes(carat, price))
p + geom_point(aes(colour = factor(cut))) + a_main_color()








library(corrplot)
M <- cor(mtcars)
corrplot(M, method = "ellipse", order="hclust", addrect=2)
corrplot.mixed(M, upper="ellipse")


Q <- cor(cdata)
corrplot(Q, method = "circle")





#
trees <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T428_Microhylidae/Alignments/T428_Astral.tre")

trees <- lapply(trees, root, outgroup="Kalophrynus_interlineatus", resolve.root = T)
trees <- lapply(trees, ladderize)
class(trees) <- "multiPhylo"

# for non-ultrametric trees we'll need to rescale them
for (k in 1:length(trees)) {
  tree <- chronos(trees[[k]], model="relaxed", lambda=0) # set lambda=1 for a consistent transformation
  #tree <- root(tree, outgroup="Laticauda_colubrina", resolve.root=T)
  write.tree(tree, file="/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T428_Microhylidae/Alignments/T428_Astral_RESCALED.trees", append=T)
}



curr.tree <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T428_Microhylidae/All_Alignments/T428_MCMCtree_50loci.tre")
write.tree(curr.tree, "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T428_Microhylidae/All_Alignments/T428_MCMCtree_50loci_newick.tre")


################################################################
## Distribution Maps for Microhylids
################################################################
setwd("/Users/Ian/Desktop")
library(phytools); library(mapdata)

coph.dist <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Microhylids/Cophixalus_Records_clean.csv", header=T)
    #coph.dist <- coph.dist[,c("Genus_species", "Latitude", "Longitude")]
    coph.dist <- coph.dist[,c("Name_in_Tree", "Latitude", "Longitude")]
        coph.dist <- coph.dist[complete.cases(coph.dist),]
        co.dist <- coph.dist[,1:2]; rownames(co.dist) <- coph.dist$Name_in_Tree; 
        names(co.dist)
        
dist.data <- dist.mat[,c("Tree.Clade", "DDLatitude", "DDLongitude")]; dist.data <- dist.data[complete.cases(dist.data),]
opoints <- as.matrix(dist.data[,2:3])
rownames(opoints) <- dist.data$Tree.Clade; colnames(opoints) <- c("lat", "long") 
opoints <- opoints[complete.cases(opoints),]
single.sp <- read.tree(file="/Users/Ian/Google.Drive/ANU Herp Work/Brad Projects/Single_Species.tre")

names(coph.dist)[1] <- "Name_in_Tree"
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/plot.distmaps.R")
plot.distmaps(coph.dist, new.directory="T428_DistMaps", base.map=cy_png)

co.dist <- as.matrix(coph.dist[,2:3])
rownames(co.dist) <- coph.dist$Name_in_Tree; colnames(co.dist) <- c("lat", "long")

single.sp <- read.tree("/Users/Ian/Google.Drive/R.Analyses/Microhylids/T428_Final.tre")
cols <- setNames(sample(viridis(n=Ntip(single.sp))), single.sp$tip.label)
obj <- phylo.to.map(single.sp, co.dist, rotate=TRUE, database="worldHires",
                    regions="Queensland", plot=FALSE)
plot(obj, direction="rightwards", colors=cols, ftype="off",cex.points=c(0,0.5),pts=FALSE)

cy_png <- get_googlemap(center = c(145,-13), zoom = 6, style="https://maps.googleapis.com/maps/api/staticmap?center=-24.229810153347742,-242.6166746271749&zoom=5&format=png&maptype=roadmap&style=element:labels%7Cvisibility:off&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:transit%7Cvisibility:off&size=480x360")
ggmap(cy_png)

col.pal <- colorRampPalette(brewer.pal(6,"RdYlBu")); point.colors <- col.pal(25)
names(point.colors) <- single.sp$tip.label

base.pal <- rep("black", 25); names(base.pal) <- single.sp$tip.label;
new.pal <- brewer.pal(5,"RdYlBu")
base.pal["Cophixalus_kulakula_1"] <- new.pal[1]
base.pal["Cophixalus_pakayakulangun"] <- new.pal[2]
base.pal["Cophixalus_petrophilus_1"] <- new.pal[3]
base.pal["Cophixalus_saxatilis"] <- new.pal[4]
base.pal["Cophixalus_zweifeli_1"] <- new.pal[5]

ggmap(cy_png) + geom_point(aes(x = Longitude, y = Latitude, fill=Name_in_Tree), colour = "black", 
                          data = coph.dist, size = 3, pch = 21) + theme(legend.position="bottom") + scale_fill_manual(values=base.pal)

ggmap(testo) + geom_point(aes(x = Longitude, y = Latitude, fill=Name_in_Tree), colour = "black", 
                          data = coph.dist, size = 2, pch = 21) + theme(legend.position="bottom") + scale_fill_brewer(palette="RdYlBu", "Taxon")

ggmap(testo) + geom_point(aes(x = Longitude, y = Latitude, fill=Name_in_Tree), colour = "black", 
                          data = coph.dist, size = 2, pch = 21) + theme(legend.position="bottom") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(6, "RdYlBu"))(Ntip(single.sp)), guide = guide_legend(nrow=5))


new.dist <- NULL
for (i in 1:length(unique(coph.dist$Name_in_Tree))){
  curr.dist <- filter(coph.dist, Name_in_Tree == unique(coph.dist$Name_in_Tree)[i])
  #min.dist <- filter(curr.dist, Latitude == min(Latitude))
  min.dist <- curr.dist[sample(1:nrow(curr.dist),1),]
  new.dist <- rbind(new.dist, min.dist)
}

ggmap(testo) + geom_point(aes(x = Longitude, y = Latitude, fill=Name_in_Tree), colour = "black", 
                          data = new.dist, size = 5, pch = 21) + theme(legend.position="bottom") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(6, "RdYlBu"))(Ntip(single.sp)), guide = guide_legend(nrow=5))



