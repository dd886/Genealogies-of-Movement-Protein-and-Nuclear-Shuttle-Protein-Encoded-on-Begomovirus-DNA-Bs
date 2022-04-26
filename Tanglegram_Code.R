install.packages("installr")
library(installr)
updateR()

#install.packages("Rtools", repos="https://cran.rstudio.com/bin/windows/Rtools/", quiet=TRUE)

install.packages("ape")
install.packages("ape",repos="https://cloud.r-project.org",quiet=TRUE)
install.packages("dendextend")
install.packages("viridis")
install.packages("dplyr")
install.packages("phylogram")
install.packages("phangorn",repos="https://cloud.r-project.org",quiet=TRUE)
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)

library(ape)
library(phytools)
library(dendextend)
library(viridis)
library(dplyr)
library(phylogram)

tree1 <- read.tree(file = "file 1 name")
#tree1 <- midpoint.root(tree1)
tree2 <- read.tree(file = "file 2 name")
plot(tree2)
#tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1 <- as.dendrogram(tree1)
tree1
tree2 <- as.dendrogram(tree2)
tree2
dndlist <- dendextend::dendlist(tree1, tree2)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 1.8, lab.cex = 0.3, lwd = 
                         0.5, edge.lwd = 0.5, type = "r")

entanglement(dndlist)

dev.copy(pdf, 'Tanglegram_MP_NSP_Ordered', width = 10, height = 11)
dev.off()
