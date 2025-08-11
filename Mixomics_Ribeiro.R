rm(list=ls())
setwd("~/Desktop/R")
library(mixOmics)
library(tidyverse)
library(dplyr)
library(plyr)
library(tidyr)
microbiome <- read.csv("microbiome.csv", row.names = 1, check.names = FALSE)
metabolomics <- read.csv("metabolomics.csv", row.names = 1, check.names = FALSE)

# Entry data
data <- list(
  microbiome = as.matrix(sapply(microbiome, as.numeric)),
  metabolomics = as.matrix(sapply(metabolomics, as.numeric))
)


# Unique names in each dataset
colnames(microbiome) <- make.unique(colnames(microbiome))
colnames(metabolomics) <- make.unique(colnames(metabolomics))


# Response variables
Y <- factor(c("Runoff", "Runoff", "Runoff", "Control", "Control", "Control"))

# Matrix 
design <- matrix(0.1, ncol = length(data), nrow = length(data),
                 dimnames = list(names(data), names(data)))
diag(design) <- 0

# Variables in each block
list.keepX <- list(
  microbiome = c(1,34),  # Ajuste os índices conforme necessário
  metabolomics = c(1,18)
)

# Execute DIABLO
sgccda.res <- block.splsda(X = data, Y = Y, ncomp = 2, 
                           keepX = list.keepX, design = design)

# Circos plot

  circosPlot(sgccda.res, cutoff = 0.65, block.labels.adj = -1, line = TRUE, 
             color.blocks = c('chocolate', 'purple'),
             color.cor = c("green", "red"),color.Y= c("#56B4E9","#E69F00"), 
             size.variables = 1, var.adj = .8, legend = FALSE,
             size.labels =0.001, linkWidth= 2, showIntraLinks=FALSE)
  

  microbiome<-microbiome[, -(1)]
  metabolomics<-metabolomics[, -(1:3)]
  # Dados de entrada
  # Reforçar que todos os valores nos blocos são numéricos
  data <- list(
    microbiome = as.matrix(sapply(microbiome, as.numeric)),
    metabolomics = as.matrix(sapply(metabolomics, as.numeric))
  )
  colnames(microbiome) <- make.unique(colnames(microbiome))
  colnames(metabolomics) <- make.unique(colnames(metabolomics))
  Y <- factor(c("Runoff", "Runoff", "Runoff", "Control", "Control", "Control"))
  
  X <- as.data.frame(data)
  
  phyllo.splsda <- splsda(sqrt(X), Y,ncomp=6)  # set ncomp to 10 for performance assessment later
  
  # plot the samples projected onto the first two components of the PLS-DA subspace

    plotIndiv(phyllo.splsda , comp = 1:2, 
              group = Y, ind.names = FALSE,  # colour points by class
              ellipse = TRUE, # include 95% confidence ellipse for each class
              legend = TRUE, title = '(a) PLSDA with confidence ellipses',
              size.title = rel(3), 
              size.xlabel = rel(4),
              size.ylabel = rel(4),
              size.axis = rel(4),
              size.legend = rel(2))

  # use the max.dist measure to form decision boundaries between classes based on PLS-DA data
  background = background.predict(phyllo.splsda, comp.predicted=2, dist = "max.dist")
  
  # plot the samples projected onto the first two components of the PLS-DA subspace

    plotIndiv(phyllo.splsda, comp = 1:2,
              group = Y, ind.names = FALSE,
              col.per.group= c("#56B4E9","#E69F00"),# colour points by class
              style = "ggplot2",
              cex=8,
              background = background, # include prediction background for each class
              legend = FALSE, title = "sPLS-DA",
              size.title = rel(3), 
              size.xlabel = rel(4),
              size.ylabel = rel(4),
              size.axis = rel(4),
              size.legend = rel(2))

    