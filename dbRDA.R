#######Microbial pollution disables the chemical defenses of sea fans
##Code for dbRDA analysis
library(vegan)
library(permute)
library(devtools)
microbiome <- read.csv("microbiome.csv", row.names = 1, check.names = FALSE)
metabolomics <- read.csv("metabolomics.csv", row.names = 1, check.names = FALSE)
microbiome$site <- factor(c("Runoff", "Runoff", "Runoff", "Control", "Control", "Control"))
microbiome<-microbiome%>% relocate(site)


dbRDA2 = dbrda(metabolomics~site+
                 Sphingomonas+Novosphingobium+Endozoicomonas
               , microbiome, dist="bray")
##Model inspection
plot(dbRDA2)
summary(dbRDA2)
#Terms
anova(dbRDA2)
anova(dbRDA2, by="terms",perm.max=1000)
vif.cca(dbRDA2)
#R-square and alias
RsquareAdj(dbRDA2)
alias(dbRDA2)
#Components
R.sum <- summary(dbRDA2)
R.sum$cont
R.sum$cont$importance[2, "dbRDA1"]
R.sum$cont$importance[2, "dbRDA2"]
