setwd("~/Doutorado/Novo projeto/Capitulo 1")

library(iNEXT)
library(betapart)
library(vegan)
library(ggplot2)
library(dendextend)
library(condvis)

# imprtar matriz com especies nas linhas e fitofisionomias nas colunas
data.fito <- read.csv("Dados_PE.csv", row.names = 1)
data.fito[is.na(data.fito)] <- 0

# curva de rarefação q0
q0fito <- iNEXT(data.fito, q = 0, datatype = "abundance", endpoint = 2500)
q0plot <- ggiNEXT(q0fito)
q0plot + scale_color_manual(values=c("#E69F00", "#004D40", "#897454", "#D81B60", "#1E88E5"))

# curva de rarefação q1
q1fito <- iNEXT(data.fito, q = 1, datatype = "abundance", endpoint = 2500)
q1plot <- ggiNEXT(q1fito)
q1plot + scale_color_manual(values=c("#E69F00", "#004D40", "#897454", "#D81B60", "#1E88E5"))

# curva de rarefação q2
q2fito <- iNEXT(data.fito, q = 2, datatype = "abundance", endpoint = 2500)
q2plot <- ggiNEXT(q2fito)
q2plot + scale_color_manual(values=c("#E69F00", "#004D40", "#897454", "#D81B60", "#1E88E5"))


# importar matriz com especies nas colunas e parcelas nas linhas
data_par <- read.csv("Dados_PE_parcelas.csv", row.names = 1)
data_par[is.na(data_par)] <- 0


# calcular matriz de dissimilaridade de Bray-Curtis
bird_matrices <- bray.part(data_par)
bird_matrices[[3]]

# classificar as parcelas de acordo com as fitofisionomias e atribuir cores
groups <- factor(c(rep(1,9), rep(2,11), rep(3,26), rep(4,10)), labels = c("campo","FOD", "FOM", "vassoural"))
groups2 <- factor(c(rep(1,9), rep(2,11), rep(3,16), rep(4,10), rep(5,10)), labels = c("Campo","FOD", "MN", "FOM", "Vassoural"))
col_fito2 <- factor2color(groups2, colors = c("#E69F00", "#004D40", "#D81B60", "#897454", "#1E88E5"))

bd<-betadisper(bird_matrices[[3]], group = groups2)
bd

## Perform test
anova(bd)

## Permutation test for F
permutest(bd, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(bd.HSD <- TukeyHSD(bd))
plot(bd.HSD)

## Plot  and show the groups and distances to centroids on the
## first two PCoA axes
plot(bd, col = c("#E69F00", "#004D40", "#D81B60", "#897454", "#1E88E5"))
betadistances(bd)
plot(bd, ellipse = T, hull = F, conf = 0.90, col = c("#E69F00", "#004D40", "#D81B60", "#897454", "#1E88E5"),
     main = "diversidade beta")

## Draw a boxplot of the distances to centroid for each group
boxplot(bd$distances ~ reorder(bd$group, bd$distances, median),
        col = c("#004D40", "#D81B60", "#1E88E5", "#897454", "#E69F00"),
        xlab = "Fitofisionomia", ylab = "Distância do centroide")

dendo <- as.dendrogram(hclust(bird_matrices[[3]]))
labels_colors(dendo) <- col_fito2[order.dendrogram(dendo)]
dendo <- color_branches(dendo, col = c("#897454", "#E69F00", "#E69F00", "#004D40", "#897454", "#1E88E5","#D81B60"), k = 7)
plot(dendo)
