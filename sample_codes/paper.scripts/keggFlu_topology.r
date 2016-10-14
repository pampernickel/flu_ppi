##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

# ---
# ---
# comparison of KEGG flu network as a function of PPI:
# Figure 3

library(igraph)
library(STRINGdb)
library(rentrez)
library(tm)
library(stringdist)
library(parallel)
library(MASS)

source('./sample_codes/functions/routineFuncs.r')
source('./sample_codes/functions/retrievalFuncs.R')
source('./sample_codes/functions/graphFuncs.r')

# ---
# ---
# Load pre-processed/pre-created data files data files
# Files created with ./paper.scripts/string_annotation.r
load("./r.data.files/string.v10.rda") # g
load("./r.data.files/hippie.rda.1.8") # ppi.igraph.1.8

# KEGG graphs, generated using the same methods shown in sampleScripts.r
# then merged with STRING v.10 and v.9 graphs to create object
# a list object kegg.string.graphs. kegg.string.graphs is a list
# with two slots: kegg.string.graphs[[1]] -> kegg in string v.9.05; kegg.string.graphs[[2]] -> kegg in string.v.10
load("./r.data.files/ppi_review/kegg.string.graphs.rda") #kegg.string.graphs

# Flu graph
kegg.string.graphs[[2]][[grep("05164", names(kegg.string.graphs[[2]]))]] -> kegg.flu
get.data.frame(kegg.flu) -> kegg.flu.df
V(kegg.flu)$name -> hits

# kegg graph topology
induced.subgraph(g, which(V(g)$name %in% 
                            V(kegg.flu)$name)) -> kegg.flu.string
as.numeric(apply(kegg.flu.df[,c(1,2)], 1,function(x)
  get.edge.ids(kegg.flu.string, x, directed = F))) -> true.edge.ind
E(kegg.flu.string)$weight <- rep(0.1, ecount(kegg.flu.string))
E(kegg.flu.string)$weight[which(E(kegg.flu.string)$combined_score <= 200)] <- 0.25
E(kegg.flu.string)$weight[which(E(kegg.flu.string)$combined_score >= 800)] <- 0.25
E(kegg.flu.string)$weight[true.edge.ind] <- 1.6
E(kegg.flu.string)$color <- rep("grey70", ecount(kegg.flu.string))
E(kegg.flu.string)$color[which(E(kegg.flu.string)$combined_score <= 200)] <- "blue"
E(kegg.flu.string)$color[which(E(kegg.flu.string)$combined_score >= 800)] <- "red"
E(kegg.flu.string)$color[true.edge.ind] <- "grey20"

induced.subgraph(ppi.igraph.1.8, which(V(ppi.igraph.1.8)$name %in% 
                                         V(kegg.flu)$name)) -> kegg.flu.hippie
igraph::simplify(kegg.flu.hippie, remove.loops=T) -> kegg.flu.hippie
as.numeric(apply(kegg.flu.df[,c(1,2)], 1,function(x)
  get.edge.ids(kegg.flu.hippie, x, directed = F))) -> true.edge.ind
E(kegg.flu.hippie)$weight <- rep(0.1, ecount(kegg.flu.hippie))
E(kegg.flu.hippie)$weight[which(E(kegg.flu.hippie)$score <= 0.2)] <- 0.25
E(kegg.flu.hippie)$weight[which(E(kegg.flu.hippie)$score >= 0.8)] <- 0.25
E(kegg.flu.hippie)$weight[true.edge.ind[which(true.edge.ind > 0)]] <- 1.6
E(kegg.flu.hippie)$color <- rep("grey70", ecount(kegg.flu.hippie))
E(kegg.flu.hippie)$color[which(E(kegg.flu.hippie)$score <= 0.2)] <- "blue"
E(kegg.flu.hippie)$color[which(E(kegg.flu.hippie)$score >= 0.8)] <- "red"
E(kegg.flu.hippie)$color[true.edge.ind[which(true.edge.ind > 0)]] <- "grey20"

# use the bigger graph to creat color map
colors.neg <- colorRampPalette(brewer.pal(9, "Reds"))(length(unique(V(kegg.flu.string)$z.rsa[which(as.numeric(V(kegg.flu.string)$z.rsa) < 0)])))
colors.pos <- colorRampPalette(brewer.pal(9, "Blues"))(length(unique(V(kegg.flu.string)$z.rsa[which(as.numeric(V(kegg.flu.string)$z.rsa) >= 0)])))
unique(V(kegg.flu.string)$z.rsa[which(as.numeric(V(kegg.flu.string)$z.rsa) < 0)]) -> neg.scores
unique(V(kegg.flu.string)$z.rsa[which(as.numeric(V(kegg.flu.string)$z.rsa) >= 0)]) -> pos.scores
rbind(cbind(neg.scores[order(neg.scores)], colors.neg[c(length(colors.neg):1)]),
      cbind(pos.scores[order(pos.scores)], colors.pos)) -> color.map

graph <- kegg.flu.hippie
sapply(V(graph)$z.rsa, function(x)
  color.map[which(color.map[,1] %in% x),2]) -> V(graph)$colors
as.numeric(V(graph)$z.rsa)*-1 -> inv.rsa
size=inv.rsa+abs(min(inv.rsa))+0.05

get.data.frame(graph, what="edges") -> edge.df.string
induced.subgraph(graph, 
                 which(V(graph)$name %in% c(edge.df.string$to, edge.df.string$from))) -> graph

plot(graph, layout=layout.fruchterman.reingold, 
     vertex.size=size, 
     vertex.label.dist=0.6, vertex.label.cex=0.5, 
     vertex.color=V(graph)$colors,
     vertex.label.family="sans",
     edge.color=E(graph)$color,
     vertex.label.color= "grey20",
     edge.width=E(graph)$weight) 