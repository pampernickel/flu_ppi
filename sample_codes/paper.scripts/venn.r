##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

# ---
# ---
# comparison of v9_05 and v10 networks based on host factors:
# Figures 2A and B

# check edges in common; restrict vhn.g for this comparison to human-human interactions included
# add version 10 for string
library(ggplot2)
library(igraph)
library(reshape2)
library(VennDiagram)

source('./scripts/routineFuncs.r')

# Files created with ./paper.scripts/string_annotation.r
# In the case of HIPPIE, these were created in an analogus manner from the HIPPIE v1_7 and
# HIPPIE v1_8 data dump:
# http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_v1_8.txt
# http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_v1_7.txt
load("./r.data.files/string.v10.rda") #g.10
load("./r.data.files/string.v9.rda") #g.9
load("./r.data.files/hippie.rda.1.8") # ppi.igraph.1.8
load("./r.data.files/hippie.rda.1.7") # ppi.igraph.1.7

list(g.10, g.9, ppi.igraph.1.8, ppi.igraph.1.7) -> graphs
names(graphs) <- c("STRING_v10", "STRING_v9_05", "HIPPIE_v1_8", "HIPPIE_v1_7")
lapply(graphs, function(x) ecount(x)) -> edges
lapply(graphs, function(x) vcount(x)) -> vertices
combn(4,2) -> combs.1
apply(combs.1, 2, function(x)
  graph.intersection(graphs[[x[1]]], graphs[[x[2]]], keep.all.vertices=F)) -> 
  double.int

# triple and quad ints
combn(4,3) -> combs.2
apply(combs.2, 2, function(x)
  graph.intersection(graphs[[x[1]]], graphs[[x[2]]], graphs[[x[3]]], keep.all.vertices=F)) -> 
  triple.int
graph.intersection(graphs[[1]], graphs[[2]], graphs[[3]],graphs[[4]],
                   keep.all.vertices=F) -> quad.int
lappend(double.int, triple.int) -> graph.intersection
graph.intersection[[11]] <- quad.int

sapply(graph.intersection, function(x) vcount(x)) -> common.vertices
sapply(graph.intersection, function(x) ecount(x)) -> common.edges

# create table
t.1a <- rbind(unlist(edges), 
              unlist(vertices))
colnames(t.1a) <- c("STRING_v10", "STRING_v9_05", "HIPPIE_v1_8", "HIPPIE_v1_7")

#code for venn_edges, venn_vertices
draw.quad.venn(vertices[[1]], vertices[[2]], vertices[[3]], vertices[[4]],
               common.vertices[1], common.vertices[2], common.vertices[3], common.vertices[4],
               common.vertices[5], common.vertices[6], common.vertices[7], common.vertices[8],
               common.vertices[9], common.vertices[10], common.vertices[11], 
               category = c("STRING,\nv10", "STRING,\nv9.05", "HIPPIE,\nv1.8", "HIPPIE,\nv1.7"),
               reverse = FALSE, lwd = rep(2, 4), 
               lty = rep("solid", 4), col = rep("black", 4), 
               fill = c("khaki4","lemonchiffon2", "gold4","forestgreen"), alpha = 0.4,
               cex = rep(1.25, 15),
               fontfamily = rep("sans", 15),
               cat.cex = rep(1.3, 4),
               cat.fontfamily = rep("sans", 4),
               cat.fontface = rep("bold", 4))
draw.quad.venn(edges[[1]], edges[[2]], edges[[3]], edges[[4]],
               common.edges[1], common.edges[2], common.edges[3], common.edges[4],
               common.edges[5], common.edges[6], common.edges[7], common.edges[8],
               common.edges[9], common.edges[10], common.edges[11], 
               category = c("STRING,\nv10", "STRING,\nv9.05", "HIPPIE,\nv1.8", "HIPPIE,\nv1.7"),
               reverse = FALSE, lwd = rep(2, 4), 
               lty = rep("solid", 4), col = rep("black", 4), 
               fill = c("khaki4","lemonchiffon2", "gold4","forestgreen"), alpha = 0.4,
               cex = rep(1.25, 15),
               fontfamily = rep("sans", 15),
               cat.cex = rep(1.3, 4),
               cat.fontfamily = rep("sans", 4),
               cat.fontface = rep("bold", 4))

# other checks:
# edges in v.9.05 NOT in v.10
graph.intersection(g.9, g.10, keep.all.vertices=F) -> string.i
graph.difference(g.10, g.9) -> string.new
graph.difference(g.9, string.i) -> ni.9 #percentage of g.9 ni intersection: removed edges
graph.difference(g.10, string.i) -> ni.10 #percentage of g.10 ni intersection: new edges


ecount(graph.difference(g.9, string.ip))
ecount(string.ip)/min(c(ecount(g.10.p), ecount(g.9)))