##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

library(igraph)
library(STRINGdb)

source('./scripts/stringFuncs.r')
source('./scripts/GOFuncs.r')
source('./scripts/graphAnnotationFuncs.r')

read.csv("./data/annotations/hgnc.csv")

# ---
# ---
# tested: version="10", version="9_05"
string.db <- STRINGdb$new( version="10", species=9606,score_threshold=0, input_directory="" )
string.db$get_graph() -> string.graph

# Map IAV metaanalysis results to STRING identifiers
all.genes <- read.csv("./data/metaanalysis/all.genes.csv", stringsAsFactor=F)
all.genes.map <- string.db$map(all.genes, "Symbol", removeUnmappedRows = TRUE )

# subset string.graph to contain only genes that are mapped from all.genes
induced_subgraph(string.graph, which(V(string.graph)$name %in% all.genes.map$STRING_id)) -> g

# for each edge, map back to protein ID;
# restrict aliases to relevant ones for retrieval speedup
loadAliases() -> aliases
aliases[which(aliases$protein_id %in% V(g)$name),] -> aliases
aliases[which(aliases$alias %in% hgnc$hgnc_symbol),] -> aliases

sapply(V(g)$name, function(x) 
  return(ifelse(x %in% aliases$protein_id,
                as.character(aliases$alias[which(aliases$protein_id %in% x)]), x))) -> V(g)$name
induced_subgraph(g, setdiff(c(1:vcount(g)), 
                            grep("9606", V(g)$name))) -> g

# associate z.rsa score with each node
sapply(V(g)$name, function(x)
  ifelse(x %in% all.genes$Symbol, 
         all.genes$Z_RSA[which(all.genes$Symbol %in% x)], NA)) -> t
as.numeric(t) -> t
set.vertex.attribute(g, "z.rsa", value=as.numeric(t)) -> g

getgos() -> gos
nodeToGO(g, gos$symbol) -> g
nodeToGO(g, gos.cc$symbol, "gos.cc") -> g

# save(g, file="./r.data.files/string.v10.rda")