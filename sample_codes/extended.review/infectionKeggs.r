##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

library(STRINGdb)
library(igraph)
library(rentrez)
library(XML)
library(stringr)
library(plyr)

# custom functions
source('./sample_codes/functions/retrievalFuncs.r')

# ---
# ---
# download all KEGG references linked to infection
readKEGGbrite("./data/annotations/br08901.keg") -> all.pathways
lapply(all.pathways, function(x){
  print(paste("Processing hsa", x, " ...", sep=""))
  html <- getURL(paste("http://www.genome.jp/dbget-bin/www_bget?hsa", x, sep=""))
  doc = htmlParse(html, asText=TRUE)
  links <- as.character(xpathSApply(doc, "//a/@href"))[grep("pubmed", 
                                                            as.character(xpathSApply(doc, "//a/@href")))]
  if (length(links) > 0){
    gsub("http://www.ncbi.nlm.nih.gov/pubmed/", "", links) -> links
    paste("PMID:", links, sep="") -> links
  } else {
    # try map
    html <- getURL(paste("http://www.genome.jp/dbget-bin/www_bget?map", x, sep=""))
    doc = htmlParse(html, asText=TRUE)
    links <- as.character(xpathSApply(doc, "//a/@href"))[grep("pubmed", 
                                                              as.character(xpathSApply(doc, "//a/@href")))]
    links[grep("ncbi", links)] -> links
    paste("PMID:", links, sep="") -> links
  }
  
  # consider changing the delay, or downloading references by batch
  refs <- NA
  if (length(links) > 0){
    xml.trees <- sapply(links, function(x) entrez_fetch(db="pubmed", id=x, rettype="xml", delay=3))
    xml.trees <- lapply(xml.trees, function(x) xmlTreeParse(x, useInternalNodes=T))
    abstract <- unlist(getAbstract(list(xml.trees)))
    list(links, abstract) -> refs
    names(refs) <- c("pids", "abstract")
  }
  return(refs)
}) -> pathway.refs
names(pathway.refs) <- all.pathways
save(pathway.refs, file="./data/r.data.files/kegg_refs/infectiousDiseases.rda")
# ---
# ---

# ---
# ---
# download KEGG networks linked to all infectious disease pathways available
lapply(all.pathways, function(x){
  getKEGGgraph(x)
}) -> idnets
names(idnets) <- names(all.pathways)
idnets[which(sapply(idnets, function(x) class(x)) %in% "igraph")] -> idnets
save(idnets, file="./data/r.data.files/kegg_refs/infectiousDiseaseGraphs.rda")
# ---
# ---


# ---
# ---
# download KEGG networks linked to all infectious disease pathways available
# map kegg graphs to string graphs/identifiers
# note that in some cases, the session might be disconnected (on STRING load)
# to be on the safe side, end with save.image() after lapply statement
load("./data/r.data.files/kegg_refs/infectiousDiseaseGraphs.rda")
string.db <- STRINGdb$new(version="10", species=9606,score_threshold=0, input_directory="" )
string <- string.db$get_graph()
lapply(idnets, function(x){
  vertex_map <- string.db$map(prepareMap(V(x)$name), "vertex", removeUnmappedRows = TRUE)
  mapToSTRING(vertex_map, x, mode="undirected") -> x
}) -> idnets_string
save(idnets, idnets_string, file="./data/r.data.files/kegg_refs/infectiousDiseaseGraphs.rda")

# then select edges that are in STRING
lapply(idnets_string, function(x) get.data.frame(x, what="edges")) -> df
lapply(df, function(x){
  intersect(which(x$to %in% V(string)$name), 
            which(x$from %in% V(string)$name)) -> ind
  x[ind,] -> df
}) -> df

# collapse to df with no duplicate edges, then extract references per edge from string
ref.df <- ldply (df, data.frame)
if (length(which(duplicated(ref.df))) > 0){
  ref.df[-which(duplicated(ref.df)),] -> ref.df
}

# fetch with for loop instead of an apply, as connection might time out, 
# and a partial save is advisable
pmids <- list()
for (i in 1:nrow(ref.df)){
  pmids[[i]] <- string.db$get_pubmed_interaction(ref.df$from[i], ref.df$to[i])
  if (i %in% seq(1,nrow(ref.df),by=10)){
    print(paste("Saving reference set ", i, "...", sep=""))
    save(pmids, file="./data/r.data.files/kegg_refs/infectiousDiseaseSTRING_refs.rda")
  }
}
# ---
# ---


# ---
# ---
# get network neighborhood of each kegg network in STRING; use version that is
# already mapped to gene symbols; filter graph first using default 400 filter
load("./data/r.data.files/string.v10.rda")
subgraph.edges(g, which(E(g)$combined_score >= 400), delete.vertices=T) -> g
lapply(idnets, function(x) as.undirected(x)) -> idnets
lapply(idnets, function(x) 
  lapply(V(x)$name, function(y) getNeighbors(g, y, ord=1))) -> nn
lapply(1:length(nn), function(x){
  constructGraph(nn[[x]], g, V(idnets[[x]])$name)}) -> idnets_nn
save(idnets_nn, file="./data/r.data.files/kegg_refs/infectiousDiseaseGraph_netneighborhood.rda")
# ---
# ---