##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================
library(STRINGdb)

loadAliases <- function(v=""){
  # v = version; alias tables are very different in v.9.05 and v.10
  # in this version, files are stored locally, but original files can be found here:
  # http://string.uzh.ch/permanent/string/9_1/protein_aliases/9606__protein_aliases_tf.tsv.gz 
  # http://string.embl.de/newstring_download/protein.aliases.v10/9606.protein.aliases.v10.txt.gz
  if (v %in% "9_05" | v %in% ""){
    aliases <- read.delim("./string.files/9606__protein_aliases_tf.tsv", stringsAsFactor=F)
  } else if (v %in% "10"){
    # direct read takes too long; use processed file
    aliases <- read.csv("./string.files/protein.aliases.v10.sub.csv") 
  }  
  return(aliases)
}

graph2STRING <- function(g, string.db=NULL){
  # g = igraph object with HGNC gene symbols as names (other forms not yet handled)
  # converts g to form with STRING identifiers; returns, g, STRING equivalent and
  # mapping between HGNC symbol and STRING identifiers; should preclude the need for
  # loadAliases from local by returning the initial mapping
  g -> g.s
  if (is.null(string.db)){
    string.db <- STRINGdb$new(version="10", species=9606,score_threshold=0, input_directory="")
  }
  
  vertex_map <- string.db$map(prepareMap(V(g)$name), "vertex", removeUnmappedRows = TRUE)
  as.character(unlist(sapply(V(g)$name, function(x) 
    ifelse(x %in% vertex_map$vertex, vertex_map$STRING_id[which(vertex_map$vertex %in% x)], x)))) -> V(g.s)$name
  
  # then drop vertices that cannot be mapped
  induced.subgraph(g.s, grep("9606", V(g.s)$name)) -> g.s
  list(g, g.s, vertex_map) -> res
  names(res) <- c("graph", "STRING", "vmap")
  return(res)
}