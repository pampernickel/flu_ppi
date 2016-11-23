##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

library(XML)
library(igraph)
library(graph)
library(ROntoTools)
library(KEGGgraph)
library(rentrez)
library(org.Hs.eg.db)

symbol <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(symbol)
symbol.map <- as.list(symbol[mapped_genes])

# Retrieval from PubMed
getAbstract <- function(xml.trees, call.origin="abstract"){
  # --- getText should also be customized to put in NAs
  # --- when the text is not found
  # --- should be able to do with calls from abstract extraction,
  # --- full text extraction modes (pmc)
  text <- list()
  for (i in 1:length(xml.trees)){
    current.tree <- xml.trees[[i]]
    sub.text <- list()
    if (call.origin=="abstract"){
      if (length(current.tree) > 0){
        for (j in 1:length(current.tree)){
          getNodeSet(current.tree[[j]], '//Abstract') -> xml.text
          if (length(xml.text) > 0){
            sub.text[[j]] <- xmlValue(xml.text[[1]])
          } else {
            sub.text[[j]] <- NA
          }
        }
      }
    } else if (call.origin=="pmc") {
      if (length(current.tree) > 0){
        for (j in 1:length(current.tree)){
          getNodeSet(current.tree[[j]], '//abstract') -> xml.text
          if (length(xml.text) > 0){
            sub.text[[j]] <- xmlValue(xml.text[[1]])
          } else {
            sub.text[[j]] <- NA
          }
        }
      }
    }
    text[[i]] <- sub.text
  } 
  return(text)
}

checkXML <- function(file){
  file.info(file)$size
}

# Retrieval from KEGG
getKEGGgraph <- function(pid){  
  tmpXML <- "test_kegg.xml"
  # retrieveKGML works fine in Linux systems, but not in Mac
  # retrieveKGML(pid, organism='hsa', destfile=tmpXML, method="wget")
  # alternative implementation
  getKGMLurl(pathwayid = pid, organism = "hsa") -> kgml
  download.file(kgml, destfile = tmpXML)
  
  # check test_kegg.xml if it contains a valid (i.e. Homo sapiens) network
  # or if it is empty -- in which case the map is comprised purely of 
  # non-human proteins (and no links to human protein interactants)
  kegg.i <- NA
  if (checkXML(tmpXML) > 10){
    kegg.g <- parseKGML2Graph(tmpXML,expandGenes=TRUE)  
    kegg.g.n <- getKEGGnodeData(kegg.g)
    mapkGedgedata <- getKEGGedgeData(kegg.g)
    
    outs <- sapply(edges(kegg.g), length) > 0
    ins <- sapply(inEdges(kegg.g), length) > 0
    ios <- outs | ins
    
    ioGeneID <- translateKEGGID2GeneID(names(ios))
    nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
    nAttrs <- list()
    nAttrs$fillcolor <- makeAttr(kegg.g, "lightgrey", list(orange=names(ios)[ios]))
    nAttrs$label <- as.character(nodesNames)
    igraph.from.graphNEL(kegg.g, name = TRUE, weight = TRUE,
                         unlist.attrs = TRUE) -> kegg.i
    as.character(sapply(gsub("hsa:", "", V(kegg.i)$name), function(x) 
      nodesNames[which(names(nodesNames) %in% x)])) -> V(kegg.i)$name
  }
  return(kegg.i)
}

makeAttr <- function(graph, default, valNodeList) {
  tmp <- nodes(graph)
  x <- rep(default, length(tmp)); names(x) <- tmp
  
  if(!missing(valNodeList)) {
    stopifnot(is.list(valNodeList))
    allnodes <- unlist(valNodeList)
    stopifnot(all(allnodes %in% tmp))
    for(i in seq(valNodeList)) {
      x[valNodeList[[i]]] <- names(valNodeList)[i]
    }
  }
  return(x)
}

# Retrieval, mapping to and from STRING
prepareMap <- function(id){
  # prepares data for mapping to STRING ids
  # id: char string, set of gene names
  as.data.frame(id) -> df
  colnames(df) <- "vertex"
  return(df)
}

mapToSTRING2 <- function(vertex_map, id){
  sapply(id, function(x) return(vertex_map$STRING_id[which(vertex_map$vertex %in% x)])) -> id
  return(id)
}

mapToSTRING <- function(vertex_map, g, mode=c("directed", "undirected")){
  get.data.frame(g) -> df
  mapToSTRING2(vertex_map, df$from) -> df$from
  mapToSTRING2(vertex_map, df$to) -> df$to
  if (mode == "directed"){
    graph.data.frame(d=df, directed=T)  -> g
  } else {
    graph.data.frame(d=df, directed=F)  -> g
  }
  return(g)
}

readKEGGbrite <- function(path){
  readLines(path) -> txt
  txt[grep("Infectious", txt)] -> headings
  grep("Infectious", txt) -> ind
  # given that the interest is in Bacterial and Viral diseases
  ind[1] -> s
  ind[3] -> e
  unlist(lapply(strsplit(txt[s:e], " "), function(x) x[grep("0", x)])) -> all.pathways
  return(all.pathways)
}