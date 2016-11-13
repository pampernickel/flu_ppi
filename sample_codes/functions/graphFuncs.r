##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

require(plyr)

calcSimilarity <- function(graph.1, graph.2){
  # http://lists.nongnu.org/archive/html/igraph-help/2008-04/msg00017.html
  # --- pairwise comparison of graphs in terms of overlapping edges,
  # --- implicitly, edges
  # sum(get.adjacency(sp.g$subgraph) != get.adjacency(mst.g.0$subgraph))
  # works only with graphs of the same size
  # graph.intersection(sp.g$subgraph, mst.g.0$subgraph)
  # doesn't work
  
  length(E(graph.2))+length(E(graph.1)) -> edge.tot
  # --- for each graph, check overlapping edge sequences
  get.edgelist(graph.1) -> graph.1
  get.edgelist(graph.2) -> graph.2
  
  apply(graph.1, 1, function(x) checkEdge(x[1], x[2], graph.2)) -> in.2
  apply(graph.2, 1, function(x) checkEdge(x[1], x[2], graph.1)) -> in.1
  
  length(c(which(in.1 %in% TRUE), which(in.2 %in% TRUE)))/edge.tot -> similarity.score
  return(similarity.score)
}

calcVertexSimilarity <- function(graph.1, graph.2){
  # --- pairwise comparison of graphs in terms of overlapping
  # --- vertices, regardless of edge included
  length(V(graph.1)$name)+length(V(graph.2)$name) -> v.tot
  length(which(V(graph.1)$name %in% V(graph.2)$name))+
    length(which(V(graph.2)$name %in% V(graph.1)$name)) -> subtot
  subtot/v.tot -> similarity
  return(similarity)
}


getNeighbors <- function(graph, queries, ord=1){
  # wrapper for igraph.neighborhood function
  unique(unlist(c(which(V(graph)$name %in% queries), 
                  ego(graph, which(V(graph)$name %in% queries), order=ord)))) -> ind
  return(V(graph)$name[ind])
}

getNodeAttribute <- function(graph, attribute, node){
  # get attribute for a given node in a graph
  a <- get.vertex.attribute(graph, attribute, index=which(V(graph)$name %in% node))
  if (length(a) > 0 && (attribute %in% "gos" || attribute %in% "gos.cc")){
    a <- strsplit(a, " ")[[1]][-which(strsplit(a, " ")[[1]] %in% "")]
  }
  return(a)
}

subgraphFromDataFrame <- function(graph, df){
  # given a df, find corresponding edges in graph
  # and create a subgraph from these edges
  apply(df, 1, function(x)
    get.edge.ids(graph, c(as.character(x[1]), as.character(x[2])))) -> ind
  unlist(ind) -> ind
  if (length(which(ind %in% 0))>0){
    ind[-which(ind %in% 0)] -> ind
  }
  graph.sub <- subgraph.edges(graph, as.numeric(ind), delete.vertices = TRUE)
  return(graph.sub)
}

constructGraph <- function(nn, ppis, p){
  # creates a graph from a list nn; each slot of nn contains 
  # a list of neighbors of some proteins p
  # length of nn[[i]] == length(p)
  # then, for each constructed subgraph, get
  # induce a subgraph from the original (so that attributes are
  # maintained)
  if (length(p) != length(nn)){
    stop("List of network neighborhoods (nn) and proteins (p) must have equal lengths")
  }
  
  lapply(1:length(nn), function(x){
    gdf <- matrix(0, nrow=0, ncol=2)
    colnames(gdf) <- c("source", "target")
    if (length(nn[[x]]) > 0){
      nn[[x]][which(nn[[x]] %ni% p[x])] -> fin
      cbind(rep(p[x], length(fin)), fin) -> t
      colnames(t) <- colnames(gdf)
      rbind(gdf, t) -> gdf
    }
    return(gdf)
  }) -> sg.df
  sg.df[which(sapply(sg.df, function(x) nrow(x)) > 0)] -> sg.df
  sg.df <- ldply (sg.df, data.frame)
  
  apply(sg.df, 1, function(y) get.edge.ids(ppis, y, directed = F)) -> eids
  subgraph.edges(ppis, eids, delete.vertices = T) -> sg
  
  return(sg)
}