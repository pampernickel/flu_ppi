##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB (Delorenzi/Stertz/Kunz)
##
##========================================================

`%ni%` <- Negate(`%in%`)

addGO <- function(graph, ind, name, attribute.name){
  a <- get.vertex.attribute(graph, attribute.name, index=ind)
  graph <- set.vertex.attribute(graph, attribute.name,
                                index=ind, value=paste(a, name, sep =" "))
  return(graph)
}

nodeToGO <- function(graph, xx.s, attribute.name){
  graph <- set.vertex.attribute(graph, attribute.name,
                                value=rep("", length(V(graph)$name)))
  lapply(xx.s, function(x) return(which(V(graph)$name %in% x))) -> go.inds
  
  # --- lapply would otherwise return a list of graphs for each go.inds element
  for (i in 1:length(go.inds)){
    addGO(graph, go.inds[[i]], names(go.inds)[i], attribute.name) -> graph
  }
  return(graph)
}