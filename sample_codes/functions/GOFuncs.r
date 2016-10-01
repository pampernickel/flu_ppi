##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB (Delorenzi/Stertz/Kunz)
##
##========================================================

# global stuff
suppressPackageStartupMessages({library(GO.db, quietly=TRUE)})
suppressPackageStartupMessages({library(org.Hs.eg.db, quietly=TRUE)})
source('./scripts/routineFuncs.r')
xx <- as.list(org.Hs.egGO2ALLEGS)

`%ni%` <- Negate(`%in%`)

# --- retrieve gos
getgos <- function(){
  suppressPackageStartupMessages({require(igraph, quietly=TRUE)})
  # --- append GO terms to graph
  sql <- "SELECT DISTINCT go_id FROM go_bp"
  dbGetQuery(org.Hs.eg_dbconn(),sql) -> goss
  which(names(xx) %in% goss[,1]) -> ind
  xx[ind] -> xx
  entrez2symbol(xx) -> xx.s
  unique(unlist(xx.s)) -> all.ids
  list(xx, xx.s) -> gos
  names(gos) <- c("entrez", "symbol")
  return(gos)
}

getGO.cc <- function(){
  suppressPackageStartupMessages({require('org.Hs.eg.db', quietly=TRUE)})
  suppressPackageStartupMessages({require(igraph, quietly=TRUE)})
  # --- append GO terms to graph
  sql <- "SELECT DISTINCT go_id FROM go_cc"
  dbGetQuery(org.Hs.eg_dbconn(),sql) -> GO.ccs
  which(names(xx) %in% GO.ccs[,1]) -> ind
  xx[ind] -> xx
  entrez2symbol(xx) -> xx.s
  unique(unlist(xx.s)) -> all.ids
  list(xx, xx.s) -> gos
  names(gos) <- c("entrez", "symbol")
  return(gos)
}

getGOnames <- function(go.terms){
  goterms <- Term(GOTERM)
  return(goterms[which(names(goterms) %in% go.terms)])
}