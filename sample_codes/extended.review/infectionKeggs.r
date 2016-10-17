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

# custom functions
source('./sample_codes/functions/retrievalFuncs.r')

# download all KEGG topologies linked to infection
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