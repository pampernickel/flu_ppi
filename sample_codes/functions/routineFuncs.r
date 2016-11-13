##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

`%ni%` <- Negate(`%in%`)

require("org.Hs.eg.db")
u2e <- toTable(org.Hs.egUNIPROT)
e2s <- toTable(org.Hs.egSYMBOL)

lappend <- function(lst, obj){
  if (is.list(obj)){
    for (i in 1:length(obj)){
      lst[[length(lst)+1]] <- obj[[i]]
    }
  } else {
    lst[[length(lst)+1]] <- obj
  }
  return(lst)
}

flatlist <- function(mylist){
  lapply(rapply(mylist, enquote, how="unlist"), eval)
}

entrez2symbol.org <- function(gene.id){
  # --- to symbols
  symbol <- NA
  if (length(gene.id) > 1){
    sapply(gene.id, function(x) return(
      ifelse(length(which(e2s$gene_id %in% x)) > 0, 
             e2s$symbol[which(e2s$gene_id %in% x)], NA))) -> symbol
    unlist(symbol) -> symbol
  } else {
    if(length(which(e2s$gene_id %in% gene.id)) > 0){
      e2s$symbol[which(e2s$gene_id %in% gene.id)] -> symbol
    }
  }
  return(symbol)
}

symbol2entrez <- function(gene.names){
  require("org.Hs.eg.db")
  e2s = toTable(org.Hs.egSYMBOL)
  gene.IDs <- rep(NA,length(gene.names))
  for (i in 1:length(gene.names)){
    if (length(which(e2s[,2] %in% gene.names[i])) > 0){
      gene.IDs[i] <- e2s[which(e2s[,2] %in% gene.names[i]),1]
    }
  }
  return(gene.IDs)
}

uniprot2entrez <- function(prot.names){
  # --- for prot.names, ensure that the isoforms are
  # --- also mapped
  substr(prot.names, 1, 6) -> prot.names
  gene.ids <- NA
  if (length(prot.names) > 1){
    sapply(prot.names, function(x) return(
      ifelse(length(which(u2e$uniprot_id %in% x)) > 0, u2e$gene_id[which(u2e$uniprot_id %in% x)], NA))) -> gene.ids
    unlist(gene.ids) -> gene.ids
  } else {
    if (length(which(u2e$uniprot_id %in% prot.names)) > 0){
      u2e$gene_id[which(u2e$uniprot_id %in% prot.names)] -> gene.ids
    }
  }
  
  return(gene.ids)
}

listToDf <- function(l){
  # convert list to data frame, with the names of the list as one column
  # and the content as the second column 
  # cumsum(sapply(l, function(x) length(x))) -> end.ind
  fin.df <- matrix(NA, nrow=0, ncol=2)
  colnames(fin.df) <- c("x","y")
  for (i in 1:length(l)){
    cbind(rep(names(l)[i], length(l[[i]])),l[[i]]) -> temp
    colnames(temp) <- c("x", "y")
    rbind(fin.df, temp) -> fin.df
  }
  return(fin.df)
}

source_https <- function(url, ...) {
  # load package
  # parse and evaluate each .R script
  # source: https://www.r-bloggers.com/source_https-sourcing-an-r-script-from-github-over-https/
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

loadDependencies <- function(){
  print("Loading scripts from pampernickel/flu_ppi...")
  source_https('https://raw.githubusercontent.com/pampernickel/flu_ppi/master/sample_codes/functions/graphFuncs.r')
  source_https('https://raw.githubusercontent.com/pampernickel/flu_ppi/master/sample_codes/functions/GOFuncs.r')
  source_https('https://raw.githubusercontent.com/pampernickel/flu_ppi/master/sample_codes/functions/graphAnnotationFuncs.r')
  source_https('https://raw.githubusercontent.com/pampernickel/flu_ppi/master/sample_codes/functions/retrievalFuncs.r')
  source_https('https://raw.githubusercontent.com/pampernickel/flu_ppi/master/sample_codes/functions/stringFuncs.r')
  print("Done.")
}