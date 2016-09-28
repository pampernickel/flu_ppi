# ---
# Selected code to demonstrate described results in
# `Finding functional modules in protein-protein interaction networks: highlights from 
# the virus-host interactome experience '
# c. MPDobay, 2016
# SystemsX Project, SIB Swiss Institute of Bioinformatics
# ---

library(STRINGdb)
library(igraph)
library(rentrez)
library(XML)
library(stringr)

# custom functions
source('./sample_codes/functions/retrievalFuncs.r')

# ---
# ---
# Code for retrieving original KEGG influenza references 
html <- getURL("http://www.genome.jp/dbget-bin/www_bget?hsa05164")
doc = htmlParse(html, asText=TRUE)
links <- as.character(xpathSApply(doc, "//a/@href"))[grep("pubmed", 
                                                          as.character(xpathSApply(doc, "//a/@href")))]
gsub("http://www.ncbi.nlm.nih.gov/pubmed/", "", links) -> links
paste("PMID:", links, sep="") -> links

# consider changing the delay, or downloading references by batch
xml.trees <- sapply(links, function(x) entrez_fetch(db="pubmed", id=x, rettype="xml", delay=3))
xml.trees <- lapply(xml.trees, function(x) xmlTreeParse(x, useInternalNodes=T))
abstract <- unlist(getAbstract(list(xml.trees)))
# ---
# ---


# ---
# ---
# Basic filtering example
# Example 1: Retain edges that contain a keyword with minimum frequency above the median
# this example uses a simple grep, without keyword stemming

# Retrieve IAV KEGG network
getKEGGgraph("05164") -> iav

# load STRINGdb; versions used in are 10 and 9_05
string.db <- STRINGdb$new(version="10", species=9606,score_threshold=0, input_directory="" )
string <- string.db$get_graph()
vertex_map <- string.db$map(prepareMap(V(iav)$name), "vertex", removeUnmappedRows = TRUE)
mapToSTRING(vertex_map, iav, mode="undirected") -> iavs

# note that iav contains complexes, which are NOT matched in STRING, but which nonetheless
# have a corresponding ID:
vcount(iavs) # all iav vertices
length(which(V(iavs)$name %in% V(string)$name)) # all iav vertices in STRING

# Get references linked to edges; consider putting a pause between
# call "get_pubmed_interaction" as the request might be blocked.
# Restrict query to all edges of iav where source and target vertices are
# in STRING 
get.data.frame(iavs, what="edges") -> df
intersect(which(df$to %in% V(string)$name),
          which(df$from %in% V(string)$name)) -> ind
df[ind,] -> df

# in this example, only the first 5 edges of iavs are considered
# abstract retrieval can be optimized by considering only unique
# abstracts, viz. unique(unlist(pmids)), then mapping these back to the
# edges
pmids <- apply(df[1:5,1:2], 1, function(x) 
              string.db$get_pubmed_interaction(x[1], x[2]))
xml.trees <- lapply(pmids, function(x) sapply(x[grep("PMID:",x)], function(y) 
    entrez_fetch(db="pubmed", id=y, rettype="xml", delay=3)))
xml.trees <- lapply(xml.trees, function(x) 
  sapply(x, function(y) xmlTreeParse(y, useInternalNodes=T)))
abstract <- getAbstract(xml.trees)

# single keyword filter, no filtering 
# in this example, only 2/5 edges are supported by edges
# that contain the keyword influenza
keyword <- "influenza"
sapply(abstract, function(x){
    unlist(x)-> x
    res <- F
    if (length(grep(keyword, x, ignore.case=T)) > 0){
      res <- T
    }
    return(res)
}) -> match
match

# combination filters
# in this example, 0/5 edges are supported by edges
# that contain the keyword influenza
keywords <- c("influenza", "endocytosis", "rab5")
combn(keywords, 2) -> keywords
apply(keywords, 2, function(y){
  as.character(sapply(abstract, function(x){
    unlist(x)-> x
    res <- F
    if (length(intersect(grep(y[1], x, ignore.case=T),
                  grep(y[2], x, ignore.case=T))) > 0){
      res <- T
    }
    return(res)
  })) -> match
}) -> match
match

# -------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #
# Addenda
# Requirements
# sessionInfo()
# R version 3.2.2 (2015-08-14)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu precise (12.04.5 LTS)
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] stringr_1.0.0          reshape2_1.4.1         stringdist_0.9.4.1     ggplot2_1.0.1          gridExtra_2.0.0       
# [6] tm_0.6-2               NLP_0.1-8              mygene_1.2.3           GenomicFeatures_1.20.6 GenomicRanges_1.20.6  
# [11] KEGG.db_3.1.2          org.Hs.eg.db_3.1.2     AnnotationDbi_1.30.1   GenomeInfoDb_1.4.3     IRanges_2.2.5         
# [16] S4Vectors_0.6.1        Biobase_2.28.0         BiocGenerics_0.14.0    ROntoTools_1.8.1       Rgraphviz_2.12.0      
# [21] KEGGgraph_1.26.0       KEGGREST_1.8.1         boot_1.3-17            graph_1.46.0           rentrez_1.0.0         
# [26] XML_3.98-1.3           STRINGdb_1.8.1         hash_2.2.6             gplots_2.17.0          RColorBrewer_1.1-2    
# [31] plotrix_3.6            RCurl_1.95-4.7         bitops_1.0-6           igraph_1.0.1           plyr_1.8.3            
# [36] sqldf_0.4-10           RSQLite_1.0.0          DBI_0.3.1              gsubfn_0.6-6           proto_0.3-10          
# [41] png_0.1-7             
# 
# loaded via a namespace (and not attached):
# [1] httr_1.0.0              jsonlite_0.9.19         splines_3.2.2           gtools_3.5.0           
# [5] Formula_1.2-1           latticeExtra_0.6-26     Rsamtools_1.20.5        slam_0.1-32            
# [9] lattice_0.20-33         chron_2.3-47            digest_0.6.9            XVector_0.8.0          
# [13] colorspace_1.2-6        biomaRt_2.24.1          zlibbioc_1.14.0         scales_0.3.0           
# [17] gdata_2.17.0            BiocParallel_1.2.22     nnet_7.3-11             survival_2.38-3        
# [21] magrittr_1.5            MASS_7.3-44             foreign_0.8-66          tools_3.2.2            
# [25] munsell_0.4.2           cluster_2.0.3           lambda.r_1.1.7          Biostrings_2.36.4      
# [29] caTools_1.17.1          futile.logger_1.4.1     tcltk_3.2.2             gtable_0.1.2           
# [33] curl_0.9.3              R6_2.1.1                GenomicAlignments_1.4.2 rtracklayer_1.28.10    
# [37] Hmisc_3.17-0            futile.options_1.0.0    KernSmooth_2.23-15      stringi_0.5-5          
# [41] Rcpp_0.11.6             rpart_4.1-10            acepack_1.3-3.3      
# -------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------- #