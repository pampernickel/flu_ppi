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