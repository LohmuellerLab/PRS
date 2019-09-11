
#library(data.table)
args = commandArgs(trailingOnly=T)
trait <- read.table(gzfile(args[1]), header=T)
#trait <- fread(args[1], header=T, data.table = F)
probs <- read.table("data/probs.txt", header=T, sep=",")
pp <- 0.46
trait$va <- 2*trait$minor_AF*(1-trait$minor_AF)*trait$beta^2

get_prob <- function(MAF){
  if (is.na(MAF)){
    return(5e-16)
  }
  if (MAF < 0.001){
    MAF = 0.001
  } else if (MAF == 0){
    return(5e-16)
  } 
  idx = round(MAF, 3)*1000
  pip = probs[idx,3]
  pi = probs[idx,2]
  if(is.na(pip) | is.na(pi)){
    return(5e-16)
  }
  if(pip == 2e-16 | pi == 2e-16){
    return(5e-16)
  }
  return(pip*pp/pi)
}

x <- sapply(trait$minor_AF, get_prob)

trait$prob <- x
private_va <- sum(subset(trait, trait$prob > 0.1)$va, na.rm=T)
total_va <- sum(trait$va, na.rm=T)
cat(private_va/total_va, "\t", private_va, "\t", total_va, "\n")
