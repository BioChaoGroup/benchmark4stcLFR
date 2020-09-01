# library for JAN20 data analysis
#functions

#https://www.biostars.org/p/9335/
matcher <- function(pattern, x) {
  
  ind = gregexpr(pattern, x)[[1]]
  start = as.numeric(ind)
  if(length(start)==1 && start <1){
    return(NA)
  }else{
    end = start + attr(ind, "match.length") - 2
    res <- as.numeric(apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2])))
    res[which(is.na(res))] <- 1
    return(res)
  }
}

doone <- function(c, cigar) {
  pat <- paste("\\d*", c , sep="")
  sum(as.numeric(matcher(pat, cigar)), na.rm=T)
}


## takes a cigar string and parses it, not very fast but...
cigarsums <- function(cigar, chars=c("M","N","D","I","S","H", "P", "X", "=")) {
  sapply (chars, doone, cigar)
}


n50 <- function(x){
  x1 <- rev(sort(x))
  cumsum.x1 <- cumsum(x1)
  find.n50 <- which(cumsum.x1 >= sum(x1)/2)
  return(x1[find.n50[1]])
}

md5 <- function(f,m){
  if(tools::md5sum(f)==m){
    warning("MD5 check passed.\n")
    return(f)
  }else{
    stop("MD5 check failed. Please check the version and content!\n")
  }
}

getMapPropFun <- function(d){
  uOrder <- c("ALL","SSU","ITS1","5.8S","ITS2","LSU")
  ulist <- data.frame(key=c("ALL","SSU","16S","ITS","ITS1","5.8S","ITS2","LSU","23S"),
                      ind=c(1,2,2,3,3,4,5,6,6))
  getM <- stri_match_all_regex(d,"(ALL|23S|LSU|ITS|ITS1|ITS2|5.8S|SSU|16S)\\((\\d+)%\\)")
  getN <- length(getM)
  res <- matrix(0,ncol = 6,nrow = getN)
  for(i in 1:length(getM)){
    geti <- getM[[i]]
    for(j in 1:nrow(geti)){
      res[i,ulist$ind[which(ulist$key==geti[j,2])]] <- geti[j,3]
    }
  }
  res <- matrix(as.numeric(res),ncol=6)
  colnames(res) <- uOrder
  return(res)
}