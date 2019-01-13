#Call from BWA2ALL.Rmd
stack.fun <- function(d,..){
  d$pos5 <- ifelse(d$POS<d$PNEXT,d$POS,d$PNEXT)
  d$pos3 <- d$pos5 + abs(d$TLEN)
  d$start <- ifelse(d$TLEN>0,d$pos5,d$pos3)
  d$end <- ifelse(d$TLEN<0,d$pos5,d$pos3)
  
  d <- d[order(d$pos5),]
  
  stack <- 1
  bin.start <- 1
  d$stack <- 1
  d$bin.end <- d$pos3
  end.prv <- d$pos3[1]
  d$bstack <- 1
  bcount <- 1
  if(nrow(d)>1){
    for(i in 2:nrow(d)){
      pos.now <- d$pos5[i]
      if(pos.now < end.prv){
        stack = stack + 1
        d$stack[i] <- stack
        
      }else{
        d$bin.end[bin.start:(i-1)] <- end.prv
        bin.start <- i
        stack = 1
        d$stack[i] = 1
      }
      end.prv <- ifelse(d$pos3[i]>end.prv,d$pos3[i],end.prv)
      #barcode identity
      if(grepl("0000",as.character(d$barcode[i]))){
        d$bstack[i] <- 0
      }else if(bpos$bstack[which(bpos$b==as.character(d$barcode[i]))]==0){
        bcount <- bcount +1 
        bpos$bstack[which(bpos$b==as.character(d$barcode[i]))] <- bcount
      }
      d$bstack[i] <- bpos$bstack[which(bpos$b==as.character(d$barcode[i]))]
    }
    d$bin.end[bin.start:nrow(d)] <- end.prv
    
  }
  return(d)
}

stack.pair <- function(d,..){
  d$pos5 <- ifelse(d$POS<d$PNEXT,d$POS,d$PNEXT)
  d$pos3 <- d$pos5 + abs(d$TLEN)
  d$start <- ifelse(d$TLEN>0,d$pos5,d$pos3)
  d$end <- ifelse(d$TLEN<0,d$pos5,d$pos3)
  
  d <- d[order(d$pos5),]
  
  stack <- 1
  bin.start <- 1
  d$stack <- 1
  d$bin.end <- d$pos3
  end.prv <- d$pos3[1]
  d$bstack <- 1
  bcount <- 1
  if(nrow(d)>1){
    for(i in 2:nrow(d)){
      pos.now <- d$pos5[i]
      if(pos.now < end.prv){
        stack = stack + 1
        d$stack[i] <- stack
        
      }else{
        d$bin.end[bin.start:(i-1)] <- end.prv
        bin.start <- i
        stack = 1
        d$stack[i] = 1
      }
      end.prv <- ifelse(d$pos3[i]>end.prv,d$pos3[i],end.prv)
      #barcode identity
      if(grepl("0000",as.character(d$barcode[i]))){
        d$bstack[i] <- 0
      }else if(bpos$bstack[which(bpos$b==as.character(d$barcode[i]))]==0){
        bcount <- bcount + 1
        bpos$bstack[which(bpos$b==as.character(d$barcode[i]))] <- ceiling(bcount)
      }
      d$bstack[i] <- bpos$bstack[which(bpos$b==as.character(d$barcode[i]))]
    }
    d$bin.end[bin.start:nrow(d)] <- end.prv
    
  }
  return(d)
}

bstack.fun <- function(d){
  bpos <- data.frame(b=unique(as.character(d$barcode)),bstack=0)
  d$bstack <- 1
  bcount <- 1
  if(nrow(d)>1){
    for(i in 2:nrow(d)){
      d$bstack[i] <- bpos$bstack[which(bpos$b==as.character(d$barcode[i]))]
    }
  }
  return(d)
}

codeType.fun <- function(d){
  if(grepl("0000",as.character(d$barcode[1]))){
    d$codeType <- "undetected"
  }else{
    bnum <- nrow(d)
    d$codeType <- ifelse(bnum>0,"multiple","unique")
  }
  return(d)
}
