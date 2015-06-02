# revised on 2-20-10
# - fix error in pass.structure: reverse rank.combined, so that big sig.value
#  are ranked with small numbers (1, 2, ...)
# - fix error on get.ez.tt.all: get ez.cutoff from sorted e.z

#
# modified EM procedure to compute empirical CDF more precisely - 09/2009



# this file contains the functions for  
# 1. computing the correspondence profile (upper rank intersection and derivatives)
# 2. inference of copula mixture model
#
# It also has functions for
# 1. reading peak caller results
# 2. processing and matching called peaks
# 3. plotting results


################ read peak caller results

# process narrow peak format
# some peak callers may not report q-values, p-values or fold of enrichment
# need further process before comparison
#
# stop.exclusive: Is the basepair of peak.list$stop exclusive? In narrowpeak and broadpeak format they are exclusive.
# If it is exclusive, we need subtract peak.list$stop by 1 to avoid the same basepair being both a start and a stop of two 
# adjacent peaks, which creates trouble for finding correct intersect  
process.narrowpeak <- function(narrow.file, chr.size, half.width=NULL, summit="offset", stop.exclusive=T, broadpeak=F){


  aa <- read.table(narrow.file)

  if(broadpeak){
    bb.ori <- data.frame(chr=aa$V1, start=aa$V2, stop=aa$V3, signal.value=aa$V7, p.value=aa$V8, q.value=aa$V9)
  }else{
    bb.ori <- data.frame(chr=aa$V1, start=aa$V2, stop=aa$V3, signal.value=aa$V7, p.value=aa$V8, q.value=aa$V9, summit=aa$V10)
  }
    
  if(summit=="summit"){
    bb.ori$summit <- bb.ori$summit-bb.ori$start # change summit to offset to avoid error when concatenating chromosomes
  }
 
  bb <- concatenate.chr(bb.ori, chr.size)

  #bb <- bb.ori

  # remove the peaks that has the same start and stop value
  bb <- bb[bb$start != bb$stop,]

  if(stop.exclusive==T){
    bb$stop <- bb$stop-1
  }

  if(!is.null(half.width)){
    bb$start.ori <- bb$start    #Anshul changed this
    bb$stop.ori <- bb$stop 	#Anshul changed this

    # if peak is narrower than the specified window, stay with its width
    # otherwise chop wider peaks to specified width
    width <- bb$stop-bb$start +1
    is.wider <- width > 2*half.width

    if(summit=="offset" | summit=="summit"){ # if summit is offset from start
      bb$start[is.wider] <- bb$start.ori[is.wider] + bb$summit[is.wider]-half.width
      bb$stop[is.wider] <- bb$start.ori[is.wider] + bb$summit[is.wider]+half.width
    } else { 
      if(summit=="unknown"){
        bb$start[is.wider] <- bb$start.ori[is.wider]+round(width[is.wider]/2) - half.width
        bb$stop[is.wider] <- bb$start.ori[is.wider]+round(width[is.wider]/2) + half.width
      }
    }

    bb$start.ori <- bb.ori$start    #Anshul changed this
    bb$stop.ori <- bb.ori$stop      #Anshul changed this
  }

  bb <- clean.data(bb)
  invisible(list(data.ori=bb.ori, data.cleaned=bb))
}

# clean data 
# and concatenate chromosomes if needed
clean.data <- function(adata){

  # remove the peaks that has the same start and stop value
  adata <- adata[adata$start != adata$stop,]

  # if some stops and starts are the same, need fix them
  stop.in.start <- is.element(adata$stop, adata$start)
  n.fix <- sum(stop.in.start)
  if(n.fix >0){
    print(paste("Fix", n.fix, "stops\n"))
    adata$stop[stop.in.start] <- adata$stop[stop.in.start]-1 
  }  
 
  return(adata) 
}

# concatenate peaks
# peaks: the dataframe to have all the peaks
# chr.file: the file to keep the length of each chromosome 
# chr files should come from the species that the data is from
#concatenate.chr <- function(peaks, chr.size){

 # chr.size <- read.table(chr.file)
#  chr.o <- order(chr.size[,1])
#  chr.size <- chr.size[chr.o,]
#
#  chr.shift <- cumsum(c(0, chr.size[-nrow(chr.size),2]))
#  chr.size.cum <- data.frame(chr=chr.size[,1], shift=chr.shift)  
#
#  for(i in 1:nrow(chr.size)){
#    is.in <- as.character(peaks$chr) == as.character(chr.size.cum$chr[i])
#    if(sum(is.in)>0){
#      peaks[is.in,]$start <- peaks[is.in,]$start + chr.size.cum$shift[i]
#      peaks[is.in,]$stop <- peaks[is.in,]$stop + chr.size.cum$shift[i]
#    }
#  }
#
#  invisible(peaks)
#}




# concatenate peaks
# peaks: the dataframe to have all the peaks
# chr.file: the file to keep the length of each chromosome 
# chr files should come from the species that the data is from
concatenate.chr <- function(peaks, chr.size){

 # chr.size <- read.table(chr.file)
  chr.o <- order(chr.size[,1])
  chr.size <- chr.size[chr.o,]

  chr.shift <- cumsum(c(0, chr.size[-nrow(chr.size),2]))
  chr.size.cum <- data.frame(chr=chr.size[,1], shift=chr.shift)  

  peaks$start.ori <- peaks$start
  peaks$stop.ori <- peaks$stop
  
  for(i in 1:nrow(chr.size)){
    is.in <- as.character(peaks$chr) == as.character(chr.size.cum$chr[i])
    if(sum(is.in)>0){
      peaks[is.in,]$start <- peaks[is.in,]$start + chr.size.cum$shift[i]
      peaks[is.in,]$stop <- peaks[is.in,]$stop + chr.size.cum$shift[i]
    }
  }

  invisible(peaks)
}


deconcatenate.chr <- function(peaks, chr.size){

  chr.o <- order(chr.size[,1])
  chr.size <- chr.size[chr.o,]

  chr.shift <- cumsum(c(0, chr.size[-nrow(chr.size),2]))
  chr.size.cum <- data.frame(chr=chr.size[,1], shift=chr.shift)  

  peaks$chr <- rep(NA, nrow(peaks))
  
  for(i in 1:(nrow(chr.size.cum)-1)){
    is.in <- peaks$start > chr.size.cum[i,2] & peaks$start <= chr.size.cum[i+1, 2]
    if(sum(is.in)>0){
      peaks[is.in,]$start <- peaks[is.in,]$start - chr.size.cum[i,2]
      peaks[is.in,]$stop <- peaks[is.in,]$stop - chr.size.cum[i,2]+1    
      peaks[is.in,]$chr <- chr.size[i,1]
    }
  }

  if(i == nrow(chr.size.cum)){
    is.in <- peaks$start > chr.size.cum[i, 2]
    if(sum(is.in)>0){
      peaks[is.in,]$start <- peaks[is.in,]$start - chr.size.cum[i,2]
      peaks[is.in,]$stop <- peaks[is.in,]$stop - chr.size.cum[i,2]+1    
      peaks[is.in,]$chr <- chr.size[i,1]
    }
  }
  
  invisible(peaks)
}

################ preprocessing peak calling output


# 
# read two calling results and sort by peak starting locations, 
# then find overlap between peaks
# INPUT:
#   rep1: the 1st replicate
#   rep2: the 2nd replicate
# OUTPUT:
#   id1, id2: the labels for the identified peaks on the replicates
find.overlap <- function(rep1, rep2){

  o1 <- order(rep1$start)
  rep1 <- rep1[o1,]
    
  o2 <- order(rep2$start)
  rep2 <- rep2[o2,]

  n1 <- length(o1)
  n2 <- length(o2)
  
  # assign common ID to peaks
  id1 <- rep(0, n1) # ID assigned on rep1
  id2 <- rep(0, n2) # ID assigned on rep2
  id <- 1 # keep track common id's
  
  # check if two replicates overlap with each other
  i <- 1
  j <- 1

  while(i <= n1|| j <= n2){

    # && (id1[n1] ==0 || id2[n2] ==0)
    
    # if one list runs out
    if(i > n1 && j < n2){
      
      j <- j+1
      id2[j] <- id
      id <- id +1
      next
    } else{
      if(j > n2 && i < n1){
        i <- i+1        
        id1[i] <- id
        id <- id +1
        next
      } else {
        if(i >= n1 && j >=n2)
          break
      }
    }

    # if not overlap

    if(!(rep1$start[i] <= rep2$stop[j] && rep2$start[j] <= rep1$stop[i])){

      # at the start of loop, when both are not assigned an ID
      # the one locates in front is assigned first
      if(id1[i] ==0 && id2[j]==0){
        if(rep1$stop[i] < rep2$stop[j]){
          id1[i] <- id
        } else {
          id2[j] <- id
        }
      } else { # in the middle of the loop, when one is already assigned
      # The one that has not assigned gets assigned
      #  if(id1[i] ==0){ # id1[i] is not assigned
      #    id1[i] <- id
      #  } else { # id2[i] is not assigned
      #    id2[j] <- id 
      #  }

        # order the id according to location
        if(rep1$stop[i] <= rep2$stop[j]){
          id1[i] <- max(id2[j], id1[i])
          id2[j] <- id  
        } else {
          if(rep1$stop[i] > rep2$stop[j]){
            id2[j] <- max(id1[i], id2[j])
            id1[i] <- id
          }
        }
          
      }
      
      id <- id +1
      
    } else { # if overlap
    
      if(id1[i] == 0 && id2[j] == 0){ # not assign label yet
        id1[i] <- id 
        id2[j] <- id
        id <- id +1
      } else { # one peak is already assigned label, the other is 0
        
        id1[i] <- max(id1[i], id2[j]) # this is a way to copy the label of the assigned peak without knowing which one is already assigned
        id2[j] <- id1[i] # syncronize the labels        
      }
      
    }
    
    if(rep1$stop[i] < rep2$stop[j]){
      i <- i+1
    } else {
      j <- j+1
    } 
    
  }

  invisible(list(id1=id1, id2=id2))
  
}

# Impute the missing significant value for the peaks called only on one replicate.
# value 
# INPUT:
#   rep1, rep2: the two peak calling output 
#   id1, id2: the IDs assigned by function find.overlap, vectors
#        If id1[i]==id2[j], peak i on rep1 overlaps with peak j on rep2
#   p.value.impute: the significant value to impute for the missing peaks 
# OUTPUT:   
#   rep1, rep2: peaks ordered by the start locations with imputed peaks
#   id1, id2: the IDs with imputed peaks
fill.missing.peaks <- function(rep1, rep2, id1, id2, p.value.impute){

#   rep1 <- data.frame(chr=rep1$chr, start=rep1$start, stop=rep1$stop, sig.value=rep1$sig.value)
#   rep2 <- data.frame(chr=rep2$chr, start=rep2$start, stop=rep2$stop, sig.value=rep2$sig.value)   
   
   o1 <- order(rep1$start)
   rep1 <- rep1[o1,]
    
   o2 <- order(rep2$start)
   rep2 <- rep2[o2,]  
     
   entry.in1.not2 <- !is.element(id1, id2)
   entry.in2.not1 <- !is.element(id2, id1)

   if(sum(entry.in1.not2) > 0){
   
     temp1 <- rep1[entry.in1.not2, ]

     # impute sig.value
     temp1$sig.value <- p.value.impute
     temp1$signal.value <- p.value.impute
     temp1$p.value <- p.value.impute
     temp1$q.value <- p.value.impute
     
     rep2.filled <- rbind(rep2, temp1)
     id2.filled <- c(id2, id1[entry.in1.not2])
   } else {
     id2.filled <- id2
     rep2.filled <- rep2
   }

   if(sum(entry.in2.not1) > 0){

     temp2 <- rep2[entry.in2.not1, ]

     # fill in p.values to 1
     temp2$sig.value <- p.value.impute
     temp2$signal.value <- p.value.impute
     temp2$p.value <- p.value.impute
     temp2$q.value <- p.value.impute
   

     # append to the end
     rep1.filled <- rbind(rep1, temp2)

     id1.filled <- c(id1, id2[entry.in2.not1])
   } else {
     id1.filled <- id1
     rep1.filled <- rep1
   }

   # sort rep1 and rep2 by the same id
   o1 <- order(id1.filled)
   rep1.ordered <- rep1.filled[o1, ]

   o2 <- order(id2.filled)
   rep2.ordered <- rep2.filled[o2, ]   
   
   invisible(list(rep1=rep1.ordered, rep2=rep2.ordered,
                  id1=id1.filled[o1], id2=id2.filled[o2]))
 }

# Merge peaks with same ID on the same replicates 
# (They are generated if two peaks on rep1 map to the same peak on rep2)
# need peak.list have 3 columns: start, stop and sig.value 
merge.peaks.best <- function(peak.list, id){

  i <- 1
  j <- 1
  dup.index <- c()
  sig.value <- c()
  start.new <- c()
  stop.new <- c()
  id.new <- c()

  # original data
  chr <- c()
  start.ori <- c()
  stop.ori <- c()
  
  signal.value <- c()
  p.value <- c()
  q.value <- c()

  while(i < length(id)){
    
    if(id[i] == id[i+1]){
      dup.index <- c(dup.index, i, i+1) # push on dup.index
    } else {
      if(length(dup.index)>0){ # pop from dup.index        
    #    sig.value[j] <- mean(peak.list$sig.value[unique(dup.index)]) # mean of -log(pvalue)
        sig.value[j] <- max(peak.list$sig.value[unique(dup.index)])
        start.new[j] <- peak.list$start[min(dup.index)]
        stop.new[j] <- peak.list$stop[max(dup.index)]
        id.new[j] <- id[max(dup.index)]
        
    #    signal.value[j] <- mean(peak.list$signal.value[unique(dup.index)])        #    p.value[j] <- mean(peak.list$p.value[unique(dup.index)]) # mean of -log(pvalue)
    #    q.value[j] <- mean(peak.list$q.value[unique(dup.index)]) # mean of -log(pvalue)     
        signal.value[j] <- max(peak.list$signal.value[unique(dup.index)]) 
        p.value[j] <- max(peak.list$p.value[unique(dup.index)]) 
        q.value[j] <- max(peak.list$q.value[unique(dup.index)]) 

        chr[j] <- as.character(peak.list$chr[min(dup.index)])
        start.ori[j] <- peak.list$start.ori[min(dup.index)]
        stop.ori[j] <- peak.list$stop.ori[max(dup.index)]
        
        dup.index <- c()
      } else { # nothing to pop
        sig.value[j] <- peak.list$sig.value[i]
        start.new[j] <- peak.list$start[i]
        stop.new[j] <- peak.list$stop[i]
        id.new[j] <- id[i]

        signal.value[j] <- peak.list$signal.value[i] 
        p.value[j] <- peak.list$p.value[i] 
        q.value[j] <- peak.list$q.value[i] 

        chr[j] <- as.character(peak.list$chr[i])
        start.ori[j] <- peak.list$start.ori[i]
        stop.ori[j] <- peak.list$stop.ori[i]
        
      }
      j <- j+1
    }
    i <- i+1
  }

  data.new <- data.frame(id=id.new, sig.value=sig.value, start=start.new, stop=stop.new, signal.value=signal.value, p.value=p.value, q.value=q.value, chr=chr, start.ori=start.ori, stop.ori=stop.ori)
  invisible(data.new)
}

# Merge peaks with same ID on the same replicates 
# (They are generated if two peaks on rep1 map to the same peak on rep2)
# need peak.list have 3 columns: start, stop and sig.value 
merge.peaks <- function(peak.list, id){

  i <- 1
  j <- 1
  dup.index <- c()
  sig.value <- c()
  start.new <- c()
  stop.new <- c()
  id.new <- c()

  # original data
  chr <- c()
  start.ori <- c()
  stop.ori <- c()
  
  signal.value <- c()
  p.value <- c()
  q.value <- c()

  while(i < length(id)){
    
    if(id[i] == id[i+1]){
      dup.index <- c(dup.index, i, i+1) # push on dup.index
    } else {
      if(length(dup.index)>0){ # pop from dup.index
        sig.value[j] <- mean(peak.list$sig.value[unique(dup.index)]) # mean of -log(pvalue)
        start.new[j] <- peak.list$start[min(dup.index)]
        stop.new[j] <- peak.list$stop[max(dup.index)]
        id.new[j] <- id[max(dup.index)]
        
        signal.value[j] <- mean(peak.list$signal.value[unique(dup.index)]) # mean of -log(pvalue)
        p.value[j] <- mean(peak.list$p.value[unique(dup.index)]) # mean of -log(pvalue)
        q.value[j] <- mean(peak.list$q.value[unique(dup.index)]) # mean of -log(pvalue)

        chr[j] <- as.character(peak.list$chr[min(dup.index)])
        start.ori[j] <- peak.list$start.ori[min(dup.index)]
        stop.ori[j] <- peak.list$stop.ori[max(dup.index)]
        
        dup.index <- c()
      } else { # nothing to pop
        sig.value[j] <- peak.list$sig.value[i]
        start.new[j] <- peak.list$start[i]
        stop.new[j] <- peak.list$stop[i]
        id.new[j] <- id[i]

        signal.value[j] <- peak.list$signal.value[i] 
        p.value[j] <- peak.list$p.value[i] 
        q.value[j] <- peak.list$q.value[i] 

        chr[j] <- as.character(peak.list$chr[i])
        start.ori[j] <- peak.list$start.ori[i]
        stop.ori[j] <- peak.list$stop.ori[i]
        
      }
      j <- j+1
    }
    i <- i+1
  }

  data.new <- data.frame(id=id.new, sig.value=sig.value, start=start.new, stop=stop.new, signal.value=signal.value, p.value=p.value, q.value=q.value, chr=chr, start.ori=start.ori, stop.ori=stop.ori)
  invisible(data.new)
}





# a wrap function to fill in missing peaks, merge peaks and impute significant values
# out1 and out2 are two peak calling outputs
pair.peaks <- function(out1, out2, p.value.impute=0){

  aa <- find.overlap(out1, out2)
  bb <- fill.missing.peaks(out1, out2, aa$id1, aa$id2, p.value.impute=0)

  cc1 <- merge.peaks(bb$rep1, bb$id1)
  cc2 <- merge.peaks(bb$rep2, bb$id2)

  invisible(list(merge1=cc1, merge2=cc2))
}



# overlap.ratio is a parameter to define the percentage of overlap
# if overlap.ratio =0, 1 basepair overlap is counted as overlap
# if overlap.ratio between 0 and 1, it is the minimum proportion of
# overlap required to be called as a match
# it is computed as the overlap part/min(peak1.length, peak2.length)
pair.peaks.filter <- function(out1, out2, p.value.impute=0, overlap.ratio=0){

  aa <- find.overlap(out1, out2)
  bb <- fill.missing.peaks(out1, out2, aa$id1, aa$id2, p.value.impute=0)

  cc1 <- merge.peaks(bb$rep1, bb$id1)
  cc2 <- merge.peaks(bb$rep2, bb$id2)

  frag12 <- cbind(cc1$start, cc1$stop, cc2$start, cc2$stop)
  
  frag.ratio <- apply(frag12, 1, overlap.middle)

  frag.ratio[cc1$sig.value==p.value.impute | cc2$sig.value==p.value.impute] <- 0

  cc1$frag.ratio <- frag.ratio
  cc2$frag.ratio <- frag.ratio

  merge1 <- cc1[cc1$frag.ratio >= overlap.ratio,]
  merge2 <- cc2[cc2$frag.ratio >= overlap.ratio,]
  
  invisible(list(merge1=merge1, merge2=merge2))
}

# x[1], x[2] are the start and end of the first fragment
# and x[3] and x[4] are the start and end of the 2nd fragment 
# If there are two fragments, we can find the overlap by ordering the
# start and stop of all the ends and find the difference between the middle two
overlap.middle  <- function(x){

  x.o <- x[order(x)]
  f1 <- x[2]-x[1]
  f2 <- x[4]-x[3]
  
  f.overlap <- abs(x.o[3]-x.o[2])
  f.overlap.ratio <- f.overlap/min(f1, f2)

  return(f.overlap.ratio)
}



#######
####### compute correspondence profile
#######

# compute upper rank intersection for one t
# tv: the upper percentile
# x is sorted by the order of paired variable
comp.uri <- function(tv, x){
  n <- length(x)
  qt <- quantile(x, prob=1-tv[1]) # tv[1] is t
#  sum(x[1:ceiling(n*tv[2])] >= qt)/n/tv[2]- tv[1]*tv[2] #tv[2] is v
  sum(x[1:ceiling(n*tv[2])] >= qt)/n

}

# compute the correspondence profile
# tt, vv: vector between (0, 1) for percentages
get.uri.2d <- function(x1, x2, tt, vv, spline.df=NULL){

  o <- order(x1, x2, decreasing=T)
  
  # sort x2 by the order of x1
  x2.ordered <- x2[o]
  
  tv <- cbind(tt, vv)
  ntotal <- length(x1) # number of peaks    

  uri <- apply(tv, 1, comp.uri, x=x2.ordered)

  # compute the derivative of URI vs t using small bins
  uri.binned <- uri[seq(1, length(uri), by=4)]
  tt.binned <- tt[seq(1, length(uri), by=4)]
  uri.slope <- (uri.binned[2:(length(uri.binned))] - uri.binned[1:(length(uri.binned)-1)])/(tt.binned[2:(length(uri.binned))] - tt.binned[1:(length(tt.binned)-1)])

  # smooth uri using spline
  # first find where the jump is and don't fit the jump
  # this is the index on the left
  # jump.left.old  <- which.max(uri[-1]-uri[-length(uri)])
  short.list.length <- min(sum(x1>0)/length(x1), sum(x2>0)/length(x2))

  if(short.list.length < max(tt)){
    jump.left <- which(tt>short.list.length)[1]-1
  } else {
    jump.left <- which.max(tt)
  }

#  reversed.index <- seq(length(tt), 1, by=-1)
#  nequal <- sum(uri[reversed.index]== tt[reversed.index])
#  temp  <- which(uri[reversed.index]== tt[reversed.index])[nequal]
#  jump.left <- length(tt)-temp
 
  if(jump.left < 6){
   jump.left <- length(tt)
  }
    
 
  if(is.null(spline.df))
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=6.4)
  else{
    uri.spl <- smooth.spline(tt[1:jump.left], uri[1:jump.left], df=spline.df)
  }
  # predict the first derivative
  uri.der <- predict(uri.spl, tt[1:jump.left], deriv=1)

  invisible(list(tv=tv, uri=uri, 
                 uri.slope=uri.slope, t.binned=tt.binned[2:length(uri.binned)], 
                 uri.spl=uri.spl, uri.der=uri.der, jump.left=jump.left,
                 ntotal=ntotal))
 }


# change the scale of uri from based on t (percentage) to n (number of peaks or basepairs)
# this is for plotting multiple pairwise URI's on the same plot 
scale.t2n <- function(uri){

  ntotal <- uri$ntotal
  tv <- uri$tv*uri$ntotal
  uri.uri <- uri$uri*uri$ntotal
  jump.left <- uri$jump.left
  uri.spl <- uri$uri.spl
  uri.spl$x <- uri$uri.spl$x*uri$ntotal 
  uri.spl$y <- uri$uri.spl$y*uri$ntotal

  t.binned <- uri$t.binned*uri$ntotal
  uri.slope <- uri$uri.slope
  uri.der <- uri$uri.der
  uri.der$x <- uri$uri.der$x*uri$ntotal
  uri.der$y <- uri$uri.der$y

  uri.n <- list(tv=tv, uri=uri.uri, t.binned=t.binned, uri.slope=uri.slope, uri.spl=uri.spl, uri.der=uri.der, ntotal=ntotal, jump.left=jump.left)
  return(uri.n)
} 




# a wrapper for running URI for peaks from peak calling results
# both data1 and data2 are calling results in narrowpeak format
compute.pair.uri <- function(data.1, data.2, sig.value1="signal.value", sig.value2="signal.value", spline.df=NULL, overlap.ratio=0){

  tt <- seq(0.01, 1, by=0.01)
  vv <- tt

  if(sig.value1=="signal.value"){
    data.1.enrich <- data.frame(chr=data.1$chr, start.ori=data.1$start.ori, stop.ori=data.1$stop.ori, start=data.1$start, stop=data.1$stop, sig.value=data.1$signal.value, signal.value=data.1$signal.value, p.value=data.1$p.value, q.value=data.1$q.value)
  } else {
    if(sig.value1=="p.value"){ 
      data.1.enrich <- data.frame(chr=data.1$chr, start.ori=data.1$start.ori, stop.ori=data.1$stop.ori, start=data.1$start, stop=data.1$stop, sig.value=data.1$p.value, signal.value=data.1$signal.value, p.value=data.1$p.value, q.value=data.1$q.value)
    } else {
      if(sig.value1=="q.value"){
        data.1.enrich <- data.frame(chr=data.1$chr, start.ori=data.1$start.ori, stop.ori=data.1$stop.ori, start=data.1$start, stop=data.1$stop, sig.value=data.1$q.value, signal.value=data.1$signal.value, p.value=data.1$p.value, q.value=data.1$q.value)
      }
    }
  }

  if(sig.value2=="signal.value"){
    data.2.enrich <- data.frame(chr=data.2$chr, start.ori=data.2$start.ori, stop.ori=data.2$stop.ori, start=data.2$start, stop=data.2$stop, sig.value=data.2$signal.value, signal.value=data.2$signal.value, p.value=data.2$p.value, q.value=data.2$q.value)
  } else {
    if(sig.value2=="p.value"){ 
      data.2.enrich <- data.frame(chr=data.2$chr, start.ori=data.2$start.ori, stop.ori=data.2$stop.ori, start=data.2$start, stop=data.2$stop, sig.value=data.2$p.value, signal.value=data.2$signal.value, p.value=data.2$p.value, q.value=data.2$q.value)
    } else {
      if(sig.value2=="q.value"){
        data.2.enrich <- data.frame(chr=data.2$chr, start.ori=data.2$start.ori, stop.ori=data.2$stop.ori, start=data.2$start, stop=data.2$stop, sig.value=data.2$q.value, signal.value=data.2$signal.value, p.value=data.2$p.value, q.value=data.2$q.value)
      }
    }
  }

  ### by peaks
  # data12.enrich <- pair.peaks(data.1.enrich, data.2.enrich)
  data12.enrich <- pair.peaks.filter(data.1.enrich, data.2.enrich, p.value.impute=0, overlap.ratio)
  uri <- get.uri.2d(as.numeric(as.character(data12.enrich$merge1$sig.value)), as.numeric(as.character(data12.enrich$merge2$sig.value)), tt, vv, spline.df=spline.df)
  uri.n <- scale.t2n(uri)

  return(list(uri=uri, uri.n=uri.n, data12.enrich=data12.enrich, sig.value1=sig.value1, sig.value2=sig.value2))


}



# compute uri for matched sample
get.uri.matched <- function(data12, df=10){

  tt <- seq(0.01, 1, by=0.01)
  vv <- tt
  uri <- get.uri.2d(data12$sample1$sig.value, data12$sample2$sig.value, tt, vv, spline.df=df)

  # change scale from t to n
  uri.n <- scale.t2n(uri)

  return(list(uri=uri, uri.n=uri.n))
  
}

# map.uv is a pair of significant values corresponding to specified consistency FDR
# assuming values in map.uv and qvalue are linearly related
# data.set is the original data set
# sig.value is the name of the significant value in map.uv, say enrichment
# nominal.value is the one we want to map to, say q-value
# 
map.sig.value <- function(data.set, map.uv, nominal.value){

  index.nominal <- which(names(data.set$merge1)==nominal.value)
  nentry <- nrow(map.uv)  
  map.nominal <- rbind(map.uv[, c("sig.value1", "sig.value2")])

  for(i in 1:nentry){

    map.nominal[i, "sig.value1"] <- data.set$merge1[unique(which.min(abs(data.set$merge1$sig.value-map.uv[i, "sig.value1"]))), index.nominal]
    map.nominal[i, "sig.value2"] <- data.set$merge2[unique(which.min(abs(data.set$merge2$sig.value-map.uv[i, "sig.value2"]))), index.nominal]
  }

  invisible(map.nominal)
}


############### plot correspondence profile

# plot multiple comparison wrt one template
# uri.list contains the total number of peaks
# plot.missing=F: not plot the missing points on the right 
plot.uri.group <- function(uri.n.list, plot.dir, file.name=NULL, legend.txt, xlab.txt="num of significant peaks", ylab.txt="num of peaks in common", col.start=0, col.txt=NULL, plot.missing=F, title.txt=NULL){

  if(is.null(col.txt))
    col.txt <- c("black", "red", "purple", "green", "blue", "cyan", "magenta", "orange", "grey")

  n <- length(uri.n.list)  

  ntotal <- c()
  for(i in 1:n)
    ntotal[i] <- uri.n.list[[i]]$ntotal

  jump.left <- c()
  jump.left.der <- c()
  ncommon <- c()
  for(i in 1:n){
#    jump.left[i]  <- which.max(uri.n.list[[i]]$uri[-1]-uri.n.list[[i]]$uri[-length(uri.n.list[[i]]$uri)])
#    if(jump.left[i] < 6)
#      jump.left[i] <- length(uri.n.list[[i]]$uri)

##  reversed.index <- seq(length(uri.n.list[[i]]$tv[,1]), 1, by=-1)
##  nequal <- sum(uri.n.list[[i]]$uri[reversed.index]== uri.n.list[[i]]$tv[reversed.index,1])
##  temp  <- which(uri.n.list[[i]]$uri[reversed.index]== uri.n.list[[i]]$tv[reversed.index,1])[nequal]
##  jump.left[i] <- length(uri.n.list[[i]]$tv[,1])-temp
##print(uri.n.list[[i]]$uri)
##print(uri.n.list[[i]]$tv[,1])
##   jump.left[i] <- uri.n.list[[i]]$jump.left

#    jump.left.der[i] <- sum(uri.n.list[[i]]$t.binned < uri.n.list[[i]]$uri.der$x[length(uri.n.list[[i]]$uri.der$x)])

    jump.left[i] <- uri.n.list[[i]]$jump.left
    jump.left.der[i] <- jump.left[i]
    ncommon[i] <- uri.n.list[[i]]$tv[jump.left[i],1]
  }


  if(plot.missing){
    max.peak <- max(ntotal)
  } else {
    max.peak <- max(ncommon)*1.05
  }

  if(!is.null(file.name)){
    postscript(paste(plot.dir, "uri.", file.name, sep=""))
    par(mfrow=c(1,1), mar=c(5,5,4,2))
  }

  plot(uri.n.list[[1]]$tv[,1], uri.n.list[[1]]$uri, type="n", xlab=xlab.txt, ylab=ylab.txt, xlim=c(0, max.peak), ylim=c(0, max.peak), cex.lab=2)

  for(i in 1:n){

    if(plot.missing){ 
      points(uri.n.list[[i]]$tv[,1], uri.n.list[[i]]$uri, col=col.txt[i+col.start], cex=0.5 )
    } else {
      points(uri.n.list[[i]]$tv[1:jump.left[i],1], uri.n.list[[i]]$uri[1:jump.left[i]], col=col.txt[i+col.start], cex=0.5)
    }
    lines(uri.n.list[[i]]$uri.spl, col=col.txt[i+col.start], lwd=4)
  }
  abline(coef=c(0,1), lty=3)
  legend(0, max.peak, legend=legend.txt, col=col.txt[(col.start+1):length(col.txt)], lty=1, lwd=3, cex=2)
  if(!is.null(title))
    title(title.txt)

  if(!is.null(file.name)){
    dev.off()
  }

  if(!is.null(file.name)){
    postscript(paste(plot.dir, "duri.", file.name, sep=""))
    par(mfrow=c(1,1), mar=c(5,5,4,2))
  } 
  plot(uri.n.list[[1]]$t.binned, uri.n.list[[1]]$uri.slope, type="n", xlab=xlab.txt, ylab="slope", xlim=c(0, max.peak), ylim=c(0, 1.5), cex.lab=2)

  for(i in 1:n){
#    if(plot.missing){ 
#      points(uri.n.list[[i]]$t.binned, uri.n.list[[i]]$uri.slope, col=col.txt[i+col.start], cex=0.5)
#    } else {
#      points(uri.n.list[[i]]$t.binned[1:jump.left.der[i]], uri.n.list[[i]]$uri.slope[1:jump.left.der[i]], col=col.txt[i+col.start], cex=0.5)
#    }
    lines(uri.n.list[[i]]$uri.der, col=col.txt[i+col.start], lwd=4)
  }
  abline(h=1, lty=3)
  legend(0.5*max.peak, 1.5, legend=legend.txt, col=col.txt[(col.start+1):length(col.txt)], lty=1, lwd=3, cex=2)

  if(!is.null(title))
    title(title.txt)

  if(!is.null(file.name)){
    dev.off()
  }
  
}



#######################
####################### copula fitting for matched peaks
#######################

# estimation from mixed copula model 

# 4-5-09
# A nonparametric estimation of mixed copula model


# updated

# c1, c2, f1, f2, g1, g2 are vectors
# c1*f1*g1 and c2*f2*g2 are copula densities for the two components
# xd1 and yd1 are the values of marginals for the first component
# xd2 and yd2 are the values of marginals for the 2nd component
#
# ez is the prob for being in the consistent group
get.ez <- function(p, c1, c2, xd1, yd1, xd2, yd2){

  return(p*c1*xd1*yd1/(p*c1*xd1*yd1 + (1-p)*c2*xd2*yd2))
}

# checked

# this is C_12 not the copula density function c=C_12 * f1* f2
# since nonparametric estimation is used here for f1 and f2, which
# are constant throughout the iterations, we don't need them for optimization
# 
# bivariate gaussian copula function
# t and s are vectors of same length, both are percentiles 
# return a vector
gaussian.cop.den <- function(t, s, rho){

  A <- qnorm(t)^2 + qnorm(s)^2
  B <- qnorm(t)*qnorm(s)

  loglik <-  -log(1-rho^2)/2 - rho/(2*(1-rho^2))*(rho*A-2*B)

  return(exp(loglik))
}

clayton.cop.den <- function(t, s, rho){

  if(rho > 0)
    return(exp(log(rho+1)-(rho+1)*(log(t)+log(s))-(2+1/rho)*log(t^(-rho) + s^(-rho)-1)))

  if(rho==0)
    return(1)

  if(rho<0)
    stop("Incorrect Clayton copula coefficient")
  
}


# checked
# estimate rho from Gaussian copula
mle.gaussian.copula <- function(t, s, e.z){

  # reparameterize to bound from rho=+-1
  l.c <- function(rho, t, s, e.z){
#    cat("rho=", rho, "\n")
    sum(e.z*log(gaussian.cop.den(t, s, rho)))}

  rho.max <- optimize(f=l.c, c(-0.998, 0.998), maximum=T, tol=0.00001, t=t, s=s, e.z=e.z)

#print(rho.max$m)

#cat("cor=", cor(qnorm(t)*e.z, qnorm(s)*e.z), "\t", "rho.max=", rho.max$m, "\n")
#  return(sign(rho.max$m)/(1+rho.max$m))
  return(rho.max$m)
}


# estimate mle from Clayton copula, 
mle.clayton.copula <- function(t, s, e.z){

  l.c <- function(rho, t, s, e.z){
    lc <- sum(e.z*log(clayton.cop.den(t, s, rho)))
#    cat("rho=", rho, "\t", "l.c=", lc, "\n")
    return(lc)
  }

  rho.max <- optimize(f=l.c, c(0.1, 20), maximum=T, tol=0.00001, t=t, s=s, e.z=e.z)

  return(rho.max$m)
}



# updated
# mixture likelihood of two gaussian copula
# nonparametric and ranked transformed
loglik.2gaussian.copula <- function(x, y, p, rho1, rho2, x.mar, y.mar){
 
  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
  c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)

  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}

loglik.2copula <- function(x, y, p, rho1, rho2, x.mar, y.mar, copula.txt){

  px.1 <- pdf.cdf$px.1
  px.2 <- pdf.cdf$px.2
  py.1 <- pdf.cdf$py.1
  py.2 <- pdf.cdf$py.2

  if(copula.txt=="gaussian"){
    c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
    c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  } else {
    if(copula.txt=="clayton"){
      c1 <- clayton.cop.den(px.1$cdf, py.1$cdf, rho1)
      c2 <- clayton.cop.den(px.2$cdf, py.2$cdf, rho2)
    }
  }  
  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}


# estimate the marginals of each component using histogram estimator in EM
# return the density, breaks, and cdf of the histogram estimator 
est.mar.hist <- function(x, e.z, breaks){

  binwidth <- c()
  nbin <- length(breaks)-1
  nx <- length(x) 

  # the histogram
  x1.pdf <- c()
  x2.pdf <- c()
  x1.cdf <- c()
  x2.cdf <- c()

  # the pdf for each point
  x1.pdf.value <- rep(NA, nx)
  x2.pdf.value <- rep(NA, nx)

  x1.cdf.value <- rep(NA, nx)
  x2.cdf.value <- rep(NA, nx) 

  for(i in 1:nbin){

    binwidth[i] <- breaks[i+1] - breaks[i]
    if(i < nbin)
      in.bin <- x>= breaks[i] & x < breaks[i+1]
    else    # last bin
      in.bin <- x>= breaks[i] & x <=breaks[i+1]

    # each bin add one observation to avoid empty bins
    # multiple (nx+nbin)/(nx+nbin+1) to avoid blowup when looking up for
    # quantiles 
    x1.pdf[i] <- (sum(e.z[in.bin])+1)/(sum(e.z)+nbin)/binwidth[i]*(nx+nbin)/(nx+nbin+1)        
    x2.pdf[i] <- (sum(1-e.z[in.bin])+1)/(sum(1-e.z)+nbin)/binwidth[i]*(nx+nbin)/(nx+nbin+1) 


#    x1.pdf[i] <- sum(e.z[in.bin])/sum(e.z)/binwidth[i]*nx/(nx+1)        
#    x2.pdf[i] <- sum(1-e.z[in.bin])/sum(1-e.z)/binwidth[i]*nx/(nx+1) 
    
# treat each bin as a value for a discrete variable    
#    x1.cdf[i] <- sum(x1.pdf[1:i]*binwidth[1:i])
#    x2.cdf[i] <- sum(x2.pdf[1:i]*binwidth[1:i])


    # cumulative density before reaching i
    if(i>1){
      x1.cdf[i] <- sum(x1.pdf[1:(i-1)]*binwidth[1:(i-1)])
      x2.cdf[i] <- sum(x2.pdf[1:(i-1)]*binwidth[1:(i-1)])    
    } else{
      x1.cdf[i] <- 0
      x2.cdf[i] <- 0
    }

    # make a vector of nx to store the values of pdf and cdf for each x
    # this will speed up the computation dramatically
    x1.pdf.value[in.bin] <- x1.pdf[i]
    x2.pdf.value[in.bin] <- x2.pdf[i]

    x1.cdf.value[in.bin] <- x1.cdf[i] + x1.pdf[i]*(x[in.bin]-breaks[i])
    x2.cdf.value[in.bin] <- x2.cdf[i] + x2.pdf[i]*(x[in.bin]-breaks[i])      
  }

#  x1.cdf <- cumsum(x1.pdf*binwidth)
#  x2.cdf <- cumsum(x2.pdf*binwidth)

  f1 <-list(breaks=breaks, density=x1.pdf, cdf=x1.cdf)
  f2 <-list(breaks=breaks, density=x2.pdf, cdf=x2.cdf)

  f1.value <- list(pdf=x1.pdf.value, cdf=x1.cdf.value)
  f2.value <- list(pdf=x2.pdf.value, cdf=x2.cdf.value)

  return(list(f1=f1, f2=f2, f1.value=f1.value, f2.value=f2.value))
}

# estimate the marginal cdf from rank
est.cdf.rank <- function(x, conf.z){

  # add 1 to prevent blow up
  x1.cdf <- rank(x[conf.z==1])/(length(x[conf.z==1])+1)

  x2.cdf <- rank(x[conf.z==0])/(length(x[conf.z==0])+1)

  return(list(cdf1=x1.cdf, cdf2=x2.cdf))
}

# df is a density function with fields: density, cdf and breaks, x is a scalar
get.pdf <- function(x, df){

  if(x < df$breaks[1])
    cat("x is out of the range of df\n")

  index <- which(df$breaks >= x)[1]

  if(index==1)
    index <- index +1
  return(df$density[index-1])  
}

# get cdf from histgram estimator for a single value
get.cdf <- function(x, df){

  index <- which(df$breaks >= x)[1]
  if(index==1)
    index <- index +1
  return(df$cdf[index-1])   
}

# df is a density function with fields: density, cdf and breaks
get.pdf.cdf <- function(x.vec, df){

  x.pdf <- sapply(x.vec, get.pdf, df=df)
  x.cdf <- sapply(x.vec, get.cdf, df=df) 
  return(list(cdf=x.cdf, pdf=x.pdf))
}

# E-step
# x and y are the original observations or ranks
# rho1 and rho2 are the parameters of each copula
# f1, f2, g1, g2 are functions, each is a histogram 
e.step.2gaussian <- function(x, y, p, rho1, rho2, x.mar, y.mar){

  # get pdf and cdf of each component from functions in the corresponding component 
  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
  c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  
  return(get.ez(p, c1, c2, px.1$pdf, py.1$pdf, px.2$pdf, py.2$pdf))
}

# E-step
# rho1 and rho2 are the parameters of each copula 
e.step.2copula <- function(x, y, p, rho1, rho2, x.mar, y.mar, copula.txt){

  # get pdf and cdf of each component from functions in the corresponding component 
  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  if(copula.txt=="gaussian"){
    c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
    c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  } else {
    if(copula.txt=="clayton"){
      c1 <- clayton.cop.den(px.1$cdf, py.1$cdf, rho1)
      c2 <- clayton.cop.den(px.2$cdf, py.2$cdf, rho2)
    } 
  }
  return(get.ez(p, c1, c2, px.1$pdf, py.1$pdf, px.2$pdf, py.2$pdf))
}




# M-step
m.step.2gaussian <- function(x, y, e.z, breaks){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  
  rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 

  p <- sum(e.z)/length(e.z) 

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar))
}

m.step.2copula <- function(x, y, e.z, breaks, copula.txt){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

  px.1 <- get.pdf.cdf(x, x.mar$f1)
  px.2 <- get.pdf.cdf(x, x.mar$f2)
  py.1 <- get.pdf.cdf(y, y.mar$f1)
  py.2 <- get.pdf.cdf(y, y.mar$f2)

  if(copula.txt=="gaussian"){
    rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  
    rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
  } else {
    if(copula.txt=="clayton"){
      rho1 <- mle.clayton.copula(px.1$cdf, py.1$cdf, e.z)  
      rho2 <- mle.clayton.copula(px.2$cdf, py.2$cdf, 1-e.z)      
    }
  }
  
  p <- sum(e.z)/length(e.z) 

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar))
}



# E-step: pass values
# x and y are the original observations or ranks
# rho1 and rho2 are the parameters of each copula
# f1, f2, g1, g2 are functions, each is a histogram 
e.step.2gaussian.value <- function(x, y, p, rho1, rho2, pdf.cdf){

  c1 <- gaussian.cop.den(pdf.cdf$px.1$cdf, pdf.cdf$py.1$cdf, rho1)
  c2 <- gaussian.cop.den(pdf.cdf$px.2$cdf, pdf.cdf$py.2$cdf, rho2)
  
  e.z <- get.ez(p, c1, c2, pdf.cdf$px.1$pdf, pdf.cdf$py.1$pdf, 
               pdf.cdf$px.2$pdf, pdf.cdf$py.2$pdf)
  return(e.z)
}


e.step.2copula.value <- function(x, y, p, rho1, rho2, pdf.cdf, copula.txt){

  if(copula.txt =="gaussian"){
    c1 <- gaussian.cop.den(pdf.cdf$px.1$cdf, pdf.cdf$py.1$cdf, rho1)
    c2 <- gaussian.cop.den(pdf.cdf$px.2$cdf, pdf.cdf$py.2$cdf, rho2)
  } else {
    if(copula.txt =="clayton"){
      c1 <- clayton.cop.den(pdf.cdf$px.1$cdf, pdf.cdf$py.1$cdf, rho1)
      c2 <- clayton.cop.den(pdf.cdf$px.2$cdf, pdf.cdf$py.2$cdf, rho2)      
    }
  }
  
  e.z <- get.ez(p, c1, c2, pdf.cdf$px.1$pdf, pdf.cdf$py.1$pdf, 
               pdf.cdf$px.2$pdf, pdf.cdf$py.2$pdf)
  return(e.z)
}


# M-step: pass values
m.step.2gaussian.value <- function(x, y, e.z, breaks, fix.rho2){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

#  px.1 <- get.pdf.cdf(x, x.mar$f1)
#  px.2 <- get.pdf.cdf(x, x.mar$f2)
#  py.1 <- get.pdf.cdf(y, y.mar$f1)
#  py.2 <- get.pdf.cdf(y, y.mar$f2)

  px.1 <- x.mar$f1.value
  px.2 <- x.mar$f2.value
  py.1 <- y.mar$f1.value
  py.2 <- y.mar$f2.value

  rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  

  if(!fix.rho2)
    rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
  else
    rho2 <- 0

  p <- sum(e.z)/length(e.z) 

  pdf.cdf <- list(px.1=px.1, px.2=px.2, py.1=py.1, py.2=py.2)

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar,
              pdf.cdf=pdf.cdf))
}

m.step.2gaussian.value2 <- function(x, y, e.z, breaks, fix.rho2, x.mar, y.mar){

  # compute f1, f2, g1 and g2
#  x.mar <- est.mar.hist(x, e.z, breaks)
#  y.mar <- est.mar.hist(y, e.z, breaks)  

#  px.1 <- get.pdf.cdf(x, x.mar$f1)
#  px.2 <- get.pdf.cdf(x, x.mar$f2)
#  py.1 <- get.pdf.cdf(y, y.mar$f1)
#  py.2 <- get.pdf.cdf(y, y.mar$f2)

  px.1 <- x.mar$f1.value
  px.2 <- x.mar$f2.value
  py.1 <- y.mar$f1.value
  py.2 <- y.mar$f2.value

  rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  

  if(!fix.rho2)
    rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
  else
    rho2 <- 0

  p <- sum(e.z)/length(e.z) 

  pdf.cdf <- list(px.1=px.1, px.2=px.2, py.1=py.1, py.2=py.2)

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar,
              pdf.cdf=pdf.cdf))
}



m.step.2copula.value <- function(x, y, e.z, breaks, fix.rho2, copula.txt){

  # compute f1, f2, g1 and g2
  x.mar <- est.mar.hist(x, e.z, breaks)
  y.mar <- est.mar.hist(y, e.z, breaks)  

#  px.1 <- get.pdf.cdf(x, x.mar$f1)
#  px.2 <- get.pdf.cdf(x, x.mar$f2)
#  py.1 <- get.pdf.cdf(y, y.mar$f1)
#  py.2 <- get.pdf.cdf(y, y.mar$f2)

  px.1 <- x.mar$f1.value
  px.2 <- x.mar$f2.value
  py.1 <- y.mar$f1.value
  py.2 <- y.mar$f2.value

  if(copula.txt=="gaussian"){
    rho1 <- mle.gaussian.copula(px.1$cdf, py.1$cdf, e.z)  
    
    if(!fix.rho2)
      rho2 <- mle.gaussian.copula(px.2$cdf, py.2$cdf, 1-e.z) 
    else
      rho2 <- 0
  } else {

    if(copula.txt=="clayton"){
      rho1 <- mle.clayton.copula(px.1$cdf, py.1$cdf, e.z)  
    
      if(!fix.rho2)
        rho2 <- mle.clayton.copula(px.2$cdf, py.2$cdf, 1-e.z) 
      else
        rho2 <- 0
    }    
  }
    
  p <- sum(e.z)/length(e.z) 

  pdf.cdf <- list(px.1=px.1, px.2=px.2, py.1=py.1, py.2=py.2)

  return(list(p=p, rho1=rho1, rho2=rho2, x.mar=x.mar, y.mar=y.mar,
              pdf.cdf=pdf.cdf))
}




# updated
# mixture likelihood of two gaussian copula
# nonparametric and ranked transformed
loglik.2gaussian.copula.value <- function(x, y, p, rho1, rho2, pdf.cdf){

  px.1 <- pdf.cdf$px.1
  px.2 <- pdf.cdf$px.2
  py.1 <- pdf.cdf$py.1
  py.2 <- pdf.cdf$py.2

  c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
  c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)

  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}



# updated
# mixture likelihood of two gaussian copula
# nonparametric and ranked transformed
loglik.2copula.value <- function(x, y, p, rho1, rho2, pdf.cdf, copula.txt){

  px.1 <- pdf.cdf$px.1
  px.2 <- pdf.cdf$px.2
  py.1 <- pdf.cdf$py.1
  py.2 <- pdf.cdf$py.2

  if(copula.txt=="gaussian"){
    c1 <- gaussian.cop.den(px.1$cdf, py.1$cdf, rho1)
    c2 <- gaussian.cop.den(px.2$cdf, py.2$cdf, rho2)
  } else {
    if(copula.txt=="clayton"){
      c1 <- clayton.cop.den(px.1$cdf, py.1$cdf, rho1)
      c2 <- clayton.cop.den(px.2$cdf, py.2$cdf, rho2)
    }
  }

  sum(log(p*c1*px.1$pdf*py.1$pdf + (1-p)*c2*px.2$pdf*py.2$pdf))
}



# EM for 2 Gaussian, speed up computation, unfinished

em.2gaussian.quick <- function(x, y, p0, rho1.0, rho2.0, eps, fix.p=F, stoc=T, fix.rho2=T){

  x <- rank(x, tie="random")
  y <- rank(y, tie="random")

#  x <- rank(x, tie="average")
#  y <- rank(y, tie="average")

  # nbin=20
  xy.min <- min(x, y)
  xy.max <- max(x, y)
  binwidth <- (xy.max-xy.min)/50
  breaks <- seq(xy.min-binwidth/100, xy.max+binwidth/100, by=(xy.max-xy.min+binwidth/50)/50)
#  breaks <- seq(xy.min, xy.max, by=binwidth)
  

  # initiate marginals 
  # initialization: first p0 data has 
#  e.z <- e.step.2gaussian(x, y, p0, rho1.0, rho2.0, x0.mar, y0.mar) # this starting point assumes two components are overlapped

  e.z <- c(rep(0.9, round(length(x)*p0)), rep(0.1, length(x)-round(length(x)*p0)))

  if(!stoc)
    para <- m.step.2gaussian.value(x, y, e.z, breaks, fix.rho2)
  else 
    para <- m.step.2gaussian.stoc.value(x, y, e.z, breaks, fix.rho2)


  if(fix.p){
    p <- p0
  } else {
    p <- para$p  
  }

  if(fix.rho2){
    rho2 <- rho2.0
  } else {
    rho2 <- para$rho2
  }

#  rho1 <- 0.8
  rho1 <- para$rho1

  l0 <- loglik.2gaussian.copula.value(x, y, p, rho1, rho2, para$pdf.cdf)

  loglik.trace <- c()
  loglik.trace[1] <- l0
#  loglik.trace[2] <- l1
  to.run <- T

  i <- 2

  # this two lines to remove
#  x.mar <- est.mar.hist(x, e.z, breaks)
#  y.mar <- est.mar.hist(y, e.z, breaks)  
  
  while(to.run){

    e.z <- e.step.2gaussian.value(x, y, p, rho1, rho2, para$pdf.cdf) 
    if(!stoc)
      para <- m.step.2gaussian.value(x, y, e.z, breaks, fix.rho2)
    else
      para <- m.step.2gaussian.stoc.value(x, y, e.z, breaks, fix.rho2)

    # fix x.mar and y.mar : to remove
#    if(!stoc)
#      para <- m.step.2gaussian.value2(x, y, e.z, breaks, fix.rho2, x.mar, y.mar)
#    else
#      para <- m.step.2gaussian.stoc.value(x, y, e.z, breaks, fix.rho2)

    
    if(fix.p){
      p <- p0
    } else {
      p <- para$p  
    }

    if(fix.rho2){
      rho2 <- rho2.0
    } else {
      rho2 <- para$rho2
    }

#    rho1 <- 0.8
    rho1 <- para$rho1

  #  l0 <- l1
    l1 <- loglik.2gaussian.copula.value(x, y, p, rho1, rho2, para$pdf.cdf)
    loglik.trace[i] <- l1

#cat("l1=", l1, "\n") 

    # Aitken acceleration criterion
    if(i > 2){
      l.inf <- loglik.trace[i-2] + (loglik.trace[i-1] - loglik.trace[i-2])/(1-(loglik.trace[i]-loglik.trace[i-1])/(loglik.trace[i-1]-loglik.trace[i-2])) 
      to.run <- abs(l.inf - loglik.trace[i]) > eps 
#cat("para=", "p=", para$p, " rho1=", rho1, " rho2=", rho2, "\n")
#cat("l.inf=", l.inf, "\n")
#cat(l.inf-loglik.trace[i], "\n")     
    }

    i <- i+1
  }

  bic <- -2*l1 + (2*(length(breaks)-1+1)+1-fix.p-fix.rho2)*log(length(x)) # parameters
  return(list(para=list(p=para$p, rho1=rho1, rho2=rho2), 
              loglik=l1, bic=bic, e.z=e.z, conf.z = para$conf.z, 
              loglik.trace=loglik.trace, x.mar=para$x.mar, y.mar=para$y.mar,
              breaks=breaks))
}



em.2copula.quick <- function(x, y, p0, rho1.0, rho2.0, eps, fix.p=F, stoc=T, fix.rho2=T, copula.txt, nbin=50){

  x <- rank(x, tie="random")
  y <- rank(y, tie="random")

#  x <- rank(x, tie="first")
#  y <- rank(y, tie="first")

  # nbin=50
  xy.min <- min(x, y)
  xy.max <- max(x, y)
  binwidth <- (xy.max-xy.min)/50
  breaks <- seq(xy.min-binwidth/100, xy.max+binwidth/100, by=(xy.max-xy.min+binwidth/50)/nbin)  
#  breaks <- seq(xy.min, xy.max, by=binwidth)
  
  # initiate marginals 
  # initialization: first p0 data has 
#  e.z <- e.step.2gaussian(x, y, p0, rho1.0, rho2.0, x0.mar, y0.mar) # this starting point assumes two components are overlapped

  e.z <- c(rep(0.9, round(length(x)*p0)), rep(0.1, length(x)-round(length(x)*p0)))


  if(!stoc)
    para <- m.step.2copula.value(x, y, e.z, breaks, fix.rho2, copula.txt)
  else 
    para <- m.step.2copula.stoc.value(x, y, e.z, breaks, fix.rho2, copula.txt)

  if(fix.p){
    p <- p0
  } else {
    p <- para$p  
  }

  if(fix.rho2){
    rho2 <- rho2.0
  } else {
    rho2 <- para$rho2
  }

  l0 <- loglik.2copula.value(x, y, p, para$rho1, rho2, para$pdf.cdf, copula.txt)

  loglik.trace <- c()
  loglik.trace[1] <- l0
#  loglik.trace[2] <- l1
  to.run <- T

  i <- 2

  while(to.run){

    e.z <- e.step.2copula.value(x, y, p, para$rho1, rho2, para$pdf.cdf, copula.txt) 
    if(!stoc)
      para <- m.step.2copula.value(x, y, e.z, breaks, fix.rho2, copula.txt)
    else
      para <- m.step.2copula.stoc.value(x, y, e.z, breaks, fix.rho2, copula.txt)

    if(fix.p){
      p <- p0
    } else {
      p <- para$p  
    }

    if(fix.rho2){
      rho2 <- rho2.0
    } else {
      rho2 <- para$rho2
    }


  #  l0 <- l1
    l1 <- loglik.2copula.value(x, y, p, para$rho1, rho2, para$pdf.cdf, copula.txt)
    loglik.trace[i] <- l1

cat("l1=", l1, "\n") 

    # Aitken acceleration criterion
    if(i > 2){
      l.inf <- loglik.trace[i-2] + (loglik.trace[i-1] - loglik.trace[i-2])/(1-(loglik.trace[i]-loglik.trace[i-1])/(loglik.trace[i-1]-loglik.trace[i-2])) 
      to.run <- abs(l.inf - loglik.trace[i]) > eps 
cat("para=", "p=", para$p, " rho1=", para$rho1, " rho2=", rho2, "\n")
#cat("l.inf=", l.inf, "\n")
#cat(l.inf-loglik.trace[i], "\n")     
    }

    i <- i+1
  }

  bic <- -2*l1 + (2*(length(breaks)-1+1)+1-fix.p-fix.rho2)*log(length(x)) # parameters
  return(list(para=list(p=para$p, rho1=para$rho1, rho2=rho2), 
              loglik=l1, bic=bic, e.z=e.z, conf.z = para$conf.z, 
              loglik.trace=loglik.trace, x.mar=para$x.mar, y.mar=para$y.mar,
              breaks=breaks))
}


#######################
####################### fit EM procedure for the matched peaks
#######################

# remove the unmatched ones
#rm.unmatch <- function(sample1, sample2, p.value.impute=0){
#
#  sample1.prune <- sample1[sample1$sig.value > p.value.impute & sample2$sig.value > p.value.impute,]
#  sample2.prune <- sample2[sample1$sig.value > p.value.impute & sample2$sig.value > p.value.impute,]
# 
#  invisible(list(sample1=sample1.prune$sig.value, sample2=sample2.prune$sig.value))
#}


# fit 2-component model
#fit.em <- function(sample12, fix.rho2=T){
#
#  prune.sample <- rm.unmatch(sample12$merge1, sample12$merge2)
#
#  em.fit <- em.2gaussian.quick(-prune.sample$sample1, -prune.sample$sample2,
# p0=0.5, rho1.0=0.7, rho2.0=0, eps=0.01, fix.p=F, stoc=F, fix.rho2)
#
#  invisible(list(em.fit=em.fit, data.pruned=prune.sample))
#}


rm.unmatch <- function(sample1, sample2, p.value.impute=0){

  sample1.prune <- sample1[sample1$sig.value > p.value.impute & sample2$sig.value > p.value.impute,]
  sample2.prune <- sample2[sample1$sig.value > p.value.impute & sample2$sig.value > p.value.impute,]
 
  invisible(list(sample1=sample1.prune, sample2=sample2.prune))
}


# fit 2-component model
fit.em <- function(sample12, fix.rho2=T){

  prune.sample <- rm.unmatch(sample12$merge1, sample12$merge2)

  em.fit <- em.2gaussian.quick(-prune.sample$sample1$sig.value, -prune.sample$sample2$sig.value,
 p0=0.5, rho1.0=0.7, rho2.0=0, eps=0.01, fix.p=F, stoc=F, fix.rho2)

  invisible(list(em.fit=em.fit, data.pruned=prune.sample))
}



fit.2copula.em <- function(sample12, fix.rho2=T, copula.txt){

  prune.sample <- rm.unmatch(sample12$merge1, sample12$merge2)

#  o <- order(prune.sample$sample1)
#  n <- length(prune.sample$sample1)
    
#  para <- init(prune.sample$sample1$sig.value, prune.sample$sample2$sig.value, c(rep(0, round(n/3)), rep(c(0,1), round(n/6)), rep(1, n-round(n/3)-round(n/6))))

#  temp <- init.dist(f0, f1)
  para <- list()
  para$rho <- 0.6
  para$p <- 0.3
  para$mu <- 2.5
  para$sigma <- 1
##  para$mu <- -temp$mu
##  para$sigma <- temp$sigma
#cat("mu=", para$mu, "sigma=", para$sigma, "\n")
  
#  em.fit <- em.transform.1loop(-prune.sample$sample1, -prune.sample$sample2,
  cat("EM is running")
  em.fit <- em.transform(prune.sample$sample1$sig.value, prune.sample$sample2$sig.value, para$mu, para$sigma, para$rho, para$p, eps=0.01)

  invisible(list(em.fit=em.fit, data.pruned=prune.sample))
}




# fit 1-component model
fit.1.component <- function(data.pruned, breaks){

#  gaussian.1 <- fit.gaussian.1(-data.pruned$sample1$sig.value, -data.pruned$sample2$sig.value, breaks)
#  clayton.1 <- fit.clayton.1(-data.pruned$sample1$sig.value, -data.pruned$sample2$sig.value, breaks)

  gaussian.1 <- fit.gaussian.1(-data.pruned$sample1, -data.pruned$sample2, breaks)
  clayton.1 <- fit.clayton.1(-data.pruned$sample1, -data.pruned$sample2, breaks)

  return(list(gaussian.1=gaussian.1, clayton.1=clayton.1))
}



#################
# Fit a single component  
#################

# a single gaussian copula
# if breaks=NULL, use empirical pdf, otherwise use histogram estimate
fit.gaussian.1 <- function(x, y, breaks=NULL){

  # rank transformed and compute the empirical cdf
  t <- emp.mar.cdf.rank(x)
  s <- emp.mar.cdf.rank(y)

  mle.rho <- mle.gaussian.copula(t, s, rep(1, length(t)))

  c1 <- gaussian.cop.den(t, s, mle.rho)
cat("c1", sum(log(c1)), "\n")

  if(is.null(breaks)){
    f1 <- emp.mar.pdf.rank(t)
    f2 <- emp.mar.pdf.rank(s)
  } else {
    x.mar <- est.mar.hist(rank(x), rep(1, length(x)), breaks)
    y.mar <- est.mar.hist(rank(y), rep(1, length(y)), breaks)

    f1 <- x.mar$f1.value$pdf  # only one component
    f2 <- y.mar$f1.value$pdf
  }


cat("f1", sum(log(f1)), "\n")
cat("f2", sum(log(f2)), "\n")

  loglik <- sum(log(c1)+log(f1)+log(f2))

  bic <- -2*loglik + log(length(t))*(1+length(breaks)-1)

  return(list(rho=mle.rho, loglik=loglik, bic=bic))
}


# a single Clayton copula
fit.clayton.1 <- function(x, y, breaks=NULL){

  # rank transformed and compute the empirical cdf
  t <- emp.mar.cdf.rank(x)
  s <- emp.mar.cdf.rank(y)

  mle.rho <- mle.clayton.copula(t, s, rep(1, length(t)))

  c1 <- clayton.cop.den(t, s, mle.rho)

  if(is.null(breaks)){
    f1 <- emp.mar.pdf.rank(t)
    f2 <- emp.mar.pdf.rank(s)
  } else {
    x.mar <- est.mar.hist(rank(x), rep(1, length(x)), breaks)
    y.mar <- est.mar.hist(rank(y), rep(1, length(y)), breaks)

    f1 <- x.mar$f1.value$pdf  # only one component
    f2 <- y.mar$f1.value$pdf
  }

  loglik <- sum(log(c1)+log(f1)+log(f2))

  bic <- -2*loglik + log(length(t))*(1+length(breaks)-1)

  return(list(rho=mle.rho, tau=rho/(rho+2), loglik=loglik, bic=bic)) 
}

## obsolete function (01-06-2010)
## compute the average posterior probability to belong to the random component
## for peaks selected at different cutoffs 
comp.uri.ez <- function(tt, u, v, e.z){

   u.t <- quantile(u, prob=(1-tt))
   v.t <- quantile(v, prob=(1-tt))

 #  ez <- mean(e.z[u >= u.t & v >=u.t]) Is this wrong?
   ez <- mean(e.z[u >= u.t & v >=v.t])

   return(ez)
}

## obsolete function (01-06-2010)
# compute the largest posterior error probability corresponding to
# the square centered at the origin and spanned top tt% on both coordinates
# so the consistent low rank ones are excluded
# boundary.txt: either "max" or "min", if it is error prob, use "max"
comp.ez.cutoff <- function(tt, u, v, e.z, boundary.txt){

   u.t <- quantile(u, prob=(1-tt))
   v.t <- quantile(v, prob=(1-tt))

   if(boundary.txt == "max"){
 #    ez.bound <- max(e.z[u >= u.t & v >=u.t])
     ez.bound <- max(e.z[u >= u.t & v >=v.t])
   } else {
 #    ez.bound <- min(e.z[u >= u.t & v >=u.t])
     ez.bound <- min(e.z[u >= u.t & v >=v.t])     
   }

   return(ez.bound)

}

# obsolete functions: 01-06-2010
# compute the error rate
# u.t and v.t are the quantiles
# this one is used for the plots generated initially in the brief writeup  
# and it was used for processing merged data in July before the IDR definition
# is formalized
# It does not implement the current definition of IDR
get.ez.tt.old  <- function(em.fit, reverse=T, fdr.level=c(0.01, 0.05, 0.1)){

  u <- em.fit$data.pruned$sample1
  v <- em.fit$data.pruned$sample2

  tt <- seq(0.01, 0.99, by=0.01)
  if(reverse){ 
    e.z <-  1-em.fit$em.fit$e.z # this is the error prob
    uri.ez <- sapply(tt, comp.uri.ez, u=u, v=v, e.z=e.z)
    ez.bound <- sapply(tt, comp.ez.cutoff, u=u, v=v, e.z=e.z, boundary.txt="max") 
  } else {
    e.z <-  em.fit$em.fit$e.z
    uri.ez <- sapply(tt, comp.uri.ez, u=u, v=v, e.z=e.z)
    ez.bound <- sapply(tt, comp.ez.cutoff, u=u, v=v, e.z=e.z, boundary.txt="min") 
  }

  u.t <- quantile(u, prob=(1-tt))
  v.t <- quantile(v, prob=(1-tt))  

  # find the levels on the two replicates
  sig.value1 <- c()
  sig.value2 <- c()
  error.prob.cutoff <- c()
  n.selected.match <- c()

  for(i in 1:length(fdr.level)){

    # find which uri.ez is closet to fdr.level
    index <- which.min(abs(uri.ez - fdr.level[i]))
    sig.value1[i] <- u.t[index]
    sig.value2[i] <- v.t[index]
    error.prob.cutoff[i] <- ez.bound[index]  
    if(reverse){
      n.selected.match[i] <- sum(e.z<=ez.bound[index])    
    } else {
      n.selected.match[i] <- sum(e.z>=ez.bound[index])    
    }
  }   

  # output the cutoff of posterior probability, signal values on two replicates
  map.uv <- cbind(error.prob.cutoff, sig.value1, sig.value2, n.selected.match)

  return(list(n=tt*length(u), uri.ez=uri.ez, u.t=u.t, v.t=v.t, tt=tt, fdr.level=fdr.level,  map.uv=map.uv, e.z=e.z, error.prob.cutoff=error.prob.cutoff))
}

# created: 01-06-2010
# Output IDR at various number of selected peaks
# Find cutoff (idr cutoff, sig.value cutoff on each replicate) for specified IDR level
# IDR definition is similar to FDR
get.ez.tt <- function(em.fit, idr.level=c(0.01, 0.05, 0.1)){

#  u <- em.fit$data.pruned$sample1$sig.value
#  v <- em.fit$data.pruned$sample2$sig.value
  u <- em.fit$data.pruned$sample1
  v <- em.fit$data.pruned$sample2
  
  e.z <-  1-em.fit$em.fit$e.z # this is the error prob
  
  o <- order(e.z)
  e.z.ordered <- e.z[o]
  n.select <- c(1:length(e.z))
  IDR <- cumsum(e.z.ordered)/n.select

  u.o <- u[o]
  v.o <- v[o]

  n.level <- length(idr.level)
#  sig.value1 <- rep(NA, n.level)
#  sig.value2 <- rep(NA, n.level)
  ez.cutoff <- rep(NA, n.level)
  n.selected <- rep(NA, n.level)
  
  for(i in 1:length(idr.level)){

    # find which uri.ez is closet to fdr.level
    index <- which.min(abs(IDR - idr.level[i]))
#    sig.value1[i] <- min(u.o[1:index])
#    sig.value2[i] <- min(v.o[1:index])
    ez.cutoff[i] <- e.z[index]      
    n.selected[i] <- sum(e.z<=ez.cutoff[i])    
  }   

  # output the cutoff of posterior probability, number of selected overlapped peaks 
#  map.uv <- cbind(ez.cutoff, sig.value1, sig.value2, n.selected)

  map.uv <- cbind(ez.cutoff, n.selected)

  return(list(n=n.select, IDR=IDR, idr.level=idr.level, map.uv=map.uv))
}   
  
#  return(list(n=tt*length(u), uri.ez=uri.ez,  fdr.level=fdr.level,  map.uv=map.uv, e.z=e.z, error.prob.cutoff=error.prob.cutoff))  
  




### compute the mean of the marginals
get.mar.mean <- function(em.out){

  x.f1 <- em.out$x.mar$f1
  x.f2 <- em.out$x.mar$f2

  y.f1 <- em.out$y.mar$f1
  y.f2 <- em.out$y.mar$f2

  x.stat1 <- get.hist.mean(x.f1)
  x.stat2 <- get.hist.mean(x.f2)
  y.stat1 <- get.hist.mean(y.f1)
  y.stat2 <- get.hist.mean(y.f2)

  return(list(x.mean1=x.stat1$mean, x.mean2=x.stat2$mean, 
              y.mean1=y.stat1$mean, y.mean2=y.stat2$mean,
              x.sd1=x.stat1$sd, x.sd2=x.stat2$sd, 
              y.sd1=y.stat1$sd, y.sd2=y.stat2$sd
              ))

}


# compute the mean of marginals
get.hist.mean  <- function(x.f){

  nbreaks <- length(x.f$breaks)
  x.bin <- x.f$breaks[-1]-x.f$breaks[-nbreaks]

  x.mid <- (x.f$breaks[-nbreaks]+x.f$breaks[-1])/2
  x.mean <- sum(x.mid*x.f$density*x.bin)
  x.sd <- sqrt(sum(x.mid*x.mid*x.f$density*x.bin)-x.mean^2)
  
  return(list(mean=x.mean, sd=x.sd))
}

get.hist.var <- function(x.f){

  nbreaks <- length(x.f$breaks)
  x.bin <- x.f$breaks[-1]-x.f$breaks[-nbreaks]

  x.mid <- (x.f$breaks[-nbreaks]+x.f$breaks[-1])/2
  x.mean <- sum(x.mid*x.f$density*x.bin)

  return(mean=x.mean)  
}

# obsolete function (01-06-2010)
# plot 
plot.ez.group.old <- function(ez.list, plot.dir, file.name=NULL, legend.txt, y.lim=NULL, xlab.txt="num of significant peaks",  ylab.txt="avg posterior prob of being random", col.txt=NULL, title.txt=NULL){

  if(is.null(col.txt))
    col.txt <- c("black", "red", "purple", "green", "blue", "cyan", "magenta", "orange", "grey")

  x <- c()
  y <- c()

  for(i in 1:length(ez.list)){
    x <- c(x, ez.list[[i]]$n)
      
    y <- c(y, ez.list[[i]]$uri.ez)
  }

  if(is.null(y.lim))
    y.lim <- c(0, max(y))

  if(!is.null(file.name)){
    postscript(paste(plot.dir, "ez.", file.name, sep=""))
    par(mfrow=c(1,1), mar=c(5,5,4,2))
  }

  plot(x, y, ylim=y.lim, type="n", xlab=xlab.txt, ylab=ylab.txt, lwd=5, cex=5, cex.axis=2, cex.lab=2)

  for(i in 1:length(ez.list)){
    lines(ez.list[[i]]$n, ez.list[[i]]$uri.ez, col=col.txt[i], cex=2, lwd=5)    
  }

#   plot(ez.list[[1]]$u.t, y, ylim=y.lim, type="l", xlab="rep-sig", ylab=ylab.txt, lwd=5, cex=5, cex.axis=2, cex.lab=2)
#   plot(ez.list[[1]]$v.t, y, ylim=y.lim, type="l", xlab="rep-sig", ylab=ylab.txt, lwd=5, cex=5, cex.axis=2, cex.lab=2)
  

  legend(0, y.lim[2], legend=legend.txt, col=col.txt[1:length(col.txt)], lty=1, lwd=5, cex=2)

  if(!is.null(title))
    title(title.txt)

  if(!is.null(file.name)){
    dev.off()
  }
  
}


plot.ez.group <- function(ez.list, plot.dir, file.name=NULL, legend.txt, y.lim=NULL, xlab.txt="num of significant peaks",  ylab.txt="IDR", col.txt=NULL, title.txt=NULL){

  if(is.null(col.txt))
    col.txt <- c("black", "red", "purple", "green", "blue", "cyan", "magenta", "orange", "grey")

  n.entry <- length(ez.list)
  x <- rep(NA, n.entry)
  y.max <- rep(NA, n.entry)

  for(i in 1:n.entry){
    x[i] <- max(ez.list[[i]]$n)
      
    y.max[i] <- max(ez.list[[i]]$IDR)
  
  }

  if(is.null(y.lim))
    y.lim <- c(0, max(y.max))

  if(!is.null(file.name)){
    postscript(paste(plot.dir, "ez.", file.name, sep=""))
    par(mfrow=c(1,1), mar=c(5,5,4,2))
  }


  
  plot(c(0, max(x)), y.lim, ylim=y.lim, type="n", xlab=xlab.txt, ylab=ylab.txt, lwd=5, cex=5, cex.axis=2, cex.lab=2)

  q <- seq(0.01, 0.99, by=0.01)
  
  for(i in 1:length(ez.list)){

    n.plot <- round(quantile(ez.list[[i]]$n, prob=q))
    IDR.plot <- ez.list[[i]]$IDR[n.plot]
    lines(n.plot, IDR.plot, col=col.txt[i], cex=2, lwd=5)    
  }


  legend(0, y.lim[2], legend=legend.txt, col=col.txt[1:length(col.txt)], lty=1, lwd=5, cex=2)

  if(!is.null(title))
    title(title.txt)

  if(!is.null(file.name)){
    dev.off()
  }
  
}



#############################################################################
#############################################################################
# statistics about peaks selected on the individual replicates
#
# idr.level: the consistency cutoff, say 0.05
# uri.output: a list of uri.output from consistency analysis generated by batch-consistency-analysis.r
# ez.list : a list of IDRs computed from get.ez.tt using the same idr.level
#
##################


# obsolete?
# compute the error rate
# u.t and v.t are the quantiles
# 
# map back to all peaks and report the number of peaks selected
get.ez.tt.all.old  <- function(em.fit, all.data1, all.data2, idr.level){

  u <- em.fit$data.pruned$sample1
  v <- em.fit$data.pruned$sample2

  tt <- seq(0.01, 0.99, by=0.01)
#  if(reverse){ 
    e.z <-  1-em.fit$em.fit$e.z # this is the error prob
    uri.ez <- sapply(tt, comp.uri.ez, u=u, v=v, e.z=e.z)
    ez.bound <- sapply(tt, comp.ez.cutoff, u=u, v=v, e.z=e.z, boundary.txt="max") 
#  } else {
#    e.z <-  em.fit$em.fit$e.z
#    uri.ez <- sapply(tt, comp.uri.ez, u=u, v=v, e.z=e.z)
#    ez.bound <- sapply(tt, comp.ez.cutoff, u=u, v=v, e.z=e.z, boundary.txt="min") 
#  }

  u.t <- quantile(u, prob=(1-tt))
  v.t <- quantile(v, prob=(1-tt))  

  # find the levels on the two replicates
  sig.value1 <- c()
  sig.value2 <- c()
  error.prob.cutoff <- c()
  n.selected.match <- c()
  npeak.rep1 <- c()
  npeak.rep2 <- c()

  for(i in 1:length(idr.level)){

    # find which uri.ez is closet to idr.level
    index <- which.min(abs(uri.ez - as.numeric(idr.level[i])))

    sig.value1[i] <- u.t[index]
    sig.value2[i] <- v.t[index]
    error.prob.cutoff[i] <- ez.bound[index]  
    n.selected.match[i] <- sum(u>= u.t[index] & v>=v.t[index])

    npeak.rep1[i] <- sum(all.data1["sig.value"] >= sig.value1[i])
    npeak.rep2[i] <- sum(all.data2["sig.value"] >= sig.value2[i])    
  }   


  # output the cutoff of posterior probability, signal values on two replicates
  map.uv <- cbind(error.prob.cutoff, sig.value1, sig.value2, n.selected.match, npeak.rep1, npeak.rep2)

  return(list(n=tt*length(u), uri.ez=uri.ez, u.t=u.t, v.t=v.t, tt=tt, idr.level=idr.level,  map.uv=map.uv, e.z=e.z, error.prob.cutoff=error.prob.cutoff))
}


get.ez.tt.all <- function(em.fit, all.data1, all.data2, idr.level=c(0.01, 0.05, 0.1)){

  u <- em.fit$data.pruned$sample1$sig.value
  v <- em.fit$data.pruned$sample2$sig.value
#  u <- em.fit$data.pruned$sample1
#  v <- em.fit$data.pruned$sample2
  
  e.z <-  1-em.fit$em.fit$e.z # this is the error prob
  
  o <- order(e.z)
  e.z.ordered <- e.z[o]
  n.select <- c(1:length(e.z))
  IDR <- cumsum(e.z.ordered)/n.select

  u.o <- u[o]
  v.o <- v[o]

  n.level <- length(idr.level)
#  sig.value1 <- rep(NA, n.level)
#  sig.value2 <- rep(NA, n.level)
  ez.cutoff <- rep(NA, n.level)
  n.selected <- rep(NA, n.level)
  npeak.rep1 <- rep(NA, n.level)
  npeak.rep2 <- rep(NA, n.level)
  
  for(i in 1:length(idr.level)){

    # find which uri.ez is closet to fdr.level
    index <- which.min(abs(IDR - idr.level[i]))
#    sig.value1[i] <- min(u.o[1:index])
#    sig.value2[i] <- min(v.o[1:index])
    ez.cutoff[i] <- e.z.ordered[index]      # fixed on 02/20/10
    n.selected[i] <- sum(e.z<=ez.cutoff[i])
#    npeak.rep1[i] <- sum(all.data1["sig.value"] >= sig.value1[i])
#    npeak.rep2[i] <- sum(all.data2["sig.value"] >= sig.value2[i])     
  }   

  # output the cutoff of posterior probability, number of selected overlapped peaks 
  map.uv <- cbind(ez.cutoff, n.selected)

  return(list(n=n.select, IDR=IDR, idr.level=idr.level, map.uv=map.uv))
}   
  
#  return(list(n=tt*length(u), uri.ez=uri.ez,  fdr.level=fdr.level,  map.uv=map.uv, e.z=e.z, error.prob.cutoff=error.prob.cutoff))  
  





####### the following is for determining thresholds for merged dataset

############# select peaks above a given threshold
#
# pass.threshold: a simple method, passing the threshold on the threshold on the individual replicate to the pooled sample 
#                 
# sig.map.list: a list of matrix to include all the cutoff values, each row corresponds to a cutoff. The first column is idr.level
#          the 2nd column is the cutoff of ez, the rest of columns are consistency analysis for other replicates 
# sig.value.name: the name of the sig.value column
# combined: combined dataset
# nrep: number of pairs of comparisons
#
# Procedure:
# 1. Find the significant threshold corresponding to the idr cutoff on the matched peaks. 
# 2. Each time we will get two or more (if >2 replicates) cutoffs and will report the most stringent and the least stringent
#    cutoff and the number of peaks selected at those two cutoffs
#############

pass.threshold <- function(sig.map.list, sig.value.name, combined, idr.level, nrep, chr.size){

  sig.map <- c()

  # choose idr.level
  idr.index <- which(rbind(sig.map.list[[1]])[,1] == idr.level)
  if(length(i) ==0){
    print("no level matches specified idr.level")
    return(-1)
  }

  for(i in 1:length(sig.map.list))
    sig.map <- c(sig.map, rbind(sig.map.list[[i]])[idr.index, c("sig.value1", "sig.value2")])
  
  
  npeak.tight <- c()  
  npeak.loose <- c()


  max.sig <- max(sig.map)
  min.sig <- min(sig.map)
  selected.sig.tight <- combined[combined[,sig.value.name]>=max.sig, ]
  selected.sig.loose <- combined[combined[,sig.value.name]>=min.sig, ]

  selected.sig.tight <- deconcatenate.chr(selected.sig.tight, chr.size)[,c("chr", "start", "stop", "signal.value", "p.value", "q.value")]
  selected.sig.loose <- deconcatenate.chr(selected.sig.loose, chr.size)[,c("chr", "start", "stop", "signal.value", "p.value", "q.value")]
  
  npeak.tight <- nrow(selected.sig.tight)
  npeak.loose <- nrow(selected.sig.loose)
  
  
  npeak.stat <- list(idr.level=idr.level, max.sig=max.sig, min.sig=min.sig, npeak.tight=npeak.tight, npeak.loose=npeak.loose)

  invisible(list(npeak.stat=npeak.stat, combined.selected.tight=selected.sig.tight, combined.selected.loose=selected.sig.loose))  
}

#################
# pass the regions selected from consistency analysis to combined data
# Threshold is determined on the replicates, the regions above the threshold are selected
# then peaks on the combined data are selected from the selected regions
#
# To avoid being too stringent, regions satisfying the following conditions are selected
# 1. regions above the significant threshold determined by consistency analysis on either replicate
# 2. regions that have consistent low peaks, i.e. posterior prob > threshold but not passing the significant threshold
#
# This method doesn't make a difference when using different thresholds
#################

pass.region <- function(sig.map.list, uri.output, ez.list, em.output, combined, idr.level, sig.value.impute=0, chr.size){
 
  combined <- combined[, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")]
  npair <- length(uri.output) # number of pairs of consistency analysis
  combined.region <- c()

  # choose idr.level
  idr.index <- which(rbind(sig.map.list[[1]])[,1] == idr.level)
  if(length(idr.index) ==0){
    print("no level matches specified idr.level")
    return(-1)
  }

    for(j in 1:npair){
      # select peaks from individual replicates using individual cutoff
      above.1 <- uri.output[[j]]$data12.enrich$merge1["sig.value"] >= ez.list[[j]]$map.uv[idr.index,"sig.value1"]
      above.2 <- uri.output[[j]]$data12.enrich$merge1["sig.value"] >= ez.list[[j]]$map.uv[idr.index,"sig.value2"]
      selected.sig.rep1 <- uri.output[[j]]$data12.enrich$merge1[above.1, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")]
      selected.sig.rep2 <- uri.output[[j]]$data12.enrich$merge2[above.2, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")] 
      
      # find the peaks that are overlapped with reliable peaks in the individual replicates
      overlap.1 <- pair.peaks(selected.sig.rep1, combined)$merge2
      overlap.2 <- pair.peaks(selected.sig.rep2, combined)$merge2

      # choose the ones with significant value > 0, which are the overlapped ones

      combined.in1 <- overlap.1[overlap.1$sig.value > sig.value.impute, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")]
      combined.in2 <- overlap.2[overlap.2$sig.value > sig.value.impute, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")]

      ## consistent low significant ones
      ## first find consistenct ones, ie. high posterior prob
      # is.consistent <- ez.list[[j]]$e.z < ez.list[[j]]$ez.cutoff 

      # data.matched <- keep.match(uri.output[[j]]$data12.enrich$merge1[!above.1, ], uri.output[[j]]$data12.enrich$merge2[!above.2, ], sig.value.impute=0)
      # data.matched$sample1 <- data.matched$sample1[, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")]
      # data.matched$sample2 <- data.matched$sample2[, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")]

      # consistent.in1 <- data.matched$sample1[is.consistent, ]
      # consistent.in2 <- data.matched$sample2[is.consistent, ]

      # overlap.consistent.1 <- pair.peaks(consistent.in1, combined)$merge2
      # overlap.consistent.2 <- pair.peaks(consistent.in2, combined)$merge2

      ## choose the ones with significant value > 0, which are the overlapped ones

      # combined.consistent.in1 <- overlap.consistent.1[overlap.consistent.1$sig.value > sig.value.impute, ]
      # combined.consistent.in2 <- overlap.consistent.2[overlap.consistent.2$sig.value > sig.value.impute, ]

      # combined.region <- rbind(combined.region, combined.in1, combined.in2, combined.consistent.in1, combined.consistent.in2)

       combined.region <- rbind(combined.region, combined.in1, combined.in2)

      is.repeated <- duplicated(combined.region$start)
      combined.region <- combined.region[!is.repeated, c("start", "stop", "sig.value", "signal.value", "p.value", "q.value")]
      
    }
    npeak <- nrow(combined.region)
  
   sig.combined <- c(min(combined.region[,"sig.value"], na.rm=T), max(combined.region[,"sig.value"], na.rm=T))

  # idr.combined <- c(min(combined.region[,"q.value"], na.rm=T), max(combined.region[,"q.value"], na.rm=T))

   npeak.stat <- list(idr.level=idr.level, npeak=npeak)

   combined.region <- deconcatenate.chr(combined.region, chr.size)[,c("chr", "start", "stop", "signal.value", "p.value", "q.value")]

  invisible(list(npeak.stat=npeak.stat, combined.selected=combined.region, sig.combined=sig.combined))
}

################
# pass structure: this method does another round of inference on the combined data
#
# To make the mixture structure comparable on the replicates and the combined data, the 2nd inference is done on the peaks
# at the reliable regions on the combined data, using rank transformed significant values. The mixture structure is estimated using my consistency analysis, which 
# estimates marginal distributions of ranks using nonparametric ways. Then the significant values are found out.
# There are several advantages to do it this way:  
# 1. The premise of passing structure is that the means and variance (i.e. distribution) of two replicates should be the same
#    The significant values on the two replicates clearly have different distributions. The structure estimated from consistency
#    analysis will generate similar rank distribution on two replicates by its setup (i.e. same number of peaks are paired up).  
# 2. Because pooled sample is a black box, the structure is more likely to be followed in the matched regions than other locations,
#    after all, we don't know what other things are. If even the structure doesn't hold on the matched regions, 
#    which is possible, let alone the other regions. Focusing on the reliable regions helps to get rid of those unknown noises.
#  
# 
# modified on 2-20-10: reverse rank.combined, make big sig.value with small
# ranks, to be consistent with f1 and f2
################ 

pass.structure <- function(uri.output, em.output, combined, idr.level, sig.value.impute, chr.size, overlap.ratio=0){

  columns.keep <- c("sig.value", "start", "stop", "signal.value", "p.value", "q.value", "chr", "start.ori", "stop.ori")
  combined <- combined[, columns.keep]
  combined.selected.all <- c()

  for(j in 1:npair){

    sample1 <- uri.output[[j]]$data12.enrich$merge1[, columns.keep]
    sample2 <- uri.output[[j]]$data12.enrich$merge2[, columns.keep]
        
    # find peaks on the matched region on the combined one
    data.matched <- keep.match(sample1, sample2, sig.value.impute=sig.value.impute)

    data.matched$sample1 <- data.matched$sample1[, columns.keep]
    data.matched$sample2 <- data.matched$sample2[, columns.keep]

    overlap.1 <- pair.peaks.filter(data.matched$sample1, combined, p.value.impute=sig.value.impute, overlap.ratio)$merge2
    overlap.2 <- pair.peaks.filter(data.matched$sample2, combined, p.value.impute=sig.value.impute, overlap.ratio)$merge2

    # choose the ones with significant value > sig.value.impute, which are the overlapped ones

    combined.in1 <- overlap.1[overlap.1$sig.value > sig.value.impute, ]
    combined.in2 <- overlap.2[overlap.2$sig.value > sig.value.impute, ]

    combined.region <- rbind(combined.in1, combined.in2)
  
    is.repeated <- duplicated(combined.region$start)
    combined.region <- combined.region[!is.repeated,]

    # now rank the peaks in matched region
    rank.combined <- rank(-combined.region$sig.value)
    
    # now transform the parameters estimated into the new scale
    npeaks.overlap <- nrow(combined.region)
    npeaks.consistent <- nrow(cbind(em.output[[j]]$data.pruned$sample1))

    
    # the breaks are the same for x and y
    f1 <- list(breaks=em.output[[j]]$em.fit$x.mar$f1$breaks*npeaks.overlap/npeaks.consistent, density=(em.output[[j]]$em.fit$x.mar$f1$density+em.output[[j]]$em.fit$y.mar$f1$density)/2)
    # the first break boundary goes up when changing scale, need set it back to be a bit smaller than 1
    f1$breaks[1] <- min(f1$breaks[1], 0.95)
    
    f2 <- list(breaks=em.output[[j]]$em.fit$x.mar$f2$breaks*npeaks.overlap/npeaks.consistent, density=(em.output[[j]]$em.fit$x.mar$f2$density+em.output[[j]]$em.fit$y.mar$f2$density)/2)
    # the first break boundary goes up when changing scale, need set it back to be a bit smaller than 1
    f2$breaks[1] <- min(f2$breaks[1], 0.95)
    
    p <- em.output[[j]]$em.fit$para$p
 
    # find the posterior probability
    errorprob.combined <- get.comp2.prob(rank.combined, p, f1, f2)

    # compute the FDR and find cutoff of posterior prob and the sig value
    o <- order(errorprob.combined)
    idr <- cumsum(errorprob.combined[o])/c(1:length(o))
    idr.index <- which(idr > idr.level)[1]
    errorprob.cutoff <- errorprob.combined[o][idr.index]

    # find the minimum significant measure among selected peaks
    sig.value <- min(combined.region$sig.value[o][1:idr.index])
  #  sig.value <- quantile(combined.region$sig.value[o][1:idr.index], prob=0.05)
#sig.value <- quantile(combined.region$sig.value[errorprob.combined<=errorprob.cutoff], prob=0.05)
    
    # apply the significant value on the whole pooled list
    combined.selected <- combined[combined$sig.value >= sig.value,]

    combined.selected.all <- rbind(combined.selected.all, combined.selected)
  }

  is.repeated <- duplicated(combined.selected.all$start)
  combined.selected.all <- combined.selected.all[!is.repeated,]

  npeak <- nrow(combined.selected.all)  
  
  npeak.stat <- list(idr.level=idr.level, npeak=npeak)

  sig.combined <- c(min(combined.selected.all[,"sig.value"], na.rm=T), max(combined.selected.all[,"sig.value"], na.rm=T))

 #  idr.combined <- c(min(combined.selected.all[,"q.value"], na.rm=T), max(combined.selected.all[,"q.value"], na.rm=T))
 # combined.selected.all <- deconcatenate.chr(combined.selected.all, chr.size)[,c("chr", "start", "stop", "signal.value", "p.value", "q.value")]

  combined.selected.all <- combined.selected.all[,  c("chr", "start.ori", "stop.ori", "signal.value", "p.value", "q.value")]
  colnames(combined.selected.all) <- c("chr", "start", "stop", "signal.value", "p.value", "q.value")
  
  invisible(list(npeak.stat=npeak.stat, combined.selected=combined.selected.all, sig.combined=sig.combined))
}



# get the posterior probability of the 2nd component
get.comp2.prob <- function(x, p, f1, f2){

  # get pdf and cdf of each component from functions in the corresponding component 
  px.1 <- sapply(x, get.pdf, df=f1)
  px.2 <- sapply(x, get.pdf, df=f2)

  comp2prob <- 1 - p*px.1/(p*px.1+(1-p)*px.2)
  
  return(comp2prob)
}

keep.match <- function(sample1, sample2, sig.value.impute=0){

  sample1.prune <- sample1[sample1$sig.value > sig.value.impute & sample2$sig.value > sig.value.impute,]
  sample2.prune <- sample2[sample1$sig.value > sig.value.impute & sample2$sig.value > sig.value.impute,]
 
  invisible(list(sample1=sample1.prune, sample2=sample2.prune))
}


##############################################
#
# The following is for simulation
#
##############################################


# simulate gaussian copula
# u is the uniform random variable and rho is correlation coefficient 
simu.gaussian.copula <- function(u, rho){

  n <- length(u)

  # simulate y given x=qnorm(u)
  y <- qnorm(u)*rho + rnorm(n)*sqrt(1-rho^2)

  v <- pnorm(y)

  invisible(v)
}

## simulate Clayton copula from its generating function
## Genest and MacKay (1986)

phi.ori <- function(t, s){

 (t^(-s) -1)/s
}


phi.inv <- function(y, s){

  exp(-log(s*y+1)/s)
}

phi.der <- function(t, s){

  -t^(-s-1)
}

phi.der.inv <- function(y, s){

  exp(log(-y)/(-s-1))
}

get.w <- function(u, t, s){

  phi.der.inv(phi.der(u, s)/t, s)
} 

get.v <- function(w, u, s){

  phi.inv(phi.ori(w, s) - phi.ori(u, s), s) 
}

# u is a uniform random variable, s is the association parameter
simu.clayton.copula <- function(u, s){

  t <- runif(length(u))

  if(s>0){
    w <- get.w(u, t, s)
    v <- get.v(w, u, s)
    return(v)    
  }

  if(s==0){
    return(t)
  }

  if(s <0){
    print("Invalid association parameters for clayton copula")
  }
  
}



###### 09-09-09

# simulate a two-component copula mixture:
# - marginal distributions for the two variables in each component are both 
#   normal and with the same parameters 
# p is the mixing proportion of component 1
# n is the total sample size
simu.copula.2mix <- function(s1, s2, p, n, mu1, mu2, sd1, sd2, copula.txt){

  n1 <- round(n*p)
  n2 <- n-n1

  u1 <- runif(n1)
  
  if(copula.txt =="clayton")
    v1 <- simu.clayton.copula(u1, s1)
  else{
    if(copula.txt =="gaussian")
      v1 <- simu.gaussian.copula(u1, s1)
  }

  u2 <- runif(n2)

  if(copula.txt =="clayton")
    v2 <- simu.clayton.copula(u2, s2)
  else{
    if(copula.txt =="gaussian")
      v2 <- simu.gaussian.copula(u2, s2)
  }

  # generate test statistics
  sample1.1 <- qnorm(u1, mu1, sd1)
  sample1.2 <- qnorm(v1, mu1, sd1)

  sample2.1 <- qnorm(u2, mu2, sd2)
  sample2.2 <- qnorm(v2, mu2, sd2)

  return(list(u=c(u1, u2), v=c(v1, v2), 
              u.inv=c(sample1.1, sample2.1), v.inv=c(sample1.2, sample2.2),
              label=c(rep(1, n1), rep(2, n2))))
}

# using inverse of the cdf to generate original observations 

simu.copula.2mix.inv <- function(s1, s2, p, n, cdf1.x, cdf1.y, cdf2.x, cdf2.y, copula.txt){

  n1 <- round(n*p)
  n2 <- n-n1

  u1 <- runif(n1)
  
  if(copula.txt =="clayton")
    v1 <- simu.clayton.copula(u1, s1)
  else{
    if(copula.txt =="gaussian")
      v1 <- simu.gaussian.copula(u1, s1)
  }

  u2 <- runif(n2)

  if(copula.txt =="clayton")
    v2 <- simu.clayton.copula(u2, s2)
  else{
    if(copula.txt =="gaussian")
      v2 <- simu.gaussian.copula(u2, s2)
  }

  # generate test statistics
#  sample1.1 <- qnorm(u1, mu1, sd1)
#  sample1.2 <- qnorm(v1, mu1, sd1)

#  sample2.1 <- qnorm(u2, mu2, sd2)
#  sample2.2 <- qnorm(v2, mu2, sd2)
  
  sample1.x <- inv.cdf.vec(u1, cdf1.x)
  sample1.y <- inv.cdf.vec(v1, cdf1.y)

  sample2.x <- inv.cdf.vec(u2, cdf2.x)
  sample2.y <- inv.cdf.vec(v2, cdf2.y)
  
  
  return(list(u=c(u1, u2), v=c(v1, v2), 
              u.inv=c(sample1.x, sample2.x), v.inv=c(sample1.y, sample2.y),
              label=c(rep(1, n1), rep(2, n2))))
}

# obtain original observation by converting cdf into quantiles
# u is one cdf
# u.cdf is a cdf (assuming it is a histogram) and has the break points (cdf$cdf and cdf$breaks)
# the smallest value of cdf=0 and the largest =1 
inv.cdf <- function(u, u.cdf){

  # which bin it falls into
  i <- which(u.cdf$cdf> u)[1]
  q.u  <- (u - u.cdf$cdf[i-1])/(u.cdf$cdf[i] - u.cdf$cdf[i-1])* (u.cdf$breaks[i]-u.cdf$breaks[i-1]) + u.cdf$breaks[i-1]

  return(q.u)
}

inv.cdf.vec <- function(u, u.cdf){

  # check if cdf has the right range (0, 1)  
  ncdf <- length(u.cdf$cdf)
  nbreaks <- length(u.cdf$breaks)
  
  if(ncdf == nbreaks-1 & u.cdf$cdf[ncdf]< 1)
    u.cdf[ncdf] <- 1
    
  q.u <- sapply(u, inv.cdf, u.cdf)

  return(q.u) 
}

# here we simulate a likely real situation
# the test statistics from two normal distributions
# according to their labels, then convert them into p-values w.r.t H0 using
# one-sided test.
# The test statistics are correlated for the signal component and independent
# for the noise component
# For the signal component, Y = X + eps, where eps ~ N(0, sigma^2)
simu.test.stat <- function(p, n, mu1, sd1, mu0, sd0, sd.e){

  # first component - signal
  n.signal <- round(n*p)
  n.noise <- n - n.signal

  # labels
  labels <- c(rep(1, n.signal), rep(0, n.noise))
  
  # test statistics for signal and noise
  mu.signal <- rnorm(n.signal, mu1, sd1)
  x.signal <- mu.signal + rnorm(n.signal, 0, sd.e)
  x.noise <- rnorm(n.noise, mu0, sd0) + rnorm(n.noise, 0, sd.e)

  y.signal <- mu.signal + rnorm(n.signal, 0, sd.e)
                                        # sd.e can be dependent on signal
  y.noise <- rnorm(n.noise, mu0, sd0) + rnorm(n.noise, 0, sd.e)

  # concatenate
  x <- c(x.signal, x.noise)
  y <- c(y.signal, y.noise)
  
  # convert to p-values based on H0
  p.x <- 1-pnorm(x, mu0, sqrt(sd0^2+sd.e^2))
  p.y <- 1-pnorm(y, mu0, sqrt(sd0^2+sd.e^2))

  return(list(p.x=p.x, p.y=p.y, x=x, y=y, labels=labels))
  
}

# compute the tradeoff and calibration
forward.decoy.tradeoff.ndecoy <- function(xx, labels, ndecoy){

  xx <- round(xx, 5) 
  o <- order(xx, decreasing=T)

  rand <- 1-labels # if rand==0, consistent
  # order the random indicator in the same order
  rand.o <- rand[o]

  if(sum(rand.o) > ndecoy){
    index.decoy <- which(cumsum(rand.o)==ndecoy)
  } else {
    index.decoy <- which(cumsum(rand.o)==sum(rand.o))
  }
  
  cutoff.decoy <- xx[o][index.decoy]
    
  # only consider the unique ones
  cutoff.unique <- unique(xx[o])

  cutoff <- cutoff.unique[cutoff.unique >= cutoff.decoy[length(cutoff.decoy)]]

   get.decoy.count <- function(cut.off){
     above <- rep(0, length(xx))
     above[xx >= cut.off] <- 1
     decoy.count <- sum(above==1 & rand==1)     
     return(decoy.count)
   }

   get.forward.count <- function(cut.off){
     above <- rep(0, length(xx))
     above[xx >= cut.off] <- 1
     forward.count <- sum(above==1 & rand==0)
     return(forward.count)
   }

   get.est.fdr <- function(cut.off){
     above <- rep(0, length(xx))
     above[xx >= cut.off] <- 1
     est.fdr <- 1-mean(xx[above==1])
     return(est.fdr)
   }

  # assuming rand=0 is right
   get.false.neg.count <- function(cut.off){
     below <- rep(0, length(xx))
     below[xx < cut.off] <- 1
     false.neg.count <- sum(below==1 & rand==0)
     return(false.neg.count)
   }

  get.false.pos.count <- function(cut.off){
     above <- rep(0, length(xx))
     above[xx >= cut.off] <- 1
     false.pos.count <- sum(above==1 & rand==1)
     return(false.pos.count)
   } 

   decoy <- sapply(cutoff, get.decoy.count)
   forward <- sapply(cutoff, get.forward.count)

   est.fdr <- sapply(cutoff, get.est.fdr)
   emp.fdr <- decoy/(decoy+forward)
  
   # compute specificity and sensitivity
   # assuming rand=1 is wrong and rand=0 is right
   false.neg <- sapply(cutoff, get.false.neg.count)
   false.pos <- sapply(cutoff, get.false.pos.count)
  
   true.pos <- sum(rand==0)-false.neg
   true.neg <- sum(rand==1)-false.pos
  
   sensitivity <- true.pos/(true.pos+false.neg)
   specificity <- true.neg/(true.neg+false.pos)
  
   return(list(decoy=decoy, forward=forward, cutoff=cutoff, est.fdr=est.fdr, emp.fdr=emp.fdr, sensitivity=sensitivity, specificity=specificity))  
}


# compute the em for jackknife and all data, and find FDR
get.emp.jack <- function(a, p0){

  nobs <- length(a$labels)
  est <- list()
  est.all <- list()

  temp.all <- em.transform(-a$p.x, -a$p.y, mu=1.5, sigma=1.4, rho=0.4, p=0.7, eps=0.01)
#  temp.all <- em.2copula.quick(a$p.x, a$p.y, p0=p0, rho1.0=0.7,
#      rho2.0=0, eps=0.01, fix.p=T, stoc=F, fix.rho2=T, "gaussian")

  est.all$p <- temp.all$para$p
  est.all$rho1 <- temp.all$para$rho1
  est.all$FDR <- get.FDR(temp.all$e.z)

  FDR <- list()
  p <- c()
  rho1 <- c()


  for(i in 1:nobs){

    temp <- em.transform(-a$p.x[-i], -a$p.y[-i], mu=1.5, sigma=1.4, rho=0.4, p=0.7, eps=0.01)    
#    temp <- em.2copula.quick(a$p.x[-i], a$p.y[-i], p0=p0, rho1.0=0.7,
#      rho2.0=0, eps=0.01, fix.p=T, stoc=F, fix.rho2=T, "gaussian")

    est[[i]] <- list(p=temp$para$p, rho1=temp$para$rho1, FDR=get.FDR(temp$e.z))

    FDR[[i]] <- est[[i]]$FDR # this is the FDR for top n peaks
    p[i] <- est[[i]]$p
    rho1[i] <- est[[i]]$rho1 
  }

  est.jack <- list(FDR=FDR, p=p, rho1=rho1) 
  return(list(est.jack=est.jack, est.all=est.all))
}


# get the npeaks corresponding to the nominal FDR estimated from the sample
# and find the corresponding FDR from the entire data
get.FDR.jack <- function(est, FDR.nominal){
  
  nobs <- length(est$est.jack$FDR)
  FDR.all <- c()  
  top.n <- c()
  
  for(i in 1:nobs){
    top.n[i] <- max(which(est$est.jack$FDR[[i]] <= FDR.nominal))
    FDR.all[i] <- est$est.all$FDR[top.n[i]]
  }
  
  invisible(list(FDR.all=FDR.all, top.n=top.n))
}

# compute Jackknife peudonumber
# a is the dataset
get.emp.IF <- function(a, p0){

  nobs <- length(a$labels)
  est <- list()
  est.all <- list()

  temp.all <- em.2copula.quick(a$p.x, a$p.y, p0=p0, rho1.0=0.7,
      rho2.0=0, eps=0.01, fix.p=T, stoc=F, fix.rho2=T, "gaussian")

  est.all$p <- temp.all$para$p
  est.all$rho1 <- temp.all$para$rho1
  est.all$FDR <- get.FDR(temp.all$e.z)

  IF.FDR <- list()
  IF.p <- c()
  IF.rho1 <- c()
  
  for(i in 1:nobs){
    
    temp <- em.2copula.quick(a$p.x[-i], a$p.y[-i], p0=p0, rho1.0=0.7,
      rho2.0=0, eps=0.01, fix.p=T, stoc=F, fix.rho2=T, "gaussian")

    est[[i]] <- list(p=temp$para$p, rho1=temp$para$rho1, FDR=get.FDR(temp$e.z))

    IF.FDR[[i]] <- (nobs-1)*(est.all$FDR[-nobs] - est[[i]]$FDR) # this is the FDR for top n peaks
    IF.p[i] <- (nobs-1)*(est.all$p - est[[i]]$p)
    IF.rho1[i] <- (nobs-1)*(est.all$rho1 - est[[i]]$rho1) 
  }

  emp.IF <- list(FDR=IF.FDR, p=IF.p, rho1=IF.rho1) 

  invisible(list(emp.IF=emp.IF, est.all=est.all, est=est))
}

# e.z is the posterior probability of being in signal component
get.FDR <- function(e.z){

  e.z.o <- order(1-e.z)
  FDR <- cumsum(1-e.z[e.z.o])/c(1:length(e.z.o))

  invisible(FDR)
}

# get the FDR of selecting the top n peaks
# IF.est is the sample influence function
# top.n
get.IF.FDR <- function(IF.est, top.n){

  nobs <- length(IF.est$emp.IF$FDR)
  FDR <- c()  
  
  # influence function of p
  for(i in 1:nobs)
    FDR[i] <- IF.est$emp.IF$FDR[[i]][top.n]
  
  invisible(FDR)
}

# get the sample influence function for FDR at a given FDR size
# 1. find the number of peaks selected at a given FDR computed from all obs
# 2. use the number to find the sample influence function for FDR
# IF.est$est.all is the FDR with all peaks
get.IF.FDR.all <- function(IF.est, FDR.size){

  top.n <- which.min(abs(IF.est$est.all$FDR -FDR.size))
  nobs <- length(IF.est$est.all$FDR)
  FDR <- c()  
  
  # influence function of p
  for(i in 1:nobs)
    FDR[i] <- IF.est$emp.IF$FDR[[i]][top.n]
  
  invisible(list(FDR=FDR, top.n=top.n))
}

plot.simu.uri <- function(x, y){

  tt <- seq(0.01, 0.99, by=0.01)
  uri <- sapply(tt, comp.uri.prob, u=x, v=y)
  uri.thin <- uri[seq(1, length(tt), by=3)]
  tt.thin <- tt[seq(1, length(tt), by=3)]
  duri <- (uri.thin[-1]-uri.thin[-length(uri.thin)])/(tt.thin[-1]-tt.thin[-length(tt.thin)])
  uri.spl <- smooth.spline(tt, uri, df=6.4)
  uri.der <- predict(uri.spl, tt, deriv=1)

  par(mfrow=c(2,2))
  plot(x[1:n0], y[1:n0])
  points(x[(n0+1):n], y[(n0+1):n], col=2)
  plot(rank(-x)[1:n0], rank(-y)[1:n0])
  points(rank(-x)[(1+n0):n], rank(-y)[(1+n0):n])
  plot(tt, uri)
  lines(c(0,1), c(0,1), lty=2)
  title(paste("rho1=", rho1, " rho2=", rho2, "p=", p, sep=""))
  plot(tt.thin[-1], duri)
  lines(uri.der)
  abline(h=1)
  invisible(list(x=x, y=y, uri=uri, tt=tt, duri=duri, tt.thin=tt.thin, uri.der=uri.der))

}


###### new fitting procedure




# 1. rank pairs

# 2. initialization
# 3. convert to pseudo-number

# 4. EM

# need plugin and test
# find the middle point between the bins
get.pseudo.mix <- function(x, mu, sigma, rho, p){

  
  # first compute cdf for points on the grid
  # generate 200 points between [-3, mu+3*sigma]
  nw <- 1000
  w <- seq(min(-3, mu-3*sigma), max(mu+3*sigma, 3), length=nw) 
  w.cdf <- p*pnorm(w, mean=mu, sd=sigma) + (1-p)*pnorm(w, mean=0, sd=1)

  i <- 1

  quan.x <- rep(NA, length(x))

  for(i in c(1:nw)){
    index <- which(x >= w.cdf[i] & x < w.cdf[i+1])
    quan.x[index] <- (x[index]-w.cdf[i])*(w[i+1]-w[i])/(w.cdf[i+1]-w.cdf[i]) +w[i]
  }

  index <- which(x < w.cdf[1])
  if(length(index)>0)
    quan.x[index] <- w[1]

  index <- which(x > w.cdf[nw])
  if(length(index)>0)
    quan.x[index] <- w[nw]  
  
#  linear.ext <- function(x, w, w.cdf){
  # linear interpolation
#    index.up <- which(w.cdf>= x)[1]
#    left.index <- which(w.cdf <=x)
#    index.down <- left.index[length(left.index)]
#    quan.x <- (w[index.up] + w[index.down])/2  
#  }
  
#  x.pseudo <- sapply(x, linear.ext, w=w, w.cdf=w.cdf)

#  invisible(x.pseudo)
  invisible(quan.x)
}


# EM to compute the latent structure
# steps:
# 1. raw values are first transformed into pseudovalues
# 2. EM is used to compute the underlining structure, which is a mixture
#    of two normals
em.transform <- function(x, y, mu, sigma, rho, p, eps){
  
  x.cdf.func <- ecdf(x)
  y.cdf.func <- ecdf(y)
  afactor <- length(x)/(length(x)+1)
  x.cdf <- x.cdf.func(x)*afactor
  y.cdf <- y.cdf.func(y)*afactor
  
  # initialization
  para <- list()
  para$mu <- mu
  para$sigma <- sigma
  para$rho <- rho
  para$p <- p  

  j <- 1
  to.run <- T
  loglik.trace <- c()
  loglik.inner.trace <- c()
  
  #to.run.inner <- T
  z.1 <- get.pseudo.mix(x.cdf, para$mu, para$sigma, para$rho, para$p)
  z.2 <- get.pseudo.mix(y.cdf, para$mu, para$sigma, para$rho, para$p)

#  cat("length(z1)", length(z.1), "\n")
  while(to.run){
    
    # get pseudo value in each cycle
#    z.1 <- get.pseudo.mix(x.cdf, para$mu, para$sigma, para$rho, para$p)
#    z.2 <- get.pseudo.mix(y.cdf, para$mu, para$sigma, para$rho, para$p)

    i <- 1
    while(to.run){
      
      # EM for latent structure
      e.z <- e.step.2normal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)
      para <- m.step.2normal(z.1, z.2, e.z)
#para$rho <- rho
#para$p <- p    
#para$mu <- mu
#para$sigma <- sigma    
      if(i > 1)
        l.old <- l.new
    
      # this is just the mixture likelihood of two-component Gaussian
      l.new <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)

      loglik.inner.trace[i] <- l.new 

      if(i > 1){
        to.run <- loglik.inner.trace[i]-loglik.inner.trace[i-1]>eps         
      }
        
    
#      if(i > 2){
#        l.inf <- loglik.inner.trace[i-2] + (loglik.inner.trace[i-1] - loglik.inner.trace[i-2])/(1-(loglik.inner.trace[i]-loglik.inner.trace[i-1])/(loglik.inner.trace[i-1]-loglik.inner.trace[i-2]))

#        if(loglik.inner.trace[i-1]!=loglik.inner.trace[i-2])
#          to.run <- abs(l.inf - loglik.inner.trace[i]) > eps
#        else
#          to.run <- F
          
#      }

      cat("loglik.inner.trace[", i, "]=", loglik.inner.trace[i], "\n")
    cat("mu=", para$mu, "sigma=", para$sigma, "p=", para$p, "rho=", para$rho, "\n\n")
      
      i <- i+1
    }
    

    # get pseudo value in each cycle
    z.1 <- get.pseudo.mix(x.cdf, para$mu, para$sigma, para$rho, para$p)
    z.2 <- get.pseudo.mix(y.cdf, para$mu, para$sigma, para$rho, para$p)

    if(j > 1)
      l.old.outer <- l.new.outer

    l.new.outer <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)

    loglik.trace[j] <- l.new.outer
    
    if(j == 1)
      to.run <- T
    else{ # stop when iteration>100
      if(j > 100)
        to.run <- F
      else
        to.run <- l.new.outer - l.old.outer > eps
    }

#    if(j %% 10==0)
      cat("loglik.trace[", j, "]=", loglik.trace[j], "\n")
    cat("mu=", para$mu, "sigma=", para$sigma, "p=", para$p, "rho=", para$rho, "\n")
    
    j <- j+1
  }

  bic <- -2*l.new + 4*log(length(z.1))
  
  return(list(para=list(p=para$p, rho=para$rho, mu=para$mu, sigma=para$sigma),
              loglik=l.new, bic=bic, e.z=e.z, loglik.trace=loglik.trace))
}  




# compute log-likelihood for mixture of two bivariate normals
loglik.2binormal <- function(z.1, z.2, mu, sigma, rho, p){

  l.m <- sum(d.binormal(z.1, z.2, 0, 1, 0)+log(p*exp(d.binormal(z.1, z.2, mu, sigma, rho)-d.binormal(z.1, z.2, 0, 1, 0))+(1-p)))
  
#  l.m <- sum((p*d.binormal(z.1, z.2, mu, sigma, rho) + (1-p)*d.binormal(z.1, z.2, 0, 1, 0)))
  return(l.m) 
}

# check this when rho=1

# density of binomial distribution with equal mean and sigma on both dimensions
d.binormal <- function(z.1, z.2, mu, sigma, rho){

  loglik <- (-log(2)-log(pi)-2*log(sigma) - log(1-rho^2)/2 - (0.5/(1-rho^2)/sigma^2)*((z.1-mu)^2 -2*rho*(z.1-mu)*(z.2-mu) + (z.2-mu)^2))

  return(loglik)
}

# E-step for computing the latent strucutre
# e.z is the prob to be in the consistent group
# e.step for estimating posterior prob
# z.1 and z.2 can be vectors or scalars
e.step.2normal <- function(z.1, z.2, mu, sigma, rho, p){

  e.z <- p/((1-p)*exp(d.binormal(z.1, z.2, 0, 1, 0)-d.binormal(z.1, z.2, mu, sigma, rho))+ p)
  
  invisible(e.z)
}

# M-step for computing the latent structure
# m.step for estimating proportion, mean, sd and correlation coefficient
m.step.2normal <- function(z.1, z.2, e.z){

  p <- mean(e.z)
  mu <- sum((z.1+z.2)*e.z)/2/sum(e.z) 
  sigma <- sqrt(sum(e.z*((z.1-mu)^2+(z.2-mu)^2))/2/sum(e.z))
  rho <- 2*sum(e.z*(z.1-mu)*(z.2-mu))/(sum(e.z*((z.1-mu)^2+(z.2-mu)^2)))

  return(list(p=p, mu=mu, sigma=sigma, rho=rho))
}


# assume top p percent of observations are true
# x and y are ranks, estimate
init <- function(x, y, x.label){
  
  x.o <- order(x)

  x.ordered <- x[x.o]
  y.ordered <- y[x.o]
  x.label.ordered <- x.label[x.o]
  
  n <- length(x)
  p <- sum(x.label)/n
  
  rho <- cor(x.ordered[1:ceiling(p*n)], y.ordered[1:ceiling(p*n)])

  temp <- find.mu.sigma(x.ordered, x.label.ordered)
  mu <- temp$mu
  sigma <- temp$sigma
  
  invisible(list(mu=mu, sigma=sigma, rho=rho, p=p))

}

# find mu and sigma if the distributions of marginal ranks are known
# take the medians of the two dist and map back to the original
init.dist <- function(f0, f1){

  # take the median in f0
  index.median.0 <- which(f0$cdf>0.5)[1]
  q.0.small <- f0$cdf[index.median.0] # because f0 and f1 have the same bins
  q.1.small <- f1$cdf[index.median.0]

  # take the median in f1
  index.median.1 <- which(f1$cdf>0.5)[1]
  q.0.big <- f0$cdf[index.median.1] # because f0 and f1 have the same bins
  q.1.big <- f1$cdf[index.median.1]

  # find pseudo value for x.middle[1] on normal(0,1) 
  pseudo.small.0 <- qnorm(q.0.small, mean=0, sd=1)
  pseudo.small.1 <- qnorm(q.1.small, mean=0, sd=1)

  # find pseudo value for x.middle[2] on normal(0,1) 
  pseudo.big.0 <- qnorm(q.0.big, mean=0, sd=1)
  pseudo.big.1 <- qnorm(q.1.big, mean=0, sd=1)

  mu <- (pseudo.small.0*pseudo.big.1 - pseudo.small.1*pseudo.big.0)/(pseudo.big.1-pseudo.small.1) 

  sigma <- (pseudo.small.0-mu)/pseudo.small.1

  return(list(mu=mu, sigma=sigma))  
}

# generate labels

# find the part of data with overlap

# find the percentile on noise and signal

# Suppose there are signal and noise components, with mean=0 and sd=1 for noise
# x and x.label are the rank of the observations and their labels,
# find the mean and sd of the other component
# x.label takes values of 0 and 1
find.mu.sigma <- function(x, x.label){

  x.0 <- x[x.label==0]
  x.1 <- x[x.label==1]

  n.x0 <- length(x.0)
  n.x1 <- length(x.1)

  x.end <- c(min(x.0), min(x.1), max(x.0), max(x.1))
  o <- order(x.end)
  x.middle <- x.end[o][c(2,3)]

  # the smaller end of the overlap
  q.1.small <- mean(x.1 <= x.middle[1])*n.x1/(n.x1+1)
  q.0.small <- mean(x.0 <= x.middle[1])*n.x0/(n.x0+1)

  # the bigger end of the overlap
  q.1.big <- mean(x.1 <= x.middle[2])*n.x1/(n.x1+1)
  q.0.big <- mean(x.0 <= x.middle[2])*n.x0/(n.x0+1)

  # find pseudo value for x.middle[1] on normal(0,1) 
  pseudo.small.0 <- qnorm(q.0.small, mean=0, sd=1)
  pseudo.small.1 <- qnorm(q.1.small, mean=0, sd=1)

  # find pseudo value for x.middle[2] on normal(0,1) 
  pseudo.big.0 <- qnorm(q.0.big, mean=0, sd=1)
  pseudo.big.1 <- qnorm(q.1.big, mean=0, sd=1)

  mu <- (pseudo.small.0*pseudo.big.1 - pseudo.small.1*pseudo.big.0)/(pseudo.big.1-pseudo.small.1) 

  sigma <- (pseudo.small.0-mu)/pseudo.small.1

  return(list(mu=mu, sigma=sigma))
}
