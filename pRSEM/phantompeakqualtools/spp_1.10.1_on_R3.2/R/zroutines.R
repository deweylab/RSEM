#library(caTools)
#dyn.load("src/bed2vector.so");
#dyn.load("src/wdl.so");
#dyn.load("src/peaks.so");
#dyn.load("src/cdensum.so");


# -------- ROUTINES FOR READING IN THE DATA FILES ------------
# fix.chromosome.names : remove ".fa" suffix from match sequence names
read.eland.tags <- function(filename,read.tag.names=F,fix.chromosome.names=T,max.eland.tag.length=-1,extended=F,multi=F) {
  if(read.tag.names) { rtn <- as.integer(1); } else { rtn <- as.integer(0); };
  storage.mode(max.eland.tag.length) <- "integer";
  callfunction <- "read_eland";
  if(extended) { callfunction <- "read_eland_extended"; };
  if(multi) { callfunction <- "read_eland_multi"; };
  tl <- lapply(.Call(callfunction,filename,rtn,max.eland.tag.length),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    if(read.tag.names) {
      d$s <- d$s[xo];
    }
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  if(read.tag.names) {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),names=lapply(tl,function(d) d$s)));
  } else {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
  }
}

read.tagalign.tags <- function(filename,fix.chromosome.names=T,fix.quality=T) {
  tl <- lapply(.Call("read_tagalign",filename),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    #if(fix.quality) {
    #  d$n <- 4-cut(d$n,breaks=c(0,250,500,750,1000),labels=F)
    #}
    if(fix.quality) { # Anshul: changed the way the quality field is processed
      if (min(d$n)<0.5){
        d$n = ceiling(1000/4^d$n);
      }
      break.vals <- unique(sort(c(0,unique(d$n))));
      d$n <- length(break.vals)-1-cut(d$n,breaks=break.vals,labels=F);
    }    
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
}


read.short.arachne.tags <- function(filename,fix.chromosome.names=F) {
  tl <- lapply(.Call("read_arachne",filename),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
}


read.arachne.tags <- function(filename,fix.chromosome.names=F) {
  tl <- lapply(.Call("read_arachne_long",filename),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    d$l <- d$l[xo];
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),length=lapply(tl,function(d) d$l)));
}

read.bowtie.tags <- function(filename,read.tag.names=F,fix.chromosome.names=F) {
  if(read.tag.names) { rtn <- as.integer(1); } else { rtn <- as.integer(0); };
  tl <- lapply(.Call("read_bowtie",filename,rtn),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    if(read.tag.names) {
      d$s <- d$s[xo];
    }
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  if(read.tag.names) {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),names=lapply(tl,function(d) d$s)));
  } else {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
  }
}

read.bam.tags <- function(filename,read.tag.names=F,fix.chromosome.names=F) {
  if(read.tag.names) { rtn <- as.integer(1); } else { rtn <- as.integer(0); };
  tl <- lapply(.Call("read_bam",filename,rtn),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    if(read.tag.names) {
      d$s <- d$s[xo];
    }
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  if(read.tag.names) {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),names=lapply(tl,function(d) d$s)));
  } else {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
  }
}


read.helicos.tags <- function(filename,read.tag.names=F,fix.chromosome.names=F,include.length.info=T) {
  if(read.tag.names) { rtn <- as.integer(1); } else { rtn <- as.integer(0); };
  tl <- lapply(.Call("read_helicostabf",filename,rtn),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    d$l <- d$l[xo];
    if(read.tag.names) {
      d$s <- d$s[xo];
    }
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  if(read.tag.names) {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),length=lapply(tl,function(d) d$l),names=lapply(tl,function(d) d$s)));
  } else {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),length=lapply(tl,function(d) d$l)));
  }
}

read.maqmap.tags <- function(filename,read.tag.names=F,fix.chromosome.names=T) {
  if(read.tag.names) { rtn <- as.integer(1); } else { rtn <- as.integer(0); };
  tl <- lapply(.Call("read_maqmap",filename,rtn),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    if(read.tag.names) {
      d$s <- d$s[xo];
    }
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  if(read.tag.names) {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),names=lapply(tl,function(d) d$s)));
  } else {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
  }
}


read.bin.maqmap.tags <- function(filename,read.tag.names=F,fix.chromosome.names=T) {
  if(read.tag.names) { rtn <- as.integer(1); } else { rtn <- as.integer(0); };
  tl <- lapply(.Call("read_binmaqmap",filename,rtn),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    if(read.tag.names) {
      d$s <- d$s[xo];
    }
    return(d);
  });
  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  if(read.tag.names) {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n),names=lapply(tl,function(d) d$s)));
  } else {
    return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
  }
}


# read in tags from an extended eland format with match length information
read.meland.tags <- function(filename,read.tag.names=F,fix.chromosome.names=T) {
  if(read.tag.names) { rtn <- as.integer(1); } else { rtn <- as.integer(0); };
  tl <- lapply(.Call("read_meland",filename,rtn),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    d$l <- d$l[xo];
    if(read.tag.names) {
      d$s <- d$s[xo];
    }
    return(d);
  });

  if(fix.chromosome.names) {
    # remove ".fa"
    names(tl) <- gsub("\\.fa","",names(tl))
  }
  # separate tags and quality
  chrl <- names(tl); names(chrl) <- chrl;
  # reformulate quality scores into monotonic integers
  ml <- max(unlist(lapply(tl,function(d) max(d$l))));
  qual <- lapply(chrl,function(chr) (ml-tl[[chr]]$l)+tl[[chr]]$n/10);
  if(read.tag.names) {
    return(list(tags=lapply(tl,function(d) d$t),quality=qual,names=lapply(tl,function(d) d$s)));
  } else {
    return(list(tags=lapply(tl,function(d) d$t),quality=qual));
  }
}

# -------- ROUTINES FOR ASSESSING BINDING PATTERN AND SELECTING INFORMATIVE TAGS  ------------

# removes tag positions that have anomalously high counts on both strands
# z - z-score used to determine anomalous bins
# zo - z used to filter out one-strand matches
# trim.fraction - fraction of top bins to discard when calculating overall background density
remove.tag.anomalies <- function(data, bin=1,trim.fraction=1e-3,z=5,zo=3*z) {
  
  t.remove.tag.anomalies <- function(tv,bin=1,trim.fraction=1e-3,z=5,zo=3*z,return.indecies=F) {
    tt <- table(floor(tv/bin));

    # trim value
    stt <- sort(as.numeric(tt));
    stt <- stt[1:(length(stt)*(1-trim.fraction))];
    mtc <- mean(stt); tcd <- sqrt(var(stt));

    thr <- max(1,ceiling(mtc+z*tcd));
    thr.o <- max(1,ceiling(mtc+zo*tcd));
    # filter tt
    tt <- tt[tt>=thr]
    # get + and - tags
    tp <- as.numeric(names(tt));
    pti <- tp>0;
    it <- intersect(tp[pti],(-1)*tp[!pti]);
    # add one-strand matches
    it <- unique(c(it,tp[tt>=thr.o]));
    sit <- c(it,(-1)*it);
    
    if(bin>1) {
      sit <- sit*bin;
      sit <- c(sit,unlist(lapply(1:bin,function(i) sit+i)))
    }
    if(return.indecies) {
      return(!tv %in% sit);
    } else {
      return(tv[!tv %in% sit]);
    }
  }

  vil <- lapply(data$tags,t.remove.tag.anomalies,return.indecies=T,bin=bin,trim.fraction=trim.fraction,z=z,zo=zo);
  chrl <- names(data$tags); names(chrl) <- chrl;
  data$tags <- lapply(chrl,function(chr) data$tags[[chr]][vil[[chr]]]);
  # count tags to remove empty chromosomes
  nt <- unlist(lapply(data$tags,length));
  if(any(nt==0)) {
    data$tags <- data$tags[nt!=0]
  }
  
  if(!is.null(data$quality)) {
    data$quality <- lapply(chrl,function(chr) data$quality[[chr]][vil[[chr]]]);
    data$quality <- data$quality[nt!=0];
  }
  if(!is.null(data$names)) {
    data$names <- lapply(chrl,function(chr) data$names[[chr]][vil[[chr]]]);
    data$names <- data$names[nt!=0];
  }
  
  return(data);
}

# caps or removes tag positions that are significantly higher than local background
remove.local.tag.anomalies <- function(tags,window.size=200,eliminate.fold=10,cap.fold=4,z.threshold=3) {
  lapply(tags,filter.singular.positions.by.local.density,window.size=2e2,eliminate.fold=10,cap.fold=4,z.threshold=3);
}



# assess strand cross-correlation, determine peak position, determine appropriate window size
# for binding detection.
get.binding.characteristics <- function(data,srange=c(50,500),bin=5,cluster=NULL,debug=F,min.tag.count=1e3,acceptance.z.score=3,remove.tag.anomalies=T,anomalies.z=5,accept.all.tags=F) {
  if(remove.tag.anomalies) {
    data <- remove.tag.anomalies(data,z=anomalies.z);
  }
  
  # take highest quality tag bin
  if(!is.null(data$quality) & !accept.all.tags) {
    min.bin <- min(unlist(lapply(data$quality,min)))
    chrl <- names(data$tags); names(chrl) <- chrl;
    otl <- lapply(chrl,function(chr) data$tags[[chr]][data$quality[[chr]]==min.bin]);
  } else {
    otl <- data$tags;
  }
  # remove empty chromosomes
  otl <- otl[unlist(lapply(otl,length))!=0];


  # calculate strand scc
  if(!is.null(cluster)) {
    cc <- clusterApplyLB(cluster,otl,tag.scc,srange=srange,bin=bin);
    names(cc) <- names(otl); 
  } else {
    cc <- lapply(otl,tag.scc,srange=srange,bin=bin);
  }
  ccl<-list(sample=cc);
  ccl.av <- lapply(names(ccl),t.plotavcc,type='l',ccl=ccl,return.ac=T,ttl=list(sample=otl),plot=F)[[1]]
  ccl.av <- data.frame(x=as.numeric(names(ccl.av)),y=as.numeric(ccl.av));
  
  # find peak
  pi <- which.max(ccl.av$y);
  
  # determine width at third-height
  th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/3+ccl.av$y[length(ccl.av$y)]
  whs <- max(ccl.av$x[ccl.av$y>=th]);
  
  if (! is.integer(whs)) { # Anshul: added this to avoid situations where whs ends up being -Inf
  	whs <- ccl.av$x[ min(c(2*pi,length(ccl.av$y))) ]
  }

  # determine acceptance of different quality bins
  
  # calculates tag scc for the best tags, and combinations of best tag category with every other category
  # for subsequent selection of acceptable categories
  scc.acceptance.calc <- function() {

    qr <- range(unlist(lapply(data$quality,range)))

    # start with best tags

    # determine half-width for scc calculations
    pi <- which.max(ccl.av$y);

    # determine width at half-height
    th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/2+ccl.av$y[length(ccl.av$y)]
    lwhs <- max(ccl.av$x[ccl.av$y>=th])-ccl.av$x[pi];
    lwhs <- max(c(20,bin*10,lwhs));
    srange <- ccl.av$x[pi]+c(-lwhs,lwhs)

    # calculate chromosome-average scc
    t.scc <- function(tags) {
      if(is.null(cluster)) {
        cc <- lapply(tags,tag.scc,srange=srange,bin=bin);
      } else {
        cc <- clusterApplyLB(cluster,tags,tag.scc,srange=srange,bin=bin); names(cc) <- names(tags);
      }
      return(t.plotavcc(1,type='l',ccl=list(cc),ttl=list(tags),plot=F,return.ac=T))
    }


    # returns info list for a given tag length (lv), mismatch count (nv)
    t.cat <- function(qual) {
      # construct tag set
      if(qual==qr[1]) {
        ts <- otl;
      } else {
        nts <- names(otl); names(nts) <- nts;
        # select tags
        at <- lapply(nts,function(chr) data$tags[[chr]][data$quality[[chr]]==qual]);
        ntags <- sum(unlist(lapply(at,length)));
        if(ntags<min.tag.count) { return(NULL); }

        # append to otl
        ts <- lapply(nts,function(nam) c(otl[[nam]],at[[nam]]));
      }

      return(t.scc(ts));
    }


    # calculate cross-correlation values for each quality bin
    ql <- sort(unique(unlist(lapply(data$quality,unique)))); names(ql) <- ql;

    qccl <- lapply(ql,t.cat);

    # acceptance tests
    ac <- c(T,unlist(lapply(qccl[-1],function(d) if(is.null(d)) { return(F) } else { t.test(d-qccl[[as.character(min.bin)]],alternative="greater")$p.value<pnorm(acceptance.z.score,lower.tail=F) }))); names(ac) <- names(qccl);
    return(list(informative.bins=ac,quality.cc=qccl))
  }

  if(accept.all.tags | is.null(data$quality)) {
    return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs))    
  } else {
    acc <- scc.acceptance.calc();
    return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs,quality.bin.acceptance=acc));
  }

}


# select a set of informative tags based on the pre-calculated binding characteristics
select.informative.tags <- function(data,binding.characteristics=NULL) {
  if(is.null(binding.characteristics)) {
    return(data$tags);
  }
  if(is.null(binding.characteristics$quality.bin.acceptance)) {
    cat("binding characteristics doesn't contain quality selection info, accepting all tags\n");
    return(data$tags);
  }

  ib <- binding.characteristics$quality.bin.acceptance$informative.bins;
  abn <- names(ib)[ib]

  chrl <- names(data$tags); names(chrl) <- chrl;
  lapply(chrl,function(chr) {
    data$tags[[chr]][as.character(data$quality[[chr]]) %in% abn]
  })
}

# -------- ROUTINES FOR CALLING BINDING POSITIONS  ------------

# determine binding positions
# signal.data - IP tag lists
# control.data - input tag lists
# e.value - desired E-value threshold (either E-value or FDR threshold must be provided)
# fdr - desired FDR threshold
# min.dist - minimal distance between detected positions
# tag.count.whs - size of the window to be used to estimate confidence interval of the peak fold enrichment ratios
# enrichmnent.z - Z-score defining the desired confidence level for enrichment interval estimates
# enrichment.background.scales - define how many tiems larger should be the window for estimating background
#                                tag density when evaluating peak enrichment confidence intervals.
#                                If multiple values are given, multiple independent interval estimates will be
#                                calculated.
# tec.filter - whether to mask out the regions that exhibit significant background enrichment
# tec.window.size, tec.z - window size and Z-score for maksing out significant background enrichment regions
#
# If the control.data is not provided, the method will assess significance of the determined binding positions
# based on the randomizations of the original data. The following paramters control such randomizations:
# n.randomizations - number of randomizations to be performed
# shuffle.window - size of the bin that defines the tags that are kept together during randomization.
#                  value of 0 means that all tags are shuffled independently
#
# Binding detection methods: 
# tag.wtd - default method.
#           must specify parameter "whs", which is the half-size of the window used to calculate binding scores
# tag.lwcc - LWCC method;
#           must specify whs - a size of the window used to calculate binding scores
#           can specify isize (default=15bp) - size of the internal window that is masked out
find.binding.positions <- function(signal.data,f=1,e.value=NULL,fdr=NULL, masked.data=NULL,control.data=NULL,whs=200,min.dist=200,window.size=4e7,cluster=NULL,debug=T,n.randomizations=3,shuffle.window=1,min.thr=2,topN=NULL, tag.count.whs=100, enrichment.z=2, method=tag.wtd, tec.filter=T,tec.window.size=1e4,tec.z=5,tec.masking.window.size=tec.window.size, tec.poisson.z=5,tec.poisson.ratio=5, tec=NULL, n.control.samples=1, enrichment.scale.down.control=F, enrichment.background.scales=c(1,5,10), use.randomized.controls=F, background.density.scaling=T, mle.filter=F, min.mle.threshold=1, ...) {

  if(f<1) {
    if(debug) { cat("subsampling signal ... "); }
    signal.data <- lapply(signal.data,function(x) sample(x,length(x)*f))
    if(debug) {  cat("done\n"); }
  }


  if(!is.null(control.data) & !use.randomized.controls) {
    # limit both control and signal data to a common set of chromosomes
    chrl <- intersect(names(signal.data),names(control.data));
    signal.data <- signal.data[chrl];
    control.data <- control.data[chrl];
    control <- list(control.data);
  } else {
    control <- NULL;
  }
  
  prd <- lwcc.prediction(signal.data,min.dist=min.dist,whs=whs,window.size=window.size,e.value=e.value,fdr=fdr,debug=debug,n.randomizations=n.randomizations,shuffle.window=shuffle.window,min.thr=min.thr,cluster=cluster,method=method,bg.tl=control.data,mask.tl=masked.data, topN=topN, control=control,tec.filter=tec.filter,tec.z=tec.z,tec.window.size=tec.window.size, tec.masking.window.size=tec.masking.window.size, tec.poisson.z=tec.poisson.z,tec.poisson.ratio=tec.poisson.ratio, background.density.scaling=background.density.scaling, ...);

  # add tag counts
  chrl <- names(prd$npl); names(chrl) <- chrl;
  prd$npl <- lapply(chrl,function(chr) {
    pd <- prd$npl[[chr]];
    pd$nt <- points.within(abs(signal.data[[chr]]),pd$x-tag.count.whs,pd$x+tag.count.whs,return.point.counts=T);
    return(pd);
  });
  prd$f <- f;
  prd$n <- sum(unlist(lapply(signal.data,length)));
  if(!is.null(control.data)) {
    prd$n.bg <- sum(unlist(lapply(control.data,length)));
  }
  
  # calculate enrichment ratios
  prd <- calculate.enrichment.estimates(prd,signal.data,control.data=control.data,fraction=1,tag.count.whs=tag.count.whs,z=enrichment.z,scale.down.control=enrichment.scale.down.control,background.scales=enrichment.background.scales);

  if(mle.filter) {
    if(!is.null(prd$npl)) {
      if(length(prd$npl)>1) {
        mle.columns <- grep("enr.mle",colnames(prd$npl[[1]]));
        if(length(mle.columns)>1) {
          prd$npl <- lapply(prd$npl,function(d) d[apply(d[,mle.columns],1,function(x) all(x>min.mle.threshold)),])
        }
      }
    }
  }

  prd$whs <- whs;

  return(prd);
}



# -------- ROUTINES FOR WRITING OUT TAG DENSITY AND ENRICHMENT PROFILES  ------------
# calculate smoothed tag density, optionally subtracting the background
get.smoothed.tag.density <- function(signal.tags,control.tags=NULL,bandwidth=150,bg.weight=NULL,tag.shift=146/2,step=round(bandwidth/3),background.density.scaling=T,rngl=NULL,scale.by.dataset.size=F) {
  chrl <- names(signal.tags); names(chrl) <- chrl;

  if(!is.null(control.tags)) {
    bg.weight <- dataset.density.ratio(signal.tags,control.tags,background.density.scaling=background.density.scaling);
  }

  if(scale.by.dataset.size) {
    den.scaling <- dataset.density.size(signal.tags,background.density.scaling=background.density.scaling)/1e6;
  } else {
    den.scaling <- 1;
  }
  
  lapply(chrl,function(chr) {
    ad <- abs(signal.tags[[chr]]+tag.shift);
    rng <- NULL;
    if(!is.null(rngl)) {
      rng <- rngl[[chr]];
    }
    if(is.null(rng)) {
      rng <- range(ad);
    }

    ds <- densum(ad,bw=bandwidth,from=rng[1],to=rng[2],return.x=T,step=step);
    if(!is.null(control.tags)) {
      if(!is.null(control.tags[[chr]])) {
        bsd <- densum(abs(control.tags[[chr]]+tag.shift),bw=bandwidth,from=rng[1],to=rng[2],return.x=F,step=step);
        ds$y <- ds$y-bsd*bg.weight;
      }
    }
    return(data.frame(x=seq(ds$x[1],ds$x[2],by=step),y=den.scaling*ds$y))
  })
}

# get smoothed maximum likelihood estimate of the log2 signal to control enrichment ratio
get.smoothed.enrichment.mle <- function(signal.tags, control.tags, tag.shift=146/2, background.density.scaling=F, pseudocount=1,bg.weight=NULL,  ... ) {
  # determine common range
  chrl <- intersect(names(signal.tags),names(control.tags)); names(chrl) <- chrl;
  rngl <- lapply(chrl,function(chr) range(c(range(abs(signal.tags[[chr]]+tag.shift)),range(abs(control.tags[[chr]]+tag.shift)))))
  ssd <- get.smoothed.tag.density(signal.tags, rngl=rngl, ..., scale.by.dataset.size=F)
  csd <- get.smoothed.tag.density(control.tags, rngl=rngl, ..., scale.by.dataset.size=F)
  if(is.null(bg.weight)) {
    bg.weight <- dataset.density.ratio(signal.tags,control.tags,background.density.scaling=background.density.scaling);
  }
  cmle <- lapply(chrl,function(chr) { d <- ssd[[chr]]; d$y <- log2(d$y+pseudocount) - log2(csd[[chr]]$y+pseudocount) - log2(bg.weight); return(d); })
}


# returns a conservative upper/lower bound profile (log2) given signal tag list, background tag list and window scales
get.conservative.fold.enrichment.profile <- function(ftl,btl,fws,bwsl=c(1,5,25,50)*fws,step=50,tag.shift=146/2,alpha=0.05,use.most.informative.scale=F,quick.calculation=T,background.density.scaling=T,bg.weight=NULL,posl=NULL,return.mle=F) {
  # include only chromosomes with more than 2 reads
  ftl <- ftl[unlist(lapply(ftl,length))>2]
  chrl <- names(ftl); names(chrl) <- chrl;
  if(!is.null(posl)) {
    chrl <- chrl[chrl %in% names(posl)];
  }
  # calculate background tag ratio
  if(is.null(bg.weight)) {
    bg.weight <- dataset.density.ratio(ftl,btl,background.density.scaling=background.density.scaling);
  }
  lapply(chrl,function(chr) {
    if(is.null(btl[[chr]])) { bt <- c(); } else { bt <- abs(btl[[chr]]+tag.shift); }
    if(is.null(posl)) {
      x <- mbs.enrichment.bounds(abs(ftl[[chr]]+tag.shift),bt,fws=fws,bwsl=bwsl,step=step,calculate.upper.bound=T,bg.weight=bg.weight,use.most.informative.scale=use.most.informative.scale,quick.calculation=quick.calculation,alpha=alpha);
    } else {
      x <- mbs.enrichment.bounds(abs(ftl[[chr]]+tag.shift),bt,fws=fws,bwsl=bwsl,step=step,calculate.upper.bound=T,bg.weight=bg.weight,use.most.informative.scale=use.most.informative.scale,quick.calculation=quick.calculation,alpha=alpha,pos=posl[[chr]]);
    }
    # compose profile showing lower bound for enriched, upper bound for depleted regions
    ps <- rep(1,length(x$mle));
    vi <- which(!is.na(x$lb) & x$lb>1);
    ps[vi] <- x$lb[vi];
    vi <- which(!is.na(x$ub) & x$ub<1);
    ps[vi] <- x$ub[vi];
    ps <- log2(ps);
    if(is.null(posl)) {
      if(return.mle) {
        return(data.frame(x=seq(x$x$s,x$x$e,by=x$x$step),y=ps,mle=log2(x$mle),lb=log2(x$lb),ub=log2(x$ub)));
      } else {
        return(data.frame(x=seq(x$x$s,x$x$e,by=x$x$step),y=ps));
      }
    } else {
      if(return.mle) {
        return(data.frame(x=posl[[chr]],y=ps,mle=log2(x$mle),lb=log2(x$lb),ub=log2(x$ub)));
      } else {
        return(data.frame(x=posl[[chr]],y=ps));
      }
    }
  })
}


# write a per-chromosome $x/$y data structure into a wig file
writewig <- function(dat,fname,feature,threshold=5,zip=F) {
  chrl <- names(dat); names(chrl) <- chrl;
  invisible(lapply(chrl,function(chr) {
    bdiff <- dat[[chr]];
    ind <- seq(1,length(bdiff$x));
    ind <- ind[!is.na(bdiff$y[ind])];
    header <- chr==chrl[1];
    write.probe.wig(chr,bdiff$x[ind],bdiff$y[ind],fname,append=!header,feature=feature,header=header);
  }))
  if(zip) {
    zf <- paste(fname,"zip",sep=".");
    system(paste("zip \"",zf,"\" \"",fname,"\"",sep=""));
    system(paste("rm \"",fname,"\"",sep=""));
    return(zf);
  } else {
    return(fname);
  }
}



# -------- ROUTINES FOR ANALYZING SATURATION PROPERTIES  ------------

# PUBLIC
# calculate minimal saturation enrichment ratios (MSER) 
get.mser <- function(signal.data,control.data,n.chains=5,step.size=1e5, chains=NULL, cluster=NULL, test.agreement=0.99, return.chains=F, enrichment.background.scales=c(1), n.steps=1, ...) {
  if(is.null(chains)) {
    ci <- c(1:n.chains); names(ci) <- ci;
    if(is.null(cluster)) {
      chains <- lapply(ci,get.subsample.chain.calls,signal.data=signal.data,control.data=control.data,n.steps=n.steps,step.size=step.size,subsample.control=F, enrichment.background.scales=enrichment.background.scales, ...);
    } else {
      chains <- clusterApplyLB(cluster,ci,get.subsample.chain.calls,signal.data=signal.data,control.data=control.data,n.steps=n.steps,step.size=step.size,subsample.control=F, enrichment.background.scales=enrichment.background.scales, ...);
      names(chains) <- ci;
    }
  }
  cvl <- mser.chain.interpolation(chains=chains,enrichment.background.scales=enrichment.background.scales,test.agreement=test.agreement,return.lists=F);
  if(n.steps>1) {
    msers <- cvl;
  } else {
    msers <- unlist(lapply(cvl,function(d) d$me))
  }
  if(return.chains) {
    return(list(mser=msers,chains=chains));
  } else {
    return(msers);
  }
}

# PUBLIC 
# interpolate MSER dependency on tag counts
get.mser.interpolation <- function(signal.data,control.data,target.fold.enrichment=5,n.chains=10,n.steps=6,step.size=1e5, chains=NULL,  test.agreement=0.99, return.chains=F, enrichment.background.scales=c(1), excluded.steps=c(seq(2,n.steps-2)), ...) {
  msers <- get.mser(signal.data,control.data,n.chains=n.chains,n.steps=n.steps,step.size=step.size,chains=chains,test.agrement=test.agreement,return.chains=T,enrichment.background.scales=enrichment.background.scales,excluded.steps=excluded.steps, ...);

  # adjust sizes in case a subset of chromosomes was used
  mser <- mser.chain.interpolation(chains=msers$chains,enrichment.background.scales=enrichment.background.scales,test.agreement=test.agreement,return.lists=T);
  sr <- sum(unlist(lapply(signal.data,length)))/mser[[1]][[1]]$n[1];

  # Subsampling each chain requires removing a fraction of each chromosome's
  # tag list.  To get the exact step.size, this often leaves chromosomes with
  # a non-integer number of tags.  The non-integer values are floored, so each
  # chr can contribute at most 0.999.. <= 1 error to the step.size.
  floor.error <- length(msers$chains[[1]][[1]]$npl)
  intpn <- lapply(mser,function(ms) {
    lmvo <- do.call(rbind,ms)
    lmvo$n <- lmvo$n*sr;
    # Don't select rows corresponding to excluded.steps
    # Keep in mind that nd values are negative.
    lmvo <- lmvo[lmvo$nd <= (lmvo$nd[1] + floor.error) & lmvo$nd >= (lmvo$nd[1] - floor.error),];
    lmvo <- na.omit(lmvo);
    if(any(lmvo$me==1)) {
      return(list(prediction=NA));
    }
    lmvo$n <- log10(lmvo$n); lmvo$me <- log10(lmvo$me-1)
    # remove non-standard steps
    emvf <- lm(me ~ n,data=lmvo);
    tfe <- (log10(target.fold.enrichment-1)-coef(emvf)[[1]])/coef(emvf)[[2]];
    tfen <- 10^tfe;
    return(list(prediction=tfen,log10.fit=emvf));
  })
  
  if(return.chains) {
    return(list(interpolation=intpn,chains=msers$chains))
  } else {
    return(intpn);
  }
  
  return(msers);
 
}


# output binding detection results to a text file
# the file will contain a table with each row corresponding
# to a detected position, with the following columns:
# chr - chromosome or target sequence
# pos - position of detected binding site on the chromosome/sequence
# score - a score reflecting magnitude of the binding
# Evalue - E-value corresponding to the peak magnitude
# FDR - FDR corresponding to the peak magnitude
# enrichment.lb - lower bound of the fold-enrichment ratio
# enrichment.mle - maximum likelihood estimate of the fold-enrichment ratio
output.binding.results <- function(results,filename) {
  write(file=filename,"chr\tpos\tscore\tEvalue\tFDR\tenrichment.lb\tenrichment.mle",append=F);
  chrl <- names(results$npl); names(chrl) <- chrl;
  x <- lapply(chrl,function(chr) {
    d <- results$npl[[chr]];
    if(dim(d)[1]>0) {
      if(results$thr$type=="topN") {
        od <- cbind(rep(chr,dim(d)[1]),subset(d,select=c(x,y,enr,enr.mle)))
      } else {
        od <- cbind(rep(chr,dim(d)[1]),subset(d,select=c(x,y,evalue,fdr,enr,enr.mle)))
      }
      write.table(od,file=filename,col.names=F,row.names=F,sep="\t",append=T,quote=F)
    }
  })
}


# -------- LOW-LEVEL ROUTINES  ------------

# calculates tag strand cross-correlation for a range of shifts (on positive strand)
tag.scc <- function(tags,srange=c(50,250),bin=1,tt=NULL,llim=10) {
  if(is.null(tt)) {
    tt <- table(sign(tags)*as.integer(floor(abs(tags)/bin+0.5)));
  }
  if(!is.null(llim)) { l <- mean(tt); tt <- tt[tt<llim*l] }
  tc <- as.integer(names(tt));
  tt <- as.numeric(tt);

  pv <- tt; pv[tc<0]<-0;
  nv <- tt; nv[tc>0]<-0;

  pti <- which(tc>0)
  nti <- which(tc<0);

  ptc <- tc[pti];
  ntc <- (-1)*tc[nti];

  ptv <- tt[pti];
  ntv <- tt[nti];

  trng <- range(c(range(ptc),range(ntc)))
  l <- diff(trng)+1;
  rm(tc,tt);

  mp <- sum(ptv)*bin/l;   mn <- sum(ntv)*bin/l;
  ptv <- ptv-mp; ntv <- ntv-mn;
  ss <- sqrt((sum(ptv*ptv)+(l-length(ptv))*mp^2) * (sum(ntv*ntv)+(l-length(ntv))*mn^2));

  t.cor <- function(s) {
    smi <- match(ptc+s,ntc);
    return((sum(ptv[!is.na(smi)]*ntv[na.omit(smi)]) -
           mn*sum(ptv[is.na(smi)]) -
           mp*sum(ntv[-na.omit(smi)]) +
           mp*mn*(l-length(ptv)-length(ntv)+length(which(!is.na(smi)))))/ss);
  }
  shifts <- floor(seq(srange[1],srange[2],by=bin)/bin+0.5);
  scc <- unlist(lapply(shifts,t.cor)); names(scc) <- shifts*bin;
  return(scc);
}


# plot tag cross-correlation
t.plotcc <- function(ac, lab=c(10,5,7), ylab="correlation", xlab="lag", pch=19, grid.i=c(-5:5), grid.s=10, type='b', plot.grid=F, cols=c(1,2,4,"orange",8,"pink"), min.peak.x=NULL, xlim=NULL, plot.147=F, plot.max=T, rmw=1, rescale=F, legendx="right", ltys=rep(1,length(ac)), ...) {
    if(is.list(ac)) {
      cols <- cols[1:length(ac)];

      if(!is.null(xlim)) {
        vx <- as.numeric(names(ac[[1]])); vx <- which(vx>=xlim[1] & vx<=xlim[2]);
        ac[[1]] <- (ac[[1]])[vx];
      } else {
        xlim <- range(as.numeric(names(ac[[1]])));
      }


      plot(as.numeric(names(ac[[1]])),runmean(ac[[1]],rmw),type=type,pch=pch,xlab=xlab,ylab=ylab,lab=lab, col=cols[1], xlim=xlim, lty=ltys[1], ...);
      if(length(ac)>1) {
        for(i in seq(2,length(ac))) {
          irng <- range(ac[[i]]);
          vx <- as.numeric(names(ac[[i]])); vx <- which(vx>=xlim[1] & vx<=xlim[2]);
          if(rescale) {
            lines(as.numeric(names(ac[[i]])[vx]),runmean((ac[[i]][vx]-irng[1])/diff(irng)*diff(range(ac[[1]]))+min(ac[[1]]),rmw),col=cols[i],lty=ltys[i]);
          } else {
            lines(as.numeric(names(ac[[i]]))[vx],runmean(ac[[i]][vx],rmw),col=cols[i],lty=ltys[i]);
          }
        }
      }
      if(is.null(min.peak.x)) {
        m <- as.numeric(names(ac[[1]])[which.max(ac[[1]])]);
      } else {
        sac <- (ac[[1]])[which(as.numeric(names(ac[[1]]))>min.peak.x)]
        m <- as.numeric(names(sac)[which.max(sac)]);
      }
      legend(x="topright",bty="n",legend=c(names(ac)),col=cols,lty=ltys)
    } else {
      if(!is.null(xlim)) {
        vx <- as.numeric(names(ac));
        vx <- which(vx>=xlim[1] & vx<=xlim[2]);
        ac <- ac[vx];
      } else {
        xlim <- range(as.numeric(names(ac)));
      }
      
      plot(names(ac),runmean(ac,rmw),type=type,pch=pch,xlab=xlab,ylab=ylab,lab=lab, xlim=xlim, ...);
      if(is.null(min.peak.x)) {
        m <- as.numeric(names(ac)[which.max(ac)]);
      } else {
        sac <- ac[which(names(ac)>min.peak.x)]
        m <- as.numeric(names(sac)[which.max(sac)]);
      }
    }
    if(plot.147) {
      abline(v=147,lty=2,col=8);
    }
    if(plot.grid) {
      abline(v=m+grid.i*grid.s,lty=3,col="pink");
    }
    if(plot.max) {
      abline(v=m,lty=2,col=2);
      legend(x=legendx,bty="n",legend=c(paste("max at ",m,"bp",sep="")));
      return(m);
    }
  }
  
  # plot chromosome-acerage cross-correlation 
  t.plotavcc <- function(ci, main=paste(ci,"chromosome average"), ccl=tl.cc, return.ac=F, ttl=tl, plot=T, ... ) {
    cc <- ccl[[ci]];
    if(length(cc)==1)  { return(cc[[1]]) };
    if(length(cc)==0) { return(c()) };
    ac <- do.call(rbind,cc);
    # omit NA chromosomes
    ina <- apply(ac,1,function(d) any(is.na(d)));

    tags <- ttl[[ci]]; 
    avw <- unlist(lapply(tags,length));    avw <- avw/sum(avw);    
    ac <- ac[!ina,]; avw <- avw[!ina];
    ac <- apply(ac,2,function(x) sum(x*avw));
    if(plot) {
      m <- t.plotcc(ac, main=main, ...);
      if(!return.ac) { return(m) }
    }
    if(return.ac) { return(ac) }
  }

  t.plotchrcc <- function(ci,ncol=4, ccl=tl.cc, ... ) {
    cc <- ccl[[ci]];
    ac <- do.call(rbind,cc);
    par(mfrow = c(length(cc)/ncol,ncol), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
    lapply(names(cc),function(ch) { t.plotcc(cc[[ch]],main=paste(ci,": chr",ch,sep=""), ...) })
  }

  t.plotavccl <- function(ci, ccl=tl.ccl, main=paste(ci,"chromosome average"), rtl=tl, ... ) {
    #cc <- lapply(ccl[[ci]],function(x) { if(!is.null(x$M)) { x$M <- NULL;}; return(x); });
    cc <- ccl[[ci]];
    chrs <- names(cc[[1]]); names(chrs) <- chrs;
    acl <- lapply(cc,function(x) do.call(rbind,x));
    tags <- rtl[[ci]][chrs]; 
    avw <- unlist(lapply(tags,length));    avw <- avw/sum(avw);
    acl <- lapply(acl,function(ac) apply(ac,2,function(x) sum(x*avw)))
    t.plotcc(acl, main=main, ...);
  }
  
  t.plotchrccl <- function(ci,ccl=tl.ccl,ncol=4, ... ) {
    par(mfrow = c(length(cc[[1]])/ncol,ncol), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
    lapply(names(cc[[1]]),function(ch) { t.plotcc(lapply(cc,function(x) x[[ch]]),main=paste(ci,": chr",ch,sep=""), ...) })
  }

  

show.scc <- function(tl,srange,cluster=NULL) {
  if(!is.null(cluster)) {
    cc <- clusterApplyLB(cluster,tl,tag.scc,srange=srange);
    names(cc) <- names(tl); 
  } else {
    cc <- lapply(tl,tag.scc,srange=srange);
  }
  par(mfrow = c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
  ccl<-list(sample=cc);
  ccl.av <- lapply(names(ccl),t.plotavcc,type='l',ccl=ccl,xlim=srange,return.ac=F,ttl=list(sample=tl),main="")[[1]]
}

# find regions of significant tag enrichment
find.significantly.enriched.regions <- function(signal.data,control.data,window.size=500,multiplier=1,z.thr=3,mcs=0,debug=F,background.density.scaling=T,masking.window.size=window.size,poisson.z=0,poisson.ratio=4,either=F,tag.shift=146/2,bg.weight=NULL) {
  if(is.null(bg.weight)) {
    bg.weight <- dataset.density.ratio(signal.data,control.data,background.density.scaling=background.density.scaling);
  }

  if(debug) {
    cat("bg.weight=",bg.weight,"\n");
  }
  chrl <- names(signal.data); names(chrl) <- chrl; 
  tec <- lapply(chrl,function(chr) {
    d <- tag.enrichment.clusters(signal.data[[chr]],control.data[[chr]],bg.weight=bg.weight*multiplier,thr=z.thr,wsize=window.size,mcs=mcs,min.tag.count.z=poisson.z,min.tag.count.ratio=poisson.ratio,either=either,tag.shift=tag.shift);
    d$s <- d$s-masking.window.size/2; d$e <- d$e+masking.window.size/2;
    return(d);
  })
}


# given tag position vectors, find contigs of significant enrichment of signal over background
# thr - z score threshold
# mcs - minimal cluster size
# bg.weight - fraction by which background counts should be multipled
# min.tag.count.z will impose a poisson constraint based on randomized signal in parallel of background constaint (0 - no constraint)
tag.enrichment.clusters <- function(signal,background,wsize=200,thr=3,mcs=1,bg.weight=1,min.tag.count.z=0,tag.av.den=NULL,min.tag.count.thr=0,min.tag.count.ratio=4,either=F,tag.shift=146/2) {
  if(is.null(tag.av.den)) {
    tag.av.den <- length(signal)/diff(range(abs(signal)));
  }
  if(min.tag.count.z>0) {
    min.tag.count.thr <- qpois(pnorm(min.tag.count.z,lower.tail=F),min.tag.count.ratio*tag.av.den*wsize,lower.tail=F)
  } else {
    min.tag.count.thr <- 0;
  }
  
  #if(bg.weight!=1) {
  #  background <- sample(background,length(background)*(bg.weight),replace=T);
  #}
  # make up combined position, flag vectors
  pv <- abs(c(signal,background)+tag.shift);
  fv <- c(rep(1,length(signal)),rep(0,length(background)));
  po <- order(pv);
  pv <- pv[po];
  fv <- fv[po];

  #thr <- pnorm(thr,lower.tail=F);
  
  storage.mode(wsize) <- storage.mode(mcs) <- storage.mode(fv) <- "integer";
  storage.mode(thr) <- storage.mode(pv) <- "double";
  storage.mode(bg.weight) <- "double";
  storage.mode(min.tag.count.thr) <- "double";
  either <- as.integer(either);
  storage.mode(either) <- "integer";
  
  z <- .Call("find_poisson_enrichment_clusters",pv,fv,wsize,thr,mcs,bg.weight,min.tag.count.thr,either)
  return(z);
}





# estimates threshold, calculates predictions on complete data and randomized data
# input: tvl
# control - a list of control tag datasets
# no randomization is done if control is supplied
# return.rtp - return randomized tag peaks - do not fit thresholds or do actual predictions
# topN - use min threshold to do a run, return topN peaks from entire genome
# threshold - specify a user-defined threshold
lwcc.prediction <- function(tvl,e.value=NULL, fdr=0.01, chrl=names(tvl), min.thr=0, n.randomizations=1, shuffle.window=1, debug=T, predict.on.random=F, shuffle.both.strands=T,strand.shuffle.only=F, return.rtp=F, control=NULL, print.level=0, threshold=NULL, topN=NULL, bg.tl=NULL, tec.filter=T, tec.window.size=1e3,tec.z=3, tec.masking.window.size=tec.window.size, tec.poisson.z=3,tec.poisson.ratio=4, bg.reverse=T, return.control.predictions=F, return.core.data=F, background.density.scaling=T, ... ) {

  control.predictions <- NULL;
  core.data <- list();

  if(!is.null(bg.tl) & tec.filter) {
    if(debug) { cat("finding background exclusion regions ... "); }
    tec <- find.significantly.enriched.regions(bg.tl,tvl,window.size=tec.window.size,z.thr=tec.z,masking.window.size=tec.masking.window.size,poisson.z=tec.poisson.z,poisson.ratio=tec.poisson.ratio,background.density.scaling=background.density.scaling,either=T);
    if(return.core.data) {
      core.data <- c(core.data,list(tec=tec));
    }
    if(debug) { cat("done\n"); }
  }

  
  if(is.null(threshold) & is.null(topN)) { # threshold determination is needed
    # generate control predictions
    if(!is.null(control)) {
      if(debug) { cat("determining peaks on provided",length(control),"control datasets:\n");   }
      if(!is.null(bg.tl)) {
        if(bg.reverse) {
          if(debug) { cat("using reversed signal for FDR calculations\n"); }
          rbg.tl <- tvl;
        } else {
          if(debug) { cat("generating randomized (within chromosome) background ... "); }
          rbg.tl <- lapply(bg.tl,function(d) {
            if(length(d)<2) { return(d); }
            rng <- range(abs(d));
            rd <- round(runif(length(d),rng[1],rng[2]));
            nrd <- sample(1:length(rd),length(which(d<0)));
            rd[nrd] <- rd[nrd]*(-1);
            return(rd);
          })
          if(debug) { cat("done\n"); }
        }
      } else {
        rbg.tl <- NULL;
      }
      n.randomizations <- length(control);
      #signal.size <- sum(unlist(lapply(tvl,length)));
      rtp <- lapply(control,function(d) {
        # calculate tag.weight
        #tag.weight <- sum(unlist(lapply(tvl,length)))/sum(unlist(lapply(d,length)));
        tag.weight <- dataset.density.ratio(tvl,d,background.density.scaling=background.density.scaling);
        #cat("tag.weight=",tag.weight," ");
        return(window.call.mirror.binding(d,min.thr=min.thr, tag.weight=tag.weight,bg.tl=rbg.tl, debug=debug, round.up=T,background.density.scaling=background.density.scaling, ...));
        #return(window.call.mirror.binding(d,min.thr=min.thr, method=tag.wtd,wsize=200,bg.tl=control.data,window.size=window.size,debug=T,min.dist=min.dist,cluster=cluster))
      });
      if(return.core.data) {
        core.data <- c(core.data,list(rtp.unfiltered=rtp));
      }
      if(tec.filter) {
        if(debug) { cat("excluding systematic background anomalies ... "); }
        rtp <- lapply(rtp,filter.binding.sites,tec,exclude=T);
        if(debug) { cat("done\n"); }
      }
    } else {
      if(debug) { cat("determining peaks on ",n.randomizations,"randomized datasets:\n");   }
      rtp <- lapply(1:n.randomizations,function(i) {
        rd <- generate.randomized.data(tvl,shuffle.window=shuffle.window,shuffle.both.strands=shuffle.both.strands,strand.shuffle.only=strand.shuffle.only);
        return(window.call.mirror.binding(rd,min.thr=min.thr,bg.tl=bg.tl, debug=debug, ...));
        #return(window.call.mirror.binding(rd,min.thr=min.thr, method=tag.wtd,wsize=200,bg.tl=control.data,window.size=window.size,debug=T,min.dist=min.dist))
      });
    }
    if(return.control.predictions) {
      control.predictions <- rtp;
    } 
    rtp <- do.call(rbind,lapply(rtp,function(d) do.call(rbind,d))); # merge tables
    
    # generate real data predictions
    if(debug) { cat("determining peaks on real data:\n");   }
    npl <- window.call.mirror.binding(tvl,min.thr=min.thr,bg.tl=bg.tl, debug=debug, background.density.scaling=background.density.scaling, ...);
    #npl <- window.call.mirror.binding(tvl,min.thr=min.thr, method=tag.wtd,wsize=200,bg.tl=control.data,window.size=window.size,debug=T,min.dist=min.dist,cluster=cluster);
    if(return.core.data) {
      core.data <- c(core.data,list(npl.unfiltered=npl));
    }

    if(!is.null(bg.tl) & tec.filter) {
      if(debug) { cat("excluding systematic background anomalies ... "); }
      npl <- filter.binding.sites(npl,tec,exclude=T);
      if(debug) { cat("done\n"); }
    }

    # calculate E-value and FDRs for all of the peaks
    if(debug) { cat("calculating statistical thresholds\n"); }
    chrl <- names(npl); names(chrl) <- chrl;
    npld <- do.call(rbind,lapply(names(npl),function(chr) { k <- npl[[chr]]; if(!is.null(k) & dim(k)[1]>0) { k$chr <- rep(chr,dim(k)[1]) }; return(k) }))
    npld <- cbind(npld,get.eval.fdr.vectors(npld$y,rtp$y));
    # correct for n.randomizations
    npld$fdr <- npld$fdr/n.randomizations;
    npld$evalue <- npld$evalue/n.randomizations;

    if(return.core.data) {
      core.data <- c(core.data,list(npld=npld));
    }

    # determine actual thresholds
    if(is.null(e.value)) {
      if(is.null(fdr)) { fdr <- 0.01; }
      thr <- list(root=min(npld$y[npld$fdr<=fdr]),type="FDR",fdr=fdr)
      if(debug) { cat("FDR",fdr,"threshold=",thr$root,"\n");  }
    } else {
      # determine threshold based on e-value
      thr <- list(root=min(npld$y[npld$evalue<=e.value]),type="Evalue",e.value=e.value)
      if(debug) { cat("E-value",e.value,"threshold=",thr$root,"\n");  }
    }


    npld <- npld[npld$y>=thr$root,];
    if(dim(npld)[1]>0) {
      npl <- tapply(c(1:dim(npld)[1]),as.factor(npld$chr),function(ii) {df <- npld[ii,]; df$chr <- NULL; return(df) });
    } else {
      npl <- list();
    }
  } else {
    if(is.null(threshold)) {
      thr <- list(root=min.thr,type="minimal");
    } else {
      thr <- list(root=threshold,type="user specified");
    }

    cat("calling binding positions using",thr$type,"threshold (",thr$root,") :\n");
    npl <- window.call.mirror.binding(tvl=tvl,min.thr=thr$root,bg.tl=bg.tl, debug=debug, ...);
    if(!is.null(bg.tl) & tec.filter) {
      if(debug) { cat("excluding systematic background anomalies ... "); }
      npl <- filter.binding.sites(npl,tec,exclude=T);
      if(debug) { cat("done\n"); }
    }

    if(!is.null(topN)) {
      # determine threshold based on topN peaks
      ay <- unlist(lapply(npl,function(d) d$y));
      if(length(ay)>topN) {
        thr <- list(root=sort(ay,decreasing=T)[topN],type="topN",topN=topN);
        cat(paste("determined topN threshold :",thr$root,"\n"));      
        npl <- lapply(npl,function(d) d[d$y>thr$root,]);
      }
    }
  }

  if(return.core.data) {
    return(c(list(npl=npl,thr=thr),core.data));
  }
  if(return.control.predictions & !is.null(control.predictions)) {
    return(list(npl=npl,thr=thr,control.predictions=control.predictions));
  }
  return(list(npl=npl,thr=thr));
}

# window tag difference method
wtd <- function(x,y,s,e,whs=200,return.peaks=T,min.thr=5,min.dist=200,step=1,direct.count=F,tag.weight=1,bg.x=NULL,bg.y=NULL,bg.weight=1,mask.x=NULL,mask.y=NULL,ignore.masking=F, bg.whs=whs, round.up=F, ...) {
  ignore.masking <- ignore.masking | (is.null(mask.x) & is.null(mask.y));
  if(step>1) {
    x <- floor(x/step+0.5); y <- floor(y/step+0.5)
    
    if(!is.null(bg.x)) {
      bg.x <- floor(bg.x/step+0.5); bg.y <- floor(bg.y/step+0.5)  
    }
    
    if(!is.null(mask.x)) {
      mask.x <- floor(mask.x/step+0.5); mask.y <- floor(mask.y/step+0.5)  
    }

    
    whs <- floor(whs/step+0.5);
    bg.whs <- floor(bg.whs/step+0.5);
    min.dist <- floor(min.dist/step +0.5);
    s <- floor(s/step+0.5)
    e <- floor(e/step+0.5)
  }

  # scale bg.weight, since within calculation they are considered independent
  bg.weight <- bg.weight*tag.weight;

  rx <- c(s-whs,e+whs);

  # compile tag vectors
  xt <- table(x);
  xh <- integer(diff(rx)+1);
  xh[as.integer(names(xt))-rx[1]+1] <- as.integer(xt);

  yt <- table(y);
  yh <- integer(diff(rx)+1);
  yh[as.integer(names(yt))-rx[1]+1] <- as.integer(yt);

  # compile background vectors
  if(!is.null(bg.x) & length(bg.x)>0) {
    bg.subtract <- 1;

    bg.xt <- table(bg.x);
    bg.xh <- integer(diff(rx)+1);
    bg.xh[as.integer(names(bg.xt))-rx[1]+1] <- as.integer(bg.xt);
    rm(bg.xt);

    bg.yt <- table(bg.y);
    bg.yh <- integer(diff(rx)+1);
    bg.yh[as.integer(names(bg.yt))-rx[1]+1] <- as.integer(bg.yt);
    rm(bg.yt);

    # adjust bg.weight according to bg.whs
    if(bg.whs!=whs) {
      bg.weight <- bg.weight*whs/bg.whs;
    }
  } else {
    bg.subtract <- 0;
    bg.xh <- bg.yh <- c();
  }

  # record masked positions
  if(!ignore.masking) {
    if(!is.null(mask.x) & length(mask.x)>0) {
      mvx <- unique(mask.x); mvx <- setdiff(mvx,as.numeric(names(xt)));
      mvx <- mvx[mvx>=rx[1] & mvx<=rx[2]];
      xh[mvx-rx[1]+1] <- -1;
    }

    if(!is.null(mask.y) & length(mask.y)>0) {
      mvy <- unique(mask.y); mvy <- setdiff(mvy,as.numeric(names(yt)));
      mvy <- mvy[mvy>=rx[1] & mvy<=rx[2]];
      yh[mvy-rx[1]+1] <- -1;
    }
  }

  rm(xt,yt);

  if(round.up) { round.up <- 1; } else { round.up <- 0; }
  
  storage.mode(xh) <- storage.mode(yh) <- "integer";
  storage.mode(bg.xh) <- storage.mode(bg.yh) <- "integer";
  nx <- length(xh);   storage.mode(nx) <- storage.mode(whs) <- storage.mode(bg.whs) <- "integer";
  rp <- as.integer(return.peaks);
  dcon <- as.integer(direct.count);
  storage.mode(rp) <- storage.mode(min.dist) <- "integer";
  storage.mode(min.thr) <- "double";
  storage.mode(dcon) <- "integer";
  storage.mode(tag.weight) <- "double";
  storage.mode(bg.weight) <- "double";
  storage.mode(bg.subtract) <- "integer";
  storage.mode(round.up) <- "integer";
  im <- as.integer(ignore.masking);
  storage.mode(im) <- "integer";
  z <- .Call("wtd",xh,yh,whs,rp,min.dist,min.thr,dcon,tag.weight,im,bg.subtract,bg.xh,bg.yh,bg.whs,bg.weight,round.up);
  if(return.peaks) {
    return(data.frame(x=(z$x+rx[1])*step,y=z$v));
  } else {
    return(list(x=rx*step,y=z));
  }
}


tag.wtd <- function(ctv,s,e,return.peaks=T, bg.ctv=NULL,  mask.ctv=NULL, ...) {
  x <- ctv[ctv>=s & ctv<=e];
  y <- (-1)*ctv[ctv<=-s & ctv>=-e];

  if(!is.null(bg.ctv)) {
    bg.x <- bg.ctv[bg.ctv>=s & bg.ctv<=e];
    bg.y <- (-1)*bg.ctv[bg.ctv<=-s & bg.ctv>=-e];
  } else {
    bg.x <- bg.y <- NULL;
  }

  if(!is.null(mask.ctv)) {
    mask.x <- mask.ctv[mask.ctv>=s & mask.ctv<=e];
    mask.y <- (-1)*mask.ctv[mask.ctv<=-s & mask.ctv>=-e];
  } else {
    mask.x <- mask.y <- NULL;
  }

  if(length(x)==0 | length(y) ==0) {
    if(return.peaks) {
      return(data.frame(x=c(),y=c()));
    } else {
      rx <- range(c(x,y));
      return(list(x=rx,y=numeric(diff(rx)+1)));
    }
  } else {
    return(wtd(x,y,s,e,return.peaks=return.peaks,  bg.x=bg.x,bg.y=bg.y, mask.x=mask.x,mask.y=mask.y, ...))
  }
}

# shuffles tags in chromosome blocks of a specified size
# note: all coordinates should be positive
tag.block.shuffle <- function(tags,window.size=100) {
  if(length(tags)<3) {
    warning("too few tags for shuffling");
    return(tags);
  }
  rng <- range(tags);
  #if(rng[1]<0) { stop("negative tag coordinates found") }
  if(diff(rng)<=window.size) {
    warning(paste("tag range (",diff(rng),") is smaller than shuffle window size"));
    return(tags);
  }

  if(window.size==0) {
    return(as.integer(runif(length(tags),min=rng[1],max=rng[2])))
  } else if(window.size==1) {
    tt <- table(tags);
    return(rep(runif(length(tt),min=rng[1],max=rng[2]),as.integer(tt)))
  } else {
  # block positions
    bp <- tags %/% window.size;
  # block-relative tag positions
    rp <- tags %% window.size;

  # shuffle block positions
    bpu <- unique(bp);
    rbp <- range(bpu);
    bps <- as.integer(runif(length(bpu),min=rbp[1],max=rbp[2]));
    bpi <- match(bp,bpu);
    sbp <- bps[bpi];
    #sbp <- rbp[1]+match(bp,sample(rbp[1]:rbp[2]))
    return(sbp*window.size+rp);
  }
}


# calculate window cross-correlation
lwcc <- function(x,y,s,e,whs=100,isize=20,return.peaks=T,min.thr=1,min.dist=100,step=1,tag.weight=1,bg.x=NULL,bg.y=NULL,bg.weight=NULL,mask.x=NULL,mask.y=NULL,bg.whs=whs,round.up=F) {
  if(step>1) {
    x <- floor(x/step+0.5); y <- floor(y/step+0.5)
    
    if(!is.null(bg.x)) {
      bg.x <- floor(bg.x/step+0.5); bg.y <- floor(bg.y/step+0.5)  
    }
    
    if(!is.null(mask.x)) {
      mask.x <- floor(mask.x/step+0.5); mask.y <- floor(mask.y/step+0.5)  
    }

    whs <- floor(whs/step+0.5);
    bg.whs <- floor(bg.whs/step+0.5);
    isize <- floor(isize/step+0.5);
    min.dist <- floor(min.dist/step +0.5);
    s <- floor(s/step+0.5)
    e <- floor(e/step+0.5)
  }

  # scale bg.weight, since within calculation they are considered independent
  bg.weight <- bg.weight*tag.weight;

  
  rx <- c(s-whs,e+whs);
  xt <- table(x);
  xh <- integer(diff(rx)+1);
  xh[as.integer(names(xt))-rx[1]+1] <- as.integer(xt);

  yt <- table(y);
  
  yh <- integer(diff(rx)+1);
  yh[as.integer(names(yt))-rx[1]+1] <- as.integer(yt);

  # compile background vectors
  if(!is.null(bg.x) & length(bg.x)>0) {
    bg.subtract <- 1;

    bg.xt <- table(bg.x);
    bg.xh <- integer(diff(rx)+1);
    bg.xh[as.integer(names(bg.xt))-rx[1]+1] <- as.integer(bg.xt);
    rm(bg.xt);

    bg.yt <- table(bg.y);
    bg.yh <- integer(diff(rx)+1);
    bg.yh[as.integer(names(bg.yt))-rx[1]+1] <- as.integer(bg.yt);
    rm(bg.yt);

    # adjust bg.weight according to bg.whs
    bg.weight <- bg.weight*(whs-isize)/bg.whs;
  } else {
    bg.subtract <- 0;
    bg.xh <- bg.yh <- c();
  }

  # record masked positions
  if(!is.null(mask.x) & length(mask.x)>0) {
    mvx <- unique(mask.x); mvx <- setdiff(mvx,as.numeric(names(xt)));
    mvx <- mvx[mvx>=rx[1] & mvx<=rx[2]];
    
    xh[mvx-rx[1]+1] <- -1;
  }

  if(!is.null(mask.y) & length(mask.y)>0) {
    mvy <- unique(mask.y); mvy <- setdiff(mvy,as.numeric(names(yt)));
    mvy <- mvy[mvy>=rx[1] & mvy<=rx[2]];
    yh[mvy-rx[1]+1] <- -1;
  } 
  
  rm(xt,yt);
  if(round.up) { round.up <- 1; } else { round.up <- 0; }
  
  storage.mode(xh) <- storage.mode(yh) <- "integer";
  storage.mode(bg.xh) <- storage.mode(bg.yh) <- "integer";
  nx <- length(xh);   storage.mode(nx) <- storage.mode(whs) <- storage.mode(isize) <- storage.mode(bg.whs) <- "integer";
  rp <- as.integer(return.peaks);
  storage.mode(rp) <- storage.mode(min.dist) <- "integer";
  storage.mode(min.thr) <- "double";
  storage.mode(tag.weight) <- "double";
  storage.mode(bg.weight) <- "double";
  storage.mode(bg.subtract) <- "integer";
  storage.mode(round.up) <- "integer";

  # allocate return arrays
  #cc <- numeric(nx); storage.mode(cc) <- "double";
  z <- .Call("lwcc",xh,yh,whs,isize,rp,min.dist,min.thr,tag.weight,bg.subtract,bg.xh,bg.yh,bg.whs,bg.weight,round.up);
  if(return.peaks) {
    return(data.frame(x=(z$x+rx[1])*step,y=z$v));
  } else {
    return(list(x=rx*step,y=z));
  }
}


tag.lwcc <- function(ctv,s,e,return.peaks=T, bg.ctv=NULL, mask.ctv=NULL, ...) {
  x <- ctv[ctv>=s & ctv<=e];
  y <- (-1)*ctv[ctv<=-s & ctv>=-e];

  if(!is.null(bg.ctv)) {
    bg.x <- bg.ctv[bg.ctv>=s & bg.ctv<=e];
    bg.y <- (-1)*bg.ctv[bg.ctv<=-s & bg.ctv>=-e];
  } else {
    bg.x <- bg.y <- NULL;
  }

  if(!is.null(mask.ctv)) {
    mask.x <- mask.ctv[mask.ctv>=s & mask.ctv<=e];
    mask.y <- (-1)*mask.ctv[mask.ctv<=-s & mask.ctv>=-e];
  } else {
    mask.x <- mask.y <- NULL;
  }
  
  if(length(x)==0 | length(y) ==0) {
    if(return.peaks) {
      return(data.frame(x=c(),y=c()));
    } else {
      rx <- range(c(x,y));
      return(list(x=rx,y=numeric(diff(rx)+1)));
    }
  } else { 
    return(lwcc(x,y, s,e,return.peaks=return.peaks, bg.x=bg.x,bg.y=bg.y,  mask.x=mask.x,mask.y=mask.y, ...))
  }
}

# determine mirror-based binding positions using sliding window along each chromosome
# extra parameters are passed on to call.nucleosomes()
window.call.mirror.binding <- function(tvl,window.size=4e7, debug=T, cluster=NULL, bg.tl=NULL, mask.tl=NULL, background.density.scaling=T, ...) {
  chrl <- names(tvl);
  # determine bg.weight
  if(!is.null(bg.tl)) {
    bg.weight <- dataset.density.ratio(tvl,bg.tl,background.density.scaling=background.density.scaling);
  } else {
    bg.weight <- NULL;
  }
  if(debug) {
    cat("bg.weight=",bg.weight," ");
  }
  
  names(chrl) <- chrl;

  if(is.null(cluster)) {
    return(lapply(chrl,function(chr) {
      bg.ctv <- NULL; if(!is.null(bg.tl)) { bg.ctv <- bg.tl[[chr]]; };
      mask.ctv <- NULL; if(!is.null(mask.tl)) { mask.ctv <- mask.tl[[chr]]; };
      
      window.chr.call.mirror.binding(list(ctv=tvl[[chr]],bg.ctv=bg.ctv,mask.ctv=mask.ctv),window.size=window.size,chr=chr,debug=debug, bg.weight=bg.weight, bg.ctv=bg.ctv, mask.ctv=mask.ctv, ...);
    }));
  } else {
    # add bg.ctv and mask.ctv to parallel call
    tvll <- lapply(chrl,function(chr) {
      bg.ctv <- NULL; if(!is.null(bg.tl)) { bg.ctv <- bg.tl[[chr]]; };
      mask.ctv <- NULL; if(!is.null(mask.tl)) { mask.ctv <- mask.tl[[chr]]; };
      return(list(ctv=tvl[[chr]],bg.ctv=bg.ctv,mask.ctv=mask.ctv))
    });
    bl <- clusterApplyLB(cluster,tvll,window.chr.call.mirror.binding,window.size=window.size,debug=debug, bg.weight=bg.weight, ...);
    names(bl) <- chrl;
    return(bl);
  }
}

window.chr.call.mirror.binding <- function(ctvl,window.size,debug=T, chr="NA", cluster=NULL, method=tag.wtd, bg.ctv=NULL, mask.ctv=NULL, ...) {
  ctv <- ctvl$ctv; bg.ctv <- ctvl$bg.ctv; mask.ctv <- ctvl$mask.ctv;
  if(is.null(ctv)) { return(data.frame(x=c(),y=c())) }
  if(length(ctv)<2) { return(data.frame(x=c(),y=c())) }
  
  dr <- range(unlist(lapply(ctv,function(x) range(abs(x)))))
  n.windows <- ceiling(diff(dr)/window.size);
  
  
  pinfo <- c();
  if(debug) {
    cat(paste("processing ",chr," in ",n.windows," steps [",sep=""));
  }
  for(i in 1:n.windows) {
    s <- dr[1]+(i-1)*window.size;
    npn <- method(s=s, e=s+window.size,ctv=ctv, return.peaks=T, bg.ctv=bg.ctv, mask.ctv=mask.ctv, ... );
    if(length(npn) > 0) { pinfo <- rbind(pinfo,npn)  }
    if(debug) {
      cat(".");
    }
  }
  if(debug) {
    cat(paste("] done (",dim(pinfo)[1],"positions)\n"));
  } else {
    cat(".");
  }
  return(data.frame(x=pinfo[,1],y=pinfo[,2]));
}

generate.randomized.data <- function(data,shuffle.window=1,shuffle.both.strands=T,strand.shuffle.only=F,chrl=names(data)) {
  names(chrl) <- unlist(chrl);
  if(strand.shuffle.only) {
    # shuffle just strand assignment, not tag positions
    rt <- lapply(data[unlist(chrl)],function(tv) tv*sample(c(-1,1),length(tv),replace=T));
  } else {
    if(shuffle.both.strands) {
      rt <- lapply(data[unlist(chrl)],function(tv) {
        pti <- which(tv>0); return(c(tag.block.shuffle(tv[pti],window.size=shuffle.window),tag.block.shuffle(tv[-pti],window.size=shuffle.window)))
      });
    } else {
      rt <- lapply(data[unlist(chrl)],function(tv) { pti <- which(tv>0); return(c(tag.block.shuffle(tv[pti],window.size=shuffle.window),tv[-pti]))});
    }
  }
}

# determine threshold based on E value
# for efficiency chrl should include just one or two small chromosomes
# optional parameters are passed to call.nucleosomes()
determine.lwcc.threshold <- function(tvl,chrl=names(tvl),e.value=100, n.randomizations=1, min.thr=1, debug=F, tol=1e-2, shuffle.window=1, shuffle.both.strands=T, return.rtp=F, control=NULL, strand.shuffle=F, ...) {
  names(chrl) <- unlist(chrl);
  
  # determine fraction of total tags contained in the specified nucleosomes
  ntags <- sum(unlist(lapply(tvl,function(cv)  length(cv))));
  nctags <- sum(unlist(lapply(chrl, function(cn) length(tvl[[cn]]))));
  # calculate actual target E value
  if(!is.null(control)) {
    n.randomizations <- length(control);
  }
  eval <- e.value*n.randomizations*nctags/ntags
  if(eval<1) {
    warning("specified e.value and set of chromosomes results in target e.value of less than 1");
    eval <- 1;
  }
  
  if(debug) {
    cat(paste("randomizations =",n.randomizations," chromosomes =",length(chrl),"\n"))
    cat(paste("adjusted target eval =",eval,"\ngenerating randomized tag peaks ..."));
  }

  # get peaks on randomized tags
  if(is.null(control)) {
    rtp <- data.frame(do.call(rbind,lapply(1:n.randomizations,function(i) {
      if(strand.shuffle) {
        # shuffle just strand assignment, not tag positions
        rt <- lapply(tvl[unlist(chrl)],function(tv) tv*sample(c(-1,1),length(tv),replace=T));
      } else {
        if(shuffle.both.strands) {
          rt <- lapply(tvl[unlist(chrl)],function(tv) {
            pti <- which(tv>0); return(c(tag.block.shuffle(tv[pti],window.size=shuffle.window),tag.block.shuffle(tv[-pti],window.size=shuffle.window)))
          });
        } else {
          rt <- lapply(tvl[unlist(chrl)],function(tv) { pti <- which(tv>0); return(c(tag.block.shuffle(tv[pti],window.size=shuffle.window),tv[-pti]))});
        }
      }
      if(debug) {
        cat(".");
      }
      rl <- window.call.mirror.binding(rt,min.thr=min.thr, debug=F, ...);
      
      return(do.call(rbind,rl))
      #return(do.call(rbind,window.call.mirror.binding(rt,min.thr=min.thr, debug=F, whs=100,isize=10,window.size=3e7,min.dist=200)))
    })));

  } else {
    if(debug) {
      cat(" using provided controls ");
    }
    rtp <- data.frame(do.call(rbind,lapply(control,function(rt) do.call(rbind,window.call.mirror.binding(rt,min.thr=min.thr, debug=F, ...)))))
  }

  if(return.rtp) {
    return(rtp)
  }

  if(debug) {
    cat(" done\nfinding threshold .");
  }

  # determine range and starting value
  rng <- c(min.thr,max(na.omit(rtp$y)))
  
    # find E value threshold
  count.nucs.f <- function(nthr) {
    return(eval-length(which(rtp$y>=nthr)));
  }
  
  # estimate position of the root by downward bisection iterations
  mv <- c(eval); mvp <- c(rng[2]); ni <- 1;
  max.it <- 2*as.integer(log2(rng[2]/rng[1])+0.5);
  while((ni<=max.it) & (mv[1]>=0)) {
    np <- mvp[1]/2;
    npv <- count.nucs.f(np);
    mv <- c(npv,mv);
    mvp <- c(np,mvp);
    ni <- ni+1;
  }
  
  
  if(ni>max.it) {
    # determine lowest value
    if(debug) {
      cat(paste("exceeded max.it (",max.it,"), returning lowest point",signif(mvp[1],4)));
    }
    return(list(root=mvp[1]))
  } else {
    rng <- mvp[1:2];
    if(mv[2]==0) rng[2] <- mvp[3];
    if(debug) {
      cat(paste("bound to (",signif(rng[1],4),signif(rng[2],4),") "));
    }
  }
  
  # find root on the right side
  x <- uniroot(count.nucs.f,rng,tol=tol);
  #x$max <- o$par;
  #x$f.max <- (-1)*o$value;
  if(debug) {
    cat(paste(" done (thr=",signif(x$root,4),")\n"));
  }
  return(x);

}


# determine membership of points in fragments
points.within <- function(x,fs,fe,return.list=F,return.unique=F,sorted=F,return.point.counts=F) {
  if(is.null(x) | length(x) < 1) { return(c()) };
  if(!sorted) {
    ox <- rank(x,ties="first");
    x <- sort(x);
  }

  se <- c(fs,fe);
  fi <- seq(1:length(fs));
  fi <- c(fi,-1*fi);

  fi <- fi[order(se)];
  se <- sort(se);
  
  storage.mode(x) <- storage.mode(fi) <- storage.mode(se) <- "integer";
  if(return.unique) { iu <- 1; } else { iu <- 0; }
  if(return.list) { il <- 1; } else { il <- 0; }
  if(return.point.counts) { rpc <- 1; } else { rpc <- 0; }
  storage.mode(iu) <- storage.mode(il) <- storage.mode(rpc) <- "integer";
  result <- .Call("points_within",x,se,fi,il,iu,rpc);
  if(!sorted & !return.point.counts) {
    result <- result[ox];
  }
  return(result);  
}


# determine cooridnates of points x relative to signed
# positions pos within size range
get.relative.coordinates <- function(x,pos,size,sorted=F) {
  if(!sorted) {
    op <- order(abs(pos));
    x <- sort(x); pos <- pos[op];
  }
  #dyn.load("~/zhao/sc/peaks.so");
  storage.mode(x) <- storage.mode(pos) <- storage.mode(size) <- "integer";
  rf <- .Call("get_relative_coordinates",x,pos,size);
  if(!sorted) { 
    rf$i <- op[rf$i];
  } else {
    return(rf$i);
  }
  return(rf);
}

# given list of magnitude values for signal(x) and control (y),
# return a dataframe with $e.val and $fdr
get.eval.fdr.vectors <- function(x,y) {
  nx <- length(x); ny <- length(y);
  if(nx==0) { return(data.frame(evalue=c(),fdr=c())) }
  if(ny==0) { return(data.frame(evalue=rep(0,nx),fdr=rep(1,nx))) }
  ex <- ecdf(x); ey <- ecdf(y);

  evals <- (1-ey(x))*ny;
  yvals <- (1-ex(x))*nx;
  fdr <- (evals+0.5)/(yvals+0.5); # with pseudo-counts
  fdr[yvals==0] <- min(fdr); # correct for undercounts
  # find a min x corresponding to a minimal FDR
  mfdr <- min(fdr);
  mfdrmx <- min(x[fdr==mfdr]);
  # correct
  fdr[x>=mfdrmx] <- mfdr;
  return(data.frame(evalue=(evals+1),fdr=fdr));
}


# filter predictions to remove calls failling into the tag enrichment clusters ( chr list of $s/$e dfs)
filter.binding.sites <- function(bd,tec,exclude=F) {
  chrl <- names(bd); names(chrl) <- chrl;
  lapply(chrl,function(chr) {
    cbd <- bd[[chr]];
    if(is.null(cbd)) { return(NULL) };
    if(length(cbd)==0) { return(NULL) };
    if(dim(cbd)[1]>0) {
      ctec <- tec[[chr]];
      if(length(ctec$s)>0) {
        if(exclude) {
          pwi <- which(points.within(cbd$x,ctec$s,ctec$e)== -1);
        } else {
          pwi <- which(points.within(cbd$x,ctec$s,ctec$e)> -1);
        }
        return(cbd[pwi,]);
      } else {
        if(exclude) {
          return(cbd);
        } else {
          return(data.frame(x=c(),y=c()));
        }
      }
    } else {
      return(cbd);
    }
  });  
}


# PUBLIC
# generate predictions on sequential (chained) subsamples of data
# if step.size <1, it is intepreted as a fraciton and a  each subsequent subsample
# is of a size (1-fraction.step)*N (N - size of the signal data);
# otherwise the step.size is interpreted as a number of tags, and each subsample is of the size N-step.size
get.subsample.chain.calls <- function(signal.data,control.data,n.steps=NULL,step.size=1e6,subsample.control=F,debug=F,min.ntags=1e3, excluded.steps=c(), test.chromosomes=NULL, ... ) {

  if(!is.null(test.chromosomes)) {
    # adjust step size
    sz <- sum(unlist(lapply(signal.data,length)))
    signal.data <- signal.data[test.chromosomes];
    control.data <- control.data[test.chromosomes];
    
    if(step.size>1) {
      step.size <- step.size*sum(unlist(lapply(signal.data,length)))/sz;
        # cat("adjusted step.size=",step.size,"\n");
    }
  }

  if(is.null(n.steps)) {
    if(step.size<1) {
      # down to 10%
      n.steps <- log(0.1)/log(step.size);
    } else {
      n.steps <- floor(sum(unlist(lapply(signal.data,length)))/step.size)
    }
  }
  if(subsample.control & !is.null(control.data)) {
    # normalize control to the signal size
    if(debug) { cat("pre-subsampling control.\n"); }
    bg.weight <- sum(unlist(lapply(signal.data,length)))/sum(unlist(lapply(control.data,length)))
    control.data <- lapply(control.data,function(d) sample(d,length(d)*bg.weight,replace=(bg.weight>1)))
  }
  calls <- list();
  callnames <- c();
  for(i in 0:n.steps) {
    if(debug) { cat("chained subsample step",i,":\n"); }
    if(!i %in% excluded.steps) {
      ans <- list(find.binding.positions(signal.data=signal.data,control.data=control.data,debug=debug, skip.control.normalization=T, ...));
      names(ans) <- as.character(c(i));
      calls <- c(calls,ans);
      callnames <- c(callnames,i);
    }
    # subsample
    if(step.size<1) {
      # fraction steps
      f <- 1-step.size;
    } else {
      # bin steps
      sz <- sum(unlist(lapply(signal.data,length)));
      f <- (sz-step.size)/sz;
      if(f<=0) break;
    }
    if(debug) { cat("chained subsampling using fraction",f,".\n"); }
    signal.data <- lapply(signal.data,function(d) sample(d,length(d)*f));
    if(subsample.control & !is.null(control.data)) {
      control.data <- lapply(control.data,function(d) sample(d,length(d)*f));
    }
    sz <- sum(unlist(lapply(signal.data,length)));
    if(sz<min.ntags) break;
  }
  names(calls) <- callnames;
  return(calls);
}


# chain-subsample dataset and calculate MSER interpolation
mser.chain.interpolation <- function(signal.data=NULL,control.data=NULL,chains=NULL,n.chains=5,debug=F, enrichment.background.scales=c(1,5), test.agreement=0.99, agreement.distance=50, return.median=F, mean.trim=0.1, enr.field="enr", return.lists=F, ...) {
  if(is.null(chains)) {
    cn <- c(1:n.chains); names(cn) <- cn;
    tf <- function(i, ...) get.subsample.chain.calls(signal.data,control.data,debug=debug, enrichment.background.scales=enrichment.background.scales, ...);
    chains <- lapply(cn,tf,...);
  } 
  names(enrichment.background.scales) <- enrichment.background.scales;
  lapply(enrichment.background.scales,function(scale) {
    actual.enr.field <- enr.field;
    if(scale>1) {
      actual.enr.field <- paste(actual.enr.field,scale,sep=".");
    }
      
    cvl <- lapply(chains,function(chain) {
      nn <- sort(unlist(lapply(chain,function(d) d$n)),decreasing=T);
      nd <- diff(nn);
      nn <- nn[-length(nn)];
      me <- lapply(c(2:length(chain)),function(i) {
        sla <- t.precalculate.ref.peak.agreement(chain[[i-1]],chain[i],agreement.distance=agreement.distance,enr.field=actual.enr.field)
        me <- t.find.min.saturated.enr(sla,thr=1-test.agreement)
        menr <- max(min(na.omit(unlist(lapply(chain[[i-1]]$npl,function(d) d[actual.enr.field])))),min(na.omit(unlist(lapply(chain[[i]]$npl,function(d) d[actual.enr.field])))),1)
        if(me<=menr) { me <- 1; };
        return(me);
      })
      data.frame(n=nn,me=unlist(me),nd=nd);
    });
    if(return.lists) { return(cvl) }
    cvl <- na.omit(do.call(rbind,cvl));
    if(return.median) {
      tv <- tapply(cvl$me,as.factor(cvl$n),median)
    } else {
      tv <- tapply(cvl$me,as.factor(cvl$n),mean,trim=mean.trim);
    }
    df <- data.frame(n=as.numeric(names(tv)),me=as.numeric(tv));
    return(df[order(df$n,decreasing=T),])
  })
}



# returns agreement as a function of dataset size, possibly filtering peaks by min.enr threshold, and by max.fdr
chain.to.reference.comparison <- function(chains,min.enr=NULL,debug=F,agreement.distance=50, return.median=F, mean.trim=0.1, enr.field="enr",max.fdr=NULL) {
  cvl <- lapply(chains,function(chain) {
    # filter chain by fdr
    if(!is.null(max.fdr)) {
      chain <- lapply(chain,function(d) { d$npl <- lapply(d$npl,function(cd) cd[cd$fdr<=max.fdr,]); return(d); });
    }
    nn <- sort(unlist(lapply(chain,function(d) d$n)),decreasing=T);
    nn <- nn[-length(nn)];
    me <- lapply(c(2:length(chain)),function(i) {
      sla <- t.precalculate.ref.peak.agreement(chain[[1]],chain[i],agreement.distance=agreement.distance,enr.field=enr.field)
      # calculate overlap
      x <- lapply(sla,function(mpd) {
        if(!is.null(min.enr)) {

          me <- mpd$re >= min.enr;
          me[is.na(me)] <- F;
          mpd <- mpd[me,];
          ome <- mpd$oe < min.enr;
          ome[is.na(ome)] <- T;
          mpd$ov[ome] <- 0;
        }
        return(mean(mpd$ov));
      })
    })
    
    data.frame(n=nn,me=unlist(me));
  });

  cvl <- na.omit(do.call(rbind,cvl));
  if(return.median) {
    tv <- tapply(cvl$me,as.factor(cvl$n),median)
  } else {
    tv <- tapply(cvl$me,as.factor(cvl$n),mean,trim=mean.trim);
  }
  df <- data.frame(n=as.numeric(names(tv)),me=as.numeric(tv));
  return(df[order(df$n,decreasing=T),])
}


# estimates enrichment confidence interval based on 2*tag.count.whs window around each position, and a z-score (alpha/2)
# if(multiple.background.scales=T) the enrichment is also estimated using 5- and 10-fold increased background tag window
# adds $enr (lower bound), $enr.ub (upper bound) and $enr.mle fields
calculate.enrichment.estimates <- function(binding.positions,signal.data=NULL,control.data=NULL,fraction=1,tag.count.whs=100,z=2,effective.genome.size=3e9,scale.down.control=F,background.scales=c(1),bg.weight=NULL) {
  f <- fraction;
  qv <- pnorm(z,lower.tail=F);
  cn <- names(binding.positions$npl); names(cn) <- cn;

  if(is.null(control.data)) {
    # estimate from gamma distribution
    fg.lambda <- f*sum(unlist(lapply(signal.data,length)))*2*tag.count.whs/effective.genome.size;
    binding.positions$npl <- lapply(binding.positions$npl,function(d) {
      d$enr <- qgamma(qv,d$nt,scale=1)/fg.lambda;
      d$enr.ub <- qgamma(1-qv,d$nt,scale=1)/fg.lambda;
      d$enr.mle <- d$nt/fg.lambda;
      return(d);
    });      
  } else {
    # estimate using beta distribution
    if(is.null(bg.weight)) {
      bg.weight <- sum(unlist(lapply(signal.data,length)))/sum(unlist(lapply(control.data,length)))
    }
    
    if(scale.down.control) {
      # sample down control to be the same size as true signal.data (bg.weight*f)
      control.data <- lapply(control.data,function(d) sample(d,length(d)*bg.weight*f,replace=(f*bg.weight>1)))
      #bg.weight <- sum(unlist(lapply(signal.data,length)))/sum(unlist(lapply(control.data,length)))
      bg.weight <- 1/f;
      
    }

    binding.positions$enrichment.bg.weight <- bg.weight;
    binding.positions$enrichment.whs <- tag.count.whs;
    binding.positions$enrichment.z <- z;
    
    binding.positions$npl <- lapply(cn,function(chr) {
      d <- binding.positions$npl[[chr]];
      
      edf <- lapply(background.scales,function(background.width.multiplier) {
        sig.mult <- bg.weight*f/background.width.multiplier;
        nbg <- points.within(abs(control.data[[chr]]),d$x-tag.count.whs*background.width.multiplier,d$x+tag.count.whs*background.width.multiplier,return.point.counts=T,return.unique=F);
        
        nfg <- d$nt;
      
        
        # Poisson ratio Bayesian LB with non-informative prior (Clopper & Pearson 1934)
        nf <- ((nfg+0.5)/(nbg+0.5))*qf(1-qv,2*(nfg+0.5),2*(nbg+0.5),lower.tail=F)
        nf <- nf/sig.mult;
        
        ub <- ((nfg+0.5)/(nbg+0.5))*qf(qv,2*(nfg+0.5),2*(nbg+0.5),lower.tail=F)
        ub <- ub/sig.mult;
        
        mle <- (nfg+0.5)/(nbg+0.5);
        mle <- mle/sig.mult;
        if(is.null(nbg)) { nbg <- numeric(0) }
        if(is.null(nf)) { nf <- numeric(0) }
        if(is.null(ub)) { ub <- numeric(0) }
        if(is.null(mle)) { mle <- numeric(0) }
        return(data.frame(nbg=nbg,lb=nf,ub=ub,mle=mle))
      })

      adf <- do.call(cbind,lapply(c(1:length(background.scales)),function(i) {
        df <- edf[[i]];
        cn <- c("nbgt","enr","enr.ub","enr.mle");
        if(background.scales[i]>1) {
          cn <- paste(cn,as.character(background.scales[i]),sep=".");
        }
        names(df) <- cn;
        return(df);
      }))

      return(cbind(d,adf));
    });
  }
  
  return(binding.positions);
}


# precalculate peak agreement of a sampling list given a reference
t.precalculate.ref.peak.agreement <- function(ref,sf,agreement.distance=50,enr.field="enr") {
  ref <- ref$npl;
  cn <- names(ref); names(cn) <- cn;

  # for each sampling round
  lapply(sf,function(sd) {
    # calculate overlap
      
    ov <- data.frame(do.call(rbind,lapply(cn,function(chr) {
      if(dim(ref[[chr]])[1]<1) { return(cbind(ov=c(),re=c(),oe=c())) };
      pwi <- points.within(ref[[chr]]$x,sd$npl[[chr]]$x-agreement.distance,sd$npl[[chr]]$x+agreement.distance);
      pwi[pwi==-1] <- NA;
      renr <- ref[[chr]][,enr.field]
      oenr <- sd$npl[[chr]][,enr.field][pwi];
      if(length(oenr)==0) { oenr <- rep(NA,length(renr)); }
      return(cbind(ov=as.integer(!is.na(pwi)),re=renr,oe=oenr));
    })))
  })
}


# find minimal saturated enrichment given a list of replicate agreement matrices (for one fraction)
t.find.min.saturated.enr <- function(pal,thr=0.01,plot=F,return.number.of.peaks=F,plot.individual=T,return.median=F,return.vector=F) {
  nr <- length(pal);
  # merge replicate data frames
  mpd <- data.frame(do.call(rbind,pal));

  mpd$re[is.na(mpd$re)] <- Inf;
  mpd$oe[is.na(mpd$oe)] <- Inf;

  

  # round up values to avoid miscounting
  mpd$re <- round(mpd$re,digits=2);
  mpd$oe <- round(mpd$oe,digits=2);

  me <- pmin(mpd$re,mpd$oe);
  ome <- order(me,decreasing=T);
  df <- data.frame(me=me[ome],ov=mpd$ov[ome]);
  recdf <- ecdf(-mpd$re); ren <- length(mpd$re);

  # collapse equal peak heights
  xk <- tapply(df$ov,as.factor(df$me),sum); xk <- data.frame(ov=as.numeric(xk),me=as.numeric(names(xk))); xk <- xk[order(xk$me,decreasing=T),];

  
  cso <- cumsum(xk$ov)/(recdf(-xk$me)*ren);
  cso[is.na(cso)] <- 0;
  cso[!is.finite(cso)] <- 0;
  mv <- max(which(cso >= 1-thr))
  menr <- xk$me[mv];

  ir <- lapply(pal,function(d) {
    d$re[is.na(d$re)] <- Inf;
    d$oe[is.na(d$oe)] <- Inf;
        
    me <- pmin(d$re,d$oe);
    ome <- order(me,decreasing=T);
    df <- data.frame(me=me[ome],ov=d$ov[ome]);
    cso <- cumsum(df$ov)/c(1:length(df$ov));
    mv <- max(which(cso >= 1-thr))
    menr <- df$me[mv];
    return(list(df=df,menr=menr));
  });

  if(plot) {
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
    plot(df$me,cumsum(df$ov)/c(1:length(df$ov)),type='l',ylab="fraction of positions overlapping with reference",xlab="minimal enrichment of binding positions",xlim=c(min(df$me),2*menr));
    abline(h=1-thr,lty=2,col=4)
    if(plot.individual) {
      lapply(ir,function(d) {
        df <- d$df;
        lines(df$me,cumsum(df$ov)/c(1:length(df$ov)),col=8);
        abline(v=menr,col="pink",lty=3)
      });
      lines(df$me,cumsum(df$ov)/c(1:length(df$ov)),col=1);
    }
    abline(v=menr,col=2,lty=2)
    legend(x="bottomright",lty=c(1,2,1,3,2),col=c(1,2,8,"pink",4),legend=c("combined samples","combined sample MSER","individual samples","individual MSERs","consistency threshold"));
  }

  if(return.number.of.peaks) {
    mpd <- data.frame(do.call(rbind,pal));
    return(length(which(!is.na(mpd$re) & mpd$re >=menr))/nr);
  } else {
    if(return.vector) {
      return(unlist(lapply(ir,function(d) d$menr)));
    }
    if(return.median) {
      return(median(unlist(lapply(ir,function(d) d$menr))));
    } else {
      return(menr);
    }
  }
}



# determine d1/d2 dataset size ratio. If background.density.scaling=F, the ratio of tag counts is returned.
# if background.density.scaling=T, regions of significant tag enrichment are masked prior to ratio calculation.
dataset.density.ratio <- function(d1,d2,min.tag.count.z=4.3,wsize=1e3,mcs=0,background.density.scaling=T) {
  if(!background.density.scaling) {
    return(sum(unlist(lapply(d1,length)))/sum(unlist(lapply(d2,length))))
  }

  chrl <- intersect(names(d1),names(d2));
  ntc <- do.call(rbind,lapply(chrl,function(chr) {
    x1 <- tag.enrichment.clusters(abs(d1[[chr]]),c(),wsize=wsize,bg.weight=0,min.tag.count.z=min.tag.count.z,mcs=mcs,either=F)
    x2 <- tag.enrichment.clusters(abs(d2[[chr]]),c(),wsize=wsize,bg.weight=0,min.tag.count.z=min.tag.count.z,mcs=mcs,either=F)
    return(c(length(which(points.within(abs(d1[[chr]]),c(x1$s,x2$s)-wsize/2,c(x1$e,x2$e)+wsize/2)==-1)),length(which(points.within(abs(d2[[chr]]),c(x1$s,x2$s)-wsize/2,c(x1$e,x2$e)+wsize/2)==-1))))
  }))
  ntcs <- apply(ntc,2,sum);
  #print(ntcs/c(sum(unlist(lapply(d1,length))),sum(unlist(lapply(d2,length)))));
  return(ntcs[1]/ntcs[2])
}

# returns effective size of the dataset based on the same logic as dataset.density.ratio
dataset.density.size <- function(d1,min.tag.count.z=4.3,wsize=1e3,mcs=0,background.density.scaling=T) {
  if(!background.density.scaling) {
    return(sum(unlist(lapply(d1,length))))
  }

  chrl <- names(d1);
  ntc <- lapply(chrl,function(chr) {
    x1 <- tag.enrichment.clusters(abs(d1[[chr]]),c(),wsize=wsize,bg.weight=0,min.tag.count.z=min.tag.count.z,mcs=mcs,either=F)
    return(length(which(points.within(abs(d1[[chr]]),x1$s-wsize/2,x1$e+wsize/2)==-1)))
  })
  return(sum(unlist(ntc)))
}

old.dataset.density.ratio <- function(d1,d2,min.tag.count.z=4.3,wsize=1e3,mcs=0,background.density.scaling=T) {
  if(!background.density.scaling) {
    return(sum(unlist(lapply(d1,length)))/sum(unlist(lapply(d2,length))))
  }
  
  t.chromosome.counts <- function(tl) {
    lapply(tl,function(d) {
      x <- tag.enrichment.clusters(abs(d),c(),wsize=wsize,bg.weight=0,min.tag.count.z=min.tag.count.z,mcs=mcs,either=F)
      x$s <- x$s-wsize/2; x$e <- x$e+wsize/2;
      x <- regionset.intersection.c(list(x),do.union=T)
      return(c(n=length(which(points.within(abs(d),x$s,x$e)==-1)),s=diff(range(abs(d))),m=sum(x$e-x$s)));
    })
  }

  l1 <- t.chromosome.counts(d1);
  l2 <- t.chromosome.counts(d2);

  l2 <- data.frame(do.call(rbind,l2[names(l1)]));
  l1 <- data.frame(do.call(rbind,l1));

  # genome size
  gs <- sum(pmax(l1$s,l2$s))

  den1 <- sum(l1$n)/(gs-sum(l1$m))
  den2 <- sum(l2$n)/(gs-sum(l2$m))
  return(den1/den2);
}




# calculate cumulative density based on sum of scaled gaussian curves
# (by Michael Tolstorukov)
#
# vin - input vector; bw -- standard deviation, dw-gaussina cutoff in stdev; dout - output "density")
# output - if return.x=F vector of cumulative density values corresponding to integer positions described by range(vin)
# output - if return.x=T a data structure with $x and $y corresponding to the cumulative density
# optional match.wt.f is a function that will return weights for a tag vector
densum <- function(vin,bw=5,dw=3,match.wt.f=NULL,return.x=T,from=min(vin),to=max(vin),step=1)    {
  # construct vector of unique tags and their counts
  tc <- table(vin[vin>=from & vin<=to]);
  pos <- as.numeric(names(tc)); storage.mode(pos) <- "double";
  tc <- as.numeric(tc); storage.mode(tc) <- "double";
  n <- length(pos)
  # weight counts
  if(!is.null(match.wt.f)) {
    tc <- tc*match.wt.f(pos);
  }
  
  rng <- c(from,to);
  if(rng[1]<0) { stop("range extends into negative values") }
  if(range(pos)[1]<0) { stop("position vector contains negative values") }

  storage.mode(n) <- storage.mode(rng) <- storage.mode(bw) <- storage.mode(dw) <- storage.mode(step) <- "integer";
  
  spos <- rng[1]; storage.mode(spos) <- "double";

  dlength <- floor((rng[2] - rng[1])/step) + 1; # length of output array
  if(dlength<1) { stop("zero data range") }
  dout <- numeric(dlength); storage.mode(dout) <- "double";
  storage.mode(dlength) <- "integer";
  .C("cdensum",n,pos,tc,spos,bw,dw,dlength,step,dout,DUP=F);
  
  if(return.x) {
    return(list(x=c(rng[1],rng[1]+step*(dlength-1)),y=dout,step=step))
  } else {
    return(dout)
  }
}

# count tags within sliding window of a specified size
# vin - tag vector (postive values, pre-shifted)
# window.size/window.step - window characteristics
# tv - optional, pre-sorted, pre-trimmed tag vector
window.tag.count <- function(vin,window.size,window.step=1,return.x=T,from=min(vin)+floor(window.size/2),to=max(vin)-floor(window.size/2),tv=NULL) {
  whs <- floor(window.size/2);
  # select tags with margins
  if(is.null(tv)) {
    tv <- sort(vin[vin>=from-whs-1 & vin<=to+whs+1])
  }
  storage.mode(tv) <- "double";
  n <- length(tv)
  nsteps <- ceiling((to-from)/window.step);
  
  storage.mode(n) <- storage.mode(nsteps) <- storage.mode(window.size) <- storage.mode(window.step) <- "integer";
  
  spos <- from; storage.mode(spos) <- "double";

  if(nsteps<1) { stop("zero data range") }
  #dout <- integer(nsteps); storage.mode(dout) <- "integer";
  #.C("window_n_tags",n,tv,spos,window.size,window.step,nsteps,dout,DUP=F);
  dout <- .Call("cwindow_n_tags",tv,spos,window.size,window.step,nsteps);
  
  if(return.x) {
    return(list(x=c(from,from+(nsteps-1)*window.step),y=dout,step=window.step))
  } else {
    return(dout)
  }
}

# count tags in windows around specified positions (pos)
window.tag.count.around <- function(vin,window.size,pos,return.x=T,tc=NULL,sorted=F) {
  if(is.null(tc)) {
    tc <- table(vin);
  }
  if(!sorted) {
    op <- rank(pos);
    pos <- sort(pos);
  }
  storage.mode(pos) <- "double";
  tpos <- as.integer(names(tc)); storage.mode(tpos) <- "double";
  tc <- as.integer(tc); storage.mode(tc) <- "integer";
  
  whs <- floor(window.size/2);
  
  storage.mode(whs) <- "integer";
  twc <- .Call("cwindow_n_tags_around",tpos,tc,pos,whs);
  if(return.x) {
    if(sorted) {
      return(data.frame(x=pos,y=twc));
    } else {
      return(data.frame(x=pos[op],y=twc[op]));
    }
  } else {
    if(sorted) {
      return(twc);
    } else {
      return(twc[op]);
    }
  }
}

# given a tag vector (signed), identify and clean up (either remove or cap) singular positions that exceed local tag density
# vin - tag vector
# cap.fold - maximal fold over enrichment over local density allowed for a single tag position, at which the tag count is capped
# eliminate.fold - max fold enrichment that, when exceeded, results in exclusion of all the tags at that position (e.g. counted as anomaly)
# z.threshold - Z-score used to determine max allowed counts
filter.singular.positions.by.local.density <- function(tags,window.size=200,cap.fold=4,eliminate.fold=10,z.threshold=3) {
  # tabulate tag positions
  if(length(tags)<2) { return(tags); };
  
  tc <- table(tags);
  pos <- as.numeric(names(tc)); storage.mode(pos) <- "double";
  tc <- as.integer(tc); storage.mode(tc) <- "integer";
  n <- length(pos); 

  whs <- floor(window.size/2);
  
  storage.mode(n) <- storage.mode(whs) <- "integer";
  twc <- .Call("cwindow_n_tags_around",pos,tc,pos,whs);
  twc <- (twc-tc+1)/window.size; # local density

  pv <- pnorm(z.threshold,lower.tail=F)
  # exclude
  max.counts <- qpois(pv,twc*eliminate.fold,lower.tail=F)
  tc[tc>max.counts] <- 0;
  # cap
  max.counts <- qpois(pv,twc*cap.fold,lower.tail=F)
  ivi <- which(tc>max.counts);
  tc[ivi] <- max.counts[ivi]+1;

  # reconstruct tag vector
  tv <- rep(pos,tc);
  to <- order(abs(tv)); tv <- tv[to];
  return(tv);
}



# calculates enrichment bounds using multiple background scales
# ft - foreground tags (pre-shifted, positive)
# bt - background tags
# fws - foreground window size
# bwsl - background window size list
# step - window step
# rng - from/to coordinates (to will be adjusted according to step)
#
# returns: a list with $x ($s $e $step), $lb vector and $mle vector ($ub if calculate.upper.bound=T)
mbs.enrichment.bounds <- function(ft,bt,fws,bwsl,step=1,rng=NULL,alpha=0.05,calculate.upper.bound=F,bg.weight=length(ft)/length(bt),use.most.informative.scale=F,quick.calculation=F,pos=NULL) {
  # determine range
  if(is.null(rng)) {
    rng <- range(range(ft));
  }
  # foreground counts
  if(is.null(pos)) {
    fwc <- window.tag.count(ft,fws,window.step=step,from=rng[1],to=rng[2],return.x=T);
  } else {
    fwc <- window.tag.count.around(ft,fws,pos,return.x=T)
  }
  fwc$y <- fwc$y+0.5;

  zal <- qnorm(alpha/2,lower.tail=F);

  # background counts
  bt <- sort(bt);
  if(!is.null(pos)) {
    tc <- table(bt);
  }
  bgcm <- lapply(bwsl,function(bgws) {
    if(is.null(pos)) {
      window.tag.count(bt,bgws,window.step=step,from=rng[1],to=rng[2],return.x=F,tv=bt)+0.5;
    } else {
      window.tag.count.around(bt,bgws,pos,return.x=F,tc=tc)+0.5
    }
  })
  if(!is.null(pos)) {
    rm(tc);
  }

  # pick most informative scale
  if(use.most.informative.scale) {
    bgcm <- t(do.call(cbind,bgcm))
    isi <- max.col(t((bgcm)/(bwsl/fws))) # add pseudo-counts to select lowest scale in case of a tie

    bgc <- c(bgcm)[isi+dim(bgcm)[1]*(c(1:length(isi))-1)]

    if(quick.calculation) {
      rte <- fwc$y+bgc-0.25*zal*zal; rte[rte<0] <- 0;
      dn <- bgc - 0.25*zal*zal;
      lbm=(sqrt(fwc$y*bgc) - 0.5*zal*sqrt(rte))/dn;
      ivi <- which(lbm<0);
      lbm <- lbm*lbm*bwsl[isi]/fws/bg.weight;
      lbm[rte<=0] <- 1;
      lbm[dn<=0] <- 1;
      lbm[ivi] <- 1;
    } else {
      lbm <- (fwc$y/bgc)*qf(1-alpha/2,2*fwc$y,2*bgc,lower.tail=F)*bwsl[isi]/fws/bg.weight;
    }
    
    mle <- fwc$y/bgc*bwsl[isi]/fws/bg.weight; mle[is.nan(mle)] <- Inf; mle[is.na(mle)] <- Inf;
    
    rl <- list(x=list(s=fwc$x[1],e=fwc$x[2],step=fwc$step),lb=lbm,mle=mle);
    
    if(calculate.upper.bound) {
      isi <- max.col(t((-bgcm)/(bwsl/fws))) # add pseudo-counts to select highest scale in case of a tie
      bgc <- c(bgcm)[isi+dim(bgcm)[1]*(c(1:length(isi))-1)]

      if(quick.calculation) {
        ubm=(sqrt(fwc$y*bgc) + 0.5*zal*sqrt(rte))/dn;
        ivi <- which(ubm<0);
        ubm <- ubm*ubm*bwsl[isi]/fws/bg.weight;
        ubm[rte<=0] <- 1;
        ubm[ivi] <- 1;
        lbm[dn<=0] <- 1;
      } else {
        ubm <- (fwc$y/bgc)*qf(alpha/2,2*fwc$y,2*bgc,lower.tail=F)*bwsl[isi]/fws/bg.weight;
      }
      rl <- c(rl,list(ub=ubm));
    }
    return(rl);
    
  } else {
    # determine lower bounds
    lbm <- lapply(c(1:length(bgcm)),function(i) {
      nbg <- bgcm[[i]];
      if(quick.calculation) {
        rte <- fwc$y+nbg-0.25*zal*zal; rte[rte<0] <- 0;
        dn <- (nbg - 0.25*zal*zal);
        lbm=(sqrt(fwc$y*nbg) - 0.5*zal*sqrt(rte))/dn;
        ivi <- which(lbm<0);  
        lbm <- lbm*lbm*bwsl[i]/fws/bg.weight;
        lbm[rte<=0] <- 1;
        lbm[dn<=0] <- 1;
        lbm[ivi] <- 1;
        return(lbm);
      } else {
        return((fwc$y/nbg)*qf(1-alpha/2,2*fwc$y,2*nbg,lower.tail=F)*bwsl[i]/fws/bg.weight);
      }
    })
    lbm <- do.call(pmin,lbm);

    # calculate mle
    #mle <- do.call(pmin,lapply(bgcm,function(bgc) fwc/bgc))
    mle <- do.call(pmin,lapply(c(1:length(bgcm)),function(i) {
      bgc <- bgcm[[i]];
      x <- fwc$y/bgc*bwsl[i]/fws/bg.weight; x[is.nan(x)] <- Inf; x[is.na(x)] <- Inf; return(x);
    }))

    rl <- list(x=list(s=fwc$x[1],e=fwc$x[2],step=fwc$step),lb=lbm,mle=mle);
    
    if(calculate.upper.bound) {
      # determine upper bound
      ubm <- lapply(c(1:length(bgcm)),function(i) {
        nbg <- bgcm[[i]];
        if(quick.calculation) {
          rte <- fwc$y+nbg-0.25*zal*zal; rte[rte<0] <- 0;
          dn <- (nbg - 0.25*zal*zal);
          ubm=(sqrt(fwc$y*nbg) + 0.5*zal*sqrt(rte))/dn;
          ivi <- which(ubm<0);  
          ubm <- ubm*ubm*bwsl[i]/fws/bg.weight;
          ubm[rte<=0] <- 1;
          ubm[dn<=0] <- 1;
          ubm[ivi] <- 1;
          return(ubm);
        } else {
          return((fwc$y/nbg)*qf(alpha/2,2*fwc$y,2*nbg,lower.tail=F)*bwsl[i]/fws/bg.weight);
        }
      })
      ubm <- do.call(pmax,ubm);
      rl <- c(rl,list(ub=ubm));
    }

    return(rl);
  }
}

write.probe.wig <- function(chr,pos,val,fname,append=F,feature="M",probe.length=35,header=T) {
  min.dist <- min(diff(pos));
  if(probe.length>=min.dist) {
    probe.length <- min.dist-1;
    cat("warning: adjusted down wig segment length to",probe.length,"\n");
  }
  mdat <- data.frame(chr,as.integer(pos),as.integer(pos+probe.length),val)

  if(header) {
    write(paste("track type=wiggle_0 name=\"Bed Format\" description=\"",feature,"\" visibility=dense color=200,100,0 altColor=0,100,200 priority=20",sep=""),file=fname,append=append)
    write.table(mdat,file=fname,col.names=F,row.names=F,quote=F,sep=" ",append=T);
  } else {
    write.table(mdat,file=fname,col.names=F,row.names=F,quote=F,sep=" ",append=append);
  }
  
}

# returns intersection of multiple region sets
# each regionset needs to contain $s, $e and optional $v column
regionset.intersection.c <- function(rsl,max.val=-1,do.union=F) {
  # translate into position/flag form
  rfl <- lapply(rsl,function(rs) {
    rp <- c(rs$s,rs$e); rf <- c(rep(c(1,-1),each=length(rs$s)));
    
    ro <- order(rp);
    rp <- rp[ro]; rf <- rf[ro];
    if(!is.null(rs$v)) {
      rv <- c(rs$v,rs$v)[ro];
      return(data.frame(p=as.numeric(rp),f=as.integer(rf),v=as.numeric(rv)));
    } else {
      return(data.frame(p=as.numeric(rp),f=as.integer(rf)));
    }
  })
  rfd <- data.frame(do.call(rbind,lapply(1:length(rfl),function(i) {
    d <- rfl[[i]]; d$f <- d$f*i; return(d);
  })))
  rfd <- rfd[order(rfd$p),];
  if(is.null(rfd$v)) { max.val <- 0; }
  if(do.union) { ur <- 1; } else { ur <- 0; }; 
  rl <- .Call("region_intersection",as.integer(length(rfl)),as.numeric(rfd$p),as.integer(rfd$f),as.numeric(rfd$v),as.integer(max.val),as.integer(ur));
  return(data.frame(do.call(cbind,rl)));
}


# idenfity if binding peak falls within a larger region of significant tag enrichment, and if so record its booundaries
add.broad.peak.regions <- function(chip.tags,input.tags,bp,window.size=500,z.thr=2) {
  se <- find.significantly.enriched.regions(chip.tags,input.tags,window.size=window.size,z.thr=z.thr,poisson.z=0,poisson.ratio=0,either=F)
  chrl <- names(bp$npl); names(chrl) <- chrl;
  bnpl <- lapply(chrl,function(chr) {
    npl <- bp$npl[[chr]];
    if(is.null(npl) | dim(npl)[1]<1) {
      return(npl);
    }
    pi <- points.within(npl$x,se[[chr]]$s,se[[chr]]$e,return.list=T);
    
    pm <- do.call(rbind,lapply(pi,function(rl) {
      if(length(rl)>0) {
        return(range(c(se[[chr]]$s[rl],se[[chr]]$e[rl])))
      } else {
        return(c(NA,NA));
      }
    }))

    npl$rs <- pm[,1];
    npl$re <- pm[,2];
    return(npl);
  })
  bp$npl <- bnpl;
  return(bp);
}

# writing out binding results in a narrowpeak format, incorporating broad region boundaries if they are present
# if broad region info is not present, margin is used to determine region width. The default margin is equal
# to the window half size used to call the binding peaks
write.narrowpeak.binding <- function(bd,fname,margin=bd$whs,npeaks=NA) { # Anshul: added npeaks option
  if(is.null(margin)) { margin <- 50; }
  chrl <- names(bd$npl); names(chrl) <- chrl;
  md <- do.call(rbind,lapply(chrl,function(chr) {
    df <- bd$npl[[chr]];
    x <- df$x;
    rs <- df$rs; if(is.null(rs)) { rs <- rep(NA,length(x)) }
    re <- df$re; if(is.null(re)) { re <- rep(NA,length(x)) }
    #ivi <- which(is.na(rs)); if(any(ivi)) {rs[ivi] <- x[ivi]-margin;}
    ivi <- which(is.na(rs)); if(any(ivi)) {rs[ivi] <- pmax(0,x[ivi]-margin);} # Anshul: added the pmax (0, ...) to avoid negative peak starts
    ivi <- which(is.na(re)); if(any(ivi)) {re[ivi] <- x[ivi]+margin;}
    #cbind(chr,rs,re,".","0",".",df$y,-1,format(df$fdr,scientific=T,digits=3),x-rs)
    cbind(chr,rs,re,".","0",".",df$y,-1,-log10(df$fdr),x-rs) # Anshul: converted fdr to -log10    
  }))
  md <- md[order(as.numeric(md[,7]),decreasing=T),]
  if (!is.na(npeaks)) { # Anshul: added this option to print a limited number of peaks
    npeaks <- min(nrow(md),npeaks)
  	md <- md[1:npeaks,]
  }  
  write.table(md,file=fname,col.names=F,row.names=F,quote=F,sep="\t",append=F);
}


get.broad.enrichment.clusters <- function(signal.data,control.data,window.size=1e3,z.thr=3, tag.shift=146/2,background.density.scaling=F, ... ) {
  # find significantly enriched clusters
  bg.weight <- dataset.density.ratio(signal.data,control.data,background.density.scaling=background.density.scaling);
  se <- find.significantly.enriched.regions(signal.data,control.data,window.size=window.size,z.thr=z.thr,tag.shift=tag.shift, bg.weight=bg.weight, ...)
  chrl <- names(se); names(chrl) <- chrl;
  se <- lapply(chrl,function(chr) {
    d <- se[[chr]];
    if(length(d$s>1)) {
      d <- regionset.intersection.c(list(d,d),do.union=T);
      sc <- points.within(abs(signal.data[[chr]]+tag.shift),d$s,d$e,return.point.counts=T);
      cc <- points.within(abs(control.data[[chr]]+tag.shift),d$s,d$e,return.point.counts=T);
      d$rv <- log2((sc+1)/(cc+1)/bg.weight);
      return(d);
    } else {
      return(d)
    }
  })
}

write.broadpeak.info <- function(bp,fname) {
  chrl <- names(bp); names(chrl) <- chrl;
  chrl <- chrl[unlist(lapply(bp,function(d) length(d$s)))>0]
  md <- do.call(rbind,lapply(chrl,function(chr) {
    df <- bp[[chr]];
    cbind(chr,df$s,df$e,".","0",".",df$rv,-1,-1)
  }))
  md <- md[order(as.numeric(md[,7]),decreasing=T),]
  write.table(md,file=fname,col.names=F,row.names=F,quote=F,sep="\t",append=F);
}


get.clusters2 <- function(x,CL)  {
  temp <- which(diff(x) != 0)
  begin <- c(1, temp + 1)
  end <- c(temp, length(x))
  size <- end - begin + 1

  begin <- begin[size >= CL]
  end <- end[size >= CL]
  size <- size[size >= CL]

  size <- size[x[end] != 0]
  begin <- begin[x[end] != 0]
  end <- end[x[end] != 0]

  return (list(size=size,begin=begin,end=end))
}
