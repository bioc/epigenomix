.calculateCrossCorrelation <- function(object, srange=c(100, 500), bin=10, minReads=100000, mc.cores=1) {

  object <- resize(object, 1)
  object <- split(object, seqnames(object))
  reads <- lapply(object, function(o) {
               tags <- start(o)
               mInd <- as.vector(strand(o) == "-")
               tags[mInd] <- tags[mInd] * -1
               return(tags)
           })
  numReads <- sapply(reads, length)
  reads <- reads[numReads >= minReads]
  numReads <- numReads[numReads >= minReads]
 

  # internal function taken from package spp (ver. ), which is not at CRAN unfortunately :(
  # see http://compbio.med.harvard.edu/Supplements/ChIP-seq/ (Peter Kharchenko)
  tag.scc <- function(tags, srange, bin, llim=10) {
    tt <- table(sign(tags)*as.integer(floor(abs(tags)/bin+0.5)))
  
    if(!is.null(llim)) { l <- mean(tt); tt <- tt[tt<llim*l] }
    tc <- as.integer(names(tt))
    tt <- as.numeric(tt)

    pv <- tt; pv[tc<0] <- 0
    nv <- tt; nv[tc>0] <- 0

    pti <- which(tc>0)
    nti <- which(tc<0)

    ptc <- tc[pti]
    ntc <- (-1)*tc[nti]

    ptv <- tt[pti]
    ntv <- tt[nti]

    trng <- range(c(range(ptc),range(ntc)))
    l <- diff(trng)+1
    rm(tc,tt)

    mp <- sum(ptv)*bin/l;   mn <- sum(ntv)*bin/l
    ptv <- ptv-mp; ntv <- ntv-mn
    ss <- sqrt((sum(ptv*ptv)+(l-length(ptv))*mp^2) * (sum(ntv*ntv)+(l-length(ntv))*mn^2))

    t.cor <- function(s) {
      smi <- match(ptc+s,ntc)
      return((sum(ptv[!is.na(smi)]*ntv[na.omit(smi)]) -
             mn*sum(ptv[is.na(smi)]) -
             mp*sum(ntv[-na.omit(smi)]) +
             mp*mn*(l-length(ptv)-length(ntv)+length(which(!is.na(smi)))))/ss)
    }
    shifts <- floor(seq(srange[1],srange[2],by=bin)/bin+0.5)
    scc <- unlist(lapply(shifts,t.cor)); names(scc) <- shifts*bin
    return(scc)
  }

  
  res <- mclapply(reads, tag.scc, srange=srange, bin=bin, mc.cores=mc.cores)
  weights <- numReads / sum(numReads)
  resMat <- t(matrix(unlist(res), ncol=length(res)))
  colnames(resMat) <- names(res[[1]])
  res <- apply(resMat * weights, 2, sum)
  return(res)
}
  








#.calculateCrossCorrelation <- function(object, shift=100, chrs=NA, bin=1, minReads=0, mc.cores=1) {
#
#  if (any(shift < 0)) {stop("shift musst be larger or equal 0")} 
#
#  object <- resize(object, width=1)
#  object <- split(object, seqnames(object))
#  if (!is.na(chrs[1])) {
#    object <- object[is.element(names(object), chrs)]
#  }
#  numReads <- sapply(object, length)
#  object <- object[numReads >= minReads]
#  numReads <- numReads[numReads >= minReads]
#
#  chrCor <- mclapply(object, function(r) {
#      r <- keepSeqlevels(r, value=as.character(seqnames(r)[1]))
#      maxPos <- max(start(r))
#      indPlus <- strand(r) == "+"
#      
#      positive <- coverage(r[indPlus], width=maxPos)[[1]]
#      negative <- coverage(r[!indPlus], width=maxPos+max(shift))[[1]]
#
#      if (bin > 1) {
#        positive <- positive[seq(1, length(positive), bin)]
#        negative <- negative[seq(1, length(negative), bin)]
#        shift <- shift/bin
#        if (!all.equal(shift, as.integer(shift))) {
#          error("Argument shift must be a multiple of argument bin (including 00).")
#        }
#      }
#      
#      if (!all(negative[1:max(shift)] == 0)) {
#	missNulls <- max(shift) - which(negative[1:max(shift)] != 0)[1] + 1
#        positive <- c(Rle(as.integer(rep(0, missNulls))), positive)
#	negative <- c(Rle(as.integer(rep(0, missNulls))), negative)
#        maxPos <- maxPos + missNulls
#      }
#
#      meanPos <- mean(positive)
#      meanNeg <- mean(negative[1:maxPos])
#      positive <- positive - meanPos
#      negative <- negative - meanNeg
#      denom <- sqrt(sum(positive^2)) * sqrt(sum(negative[1:maxPos]^2))
#
#      cors <- numeric()
#      for (i in (shift+1)) {
#        cors <- c(cors, sum(positive * negative[i:(maxPos+i-1)]))
#      }
#      cors <- cors / denom
#      names(cors) <- as.character(shift)
#      return(cors)
#  }, mc.cores=mc.cores)
#
#  cors <- matrix(unlist(chrCor), nrow=length(chrCor), ncol=length(chrCor[[1]]))
#  weights <- numReads / sum(numReads)
#  cors <- apply(cors * weights, 2, sum)
#  names(cors) <- names(chrCor[[1]])
#  
#  return(cors)
#}

setMethod("calculateCrossCorrelation",
    signature=c(object="GRanges"),
    .calculateCrossCorrelation)
