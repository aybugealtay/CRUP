

##################################################################
#definition: headline
##################################################################

headline <- "\n*****************************************************************************************

  C                               *** This is CRUP ***                                C
  *             (C)ondition-specific (R)egulatory (U)nits (P)rediction                *	
  R                                                                                   R	
  *                                                                                   *	
  U                           (version 1.0 - 11/10/2018)                              U	
  *                                                                                   * 
  P                                     contact:                                      P 
  *                 ramisch@molgen.mpg.de,heinrich@molgen.mpg.de                      * 

*****************************************************************************************\n\n"

##################################################################
# definition: parameter spectrum:
##################################################################
        
spectrum <- matrix(c( 'norm',           'N', 0, "logical",       "computes normalized count values for .bam files",
                      'prediction',     'P', 0, "logical",       "runs CRUP - EP: (E)nhancer (P)rediction from histone modification",
                      'dynamics',       'D', 0, "logical",       "runs CRUP - ED: assigns (E)nhancer to (D)ynamic conditions",
                      'targets',        'T', 0, "logical",       "runs CRUP - ET: correlates (E)nhancer to (T)arget genes",
                      'cores',          'x', 1, "integer",       "number of cores to use (DEFAULT:1)",
                      'file',           'f', 1, "character",     "summary text file for ChIP-seq experiments",
                      'matrix',         'm', 1, "character",     "normalized data matrix (rds file format)",
                      'classifier',     'c', 1, "character",     "directory of enhancer classifier",
                      'cutoff',         'u', 1, "double",        "cutoff for probabilities [0,1] (DEFAULT: 0.5)",
                      'distance',       'd', 1, "integer",       "maximum distance (bp) for peak clustering (DEFAULT: 12500)",
                      'genome',         'g', 1, "character",     "genome used in the .bam files ('hg19', 'mm10' or 'mm9')",
                      'sequencing',     's', 1, "character",     "type of sequencing ('paired' or 'single')",
                      'outdir',         'o', 1, "character",     "output directory (DEFAULT: same as 'file' directory)",
                      'probabilities',  'p', 1, "character",     "probabilities in rds format. list: delimiter samples: ':', delimiter conditions: ','",
                      'names',          'n', 1, "character",     "aternative labels for conditions (DEFAULT: cond1,cond2, ..)",
                      'w_0',            'w', 1, "double",        "minimum difference between group means [0,1]. (DEFAULT: 0.5)",
                      'threshold',      't', 1, "double",        "threshold for p-values in [0,1]. (DEFAULT: 0.01)",
                      'threshold_c',    'C', 1, "double",        "threshold for correlation in [0.5,1]. (DEFAULT: 0.9)",
                      'len',            'l', 1, "integer",       "length of flanking region for summarizing. (DEFAULT: 1000)",
                      'regions',        'r', 1, "character",     "text file with condition-specific regions in txt format",
                      'RNA',            'E', 1, "character",     "RNA-seq experiments in bam format. list: delimiter samples: ':', delimiter conditions: ','",
                      'expression',     'e', 1, "character",     "gene expression counts for all samples and conditions",
                      'TAD',            'b', 1, "character",     ".bed file with TADs (DEFAULT: DATA/mESC_mapq30_KR_all_TADs.bed)",
                      'mapq',           'q', 1, "integer",       "minimum mapping quality (DEFAULT:10)",
                      'help',           'h', 0, "logical",       "this help message"
),ncol = 5,byrow = T)


##################################################################
# definition: valid parameter values:
##################################################################

sequencing_values <- c('paired', 'single')
genome_values <- c('hg19', 'mm10', 'mm9')

##################################################################
#definition: output messages
##################################################################

done <- function() {cat(".. done.\n")}
skip <- function() {cat("\t ..")}
startPart <- function(m) {cat(paste0("\n--- ",m," ---\n\n"))}
endPart <- function() {cat("\n\t>>> All done!\n")}

##################################################################
#definition: standard GRanges/DataFrame header:
##################################################################

GR_header <- c("seqnames", "start","end","width")
GR_header_short <- c("seqnames", "start","end")
DF_header <- c("chr", "start","end")

##################################################################
# function: re-header
##################################################################

reheader_DF <- function(DF, header){
  colnames(DF)[1:length(header)] <- header
  return(DF)
}

##################################################################
# function: check if file exists
##################################################################

check_file <- function(f){
  if (!(file.exists(f))) { 
    cat("File ",f," does not exist.\n");
    q();
  }
}

##################################################################
# function: check if outdir exists
#
#(d = outdir, alt_d = alternative dir)
##################################################################

check_outdir <- function(d, alt_d){
  if (is.null(d)) d <- dirname(alt_d)
  if (!dir.exists(d)) {
    cat(paste0("Output directory '",d,"'is not a valid directory. \n
              Output directory is set to ",dirname(alt_d)));
    d = paste0(dirname(alt_d),"/")
  }
  return(d)
}

##################################################################
# function: stop_script
##################################################################

stop_script <- function(parameter){
  cat(paste(getopt(spectrum[(spectrum[,1] %in% c(parameter, "help")),], usage = T),"\n"));
  q();
}

##################################################################
# function: check if package is installed
##################################################################

pkgTest <- function(pkg){
  if (!pkg %in% installed.packages()) {
    stop(paste0("package ",pkg, " not found."),call. = FALSE)
  }
} 

##################################################################
# function: load package
##################################################################

pkgLoad <- function(pkg) {
  pkgTest(pkg)
  cat(paste0(skip(), "load package ", pkg))
  suppressMessages(library(pkg, character.only = TRUE))
  done()
}

##################################################################
# function: partition genome into bins
##################################################################

get_binned_genome <- function(txdb, width){
  
  # get binned genome from txdb object
  binned <- tileGenome(seqinfo(txdb),
                       tilewidth = width,
                       cut.last.tile.in.chrom = TRUE)
  
  # only take autosomes and X chromosome
  seqlevels(binned, pruning.mode = "coarse") <- seqlevels(txdb)[grep("^chr[0-9]{,2}$|chrX$", seqlevels(txdb))]
  
  # delete last region in each chromosome that is smaller than 100
  return(binned[-which(width(binned) != 100)])
}


##################################################################
# function: get summarized counts in defined bins
##################################################################

get_bamProfile <- function(bampath, gr, mapq, sequencing){
  
  # gives information about path of bam and bai file (Object)
  bf = BamFile(bampath)
  
  # gives length of different chromosomes (SeqInfo class)
  si = seqinfo(bf)
  
  # fix chromosome prefix
  seqlevels(gr) = gsub("chr","",seqlevels(gr))
  if (grepl("chr",seqlevels(si)[1])) {
    seqlevels(gr) = paste0("chr", seqlevels(gr))
  }
  
  # filter chromosomes
  gr = gr[seqnames(gr) %in% seqnames(si)]
  
  # count and summarize
  if (sequencing == "paired") {
    sapply( bamProfile( bampath = bampath, 
                        gr = gr, 
                        binsize = 100, 
                        mapqual = mapq, 
                        ss = FALSE, 
                        paired.end = "midpoint", 
                        filteredFlag = 1024,
                        verbose = FALSE)@signals, 
            function(x) x)
    
  } else if (sequencing == "single") {
    sapply( bamProfile( bampath = bampath, 
                        gr = gr, 
                        binsize = 100, 
                        mapqual = mapq, 
                        ss = FALSE, 
                        shift = 100,
                        filteredFlag = 1024,
                        verbose = FALSE)@signals, 
            function(x) x)
  }
}

##################################################################
# function: quantile normalize with target
##################################################################

get_targetQuantileNorm <- function(ecdf){
  
  ### recreate HM count vector from ecdf of mESC
  x.ref <- knots(ecdf)
  y.ref <- ecdf(x.ref)
  
  # 26337756 = nb of 100 bp bins in mm10
  temp <- c(0, round(y.ref*26337756, digits = 0))
  ret <- unlist(lapply(1:length(x.ref), function(x) rep(x.ref[x], temp[x + 1] - temp[x])))
  
  return(ret)
}


##################################################################
# function: create extended data matrix (plus/minus) bins
##################################################################

extend_dataMatrix <- function(N, df, f){
  
  N_regions <- nrow(df)
  
  # make extended feature vector
  f_ext <- NULL
  for (i in 1:N) {
    f_ext <- c(f_ext, paste0(f, "_left", i ))
    f_ext <- c(f_ext, paste0(f, "_right", i ))
  }
  
  # prepare new matrix
  df_ext <- cbind(df[,c(GR_header_short, f)], matrix(0, nrow = N_regions, ncol = length(f_ext)))
  colnames(df_ext) <- c(DF_header, f, f_ext)
  
  # make extended data matrix
  region_ext <- (N + 1):(N_regions - N)
  
  for (i in 1:N) {
    df_ext[region_ext, f_ext[((2*i - 2)*length(f) + 1):((2*i - 1)*length(f))]] <- df[region_ext - i, f]
    df_ext[region_ext, f_ext[((2*i - 1)*length(f) + 1):(2*i*length(f))]] <-  df[region_ext + i, f]
  }
  
  return(df_ext)
}


#############################################################################
# function: sort peak candidates and remove overlapping peaks
# used in peak_calling function
#############################################################################

sort_peaks <- function(peaks){
  
  # sort peaks according to probability score
  peaks <- peaks[sort(score(peaks), decreasing = T, index.return = T)$ix]
  
  count <- 0
  while (length(peaks) > (count + 1)) {
    
    count <- count + 1
    overlap.to <- findOverlaps(query = peaks[count], subject = peaks)@to
    
    if (length(overlap.to) == 1) next
    
    delete.index <- sort(overlap.to, decreasing = F)[-1]
    peaks <- peaks[-delete.index]
  }
  return(peaks)
}


#############################################################################
# function: call peaks (from probabilities)
#############################################################################

get_enhancerPeaks <- function(gr, cutoff, cores){
  
  # define and merge all regions under background probability 
  candidates <- reduce(gr[which(score(gr) > cutoff)])
  
  # possible 100 bp peaks
  peaks <- gr[findOverlaps(gr, candidates)@from]
  
  # create 1100 bp peak regions
  start(peaks) <- start(peaks) - 500
  width(peaks) <- 1100
  
  # call peaks for each chromosome separately
  out <- mclapply(split(peaks, seqnames(peaks)), 
                  sort_peaks, 
                  mc.cores = cores)
  
  return(do.call("c", unname(out)))
}


##################################################################
# function: call super enhancer candidates from called peaks
##################################################################

get_enhancerCluster <- function(peaks, peak.gap, cores){
  
  # cluster candidates (neighboring enhancers within 12.5 kb)
  red.peaks <- reduce(peaks, min.gapwidth = peak.gap)
  cluster <- red.peaks[which(width(red.peaks) > 1100)]
  
  # sort cluster accroding to maximum probability
  if (length(cluster) > 0) {
    sort.max <- unlist(mclapply(1:length(cluster), 
                                function(x) max(peaks$score[findOverlaps(cluster[x], peaks)@to]),
                                mc.cores = cores))
    cluster <- cluster[sort(sort.max, index.return = T, decreasing = T)$ix]
  }  
  return(cluster)
}


##################################################################
# function: create probabilty matrix for all given files:
##################################################################
get_probabiltyMatrix <- function(files, IDs){
  
  probs <- GRanges()
  for (i in 1:length(files)) {
    for (j in 1:length(files[[i]])) {
      
      # read rds file:
      this <- readRDS(files[[i]][j])
      colnames(elementMetadata(this)) <- IDs[[i]][j]
      
      # combine:
      if(length(probs) > 0 ){
      	probs <- merge(probs, this, all = TRUE)
      }else{
	probs <- this
	}
    }
  }
  return(probs)
}


##################################################################
# function: get uninformative row indices
##################################################################

get_invalidIdx <- function(gr, w_0){
  idx = which(rowVars(as.matrix(mcols(gr))) == 0 | rowMaxs(as.matrix(mcols(gr))) < w_0)
  return(idx)
}


##################################################################
# function: generate unique key (for numbers in a matrix)
##################################################################

get_uniqueKey <- function(m, IDs){

  # length of (max) digits (after '.')
  len_digits <- max(nchar(formatC(as.matrix(mcols(m)[1])))) + 1
  
  # round
  m <- as.matrix(round(as.matrix(mcols(m)[,IDs])*10^(len_digits - 1),0))
  
  step <- 1/(10^(len_digits))
  digits <- dim(m)[2]*(len_digits)
  max.digits <- 15 #maximal number of digits...
  
  columns <- round(max.digits/(len_digits),0)
  res.final <- c()
  
  count <- 1
  repeat {
    count.max <- (count - 1) + columns
    
    if (count.max >= dim(m)[2]) {
      columns <- dim(m)[2] - count + 1
    }
    this.digits <- columns*(len_digits)
    options(digits = this.digits)
    
    j <- 1
    res <- 0
    for (i in count:(count + columns - 1)) {
      if (length(res) == 0) {
        res <- as.matrix(m[,i])*step^j
      }else{
        res <- res + as.matrix(m[,i])*step^j
      }
      j <- j + 1
    }
    count <- count + this.digits/(len_digits)
    
    # create character:
    res_char <- formatC(res*10^this.digits,
                        width = this.digits,
                        digits = 0, 
                        format = "f",
                        flag = "0",
                        mode = "double"
                      )
    
    # combine results:
    if (length(res.final) == 0) {
      res.final <- res_char
    }else{
      res.final <- paste0(res.final, res_char)
    }
    
    if (nchar(res.final[1]) == (dim(m)[2]*len_digits)) break
  }
  return(res.final)
}


##################################################################
# function: compute t test statistic
#
# #mean_a -mean_b <= w_0 == differenz is smaller or equal than w_0
##################################################################

get_tStatistic <- function(a, b, w_0){
  
  a <- as.matrix(a)
  b <- as.matrix(b)
  
  len <- dim(a)[1] 
  
  n_a <- dim(a)[2]
  n_b <- dim(b)[2]
  
  mean_a <- rowMeans(a,na.rm = T)
  mean_b <- rowMeans(b,na.rm = T)
  
  if (n_a == 1) {
    var_a <- rep(0,len)
  }else{
    var_a  <- rowVars(as.matrix(a),na.rm = T)
  }
  
  if (n_b == 1) {
    var_b <- rep(0,len)
  }else{
    var_b <- rowVars(as.matrix(b),na.rm = T)
  }
  
  if (n_a == 1 & n_b == 1) {
    S <- rep(1/100000, len) # a small number
  }else{
    S_2 <- ((n_a - 1)*var_a + (n_b - 1)*var_b)/(n_a + n_b - 2)
    S   <-  sqrt(S_2)*sqrt(1/n_a + 1/n_b)
  }
  
  S[which(S == 0)] <- 1/100000 # a small number
  diff <- (mean_a - mean_b)
  res  <- sqrt((n_a*n_b)/(n_a + n_b))*(abs(diff) - w_0)/S 
  
  #generate complete cases:
  res[is.infinite(res)]  <-  NA
  res[res < 0]  <-  NA
  
  # assign direction to pairwise comparison:
  res <- as.list(res*sign(diff))
  
  return(res)
}  


##################################################################
# function: get test statistic distribution
##################################################################

get_distribution <- function(values, keys, IDs.1, IDs.2, key.1, key.2, w_0){
  
  keys.comb <- paste0(keys[[key.1]],keys[[key.2]])
  distr <- get_tStatistic( mcols(values[duplicated(keys.comb) == F])[, IDs.1], 
                           mcols(values[duplicated(keys.comb) == F])[, IDs.2], 
                           w_0
  )
  names(distr) <- keys.comb[duplicated(keys.comb) == F]
  distr <- distr[keys.comb]
  
  return(distr)
}


##################################################################
# function: compute empirical p values
##################################################################

compute_empPvalue <- function(l, dim, x){
  
  tab <- table(abs(unlist(l)))
  res <- 1
  if (!is.na(x)) {
    idx <- which(as.numeric(names(tab)) >= abs(x))
    res <- (sum(tab[idx]) + 1)/(dim + 1)
    res <- (res*sign(x))
  }
  return(res)
}


##################################################################
# function: get empirical p values
##################################################################

get_empPvalues <- function(i, w_0, comb, IDs, dim, probs.valid, invalid.idx, probs.shuffle.valid, keys, keys.shuffle){
  
  # get IDs of i th condition
  IDs.1 = IDs[[comb[,i][1]]]
  IDs.2 = IDs[[comb[,i][2]]]
  
  #calculate t test statistic for shuffled columns
  null.distr <- get_distribution(probs.shuffle.valid, keys.shuffle, IDs.1, IDs.2, comb[,i][1], comb[,i][2], w_0 )
  
  #calculate t test statistic for original values
  distr <- get_distribution(probs.valid, keys, IDs.1, IDs.2, comb[,i][1], comb[,i][2], w_0 )
  
  #calculate empirical p values for unique values:
  p.values <- unlist(lapply( unique(distr), 
                      function(x) compute_empPvalue(null.distr, dim, x)))

  df=t(data.frame(p.values))
  colnames(df) = formatC(as.character(unique(unlist(distr))))
  tmp=df[,formatC(as.character(unlist(distr)))]

  res <- rep(1, dim)
  res[-invalid.idx] <- tmp

  #names(p.values) <- formatC(as.numeric(unique(distr)), format = "f")
  #p.values <- p.values[formatC(as.numeric(distr), format = "f")]
  #
  #res <- rep(1, dim)
  #res[-invalid.idx] <- as.numeric(unlist(p.values))
  
  return(res)
}


##################################################################
# function: get pairwise p values
##################################################################

get_pairwisePvalues <- function(probs, IDs, w_0, cores){

  # shuffle each column separately
  probs.shuffle <- probs
  elementMetadata(probs.shuffle) <-  apply(mcols(probs.shuffle), 2, sample)
 
  # define rows with uninformative entries
  invalid.idx = get_invalidIdx(probs, w_0)
  probs.valid <- probs[-invalid.idx]
  
  invalid.shuffle.idx = get_invalidIdx(probs.shuffle, w_0)
  probs.shuffle.valid <- probs.shuffle[-invalid.shuffle.idx]
  
  # create unique keys for all samples in one condition
  keys <- mclapply( as.list(seq(length(IDs))), mc.cores = cores,
                    function(i) get_uniqueKey(probs.valid, IDs[[i]]))
  
  keys.shuffle <- mclapply( as.list(seq(length(IDs))), mc.cores = cores,
                            function(i) get_uniqueKey(probs.shuffle.valid, IDs[[i]]))
  
  # permutation test for all pairwise conditions
  comb <- combn(seq(length(IDs)),2)
  p.values <- mclapply( as.list(seq(dim(comb)[2])), mc.cores = cores,
                        function(i) get_empPvalues(i,
						   w_0, 
                                                   comb,
                                                   IDs,
                                                   length(probs),
                                                   probs.valid,
                                                   invalid.idx,
                                                   probs.shuffle.valid,
                                                   keys,
						   keys.shuffle)
                        )
  return(p.values)
}


##################################################################
# function: function to translate significant peak pattern
#
# output: vector \in [0,1]
##################################################################

get_positionPattern <- function(p, i, t, pattern_empty, comb){
  
  # define region with flanking positions
  pos <- c(rev(i - 1:2), i, (i + 1:2))
  
  # all pairwise p-values at position i and +/- 2 positions:
  this <- do.call(rbind, lapply(p, function(x) x[pos]))
  idx <- (abs(this) <= t) == T
  idx.prod <- rowProds(idx)
  
  # all flanking positions have to be significant as well:
  if (sum(idx.prod) > 0) {
    this.comb <- as.matrix(comb[,(idx.prod == 1)])
    
    # check direction:
    idx.invert <- this[(idx.prod == 1),3] < 0
    if (sum(idx.invert) > 0) {
      this.comb[,idx.invert] <- this.comb[c(2,1),idx.invert]
    }
    label <- apply(this.comb,2, function(x) paste(x, collapse = ","))
    pattern <- pattern_empty
    pattern[,label] <- 1
    
    return(data.frame(pos, pattern))
  }
}

##################################################################
# function: search for significant peaks (in pairwise comparisons)
##################################################################

get_peakPattern <- function(p, threshold, IDs, cores){
  
  # specifiy all possible pairwise combinations:
  comb <- combn(seq(1:length(IDs)),2)
  comb_full <- cbind(comb, comb[c(2,1),])
  comb_full <- comb_full[,order(comb_full[1,], comb_full[2,])]
  
  # emtpy pattern vector:
  pattern_empty <- data.frame(t(rep(0,dim(comb_full)[2])))
  colnames(pattern_empty) <- apply(comb_full,2, function(x) paste(x,collapse = ","))
  
  # define positions that are significant between any two conditions
  idx.sign <- which((rowMins(abs(do.call(cbind, p))) <= threshold) == T)
    
  if (length(idx.sign) > 0) {
    
    # generate list of p value pattern for each bin
    list <- mclapply(idx.sign,  mc.cores = cores, 
                     function(i) get_positionPattern(p, i, threshold, pattern_empty, comb))
  
    # remove zero entries
    null.idx <- which(sapply(list, is.null))
    if (length(null.idx) > 0) list <- list[-null.idx]
  
    # create matrix from list
    mat <- do.call(rbind, list)
  
    # define duplicated rows
    idx.duplicated <- which(table(mat[,1]) > 1 )
    
    # create final matrix
    mat1 <- mat[!(mat[,1] %in% names(idx.duplicated)),]
    rownames(mat1) <- mat1[,1]
    colnames(mat1)[2:dim(mat1)[2]] <- seq((dim(mat1)[2] - 1))
    
    # combine vectors of duplictaed rows
    # (e.g.: v1 = (0,1,1,0), v2 = (0,1,0,0). v_combined=(0,1,0,0))
    mat_dup <- mat[(mat[,1] %in% names(idx.duplicated)),]
    mat2 <- do.call(rbind, lapply(split(mat_dup, mat_dup[,1]), function(x) as.integer(apply(x[,-1],2 , sum) > 0)))
    colnames(mat2) <- seq(dim(mat2)[2])
  
    # combine and remove zero lines
    mat_final <- rbind(mat1[,-1], mat2)
    mat_final <- mat_final[which(rowSums(mat_final) > 0),]
    
    return(mat_final)
  }
}


##################################################################
# function: combine regions
##################################################################

combine_peaks <- function(x, flank){
  
  if (length(x) > 1) {
    
    x.flank = x
    start(x.flank) = start(x.flank) - flank 
    end(x.flank) = end(x.flank) + flank
    
    red <- reduce(x.flank,with.revmap = T)
    best.idx <- unlist(lapply(mcols(red)$revmap, function(y) y[which.min(mcols(x)$score[y])]))
    
    red$score <- x$score[best.idx]
    red$index <- x$index[best.idx]
    red$pattern <- x$pattern[best.idx]
    mcols(red) <- mcols(red)[,-1]
    
    start(red) <- start(red) + flank 
    end(red) <- end(red) - flank
    
    return(red)
  }else{
    return(x)
  }
}

##################################################################
# function: combine peaks by significance pattern
##################################################################

get_combinedDiffPeaks <- function(gr, p, pattern, IDs, len){
 
  # reduce GRanges()
  idx <- as.numeric(rownames(pattern))
  gr <- gr[idx,]
  
  # get maximum pvalue from all pairwise comparisons
  score <- lapply(p, function(x) x[idx])
  score <- do.call(cbind, score)
  
  # add score, bin index and pattern to each bin:
  mcols(gr)$score <- rowMins(abs(score))
  mcols(gr)$index <- idx
  mcols(gr)$pattern <- apply(pattern, 1, function(x) paste0(x, collapse = ""))
  
  # split by p value pattern
  peaks.split = split(gr, mcols(gr)$pattern)
  
  # combine regions
  peaks = mclapply(	peaks.split, 
                    function(x) combine_peaks(x, flank = len), 
                    mc.cores = cores
  )
  peaks <- c(do.call("c", unname(peaks)))
  peaks <- combine_peaks(peaks, flank = 0)
  
  # count pattern occurence:	
  tab <- sort(table(mcols(peaks)$pattern), decreasing = T)
  
  # re-sort (by unique numbers):
  pattern <- do.call(rbind,lapply(strsplit(names(tab),""), as.numeric))
  
  # order by pattern:
  pattern.summarized <- c()
  for (i in 1:length(IDs)) {
    this <- as.matrix(pattern[,(1 + (i - 1)*(length(IDs) - 1)):((i - 1)*(length(IDs) - 1) + (length(IDs) - 1))])
    pattern.summarized <- cbind(pattern.summarized, rowSums(this))
  }
  pattern.order <- rev(do.call(order, as.list(data.frame(pattern.summarized))))
  
  mcols(peaks)$cluster = 0
  mcols(peaks)$cluster.size = 0
  
  cluster.count <- 1
  for (i in pattern.order) {
    mcols(peaks)$cluster[which(mcols(peaks)$pattern == names(tab[i]))]  <-  cluster.count
    mcols(peaks)$cluster.size[which(mcols(peaks)$pattern == names(tab[i]))]  <-  as.numeric(tab[i])
    cluster.count  <-  cluster.count + 1
  }
  peaks <- peaks[order(mcols(peaks)$cluster)]
  
  return(peaks)
}


##################################################################
# function: adjust y positioning for heatmap (y axis):
##################################################################

adjust <- function(v){
  
  mid = v[1]
  res = c(mid)
  for (i in 2:length(v)) {
    mid = (mid + v[i - 1]/2 + v[i]/2)
    res = c(res,mid)
  }
  return(res)
}

##################################################################
# function: plot heatmap of all differential enhancer regions
##################################################################

plot_heatmap <- function(mat, IDs, color_low, color_mid, color_high, x_axis_label, y_axis_label, legend_label, out.file){
  
  width = length((IDs))*2
  height = 10
  
  plot <- ggplot( mat, aes(x = X,y = Y)) + 
                  geom_tile(aes(fill = FILL, height = HEIGHT), size = 0) + 
                  facet_grid(	. ~ GRID.Y ,
                              scales = "free",
                              space = "free"
                  ) +
                  theme(	axis.ticks = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         legend.position = "bottom"
                  ) +
                  scale_fill_gradient2(	low = color_low,
                                        mid = color_mid,
                                        high = color_high,
                                        midpoint = 0.5, 
                                        name = legend_label
                  ) +
                  scale_y_continuous(	breaks = c(),
                                      expand = c(0, 0),
                                      name = y_axis_label
                  ) +
                  scale_x_discrete(	position = "top",
                                    name = x_axis_label
                  )
    ggsave(filename = out.file, plot = plot, width = width, height = height, device = "pdf", units = 'in')
}


##################################################################
# function: get summarized counts in defined bins
##################################################################

get_bamOverlaps <- function(files, IDs, txdb, singleEnd){
  
  # gives information about path of bam and bai file (Object)
  bf = BamFileList(unlist(files))
  
  # gives length of different chromosomes (SeqInfo class)
  si = seqinfo(bf)
  
  # get exons and genes
  expr.gr <- genes(txdb)
  seqlevels(expr.gr) = paste0("chr", gsub("chr","",seqlevels(expr.gr)))
  
  exons <- exonsBy(txdb, by = "gene")
  exons <- exons[which((names(exons) %in% mcols(expr.gr)$gene_id) == T)]
  
  # fix chromosome prefix:
  seqlevels(exons) <- gsub("chr","",seqlevels(exons))
  if (grepl("chr",seqlevels(si)[1])) {
    seqlevels(exons) <- paste0("chr", seqlevels(exons))
  }
  
  se <- summarizeOverlaps(exons,
                          bf,
                          singleEnd = singleEnd,
                          fragments = setdiff(c(FALSE,TRUE), singleEnd))
  
  # get counts and summarize per gene:
  counts.per.exon <- data.frame(assays(se)$counts,
                                gene = rownames(assays(se)$counts))
  counts.split <- split(counts.per.exon,
                        counts.per.exon$gene)
  counts.per.gene <- (do.call(rbind,lapply(counts.split, function(x) colSums(x[,-dim(x)[2]]))))

  # create new symmarized experiment object:
  se0 <- SummarizedExperiment( assays = SimpleList(counts = counts.per.gene),
                               colData = names(bf)
  )
  
  # stabilize the variance across the mean:
  # (The transformed data should be approximated variance stabilized
  # and also includes correction for size factors or normalization factors)
  vsd <- varianceStabilizingTransformation(suppressMessages(DESeqDataSet(se0, ~ 1)), blind = FALSE)
  counts <- assay(vsd)

  for (i in 1:dim(counts)[2]) {
        mcols(expr.gr)[,unlist(IDs)[i]] <- counts[,i]
  }
  
  expr.gr <- expr.gr[which(rowVars(counts) > 0)]
  
  return(expr.gr)
}


##################################################################
# function: create TAD if not defined
##################################################################

check_TAD <- function(t, TAD, this.region, regions){
  
  if (length(t) == 0) {
    
    precedeT <- precede(TAD, this.region)
    followT  <- follow(TAD, this.region)
    
    start <- end(TAD[rev(which(!is.na(precedeT)))[1]])
    if (length(start) == 0) start <- 1
    
    end <- start(TAD[(which(!is.na(followT)))[1]])
    if (length(end) == 0) end <- max(end(region[which(seqnames(regions) == seqnames(this.region))]))
    
    t <- GRanges(seqnames(this.region), IRanges(start, width = (end - start + 1)))
  }
  return(t)
}

##################################################################
# function: correlate probabilities and gene exression counts 
##################################################################

get_correlation <- function(i, threshold, regions.gr, TAD.gr, IDs){
  
  interactions <- data.frame(stringsAsFactors = F)
  
  # get region, associated TAD and genes within TAD
  this.region <- regions.gr[i]
  this.TAD <- subsetByOverlaps(TAD.gr, this.region)
  this.TAD <- check_TAD(this.TAD, TAD.gr, this.region, regions.gr)
  this.genes.idx <- expr.gr %within% this.TAD
  
  if (sum(this.genes.idx) > 0) {
    
    this.genes <- expr.gr[this.genes.idx,]
    cor <- apply(as.matrix(mcols(this.genes)[,IDs]), 1, 
                 function(x) cor(x, as.numeric(unlist(mcols(this.region)[,IDs]))))
    
    for (c in 1:length(cor)) {
      if (!is.na(cor[c]) && cor[c] >= threshold) {
        interactions <- rbind(	interactions, 
                               data.frame(	data.frame(this.region)[,c(GR_header_short, "cluster", IDs)],
                                           TAD_COORDINATES = paste0(this.TAD),
                                           CORRELATED_GENE = paste(mcols(this.genes)[c,1]),
                                           CORRELATION = cor[c] ))
      }
    }
    if (length(interactions) > 0) return(makeGRangesFromDataFrame(interactions, keep.extra.columns = T))
  }
}

##################################################################
# function: correlate probabilities with gene expression values
#(in same TAD; per cluster)
##################################################################

get_units <- function(regions.gr, expr.gr, TAD.gr, IDs, cores, threshold){
  
  # get correlation for each differential region:
  list <- mclapply( seq(length(regions.gr)), 
                    function(x) get_correlation(x, threshold, regions.gr, TAD.gr, IDs), 
                    mc.cores = cores
  )
  units <-  do.call("c", unname(unlist(list)))
  
  return(units)
}


