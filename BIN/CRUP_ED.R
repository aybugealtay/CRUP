##################################################################
# start run time
##################################################################

start_time <- Sys.time()

##################################################################
# get input paramter
##################################################################

if (!exists('opt')) {
  
  # path to functions
  path <- dirname(getwd())
  source(paste0(path,'/BIN/functions.R'))
  
  # library - getopt
  pkgTest("getopt")
  suppressMessages(library(getopt))
  
  # get input parameter
  opt <- getopt(spectrum)
}

##################################################################
# define fixed parameters;
##################################################################

# colors for heatmap spectrum:
color_low  <- "#D0E5F5"
color_mid  <- "#5A92BD"
color_high <- "#003964"

# axis and legend labels:
legend_label <- "Enhancer Probabilities"
y_axis_label <- "Differential Enhancer"
x_axis_label <- ""

# prefixes:
condition_prefix  <- "Condition "
ID_prefix         <- "cond"
sample_prefix     <- "Rep"

##################################################################
# check input parameter
##################################################################

parameter <- c("probabilities")
parameter.bool <- sapply(parameter, function(x) !is.null(opt[[x]]))
parameter.extra <- c("outdir", "cores", "w_0", "len", "threshold", "names")

if (!is.null(opt$help)) stop_script(c(parameter, parameter.extra))
if (sum(parameter.bool) == 0) stop_script(c(parameter, parameter.extra))

# check probability files
files <- as.list(unlist(strsplit(opt$probabilities,',')))
files <- lapply(files, function(x) unlist(strsplit(x,':')))
for (f in unlist(files)) check_file(f)

# assign IDs: 
if (is.null(opt$names)) {
  IDs <- lapply(1:length(files), function(x) paste0(ID_prefix, x, "_",seq(1:length(files[[x]]))))
}else{
  labels <- unlist(strsplit(opt$names,','))

  if (length(labels) != length(files)) {
    cat("\n\tWARNING: Number of alternative condition names ('n') is not valid.
        Names are set to default labels.\n");
    IDs <- lapply(1:length(files), function(x) paste0(ID_prefix, x, "_", seq(1:length(files[[x]]))))
  } else{
    IDs <- files
    for (i in 1:length(IDs)) {
      for (j in 1:length(IDs[[1]])) {
        IDs[[i]][j] <- paste0(labels[i], "_", j)
      }
    }
  }
}

# check output directory
opt$outdir <- check_outdir(opt$outdir, files[1])

# set default values
if (is.null(opt$cores))       opt$cores <- 1
if (is.null(opt$len))         opt$len <- 1000
if (is.null(opt$threshold)) {
  opt$threshold <- 0.01 
} else if (!is.null(opt$threshold) & (opt$threshold < 0 | opt$threshold > 1) ) {
  cat(paste0(opt$threshold, " is not in range [0,1].\n"))
  q();
}
if (is.null(opt$w_0)) {
  opt$w_0 <- 0.5
} else if (!is.null(opt$w_0) & (opt$w_0 < 0 | opt$w_0 > 1) ) {
  cat(paste0(opt$w_0, " is not in range [0,1].\n"))
  q();
}

##################################################################
# define input parameter
##################################################################

startPart("List input parameter")

files     <- lapply(files, function(x) normalizePath(x))
outdir    <- paste0(normalizePath(opt$outdir),"/")
cores     <- opt$cores
w_0       <- opt$w_0
len       <- opt$len
threshold <- opt$threshold

cat(skip(), "files: \n", 
            paste0("\t\t",unlist(lapply(1:length(files), 
                                        function(x) paste(IDs[[x]], files[[x]], sep = "\t-> "))), "\n"))
cat(skip(), "w_0: ",w_0, "\n")
cat(skip(), "len: ",len, "\n")
cat(skip(), "threshold: ",threshold, "\n")
cat(skip(), "cores: ",cores, "\n")
cat(skip(), "outdir: ",outdir, "\n")

endPart()

##################################################################
# libraries
##################################################################

startPart("Load packages")

pkgLoad("ggplot2")        # for ggplot()
pkgLoad("GenomicRanges")  # for GRanges()
pkgLoad("dplyr")          # for %>%
pkgLoad("matrixStats")    # for rowVars()
pkgLoad("parallel")       # for mclapply()

endPart()

##################################################################
# read and combine probabilities:
##################################################################

startPart("Read enhancer probabilites for all samples and conditions")
probs <- get_probabiltyMatrix(files, IDs)
endPart()

##################################################################
# get p-values
# (empirical p-values per pair comparison with sign (+/-)
# -> sign gives direction of group mean difference)
##################################################################

startPart("Calculate (pairwise) empirical p-values")
p <- get_pairwisePvalues(probs, IDs, w_0, cores)
endPart()

##################################################################
# Call dynamic enhancer peaks
##################################################################

startPart("Get condition-specific enhancer peaks")

cat(paste0(skip(), "build significance peak pattern"))
pattern <- get_peakPattern(p, threshold, IDs, cores)
done()

if (is.null(pattern)) {
  cat(skip(), "no significant peaks found between any two conditions.\n")
  q()
}

cat(paste0(skip(), "combine peaks by significance pattern"))
peaks <- get_combinedDiffPeaks(probs, p, pattern, IDs, len)
colnames(elementMetadata(peaks)) <- c("best.p.value", "index", "significance.pattern", "cluster", "cluster.size")
elementMetadata(peaks)[, unlist(IDs)] <-  mcols(probs[mcols(peaks)$index])
done()

out.txt <- paste0(outdir, paste0("dynamicEnh__w0_",w_0,"__threshold_", threshold,".txt"))
cat(paste0(skip(), "save condition-specific enhancer regions to:  ", out.txt))
write.table(data.frame(peaks)[, c(GR_header_short, "best.p.value", "cluster", "significance.pattern", unlist(IDs))], file = out.txt, quote = F, row.names = F, sep = "\t")
done()

out.bed <- paste0(outdir, paste0("dynamicEnh__w0_",w_0,"__threshold_", threshold,".bed"))
cat(paste0(skip(), "a reduced version is saved in bed file format:  ", out.bed))
write.table(data.frame(peaks)[, GR_header_short], file = out.bed, quote = F, row.names = F, col.names = F, sep = "\t")
done()

# prepare matrix for heatmap:
LABEL_COND <- gsub("_.*","", IDs)
LABEL_REP <- gsub(".*_", "", IDs)
CLUSTER <- as.numeric(mcols(peaks)$cluster)
VALUES  <- mcols(probs[mcols(peaks)$index])

mat <- data.frame(  X = unlist(lapply(LABEL_REP, function(x) rep(x, length(peaks)))),
                    Y = rep(adjust(width(peaks)), length(unlist(files))), 
                    HEIGHT = rep(width(peaks), length(unlist(files))),
                    FILL = matrix(unlist(VALUES),  ncol = 1), 
                    GRID.X = factor(rep(CLUSTER, length(unlist(files))), levels = unique(CLUSTER)),
                    GRID.Y = unlist(lapply(LABEL_COND, function(x) rep(x,  length(peaks))))
)

# plot 
out.pdf <- paste0(outdir, paste0("dynamicEnh__w0_",w_0,"__threshold_", threshold,".pdf"))
cat(paste0(skip(), "results are visualized as a heatmap:  ", out.pdf))

plot_heatmap(	mat[order(mat$GRID.X),],
              IDs,
              color_low,
              color_mid,
              color_high,
              x_axis_label,
              y_axis_label,
              legend_label,
              out.pdf
)
done()

endPart()

##################################################################
# print run time
##################################################################

run_time <- Sys.time() - start_time

startPart("Run time")
cat(paste0(skip(), format(run_time), "\n"))
endPart()



