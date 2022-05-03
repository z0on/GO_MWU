
# GO_MWU uses continuous measure of significance (such as fold-change or
# -log(p-value) ) to identify GO categories that are significantly enriches
# with either up- or down-regulated genes. The advantage - no need to impose
# arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO
# enrichment" analysis based Fisher's exact test: it will show GO categories
# over-represented among the genes that have 1 as their measure.

# On the plot, different fonts are used to indicate significance and color
# indicates enrichment with either up (red) or down (blue) regulated genes.
# No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on
# shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes
# in it; "good" genes being the ones exceeding the arbitrary absValue cutoff
# (option in gomwuPlot). For Fisher's based test, specify absValue = 0.5.
# This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu


# setup R error handling to go to stderr
options(show.error.messages = F, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, F)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
library("tools")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    "help", "h", 0, "logical",
    "scriptdir", "s", 1, "character",
    "input", "i", 1, "character",
    "goAnnotations", "a", 1, "character",
    "goDatabase", "g", 1, "character",
    "goDivision", "d", 1, "character",
    "largest", "o", 1, "numeric",
    "smallest", "m", 1, "numeric",
    "clusterheight", "c", 1, "numeric",
    "textsize", "e", 1, "numeric",
    "pcut", "p", 1, "numeric",
    "hcut", "t", 1, "numeric",
    "l1", "1", 1, "numeric",
    "l2", "2", 1, "numeric",
    "l3", "3", 1, "numeric",
    "version", "v", 0, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

if (!is.null(opt$version)) {
    cat("0.3.0\n")
    q(status = 1)
}

# enforce the following required arguments
if (is.null(opt$scriptdir)) {
    cat("'scriptdir' is required\n")
    q(status = 1)
}
if (is.null(opt$input)) {
    cat("'input' is required\n")
    q(status = 1)
}
if (is.null(opt$goDatabase)) {
    cat("'goDatabase' is required\n")
    q(status = 1)
}
if (is.null(opt$goAnnotations)) {
    cat("'goAnnotations' is required\n")
    q(status = 1)
}
if (is.null(opt$goDivision)) {
    cat("'goDivision' is required\n")
    q(status = 1)
}
if (is.null(opt$clusterheight)) {
    opt$clusterheight <- 0.25
}
if (is.null(opt$textsize)) {
    opt$textsize <- 1.2
}
if (is.null(opt$largest)) {
    opt$largest <- 0.1
}
if (is.null(opt$smallest)) {
    opt$smallest <- 5
}
if (is.null(opt$l1)) {
    opt$l1 <- 0.1
}
if (is.null(opt$l2)) {
    opt$l2 <- 0.05
}
if (is.null(opt$l3)) {
    opt$l3 <- 0.01
}
if (is.null(opt$pcut)) {
    opt$pcut <- 1e-2
}
if (is.null(opt$hcut)) {
    opt$hcut <- 0.9
}

# for testing
# #opt = list()
# opt$scriptdir = "./"
# opt$input = "toto/heats.csv"
# opt$goAnnotations = "toto/amil_defog_iso2go.tab"
# opt$goDatabase = "go.obo"
# opt$goDivision = "BP"

source_local <- function(fname) {
    # argv <- commandArgs(trailingOnly = FALSE)
    # base_dir <- dirname(substring(argv[grep("--file = ", argv)], 8))
    base_dir <- opt$scriptdir
    source(paste(base_dir, fname, sep = "/"))
}

source_local("gomwu.functions.R")


nn <- strsplit(opt$input, "[/.]")
if (length(nn[[1]]) == 3) {
    dir <- nn[[1]][1]
    name <- nn[[1]][2]
    ext <- nn[[1]][3]
} else if (length(nn[[1]]) == 2) {
    dir <- "."
    name <- nn[[1]][1]
    ext <- nn[[1]][2]
}
# It might take a few minutes for MF and BP. Do not rerun it if you just want
# to replot the data with different cutoffs, go straight to gomwuPlot. If you
# change any of the numeric values below, delete the files that were generated
# in previos runs first.

gomwuStats(opt$input, opt$goDatabase, opt$goAnnotations, opt$goDivision, opt$scriptdir,
    # replace with full path to perl executable if it is not in your system's PATH already
    perlPath = "perl",
    # a GO category will not be considered if it contains more than this fraction of the total number of genes
    largest = opt$largest,
    smallest = opt$smallest, # a GO category should contain at least this many genes to be considered.
    clusterCutHeight = opt$clusterheight, # threshold for merging similar (gene-sharing) terms. See README for details.
    # Alternative = "g" # by default the MWU test is two-tailed;
    # specify "g" or "l" of you want to test for "greater" or "less" instead.
    # Module = TRUE,Alternative="g" # un-remark this if you are analyzing a
    # SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes).
    # In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
    # Module = TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# ----------- Plotting results

# change this to a pdf output
pdf(paste0(dir, "/", "Rplots.pdf"), width = 7, height = 7)
# png(paste0(dir,"/","Rplots.png"),res=100)

results <- gomwuPlot(opt$input, opt$goAnnotations, opt$goDivision,
    # genes with the measure value exceeding this will be counted as "good genes".
    # This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing
    # Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
    absValue = -log(0.05, 10),
    # 	absValue = 1, # un-remark this if you are using log2-fold changes
    # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
    level1 = opt$l1,
    level2 = opt$l2, # FDR cutoff to print in regular (not italic) font.
    level3 = opt$l3, # FDR cutoff to print in large bold font.
    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text)
    # for better "word cloud" effect
    txtsize = opt$textsize,
    treeHeight = 0.5, # height of the hierarchical clustering tree
    # 	colors = c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    # these are default colors, uncomment and change if needed
)
dev.off()
# manually rescale the plot so the tree matches the text
# if there are too many categories displayed, try make it more stringent with level1 = 0.05,level2=0.01,level3=0.001.

# text representation of results, with actual adjusted p-values
write.table(results[[1]], paste0(dir, "/", "results.tsv"), sep = "\t")


# ------- extracting representative GOs

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut <- opt$pcut
hcut <- opt$hcut

# plotting the GO tree with the cut level (uncomment the next two lines to plot)
# plot(results[[2]],cex = 0.6)
# abline(h = hcut,col="red")

# cutting
ct <- cutree(results[[2]], h = hcut)
annots <- c()
ci <- 1
for (ci in unique(ct)) {
    message(ci)
    rn <- names(ct)[ct == ci]
    obs <- grep("obsolete", rn)
    if (length(obs) > 0) {
        rn <- rn[-obs]
    }
    if (length(rn) == 0) {
        next
    }
    rr <- results[[1]][rn, ]
    bestrr <- rr[which(rr$pval == min(rr$pval)), ]
    best <- 1
    if (nrow(bestrr) > 1) {
        nns <- sub(" .+", "", row.names(bestrr))
        fr <- c()
        for (i in 1:length(nns)) {
            fr <- c(fr, eval(parse(text = nns[i])))
        }
        best <- which(fr == max(fr))
    }
    if (bestrr$pval[best] <= pcut) {
        annots <- c(annots, sub("\\d+\\/\\d+ ", "", row.names(bestrr)[best]))
    }
}

mwus <- read.table(paste0(dir, "/", paste("MWU", opt$goDivision, name, sep = "_"), ".", ext), header = T)
best_GOs <- mwus[mwus$name %in% annots, ]
write.table(best_GOs, paste0(dir, "/", "best_go.tsv"), sep = "\t", row.names = FALSE)
