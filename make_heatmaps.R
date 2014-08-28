#######
## Script to make a heatmap out of a correlation matrix
## input
##      baseDir: Analysis directory, e.g, contains 'data', 'plots', 'src', etc
##      dataMat: Matrix of input data to generate correlation matrix
## output
##	Heatmap plots and correlation matrices saved as Rd files
## G Moyerbrailean, Wayne State University
#######

baseDir="/nfs/rprscratch/gmb/GxE/differential_expression/DEseq2/analysis/"
dataMat="lfcTable.clean.Rd"
cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
baseDir<-cargs[1];
if(length(cargs)>=2)
dataMat<-cargs[2];

library(gplots)

## Make a column-wise correlation matrix given an input matrix
makeCorrMatrix <- function(mtx){
	if (class(mat)=="matrix") {
		goodRows <- !is.na(rowSums(mat))
		mat.clean <- mat[goodRows,]
		corr <- sapply(colnames(mat.clean), function(p1){
			sapply(colnames(mat.clean), function(p2){
				c <- cor.test(mat.clean[,p1], mat.clean[,p2])
				c$estimate
			})
		})
		return(corr)
	} else {
		return(0)
	}
}

## Generate the lookup tables
## TODO - get/lookup plate #'s, don't hard-code
plateLst <- c(1,2,4,5,6,7,8,9,11,12)
cellTypeLst <- c("LCL", "PBMC", "HUVEC", "PBMC", "PBMC",  "LCL", "HUVEC", "SMC", "Mel", "Mel")
txTable <- read.table('etc/txID_to_name.txt', as.is=T, header=T, sep='\t')
rownames(txTable) <- txTable$Treatment_ID
plateTable <- data.frame(plateID=c(paste0("P", plateLst), 
	paste0("DP", plateLst)), cellType=rep(cellTypeLst, 2))
rownames(plateTable) <- plateTable$plateID

## Use lookup tables to convert ID-based labels into readable labels
gen_labels <- function(names) {
	sapply(names, function(c){
		a <- unlist(strsplit(c, "_"))
		cell <- plateTable[a[1],]$cellType
		tx <- txTable[a[2],]$Short_Name
		if (substr(a[1], 1, 1)=="D") {
			return(paste0(cell, " - ", tx, " (deep)"))
		} else{
			return(paste0(cell, " - ", tx))
		}
	})
}

#### Start
## Read in the data matrix
load(paste0(baseDir, '/data/', dataMat)) # in this case object is "lfcTable.clean"

## Make a summary table for each treatment with >= 1 deep sequence run
ind <- grep('DP', colnames(lfcTable.clean))
cols <- c(colnames(lfcTable.clean)[ind], gsub('DP', 'P', colnames(lfcTable.clean)[ind]))
cols <- cols %in% colnames(lfcTable.clean)
mat <- lfcTable.clean[,cols]
summaryMtx <- makeCorrMatrix(mat)
labs <- gen_labels(colnames(summaryMtx))
colnames(summaryMtx) <- rownames(summaryMtx) <- labs
save(summaryMtx, file=paste0(baseDir, '/data/corr.matrix.Rd'), compress=T)

system(paste0('mkdir -p ', baseDir, '/plots/heatmaps/'))
pdf(paste0(baseDir, '/plots/heatmaps/all.heatmap.2.pdf'), height=10, width=10)

## Color the columns by treatment and the rows by cell type
plate <- factor(sapply(strsplit(names(colnames(summaryMtx)),"_"),function(x){x[1]}))
treatment <- factor(sapply(strsplit(names(colnames(summaryMtx)),"_"),function(x){x[2]}))
cellType <- factor(plateTable[plate,]$cellType)
cc <- rainbow(length(levels(treatment)),start=0,end=1.0)
rc <- terrain.colors(length(levels(cellType)))

heatmap.2(summaryMtx, trace="none", cexCol=0.2, cexRow=0.2, margins=c(10,10), 
	main="All treatments", ColSideColors=cc[treatment], 
	RowSideColors=rc[cellType])
dev.off()

pdf(paste0(baseDir, '/plots/heatmaps/all.heatmap.legend.pdf'), height=13, width=10)
plot.new()
legend("topleft",      # location of the legend on the heatmap plot
    legend = levels(cellType), # category labels
    col = rc,  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
legend("topright",      # location of the legend on the heatmap plot
    legend = txTable[levels(treatment),]$Treatment_Name, # category labels
    col = cc,  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()


## Make a summary for all the deep plates
ind <- grep('DP', colnames(lfcTable.clean))
mat <- lfcTable.clean[,ind]
deepMtx <- makeCorrMatrix(mat)
labs <- gen_labels(colnames(deepMtx))
colnames(deepMtx) <- rownames(deepMtx) <- labs
save(deepMtx, file=paste0(baseDir, '/data/corr.matrix.deep.Rd'), compress=T)

system(paste0('mkdir -p ', baseDir, '/plots/heatmaps/'))

plate <- factor(sapply(strsplit(names(colnames(deepMtx)),"_"),function(x){x[1]}))
treatment <- factor(sapply(strsplit(names(colnames(deepMtx)),"_"),function(x){x[2]}))
cellType <- factor(plateTable[plate,]$cellType)
cc <- rainbow(length(levels(treatment)),start=0,end=1.0)
rc <- terrain.colors(length(levels(cellType)))

pdf(paste0(baseDir, '/plots/heatmaps/deep.heatmap.2.pdf'), height=10, width=10)
heatmap.2(deepMtx, trace="none", cexCol=1, cexRow=1, margins=c(11,11), 
	main="All deep sequencing", ColSideColors=cc[treatment], 
	RowSideColors=rc[cellType])
dev.off()

pdf(paste0(baseDir, '/plots/heatmaps/deep.heatmap.legend.pdf'), height=13, width=10)
plot.new()
legend("topleft",      # location of the legend on the heatmap plot
    legend = levels(cellType), # category labels
    col = rc,  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
legend("topright",      # location of the legend on the heatmap plot
    legend = txTable[levels(treatment),]$Treatment_Name, # category labels
    col = cc,  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()



## Make a summary for all the shallow plates
ind <- grep('DP', colnames(lfcTable.clean))
mat <- lfcTable.clean[,-ind]
shallowMtx <- makeCorrMatrix(mat)
labs <- gen_labels(colnames(shallowMtx))
colnames(shallowMtx) <- rownames(shallowMtx) <- labs
save(shallowMtx, file=paste0(baseDir, '/data/corr.matrix.shallow.Rd'), compress=T)

system(paste0('mkdir -p ', baseDir, '/plots/heatmaps/'))
pdf(paste0(baseDir, '/plots/heatmaps/shallow.heatmap.2.pdf'), height=10, width=10)

plate <- factor(sapply(strsplit(names(colnames(shallowMtx)),"_"),function(x){x[1]}))
treatment <- factor(sapply(strsplit(names(colnames(shallowMtx)),"_"),function(x){x[2]}))
cellType <- factor(plateTable[plate,]$cellType)
cc <- rainbow(length(levels(treatment)),start=0,end=1.0)
rc <- terrain.colors(length(levels(cellType)))

heatmap.2(shallowMtx, trace="none", cexCol=0.3, cexRow=0.3, margins=c(10,10), 
	main="All shallow sequencing", ColSideColors=cc[treatment], 
	RowSideColors=rc[cellType])
dev.off()

pdf(paste0(baseDir, '/plots/heatmaps/shallow.heatmap.legend.pdf'), height=13, width=10)
plot.new()
legend("topleft",      # location of the legend on the heatmap plot
    legend = levels(cellType), # category labels
    col = rc,  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
legend("topright",      # location of the legend on the heatmap plot
    legend = txTable[levels(treatment),]$Treatment_Name, # category labels
    col = cc,  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()
