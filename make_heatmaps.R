
dataDir="/wsu/home/groups/piquelab/charvey/GxE/differential_expression/DEseq2_results/"
baseDir="/nfs/rprscratch/gmb/GxE/differential_expression/DEseq2/analysis/"
cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
dataDir<-cargs[1];
if(length(cargs)>=2)
baseDir<-cargs[2];

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
txTable <- read.table('etc/txID_to_name.txt', as.is=T, header=T, sep='\t')
rownames(txTable) <- txTable$Treatment_ID
plateTable <- data.frame(plateID=c(paste0("P", c(1,2,4,5,7,8,11,12)), 
	paste0("DP", c(1,2,4,5,7,8,11,12))), cellType=rep(c("LCL", "PBMC", "HUVEC", 
		"PBMC", "LCL", "HUVEC", "Mel", "Mel"),2))
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
load(paste0(baseDir, '/data/lfcTable.clean.Rd')) # object is "lfcTable.clean"

## Make a summary table for each treatment with >= 1 deep sequence run
ind <- grep('DP', colnames(lfcTable.clean))
cols <- c(colnames(lfcTable.clean)[ind], gsub('DP', 'P', colnames(lfcTable.clean)[ind]))
mat <- lfcTable.clean[,cols]
summaryMtx <- makeCorrMatrix(mat)
labs <- gen_labels(colnames(summaryMtx))
colnames(summaryMtx) <- rownames(summaryMtx) <- labs
save(summaryMtx, file=paste0(baseDir, '/data/corr.matrix.Rd'), compress=T)

system(paste0('mkdir -p ', baseDir, '/plots/heatmaps/'))
pdf(paste0(baseDir, '/plots/heatmaps/all.heatmap.pdf'), height=10, width=10)
heatmap.2(summaryMtx, trace="none", cexCol=0.8, cexRow=0.8, margins=c(10,10), 
	main="All treatments w/ >= 1 deep sequencing run")
dev.off()


## Make a summary for all the deep plates
ind <- grep('DP', colnames(lfcTable.clean))
mat <- lfcTable.clean[,ind]
deepMtx <- makeCorrMatrix(mat)
labs <- gen_labels(colnames(deepMtx))
colnames(deepMtx) <- rownames(deepMtx) <- labs
save(deepMtx, file=paste0(baseDir, '/data/corr.matrix.deep.Rd'), compress=T)

system(paste0('mkdir -p ', baseDir, '/plots/heatmaps/'))
pdf(paste0(baseDir, '/plots/heatmaps/deep.heatmap.pdf'), height=10, width=10)
heatmap.2(deepMtx, trace="none", cexCol=1, cexRow=1, margins=c(11,11), 
	main="All deep sequencing")
dev.off()


## Make a summary for all the shallow plates
ind <- grep('DP', colnames(lfcTable.clean))
mat <- lfcTable.clean[,-ind]
shallowMtx <- makeCorrMatrix(mat)
labs <- gen_labels(colnames(shallowMtx))
colnames(shallowMtx) <- rownames(shallowMtx) <- labs
save(shallowMtx, file=paste0(baseDir, '/data/corr.matrix.shallow.Rd'), compress=T)

system(paste0('mkdir -p ', baseDir, '/plots/heatmaps/'))
pdf(paste0(baseDir, '/plots/heatmaps/shallow.heatmap.pdf'), height=10, width=10)
heatmap.2(shallowMtx, trace="none", cexCol=0.5, cexRow=0.5, margins=c(10,10), 
	main="All shallow sequencing")
dev.off()