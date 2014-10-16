
library(parallel)

## Get command-line arguments ##
# platePrefix="DP1"
# cores=6
cargs<-commandArgs(trail=T)
if (length(cargs)>0){
	platePrefix=cargs[1]}
if (length(cargs)>1){
	cores=cargs[2]}

if(cores<1){cores <- 1}
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}
LPG <- Sys.getenv("LPG")

## Get the covariate file for the indicated plate
## Loaded objects:
	## treatmentKey: mapping of treatment ID to treatment name
source(paste0(LPG, '/gmb/GxE/differential_expression/DEseq2/',
	'analysis/src/load_cv_tables.R'))

## Get the counts data for analysis
## Loaded objects:
	## data: read-counts per each gene transcript
	## anno: gene-transcript annotations
readCounts <- paste(LPG, '/scratch/charvey/GxE/derived_data/', platePrefix, 
	'/counts/GC/', platePrefix, '.data.Rd', sep='')
load(readCounts)
countsData <- data # read-counts per each gene transcript
transcriptAnno <- anno # gene-transcript annotations
rm(data); rm(anno)
n.barcodes <- dim(countsData)[2]
rownames(transcriptAnno) <- transcriptAnno$t.id

## Load the reference transcriptome
bedtranscript <- paste0(LPG, "/data/RefTranscriptome/ensGene.hg19.v2.bed.gz")

##################################################################
## Annotate transcripts with coding length, transcript length,
## and average GC content w/i coding region
gcContentFile <- paste0(LPG, "/data/RefTranscriptome/ensGene.hg19.faCount.gz")
anno2 <- read.table(gcContentFile,as.is=T, sep="\t", header=T, comment="")
rownames(anno2) <- gsub("hg19_ensGene_", "", anno2$X.seq)
anno2$avg.cg <- (anno2$C+anno2$G) / anno2$len
transcriptAnno$codLen <- anno2[transcriptAnno$t.id, "len"]
transcriptAnno$txLen <- ((transcriptAnno$c.start - transcriptAnno$start) + 
	(transcriptAnno$stop - transcriptAnno$c.stop) + transcriptAnno$codLen)
anno2$avg.cg <- (anno2$C + anno2$G) / anno2$len
transcriptAnno$avg.cg <- anno2[transcriptAnno$t.id, "avg.cg"]
rm(anno2)

## Manual conversion to R factor objects:
cv$Treatment.ID <- factor(cv$Treatment.ID)
TreatmentLevels <- levels(cv$Treatment.ID)
cv$Control.ID <- factor(cv$Control.ID)
ControlLevels <- levels(cv$Control.ID)
TreatmentOnlyLevels <- TreatmentLevels[!(TreatmentLevels %in% ControlLevels)]
ControlOnlyLevels <- TreatmentLevels[(TreatmentLevels %in% ControlLevels)]
cv$CellType <- factor(cv$CellType)
cv$CellLine <- factor(cv$CellLine)
CellLineLevels <- levels(cv$CellLine)

##  Collapsing technical replicates
allColSamples <- paste(cv$CellLine, cv$Treatment.ID, sep=".")
sp <- split( seq(along=allColSamples), allColSamples )
cdata <- ParallelSapply(sp, function(columns)
	round(rowSums( countsData[,columns,drop=FALSE] )))
cv2 <- cv[sapply(sp, `[` , 1),]


##################################################################
## CONVERSION  OF cdata TO RPKM
cs <- colSums(cdata)/1E6
cdata2 <- ParallelSapply(1:ncol(cdata), function(jj){
	cdata[,jj] / cs[jj] * 1000 / transcriptAnno$codLen
})
colnames(cdata2) <-  colnames(cdata)

## Remove very  lowly expressed transcripts
naz <- rowSums(cdata2 > 0.1)
indNaz50 <- naz > (0.9 * ncol(cdata2))
cdata2 <- cdata2[indNaz50,]
geneAnno2 <- transcriptAnno[indNaz50,]
colnames(cdata2) <- colnames(cdata2)

geneAnno2$utrLen <- geneAnno2$txLen - geneAnno2$codLen
clipT <- quantile(geneAnno2$utrLen,0.99)
geneAnno2$utrLen[geneAnno2$utrLen>clipT]<-clipT


##################################################################
## REGRESS OUT  GC EFFECT AND TRANSCRIPT LENTGHT FFECT
fit.lm <- lm(log10(cdata2[,1] + 0.001) ~ geneAnno2$avg.cg + geneAnno2$codLen + 
	geneAnno2$avg.cg * geneAnno2$codLen + geneAnno2$utrLen)
summary(fit.lm)

cdata3 <- sapply(1:ncol(cdata2),function(ii){
	fit.lm <- lm(log10(cdata2[,ii] + 0.001) ~ geneAnno2$avg.cg + 
		geneAnno2$codLen + geneAnno2$avg.cg * geneAnno2$codLen + geneAnno2$utrLen)
	x <- residuals(fit.lm)
})
colnames(cdata3) <- colnames(cdata2)
rownames(cdata3) <- rownames(cdata2)

## Save the baseline rpkm data
save(list=c('cdata3', 'cv2'), file=paste0("data/", platePrefix, ".rpkm.Rd"))


##################################################################
## Quantile normalization (by sample)
cdata4 <- apply(cdata3,2,function(x){
	qqnorm(rank(x, ties.method = "random"), plot = F)$x
})
colnames(cdata4) <- colnames(cdata3)
rownames(cdata4) <- rownames(cdata3)


##################################################################
## ADJUSTMENT
## Try two methods and compare

## Remove means
cdataAdj <- cdata4
for (cl in 1:length(unique(as.numeric(cv2$CellLine)))) {
	avg <- apply(cdataAdj[,as.numeric(cv2$CellLine)==cl], 1, mean)
	cdataAdj[,as.numeric(cv2$CellLine)==cl] <- cdataAdj[,as.numeric(cv2$CellLine)==cl] - avg
}
save(cdataAdj, file=paste0('data/', platePrefix, '.rpkm.adj.Rd'))

## Remove controls
cdataCntRm <- cdata4
drop = c()
for (cl in levels(cv2$CellLine)) {
	for (cnt in levels(cv2$Control.ID)) {
		cdataCntRm[,cv2$CellLine==cl & cv2$Control.ID==cnt] <- 
			(cdata4[,cv2$CellLine==cl & cv2$Control.ID==cnt] - 
					cdata4[,cv2$CellLine==cl & cv2$Treatment.ID==cnt])
		drop <- c(drop, which(cv2$CellLine==cl & cv2$Treatment.ID==cnt))
	}
}
# Drop the "empty" control column
cdataAdj <- cdataCntRm[,-drop] # cdataCntRm <- cdataCntRm[,-drop]
save(cdataAdj, file=paste0('data/', platePrefix, '.rpkm.cntRm.Rd'))

## THE END