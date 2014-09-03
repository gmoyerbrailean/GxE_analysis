#######
## Script to generate a matrix of transcript log fold change between Tx and control
## input
##      dataDir: Directory where DEseq2 data is stored with stats/ directory
##      baseDir: Analysis directory, e.g, contains 'data', 'plots', 'src', etc
## output
##		Matrix of plate/Tx combos by transcript IDs, with logFC as values
##		Clean matrix with all-NA rows removed
## G Moyerbrailean, Wayne State University 2014
#######

## Script to generate a matrix of transcript log fold change between Tx and control
## G Moyerbrailean, Wayne State University 2014

dataDir="/wsu/home/groups/piquelab/charvey/GxE/differential_expression/DEseq2_results/"
baseDir="/nfs/rprscratch/gmb/GxE/differential_expression/DEseq2/analysis/"
cargsdataDir="/wsu/home/groups/piquelab/charvey/GxE/differential_expression/DEseq2_results/"
baseDir="/nfs/rprscratch/gmb/GxE/differential_expression/DEseq2/analysis/"
<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
dataDir<-cargs[1];
if(length(cargs)>=2)
baseDir<-cargs[2];

## Get treatment IDs and treatment names
txTable <- read.table('etc/txID_to_name.txt', as.is=T, header=T, sep='\t')
rownames(txTable) <- txTable$Treatment_ID

## Generate the lookup tables
## TODO - get/lookup plate #'s, don't hard-code
plateList <- c(1,2,4,5,6,7,8,9,11,12)
cellTypeList <- c("LCL", "PBMC", "HUVEC", "PBMC", "PBMC",  "LCL", "HUVEC", "SMC", "Mel", "Mel")
txTable <- read.table('etc/txID_to_name.txt', as.is=T, header=T, sep='\t')
rownames(txTable) <- txTable$Treatment_ID
plateTable <- data.frame(plateID=c(paste0("P", plateList), 
	paste0("DP", plateList)), cellType=rep(cellTypeList, 2))
rownames(plateTable) <- plateTable$plateID

## TODO: make plateList user input
shalPlates <- paste0('P', plateList)
deepPlates <- paste0('DP', plateList)
allPlates <- c(shalPlates, deepPlates)

## Get master list of transcripts
cmd <- 'less ~/piquelab/data/RefTranscriptome/ensGene.hg19.v2.bed.gz | cut -f 4 | sort'
refTxome <- scan(pipe(cmd), what="")

## Find all the data for each plate-condition pair
## NOTE: this is important as some of the possible combinations may not exist
txFiles <- list.files(paste0(dataDir, 'out_data_', allPlates, '/stats/'), full.names=TRUE)
txIDs <- gsub('.*/', '', gsub('.txt', '', gsub('_DEG_stats', '', txFiles)))

## Initialize a table for logFC data
lfcTable <- matrix(data=NA, nrow=length(refTxome), ncol=length(txIDs))
rownames(lfcTable) <- refTxome
colnames(lfcTable) <- txIDs

for (p in 1:length(txFiles)){
	print(paste0('processing ', txIDs[p], '...'))
	aux <- read.table(txFiles[p], header=T, sep=' ', as.is=T)
	lfcTable[aux$t.id, txIDs[p]] <- aux$logFC
}
save(lfcTable, file=paste0(baseDir, '/data/lfcTable.Rd'), compress=T)

## Remove transcripts that are NA for all combinations
goodRows <- sapply(1:dim(lfcTable)[1], function(ii){
	sum(is.na(lfcTable[ii,])) != length(lfcTable[ii,])
})
lfcTable.clean <- lfcTable[goodRows,]
save(lfcTable.clean, file=paste0(baseDir, '/data/lfcTable.clean.Rd'), compress=T)
