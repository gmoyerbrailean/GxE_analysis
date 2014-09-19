#######
## Script to compare FC of DE genes between shallow and deep
## input
##	plate: plate number, e.g., P5
##	fdrThresh: FDR cutoff for significance, e.g., 0.2
## G Moyerbrailean, Wayne State University
#######

library(ggplot2)
require(parallel)

cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
platePrefix<-cargs[1];
if(length(cargs)>=2)
fdrThresh<-cargs[2]
if(length(cargs)>=3)
baseDir<-cargs[3]
if(length(cargs)>=4)
cores<-cargs[4]

fdrThresh <- as.numeric(fdrThresh)
cores <- as.numeric(cores)

if(cores<1){cores <- 1}
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## Get treatment IDs and treatment names
# txTable <- read.table('etc/txID_to_name.txt', as.is=T, header=T, sep='\t')
# rownames(txTable) <- txTable$Treatment_ID
source('/nfs/rprscratch/gmb/GxE/differential_expression/DEseq2/analysis/src/load_cv_tables.R')

# Data Structure - 
# plots/
#	plate/
#		Tx1.pdf
#		Tx2.pdf
# system(paste0('mkdir -p plots/fold_change/', plate))

## Get the fold change data for treatments sequenced twice
deepFiles = list.files(paste0(baseDir, "out_data_D", platePrefix, "/stats/"))
ParallelSapply(deepFiles, function(dFile){

	txID = gsub('.txt', '', gsub('.*_', '', dFile))
	sFile=gsub('DP', 'P', dFile)
	cmd = paste0("less ", baseDir, "out_data_", platePrefix, "/stats/", sFile)
	shal <- read.table(pipe(cmd), as.is=T, sep=' ', header=T)
	cmd = paste0("less ", baseDir, "out_data_D", platePrefix, "/stats/", dFile)
	deep <- read.table(pipe(cmd), as.is=T, sep=' ', header=T)
	names(shal) <- c("t.id", "qv", "pval", "logFC", "ensg", "g.id")
	names(deep) <- names(shal)

	## Get the transcripts in both data sets
	both <- shal$t.id[shal$t.id %in% deep$t.id]
	shal <- shal[shal$t.id %in% both,]
	deep <- deep[deep$t.id %in% both,]

	## sanity check - ensure sorted/matching
	shal <- shal[order(shal$t.id),]
	deep <- deep[order(deep$t.id),]
	stopifnot(shal$t.id == deep$t.id)

	## Combine and format the data for easy plotting
	dat <- data.frame(t.id=shal$t.id, ensg=shal$ensg, g.id=shal$g.id, 
		sFC=shal$logFC, sQ=shal$qv, sP=shal$pval, dFC=deep$logFC, 
		dQ=deep$qv, dP=deep$pval)
	
	## Check the correlation between shallow and deep data
	corr <- cor.test(dat$sFC, dat$dFC, method="spearman")
	# cat('#', 'cor:', corr$estimate, '\n')
	# cat('#', 'pval:', corr$p.value, '\n')

	## Assign color by significance threshold
	dat$sig <- "All"
	dat$sig[dat$sQ<fdrThresh] <- "Shallow"
	dat$sig[dat$dQ<fdrThresh] <- "Deep"
	dat$sig[dat$sQ<fdrThresh & dat$dQ<fdrThresh] <- "Both"

	## Plot!
	# png(paste0('plots/fold_change/', platePrefix, '_', txID, '.png'))
	pdf(paste0('plots/fold_change/', platePrefix, '_', txID, '.pdf'))
	plot(dat$sFC, dat$dFC, pch=20, col="black", xlab="logFC - Shallow", ylab="logFC - Deep",
		main=paste0("Plates ", platePrefix, "/D", platePrefix, ": ", treatmentKey[txID,]$Short_Name, 
		 	    " (", fdrThresh*100, "% FDR)"))
	points(dat$sFC[dat$sig=='Deep'], dat$dFC[dat$sig=='Deep'], col='red', pch=20)
	points(dat$sFC[dat$sig=='Shallow'], dat$dFC[dat$sig=='Shallow'], col='blue', pch=20)
	points(dat$sFC[dat$sig=='Both'], dat$dFC[dat$sig=='Both'], col='purple', pch=20)
	legend("topleft", 
		legend=c(paste0("All (", dim(dat)[1],")"),
				 paste0("Both (", dim(dat[dat$sig=="Both",])[1], ")"),
				 paste0("Deep only (", dim(dat[dat$sig=="Deep",])[1], ")"),
				 paste0("Shallow only (", dim(dat[dat$sig=="Shallow",])[1], ")")), 
		col=c("black", "purple", "red", "blue"), pch=20, 
		title=paste0("Spearman: ", round(corr$estimate,2), ", p=", corr$p.value))
	abline(v=0, lty=3)
	abline(h=0, lty=3)
	dev.off()

	# # Try with just the significantly DE transcripts
	# dat2 <- dat[dat$sig != "All",]
	# p1 <- ggplot(dat2, aes(x=sFC, y=dFC)) 
	# p1 + 
	# 	geom_point(alpha=0.75, aes(color=sig)) +
	# 	scale_colour_manual(values = c("purple", "red", "blue")) +
	# 	xlab("logFC - Shallow") +
	# 	ylab("logFC - Deep")
	# ggsave(paste0('plots/fold_change/', platePrefix, '_', txID, '_sig.png'))
})
