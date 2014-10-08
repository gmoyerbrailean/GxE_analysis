## Script to source in covariate table for a given plate, and the master treatment table
## Usage: source(/path/to/script/load_cv_table.R)
## Usage note: the environment you are sourcing into must have a "platePrefix" variable set, e.g., P1, DP2, etc
## G Moyerbrailean, Wayne State University 2014

LPG <- Sys.getenv("LPG")

## Load treatment and covariate data
treatmentKey <- read.table(paste(LPG, '/charvey/GxE/derived_data/covariates/GxE_treatment_key.txt', sep=''), sep='\t', as.is=TRUE, header=TRUE)
row.names(treatmentKey) <- treatmentKey$Treatment_ID
cov.file <- paste(LPG, '/charvey/GxE/derived_data/covariates/GxE_', platePrefix, '_covariates.txt', sep='')
cv <- read.table(cov.file , as.is=T, sep="\t", header=T, comment="")
cv <- cv[order(cv$Barcode.ID),]
rownames(cv) <- cv$Barcode.ID

cv$Treatment.ID <- factor(cv$Treatment.ID)
TreatmentLevels <- levels(cv$Treatment.ID)

cv$Control.ID <- factor(cv$Control.ID)
ControlLevels <- levels(cv$Control.ID)
SampleLevels <- TreatmentLevels[!(TreatmentLevels %in% ControlLevels)]

cv$CellLine <- factor(cv$CellLine)
CellLineLevels <- levels(cv$CellLine)

