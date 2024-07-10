# 0. load package
#BiocManager::install("Rtpca")
#BiocManager::install("TPP")
library(Rtpca)
library(TPP)

# 1. Import data and set dataset for analysis
load(neoantigen_objects.RData)

#load peptide data
expression_data <- as.data.frame(exprs(peptides_novel_normalized) )#contains quantitive data. relative abundance
phenotype_data  <- pData(peptides_novel_normalized)  #contains temperature
feature_data<- featureData(peptides_novel_normalized)@data #contains Novel (T/F)
#expression_data_t_pep<-as.data.frame(t(expression_data_pep))
#add peptide no as a row in expression data
expression_data$peptide_no<-rownames(expression_data)

##set config table
config_data <- phenotype_data[, c("channel", "temperature")]
# Select columns, remove row names, and drop duplicates
config_data <- unique(phenotype_data[, c("channel", "temperature")])
rownames(config_data) <- NULL
head(config_data)
# Transform config_data to use 'channel' as column names and 'temperature' as values
config_data <- config_data %>%
  distinct() %>%  # Remove any duplicates to ensure unique 'channel' values
  spread(key = channel, value = temperature)
config_data <- config_data %>%
  mutate(Experiment = "experimentData")
# Move 'Experiment' to the first column 
config_data <- config_data[c("Experiment", setdiff(names(config_data), "Experiment"))]

set_numbers <- (paste0("Set", 1:11, "_tmt16plex_"))

# import data (peptides)
trData <- tpptrImport(configTable = config_data, data = expression_data, idVar ='peptide_no',fcStr = set_numbers)
head(trData)#eset data

# 2. Run TPCA (single condition, PPI/co-agg)

##  TPCA for co-aggregation

#load pf annotation data
pf_anno<-read.delim('ProteoformAnnotation.txt')
head(pf_anno)
#remove column
pf_anno <- pf_anno %>%
  select(-peptide, Peptide_AA, peptide_no, proteoform_id)
head(pf_anno)

#change colnames
colnames(pf_anno)<-c('ensembl_id', 'protein','id')
head(pf_anno)

#run TPCA for co-aggregation (full)
vehComplexTPCA <- runTPCA(
  objList = trData,
  complexAnno = pf_anno,
  minCount = 2
)

#save result
saveRDS(vehComplexTPCA,'vehComplexTPCA.RDS')

# 3. Visualize results (PPI)

#a ROC curve for how well our data captures protein complexes:
plotComplexRoc(vehComplexTPCA_mini, computeAUC = TRUE)

# 4.  inspect significantly co-melting proteofroms
head(tpcaResultTable(vehComplexTPCA_mini) )
tpcaResultTable(vehComplexTPCA_mini) %>% filter(p_adj < 0.1)
  

