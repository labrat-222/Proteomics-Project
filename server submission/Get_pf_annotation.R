library(igraph)

#import graph object
graphs_comm<- readRDS("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/Input/graphs_comms 1.RDS")

##########
#make minimal hits_everything by processing the graphs object
id_vector <- c()
gene_vector <- c()
for(i in names(graphs_comms)){
  id_vector <- c(id_vector, unique(paste0(graphs_comms[[i]]$ioi, "_", unique(vertex_attr(graphs_comms[[i]])[["membership"]]))))
  gene_vector <- c(gene_vector, rep(graphs_comms[[i]]$ioi, length.out = length(unique(vertex_attr(graphs_comms[[i]])[["membership"]]))))
  
}

create_hits_everything <- function(graphs_comms) {
  data_list <- lapply(names(graphs_comms), function(i) {
    graph <- graphs_comms[[i]]
    data.frame(
      name = vertex_attr(graph, "name"),
      gene_symbol = vertex_attr(graph, "id"),
      protein_ids = vertex_attr(graph, "protein_ids"),
      first_protein_id = vertex_attr(graph, "first_protein_id"),
      membership = vertex_attr(graph, "membership"),
      Peptide = vertex_attr(graph, "name"),  # Assuming 'name' contains the peptide sequence
      pAdj = "N.D",
      test = "meltome"
    )
  })
  
  hits_everything <- do.call(rbind, data_list)
  return(hits_everything)
}

# Create hits_everything dataframe from graphs_comms
hits_everything <- create_hits_everything(graphs_comm)
head(hits_everything)

#####################
#assign pf column
#drop columns with 
# Create a new column 'proteoform_id'
hits_everything$proteoform_id <- paste(hits_everything$gene_symbol, hits_everything$membership, sep = "_")

# Display the first few rows to verify
head(hits_everything)

write.table(hits_everything,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/outputs/hits_everything_pf.txt',sep = "\t",col.names = TRUE)
################
#make it into annotation format, sort out unnecessary columns
hits_everything_pf<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/outputs/hits_everything_pf.txt')
head(hits_everything_pf)

###load neoantigen####
load("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/neoantigen_objects.RData")
#create datasets for different data
#expression_data <- as.data.frame(exprs(peptides_novel_normalized) )#contains quantitive data. relative abundance
#phenotype_data  <- pData(peptides_novel_normalized)  #contains temperature
feature_data<- featureData(peptides_novel_normalized)@data #contains Novel (T/F)
#expression_data_t<-as.data.frame(t(expression_data))
feature_data$peptide_no<-rownames(feature_data)#extract peptide no
feature_data_n<-feature_data

# Merge datasets by matching the peptide columns
hit_everything_merged <- merge(
  hits_everything_pf, 
  feature_data_n[, c("peptide_no", "Peptide_AA", "peptide")], 
  by.x = "Peptide", 
  by.y = "peptide", 
  all.x = TRUE
)

# Print the resulting data frame
head(hit_everything_merged,n=2)
write.table(hits_everything,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/outputs/hits_everything_merged.txt',sep = "\t",col.names = TRUE)

#drop some columns
# Keep only the specified columns and rename them
pf_anno <- hit_everything_merged[, c("Peptide", "Peptide_AA", "peptide_no", "proteoform_id")]
colnames(pf_anno) <- c("old_seq", "Peptide_AA", "peptide_no", "proteoform_id")

# Print the resulting data frame
head(pf_anno)


#try get membership
get_peptide_memberships <- function (graphs) {
  processed_graphs <- lapply(names(graphs), function (ioi) {
    g <- graphs[[ioi]]
    if (!inherits(g, "igraph") || length(V(g)) == 0) return(NULL)
    data.frame(ioi = ioi,
               peptide = V(g)$name,
               membership = V(g)$membership,
               stringsAsFactors = FALSE)
  })
  
  # Filter out NULL values and empty data frames
  valid_graphs <- Filter(Negate(is.null), processed_graphs)
  valid_graphs <- valid_graphs[sapply(valid_graphs, nrow) > 0]
  
  # Merge data frames and assign to get_mem_result
  get_mem_result <- do.call(rbind, valid_graphs)
  
  return(get_mem_result)
  
}
library(dplyr)

memberships<-get_peptide_memberships(graphs_comm)


# Ensure columns are correctly typed
memberships <- memberships %>%
  mutate(across(c(ioi, id, peptide), as.character))

data<-feature_data_n

data_memberships <- data %>%
  left_join(memberships, by = "peptide") %>%
  # filter out rows where ioi is NA and the ID contains multiple entries
  filter(!(is.na(ioi) & grepl(";", id))) %>%
  # assign ioi, if there is none
  mutate(ioi = if_else(is.na(ioi), id, ioi)) %>%
  # if the id is not in the graph names, then assign membership to 0
  mutate(membership = if_else(ioi %in% names(graphs_comm), membership, 0)) %>%
  # exclude rows with no id and no membership information
  filter(!is.na(id)) %>%
  filter(!is.na(membership)) %>%
  # create proteoform ID
  mutate(proteoform_id = paste0(ioi, "_", membership))

head(data_memberships,n=2)
write.table(hits_everything,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/outputs/data_memberships_full.txt',sep = "\t",col.names = TRUE)

#filter out unuseful columns for annottaion
# Create the new dataframe 'pf_anno' by selecting specified columns
pf_anno <- data_memberships %>%
  select(peptide, Peptide_AA, peptide_no, proteoform_id)

# Print the resulting dataframe
print(head(pf_anno))

write.table(pf_anno,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/Input/ProteoformAnnotation',sep = "\t",col.names = TRUE)


