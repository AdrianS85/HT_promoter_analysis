##############################
### LOAD DATA AND METADATA ###
##############################
source('../GRS - GJt Review Stress/FULL_DATASET/functions_for_genename_conversions.R')
source(file = 'https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# rm(list = ls(pattern = '(.*)(test)|temp(.*)'))



############################
#### PREPARE INPUT DATA ####
############################
input_genes <- readr::read_tsv('input_genes.tsv')

# Annotate gene names with entrez geneID
create_query_for_ncbi <- get_query_for_ncbi_geneID_annotation(char_vec_gene_id_to_query_with = input_genes$Gene_name, char_vec_organism = input_genes$Species, chr_gene_identifier = 'Gene name')[[1]]
entrez_geneIDs <- search_for_ids_in_ncbi(str_vector_of_ids = create_query_for_ncbi)

input_genes$entrez_geneID <-
  as.character(lapply(
    X = entrez_geneIDs,
    FUN = function(x) {
      ifelse(test = length(x$ids) != 0 ,
             yes = x$ids[[1]],
             no = '')
    }
  ))

# These two LOCs are for some reason not found when querying for [Gene name]. it gives result for quering this LOC without [Gene name]. This is wierd... Better note that
input_genes$entrez_geneID[[59]] <- '51534'
input_genes$entrez_geneID[[80]] <- '80054'
input_genes$entrez_geneID[[67]] <- '3172' # HNF4A (postion 67) is annotated in NCBI as HNF1A. Entrez automatically annotates HNF4A as HNF1A. Here, I manually attach proper number

readr::write_tsv(input_genes, 'input_genes_with_entrezID.tsv')
readr::write_file(paste(input_genes$Gene_name, collapse = '\n'), 'enrichr_input.txt') 
############################
#### PREPARE INPUT DATA ####
############################



##############################
### LOAD DATA AND METADATA ###
##############################
input_genes <- readr::read_tsv('input_genes_with_entrezID.tsv')

# BiocManager::install("BSgenome")
BS_homo <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 
TxDb_homo <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
##############################
### LOAD DATA AND METADATA ###
##############################



### !!! Read this well: https://bioconductor.org/packages/release/workflows/vignettes/generegulation/inst/doc/generegulation.html
########################################################
### GETTING PROMOTER SEQUENCES FOR GENES OF INTEREST ###
########################################################
# This returns IRangnes for every transcript of the gene. It does it correctly. Still, there are perfectly duplicated coordinates for single gene here - such that there are multiple transcripts within the same start-end sequences?
transcript_coordinates <- GenomicFeatures::transcriptsBy(x = TxDb_homo, by= "gene")[as.character(input_genes$entrez_geneID)]
save(transcript_coordinates, file = 'transcript_coordinates')



# This returns promoter sequences for every transcript
promoter_sequence <- GenomicFeatures::getPromoterSeq(query = transcript_coordinates, subject = BS_homo, upstream = 1000, downstream = 2000)
save(promoter_sequence, file = 'promoter_sequence')



# This returns unlisted sequences for every transcript
unlist_promoter_sequence <- unlist(promoter_sequence)
names(unlist_promoter_sequence) <- stringr::str_remove(string = names(unlist_promoter_sequence), pattern = '\\..*')
save(unlist_promoter_sequence, file = 'unlist_promoter_sequence')



# TRANSFAC and JASPAR are databases of transcription factor binding models.
# BiocManager::install("JASPAR2018"), BiocManager::install("TFBSTools") - package for promoter analysis;  Other packages that may be helpful: motifStack
gff_of_unlist_promoter_sequence <- get_and_write_gff_from_DNAStringSet_via_JASPAR(DNAStringSet_ = unlist_promoter_sequence, input_genes_ = input_genes, opts_species = 9606, opts_name = "HNF1A")
save(gff_of_unlist_promoter_sequence, file = 'gff_of_unlist_promoter_sequence')



# This is because there are perfectly duplicated coordinates for single gene here - such that there are multiple transcripts within the same start-end sequences. Thus, we get duplicated promoter regions. We still have multiple promoters for single gene, because different transcripts have different start sites.
unique_gff_of_unlist_promoter_sequence <- unique(gff_of_unlist_promoter_sequence)

dataset_length <- length(unique_gff_of_unlist_promoter_sequence[,1])
dataset_mean_score <- mean(unique_gff_of_unlist_promoter_sequence$score)
dataset_mean_start <- mean(unique_gff_of_unlist_promoter_sequence$start)
########################################################
### GETTING PROMOTER SEQUENCES FOR GENES OF INTEREST ###
########################################################



#############################################################
### COMPARE GENES OF INTEREST WITH RANDOM SAMPLE OF GENES ###
#############################################################
no_of_samples_for_sampling_distribution <- vector("list", length = 100)
sample_list_for_sampling_distribution <-
  lapply(
    X = no_of_samples_for_sampling_distribution,
    FUN = function(x) {
      print(paste0('...'))
      temp <- get_random_genes(sample_size_int = 250, species_str = 'Human')
      temp <- subset(x = temp, subset = !is.na(temp$`NCBI Gene ID`))
    }  )
save(sample_list_for_sampling_distribution, file = 'sample_list_for_sampling_distribution')



# works with hg19, not with hg38 - hg38 have shorter chromosomes and sometimes downstream 2000 goes beyond chromosome fucking up entire GRanges object. I dont think its worth the effort to remedy this, just use hg19
unlisted_random_promoter_sequences_for_sampling_distribution <-
  purrr::map(
    .x = sample_list_for_sampling_distribution,
    .f = function(x) {
      get_promoter_DNAStringSet_with_corresponding_geneNames_after_removal_of_unknown_genes(geneids_chrvec_ = x$`NCBI Gene ID`, txdb__ = TxDb_homo, sample_size_int_ = 80, parent_genes_table_to_geneids_chrvec_ = x, BS_genome = BS_homo, upstream_ = 1000, downstream_ = 2000)
    })
save(unlisted_random_promoter_sequences_for_sampling_distribution, file = 'unlisted_random_promoter_sequences_for_sampling_distribution')
# load('unlisted_random_promoter_sequences_for_sampling_distribution')

randoms_unlisted_lengths_vector <- purrr::map_dbl(.x = unlisted_random_promoter_sequences_for_sampling_distribution, .f = function(x) { sum(width(x[[1]])) })



# I have lost genenames for sequences at some point. It doesnt matter for my analysis though, so I will just leave it.
gff_of_unlisted_random_promoter_sequences_for_sampling_distribution <-
purrr::map(
  .x = unlisted_random_promoter_sequences_for_sampling_distribution,
  .f = function(x) {
    print(paste0('...'))
    get_and_write_gff_from_DNAStringSet_via_JASPAR(DNAStringSet_ = x[[1]], input_genes_ = x[[2]], opts_species = 9606, opts_name = "HNF1A", input_col_with_entrezID = 'NCBI Gene ID', input_col_with_geneName = 'Approved symbol')
  })
save(gff_of_unlisted_random_promoter_sequences_for_sampling_distribution, file = 'gff_of_unlisted_random_promoter_sequences_for_sampling_distribution')
# load('gff_of_unlisted_random_promoter_sequences_for_sampling_distribution')


unique_gff_of_unlisted_random_promoter_sequences_for_sampling_distribution <-
  purrr::map(
    .x = gff_of_unlisted_random_promoter_sequences_for_sampling_distribution,
    .f = function(x) {
      print(paste0('...'))
      unique(x)
    })

randoms_lengths_vector <- purrr::map_int(.x = unique_gff_of_unlisted_random_promoter_sequences_for_sampling_distribution, .f = function(x) { length(x[[1]]) })
mean(randoms_lengths_vector)
sd(randoms_lengths_vector)
density(randoms_lengths_vector)

randoms_mean_scores_vector <- purrr::map_dbl(.x = unique_gff_of_unlisted_random_promoter_sequences_for_sampling_distribution, .f = function(x) { mean(x$score) })
mean(randoms_mean_scores_vector)
sd(randoms_mean_scores_vector)
density(randoms_mean_scores_vector)

randoms_mean_start_vector <- purrr::map_dbl(.x = unique_gff_of_unlisted_random_promoter_sequences_for_sampling_distribution, .f = function(x) { mean(x$start) })
mean(randoms_mean_start_vector)
sd(randoms_mean_start_vector)
density(randoms_mean_start_vector)

#############################################################
### COMPARE GENES OF INTEREST WITH RANDOM SAMPLE OF GENES ###
#############################################################




BiocManager::install("motifStack")
rm(opts)

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- 'HNF1A'
opts[["all_versions"]] <- T
opts[["matrixtype"]] <- "PFM"

temp_resulting_PFMatrixList_ <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)

temp_PFMatrix <- temp_resulting_PFMatrixList_[[1]]



tags(temp_resulting_PFMatrixList_)
seqLogo(temp_PFMatrix)

TFBSTools::view_motifs(temp_resulting_PFMatrixList_)
# Databases like ORegAnno have attempted to collect and curate individual experimentally proven TFBSs from the literature and sets from ChIP-Seq and other genome-wide experiments
# Many of these TF binding sites (TFBSs) are being deposited into species-specific databases (FlyBase, SGD, etc) and as tracks in the UCSC browser. tracks under 'Regulation' in the UCSC Genome Browser for your species of interest. UCSC has oreganno data: Data are also available through UCSC Genome Browser (e.g., hg38 -> Regulation -> ORegAnno). http://gtrd.biouml.org/ 


# try jeszcze to: 
# https://bioconductor.org/packages/release/workflows/vignettes/generegulation/inst/doc/generegulation.html

library(MotifDb)

test <- query(MotifDb, c("DAL80"))
pcm.dal80.jaspar <- round(100 * pfm.dal80.jaspar)
matchPWM(pcm.dal80.jaspar, unlist(promoter.dal1)[[1]], "90%")

###################################################
### TRANSCRIPTION FACTOR BINDING SITES ANALYSIS ###
###################################################






###################################################
### TRANSCRIPTION FACTOR BINDING SITES ANALYSIS ###
###################################################
########################
### FURTHER ANALYSIS ###
########################

# test <- rtracklayer::readGFF(filepath = 'test.gff3')

# session <- rtracklayer::browserSession("UCSC")
# 
# rtracklayer::track(session, "targets") <- resulting_gff3_of_SiteSetList

#######################
### SIMPLE ANALYSIS ###
#######################
# genes_with_bindings_sites_detected <- unique(as.character(resulting_gff3_of_SiteSetList$seqname))
# try SiteSet methods. They have pvalue, perhaps it will help.
# try visualizing hnf1a binding sites?
#######################
### SIMPLE ANALYSIS ###
#######################

# database_str <- names(json_file)
# 
# json_file <- jsonlite::fromJSON(test)
# 
# json_file2 <- lapply(json_file[[1]], function(x) {
#   
#   temp_2 <- list()
#   
#   for (n in seq_along(x)) {
#     if (length(x[[n]]) == 1) {
#       temp_2[[n]] <- x[[n]]
#     }
#     else if (length(x[[n]]) > 1) {
#       temp_2[[n]] <- paste0(x[[n]], collapse = ', ')
#     }
#     else {
#       stop('error', call. = T)
#     }
#   }
#   
#   return(temp_2)
#   
# })
# 
# unlist_json <- lapply(X = json_file2, FUN = function(x) { unlist(x)} )
# 
# unlist_json2 <- as.data.frame(rlist::list.rbind(unlist_json))
# 
# colnames(unlist_json2) <- c('Database', 'Term', 'P-value', 'Odds_Ratio', 'Combined_Score','Genes', 'Adj_P-value', 'no_1', 'no_2')
# 
# unlist_json2$no_1 <- NULL
# unlist_json2$no_2 <- NULL
# unlist_json2$Database <- database_str


#tutaj spróbujmy wrzucić for

temp_1 <- json_file[[1]][[1]]
class(unlist(temp_1))
unlist(temp_1)



temp_2 <- list()
for (n in seq_along(temp_1)) {
  if (length(temp_1[[n]]) == 1) {
    temp_2[[n]] <- temp_1[[n]]
  }
  else if (length(temp_1[[n]]) > 1) {
    temp_2[[n]] <- paste0(temp_1[[n]], collapse = ', ')
  }
  else {
    stop('error', call. = T)
  }
}
return(temp_2)



json_file2 <- as.data.frame(do.call("cbind", json_file))
as.integer()




### MULTIPLE ALIGNMENT ###
# Muscle algoritm is easy-to-use due to being in msa package, it is very fast and reasonably accurate (Evaluating the Accuracy and Efficiency of Multiple Sequence Alignment Methods)
# We can gry using another R-compatible algorithm: MAFFT(L-INS-i) which perforomed better than Muscle here: Evaluating the Accuracy and Efficiency of Multiple Sequence Alignment Methods, though it seems to take a long time. Love you long time soldierboy. And msa takes reaal long to work either way. May need to paralelize it. This publication also suggest MAFFT: Assessing the efficiency of multiple sequence alignment programs. So lets try parallelize the shit.
# https://www.rdocumentation.org/packages/ape/versions/5.3/topics/as.alignment
# https://www.rdocumentation.org/packages/ips/versions/0.0.11/topics/mafft

# r packages for multiple alignment muscle, odseq, DECIPHER
# DECIPHER claims, that it is better to first translate nucleotide sequence to aminoacids and than align, but am not sure, cause for promoters, aminoacids do not make sense

# msa running on all 80 sequences givec below error:
# *** ERROR ***  MSA::GetLetter(0/1, 24/66577)=''/4294967295
# 
# Fatal error, exception caught.
# Error in msaFun(inputSeqs = inputSeqs, cluster = cluster, gapOpening = gapOpening,  : 
#                   MUSCLE finished by an unknown reason

# must try on shortened dataset

# Look for promoter alignment specifically

# Need to actually load the library for function to work
# library(msa)
# devtools::install_github("olafmersmann/microbenchmark")
# 
# mbm <- rbenchmark::benchmark("msa_test" = { msa_test <- msa::msa(inputSeqs = unlist_promoter_sequence_test, method = "Muscle") })






### MULTIPLE ALIGNMENT ###















# # ENSEMBL_MART_FUNCGEN The Ensembl Regulatory Build (Zerbino et al. 2015) contains a genome-wide set of regions that are likely to be involved in gene regulation.
# pa_biomart_ <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
# pa_biomart_funcgen <- biomaRt::useMart("ENSEMBL_MART_FUNCGEN")
# 
# 
# pa_test <- biomaRt::getSequence(id = 'Gapdh', 
#                                 type = "external_gene_name",
#                                 seqType="coding_gene_flank",
#                                 upstream=100, 
#                                 mart=pa_biomart_)
# 
# pa_test2 <- biomaRt::getSequence(id = 'Gapdh', 
#                                  type = "external_gene_name",
#                                  seqType="coding_gene_flank",
#                                  downstream = 100, 
#                                  mart=pa_biomart_)
# 
# ### BIOMART ###
# 
# biomartr::getMarts()
# 
# 
# pa_filters <- biomaRt::listFilters(mart = pa_biomart_)
# 
# 
# pa_mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
# 
# pa_seq <- biomaRt::getSequence(id = "BRCA1",type = "hgnc_symbol",seqType = "peptide",mart = pa_mart)
# 
# show(pa_seq)
# 
# 
# ### BIOMART ###
# 
# pa_mart_test = biomaRt::useMart('ensembl')
# 
# biomaRt::listDatasets(pa_mart_test)
# 
# biomaRt::listMarts()
# listDatasets_pa_biomart_funcgen <- biomaRt::listDatasets(mart = pa_biomart_funcgen)




# Unlist the list generated by getPromoterSeq to produce input for msa
# unlist_promoter_sequence_test <- Biostrings::unstrsplit(promoter_sequence_test)
# unlist_promoter_sequence <- Biostrings::unstrsplit(promoter_sequence)
# save(unlist_promoter_sequence, file = 'unlist_promoter_sequence')
# sum(width(unlist_promoter_sequence))
### !!! This returns composite sequence generated by adding promoters for every transcript. Perhaps we should unlist it in different way?
