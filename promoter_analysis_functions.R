# INPUT: geneids_chrvec - works with entrezIDs. many entrezIDs are not included in txdb object, thus this vector needs to be about twice as long as needed sample size OUTPUT: [[1]] - actual DNAStringSet, [[2]] - genes used in this DNAStringSet
get_random_GRangesList_from_txdb <- function(geneids_chrvec, txdb_, sample_size_int, genes_table)
{
  # This kinna works against 'subscript contains invalid names', I think I can create DNAStringSetList out of notmal list, but perhaps it would be easier to get gene ids from TxDb object?
  temp <- list()
  found <- rep(F, length(geneids_chrvec))
  for (n in seq_along(geneids_chrvec)) {
    tryCatch({
      temp[[n]] <-
        GenomicFeatures::transcriptsBy(x = txdb_, by = "gene")[as.character(geneids_chrvec[[n]])]
      found[n] <- T
    },
    error = function(ex) {
      print(paste0(n, ':  ', ex))
    })
  }

  temp2 <- list()
  for (m in seq_along(temp)) {
    if (!is.null(temp[[m]])) {
      temp2 <- rlist::list.append(temp2, unlist(temp[[m]]))
    }
  }
  
  temp_found_genes <- subset(genes_table, found)

  return(list(GenomicRanges::GRangesList(temp2)[1:sample_size_int], temp_found_genes[1:sample_size_int,]))
}



get_and_write_gff_from_DNAStringSet_via_JASPAR <- function(DNAStringSet_, input_genes_, opts_species, opts_name = "default", input_col_with_entrezID = 'entrez_geneID', input_col_with_geneName = 'Gene_name', output_file_name_ = NA, opts_all_versions = T, opts_matrixtype = "PWM")
{
  opts <- list()
  opts[["species"]] <- opts_species
  opts[["name"]] <- opts_name
  # opts[["type"]] <- "SELEX" # ChIP, SELEX, PBM
  opts[["all_versions"]] <- opts_all_versions
  opts[["matrixtype"]] <- opts_matrixtype
  
  resulting_MatrixList_ <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts) # This sould fetch proper data, checked with the JASPAR website
  
  # The SiteSet class is a container for storing a set of putative transcription factor binding sites on a nucleotide sequence (start, end, strand, score, pattern as a PWMatrix, etc.) from scaning a nucleotide sequence with the corresponding PWMatrix
  resulting_SiteSetList_ <- TFBSTools::searchSeq(x = resulting_MatrixList_, subject = DNAStringSet_)
  
  resulting_gff3_of_SiteSetList_ <- TFBSTools::writeGFF3(x = resulting_SiteSetList_, scoreType="relative")
  
  resulting_gff3_of_SiteSetList_$seqname <- recode_values_based_on_key(to_recode_chrvec = as.character(resulting_gff3_of_SiteSetList_$seqname), replace_this_chrvec = input_genes_[[input_col_with_entrezID]], with_this_chrvec = input_genes_[[input_col_with_geneName]])
  
  if(!is.na(output_file_name_)){
    write.table(x = resulting_gff3_of_SiteSetList_, file = output_file_name_, sep = '\t', row.names = F, col.names = F, quote = F)
  }
  return(resulting_gff3_of_SiteSetList_)
}



get_promoter_DNAStringSet_with_corresponding_geneNames_after_removal_of_unknown_genes <- function(geneids_chrvec_, txdb__, sample_size_int_ , parent_genes_table_to_geneids_chrvec_, BS_genome, upstream_, downstream_)
{
  print(paste0('Getting promoter sequences'))
  temp <- get_random_GRangesList_from_txdb(geneids_chrvec = geneids_chrvec_, txdb_ = txdb__, sample_size_int = sample_size_int_, genes_table = parent_genes_table_to_geneids_chrvec_)
  
  tryCatch({temp2 <- GenomicFeatures::getPromoterSeq(query = temp[[1]], subject = BS_genome, upstream = upstream_, downstream = downstream_)},
           error=function(e) {
             message(e)
             message("Returning NA")
             return(NA)
           })
  
  temp3 <- unlist(temp2)
  
  return(list(temp3, temp[[2]]))
}
