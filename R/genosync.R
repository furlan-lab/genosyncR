#' @title \code{genosync}
#'
#' @description \code{genosync} To demultiplex pooled samples, this wrapper function runs two independent methods on hash data from 
#' an input Seurat object. The first method combines kmeans and apriori association analysis, [kmeansync()]. 
#' The second method runs binary logistic regression, [logisync()]. Compares both methods for all Souporcell runs,
#' and outputs matching results. If there are no matches, all results will be output. 
#'
#' @usage genosync(seu_obj, hash_csv, soup_runs)
#'
#' @param seu_obj The input Seurat object. Must contain hash assay named \code{HTO} or \code{hto} in dgCMatrix format.
#' Must contain Souporcell genotype assays named \code{GENO} or \code{geno} with the desired k value (ex: \code{GENO5}). 
#' Assumes Souporcell was run renaming all multiplet designations to 'multiplet' (see [viewmastR::add_souporcell_seurat()]).
#'
#' @param hash_csv The input hash-sample csv file path. The csv must contain \code{Hash} and \code{Sample} columns.
#' 
#' @param soup_runs A numeric vector of Souporcell runs to iterate through, indicating the number of genotypes detected.
#'
#' @return If there are consensus matches, the kmeans and logistic results for those Souporcell runs will be returned as a list.
#' If there are no matches, all kmeans and logistic results will be output for all Souporcell runs as a list.
#'
#'
#'
#' @examples
#' \dontrun{
#' 
#'  output_list <- genosync(seu_ABCD, hash_csv='/path/to/hash_sampleABCD.csv', soup_runs=c(3:8))
#'   
#' }
#'
#' @importFrom dplyr filter mutate case_when rowwise
#' @importFrom utils read.csv
#' @importFrom methods as
#' 
#' @export

genosync <- function(seu_obj, hash_csv, soup_runs){
  
  # if hash_csv is df
  if(is.data.frame(hash_csv)){
    hashtable <- hash_csv
  # if hash_csv is filepath 
  }else{
    if(file.exists(hash_csv)){
    hashtable <- read.csv(hash_csv)
    }}
  
  # check if max_soup_run numbers are valid
  patterns <- lapply(soup_runs, function(num){paste0('GENO', num)})
  matches <- sapply(patterns, function(pattern) any(grepl(pattern, names(seu_obj@assays), ignore.case = TRUE)))
  if(!all(matches)){
    stop('Input seurat object missing a genotype assay named GENO or geno with the desired k value (ex: GENO5). Check max_soup_run value.')}

  outs_kmean <- list()
  outs_log <- list()
  # iterate over souporcell, run kmeans and log reg
  for(soup_num in soup_runs){
    message(paste0("\n", "Running kmeansync for souporcell = ", soup_num, "\n"))
    outs_kmean[[paste0('Soup_', soup_num)]] <- kmeansync(seu_obj, csv=hash_csv, soup_k=soup_num, res=TRUE)
    if(is.character(outs_kmean[[paste0('Soup_', soup_num)]][[1]])){
      message(paste0('\n', 'No significant association rules found for Souporcell = ', soup_num))
    }
    message(paste0("\n", "Running logisync for souporcell = ", soup_num, "\n"))
    outs_log[[paste0('Soup_', soup_num)]] <- logisync(seu_obj, csv=hash_csv, soup_k=soup_num, res=TRUE)
    message(strrep("_", 75))
  }
  
  # ---------------------------------------------------------------------------------
  # helper function to reorder and compare assignments
  compare_samples <- function(df1, df2){
    # split and sort genotypes
    df1 <- df1 %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Soup_sorted = paste(sort(trimws(unlist(strsplit(.data$Genotype, ",")))), collapse = ","))
    
    df2 <- df2 %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Soup_sorted = paste(sort(trimws(unlist(strsplit(.data$Genotype, ",")))), collapse = ","))

    # compare genotypes
    return(all(df1$Soup_sorted == df2$Soup_sorted))}
  # ----------------------------------------------------------------------------------
  
  # see if results match between any souporcell runs
  matches <- list()
  for(soup_num in soup_runs){
    # ensure outputs exist
    if(!is.null(outs_kmean[[paste0('Soup_', soup_num)]]) && !is.character(outs_kmean[[paste0('Soup_', soup_num)]][[1]]) && 
       !is.null(outs_log[[paste0('Soup_', soup_num)]])){
      # sample dfs
      kmean_df <- outs_kmean[[paste0('Soup_', soup_num)]][[3]]  
      log_df <- outs_log[[paste0('Soup_', soup_num)]][[3]]      
      
      # check matching
      if(compare_samples(kmean_df, log_df)){
        message(paste0("\n", "Match found for souporcell = ", soup_num, "\n"))
        # add to list
        matches[[paste0('Soup_', soup_num)]] <- list(kmeans_result = outs_kmean[[paste0('Soup_', soup_num)]], log_result = outs_log[[paste0('Soup_', soup_num)]])
      }
    }
  }
  
  # in case multiple matches occur
  if(length(matches) > 0){
    return(matches)
  # return all results if no matches
  }else{ 
    message("\n", "No matching results found. Returning all results. \n")
    names(outs_kmean) <- paste0('kmeans_result_', names(outs_kmean))
    names(outs_log) <- paste0('log_result_', names(outs_log))
    
    return(c(outs_kmean, outs_log))  
  }
}
