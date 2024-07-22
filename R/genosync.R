#' @title \code{genosync}
#'
#' @description \code{genosync} To demultiplex pooled samples, this wrapper function runs two independent methods on hash data from 
#' an input Seurat object. The first method combines kmeans and apriori association analysis, [kmeansync()]. 
#' The second method runs binary logistic regression, [logisync()]. Compares both methods for all Souporcell runs,
#' and outputs matching results. If there are no matches, all results will be output. 
#'
#' @usage genosync(seu_obj, hash_csv, max_soup_run)
#'
#' @param seu_obj The input Seurat object. Must contain hash assay named \code{HTO} or \code{hto} in dgCMatrix format.
#' Must contain Souporcell genotype assays named \code{GENO} or \code{geno} with the desired k value (ex: \code{GENO5}). 
#' Assumes Souporcell was run renaming all multiplet designations to 'multiplet' (see [viewmastR::add_souporcell_seurat()]).
#'
#' @param hash_csv The input hash-sample csv file path. The csv must contain \code{Hash} and \code{Sample} columns.
#' 
#' @param max_soup_run The highest Souporcell run, indicating the number of genotypes detected. Default is 8.
#'
#' @return If there are consensus matches, the kmeans and logistic results for those Souporcell runs will be returned.
#' If there are no matches, all kmeans and logistic results will be output for all Souporcell runs.
#'
#'
#'
#' @examples
#' \dontrun{
#' 
#'  output_list <- genosync(seu, csv='/path/to/hash_sampleABCD.csv')
#'   
#' }
#'
#' @importFrom dplyr filter mutate case_when rowwise
#' @importFrom utils read.csv
#' @importFrom methods as
#' 
#' @export

genosync <- function(seu_obj, hash_csv, max_soup_run=8){
  
  # find min number for souporcell runs
  hashtable <- read.csv(hash_csv) 
  min_genos <- length(unique(hashtable$Hash))

  # check if max_soup_run numbers are valid
  patterns <- lapply(min_genos:max_soup_run, function(num){paste0('GENO', num)})
  matches <- sapply(patterns, function(pattern) any(grepl(pattern, names(seu_obj@assays), ignore.case = TRUE)))
  if(!all(matches)){
    stop('Input seurat object missing a genotype assay named GENO or geno with the desired k value (ex: GENO5). Check max_soup_run value.')}
  
  outs_kmean <- list()
  outs_log <- list()
  # iterate over souporcell, run kmeans and log reg
  for(soup_num in c(min_genos:max_soup_run)){
    message(paste0("\n", "Running kmeansync for souporcell = ", soup_num, "\n"))
    outs_kmean[[paste0('Soup_', soup_num)]] <- kmeansync(seu_obj, hash_csv=hash_csv, soup_k=soup_num, res=TRUE)
    
    message(paste0("\n", "Running logisync for souporcell = ", soup_num, "\n"))
    outs_log[[paste0('Soup_', soup_num)]] <- logisync(seu_obj, hash_csv=hash_csv, soup_k=soup_num, res=TRUE)
    message(strrep("_", 75))
  }
  
  # ---------------------------------------------------------------------------------
  # helper function to reorder and compare assignments
  compare_samples <- function(df1, df2){
    # split and sort genotypes
    df1 <- df1 %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Soup_sorted = paste(sort(trimws(unlist(strsplit(.data$Soup, ",")))), collapse = ","))
    
    df2 <- df2 %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Soup_sorted = paste(sort(trimws(unlist(strsplit(.data$Soup, ",")))), collapse = ","))

    # compare genotypes
    return(all(df1$Soup_sorted == df2$Soup_sorted))}
  # ----------------------------------------------------------------------------------
  
  # see if results match between any souporcell runs
  matches <- list()
  for(soup_num in c(min_genos:8)){
    # ensure outputs exist
    if(!is.null(outs_kmean[[paste0('Soup_', soup_num)]]) && !is.null(outs_log[[paste0('Soup_', soup_num)]])){
      # sample dfs
      kmean_df <- outs_kmean[[paste0('Soup_', soup_num)]][[3]]  
      log_df <- outs_log[[paste0('Soup_', soup_num)]][[3]]      
      
      # check matching
      if(compare_samples(kmean_df, log_df)){
        message(paste0("\n", "Match found for souporcell = ", soup_num, "\n"))
        # add to list
        matches[[paste0('Soup_', soup_num)]] <- list(soup_num = soup_num, kmeans_result = outs_kmean[[paste0('Soup_', soup_num)]], log_result = outs_log[[paste0('Soup_', soup_num)]])
      }
    }
  }
  
  # in case multiple matches occur
  if(length(matches) > 0){
    return(matches)
  # return all results if no matches
  }else{ 
    message("\n", "No matching results found. Returning all results. \n")
    return(list(outs_kmean, outs_log))  
  }
}
