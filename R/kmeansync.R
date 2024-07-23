#' @title \code{kmeansync}
#'
#' @description \code{kmeansync} To demultiplex pooled samples, runs kmeans on hash data from an input Seurat object.
#' Hashes are assigned to kmeans clusters based on the average hash enrichment in each cluster. Apriori association rules 
#' link hashes from kmeans clusters to Souporcell genotypes. Then, genotypes can be linked to samples via the input hash-sample csv. 
#'
#' @usage kmeansync(seu_obj, csv, soup_k, conf=0.8, output_col='FinalAssignment', res=FALSE)

#'
#' @param seu_obj The input Seurat object. Must contain hash assay named \code{HTO} or \code{hto} in dgCMatrix format.
#' Must contain Souporcell genotype assay(s) named \code{GENO} or \code{geno} with the desired k value (ex: \code{GENO5}). 
#' Assumes Souporcell was run renaming all multiplet designations to 'multiplet' (see [viewmastR::add_souporcell_seurat()]
#' documentation).
#'
#' @param csv The input hash-sample csv file path. The csv must contain \code{Hash} and \code{Sample} columns.
#'
#' @param soup_k The desired Souporcell run, indicating the number of genotypes detected. The appropriate 
#' number of genotypes expected for the data should be determined prior to running kmeansync. 
#'
#' @param conf The confidence for apriori association rules, a measure of how often rules are found to be true. 
#' Default is 0.8, as recommended.
#'
#' @param output_col The name for the sample assignments column that will be added to the output Seurat object.
#' Default is \code{FinalAssignment}.
#'
#' @param res Boolean indicating additional results output. Default is FALSE. If TRUE, Seurat object with sample 
#' and hash labels will be returned, along with a dataframe linking kmeans clusters to genotypes, a dataframe 
#' linking genotypes to samples, the association rules, the silhouette graph with optimal kmeans clusters, and 
#' a graph of kmeans clustering, shaded with average hash per cluster.
#'
#' @return seu_obj Returns Seurat object with sample and hash labels by default. Note that multiplets have been excluded.
#' If \code{res=TRUE}, additional dataframes and graphs are returned as well.
#'
#'
#'
#' @examples
#' \dontrun{
#' 
#'  output_list4 <- kmeansync(seu, csv='/path/to/hash_sampleABCD.csv', soup_k=4, res=TRUE)
#'  
#'  seu6 <- kmeansync(seu, csv='/path/to/hash_sampleDEF.csv', soup_k=6, output_col='Sample_Assignment') 
#'   
#' }
#'
#' @importFrom dplyr filter mutate case_when rowwise group_by summarize_all mutate_all ungroup
#' @import PNWColors
#' @import factoextra
#' @import umap
#' @import ggforce
#' @import ggplot2
#' @importFrom arules subset apriori lhs rhs inspect
#' @importFrom Seurat NormalizeData ScaleData 
#' @importFrom SeuratObject DefaultAssay
#' @importFrom utils read.csv capture.output
#' @importFrom methods as
#' @importFrom stats kmeans as.formula 
#' 
#' @export

kmeansync <- function(seu_obj, csv, soup_k, conf=0.8, output_col='FinalAssignment', res=FALSE){
  
  UMAP1 <- UMAP2 <- lift <- cluster <- NULL
  
  if(!inherits(seu_obj, 'Seurat')){
    stop('Input must be a seurat object.')}
  
  # ensure seu_obj has HTO assay
  if(!any(names(seu_obj@assays) %in% c('HTO', 'hto'))){
    stop('Input seurat object must contain a hash assay named HTO or hto.')}
  
  # check HTO matrix
  if(!inherits(seu_obj@assays$HTO$counts, "dgCMatrix")){
    stop('HTO assay formatted incorrectly in seurat object. HTO assay counts must be a dgCMatrix.')}
  
  # ensure seu_obj has GENO assay
  pattern <- paste0('GENO', soup_k)
  if(!any(grepl(pattern, names(seu_obj@assays), ignore.case = TRUE))){
    stop('Input seurat object must contain a genotype assay named GENO or geno with the desired k value (ex: GENO5).')}
  
  # set assay
  SeuratObject::DefaultAssay(seu_obj) <- "HTO"
  # remove multiplets (avoid noise when assigning hashes)
  geno_col <- paste0('geno', soup_k)
  seu_obj <- subset(seu_obj, cells=which(seu_obj@meta.data[[geno_col]] != 'Multiplet'))
  # normalize data for kmeans
  seu_obj <- NormalizeData(seu_obj, normalization.method = "CLR", margin = 1)
  # scale data for kmeans
  seu_obj <- ScaleData(seu_obj)
  
  # read in csv
  hash_table <- read.csv(csv) 
  
  # preprocess data
  data <- Seurat::FetchData(seu_obj, c(unique(hash_table$Hash), geno_col))
  # convert to numeric for kmeans
  data[[geno_col]] <- as.numeric(data[[geno_col]])
  
  # dimensionality reduction prior to kmeans
  umap_res <- umap::umap(data)
  
  # determine optimal k for kmeans
  sil <- factoextra::fviz_nbclust(umap_res$layout, kmeans, method='silhouette', k.max=15, nstart = 250, iter.max = 10000) +
    labs(subtitle = "Silhouette method")
  # get top k value
  optimal_k <- as.numeric(sil$data$clusters[which.max(sil$data$y)])
  
  # if optimal k less than number of hashes/samples, impose a minimum k
  if(optimal_k < length(unique(hash_table$Hash))){
    optimal_k <- length(unique(hash_table$Hash))
    message(paste0('Adjusted optimal kmeans clusters to match number of hashes: ', optimal_k))
  }else{
    message(paste0('Number of clusters selected for kmeans: ', optimal_k))
  }
  
  # kmeans
  set.seed(1234)
  kmeans_res <- kmeans(umap_res$layout, centers = optimal_k, iter.max=10000, nstart = 250)
  data$cluster <- as.numeric(kmeans_res$cluster)
  
  # Assign hashes to kmeans clusters:
  hto_columns <- setdiff(colnames(data), c('cluster', geno_col))
  numeric_hto_columns <- sapply(data[hto_columns], is.numeric)
  
  # calculate average HTO enrichment per cluster
  ave_hash <- data %>% dplyr::select(which(numeric_hto_columns), cluster) %>% dplyr::group_by(cluster) %>% dplyr::summarize_all(mean, na.rm=TRUE)
  # look at max in averages to define clusters
  representative_HTO <- apply(ave_hash[, -1], 1, function(x) names(x)[which.max(x)])
  
  # assign HTO to kmeans clusters
  cluster_assignments <- data.frame(
    cluster = as.factor(1:nrow(kmeans_res$centers)),
    representative_HTO = representative_HTO)
  
  # UMAP K means cluster graph with hash assignments
  cluster_palette <- PNWColors::pnw_palette("Sailboat",6,type="discrete")
  # add color col to separate cluster_assignments df
  cluster_assignments_col <- cluster_assignments %>% mutate(colors=cluster_palette[1:nrow(cluster_assignments)])
  # order by cluster
  cluster_assignments_col <- cluster_assignments_col[order(cluster_assignments_col$cluster), ]
  # extract colors for clusters
  cluster_colors <- cluster_assignments_col$colors
  # order by representative_HTO
  cluster_assignments_col <- cluster_assignments_col[order(cluster_assignments_col$representative_HTO), ]
  # extract colors for hashes
  hash_colors <- cluster_assignments_col$colors
  # umap data
  umap_data <- data.frame(UMAP1 = umap_res$layout[,1], UMAP2 = umap_res$layout[,2], cluster = data$cluster)
  umap_data <- base::merge(umap_data, cluster_assignments_col, by = "cluster")
  # umap colored by cluster
  graph <- ggplot2::ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
    geom_point(size = 3) +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_manual(values = cluster_colors) +
    labs(col = 'Cluster', title = paste0('UMAP with Kmeans Hash Assignment, k = ', optimal_k))
  
  # Add outline/background based on hash
  graph <- graph + 
    ggforce::geom_mark_ellipse(aes(fill = representative_HTO), color = NA, alpha = 0.3) + 
    scale_fill_manual(values = hash_colors) +
    labs(fill = 'Hash')
  
  # apriori association between kmeans clusters and genotypes
  apriori_data <- data %>% dplyr::select(!!sym(geno_col), cluster) %>% dplyr::mutate_all(as.factor) %>% as("transactions")
  rules <- arules::apriori(apriori_data, parameter = list(support = 0.05, confidence = conf, minlen=2), control=list(verbose = FALSE))
 
  # Link hashes to association results
  pairs <- data.frame() 
  # only keep rules with lift > 1.5 (indicates likelihood, 1 = not likely)
  sig_rules <- arules::subset(rules, subset = lift > 1.5) 
  # extract rule pairs
  for (i in 1:length(sig_rules)) {
    lhs_string <- as(lhs(sig_rules[i]), 'list')
    rhs_string <- as(rhs(sig_rules[i]), 'list')
    
    # merge lhs and rhs 
    pair <- paste(lhs_string, rhs_string, sep = " ")
    # add to df
    pairs <- rbind(pairs, data.frame(pair = pair))
  }
  
  # reorder pairs for consistency (separate kmeans clusters and genotypes)
  pairs$pair <- sapply(pairs$pair, function(pair){
    elements <- unlist(strsplit(pair, " "))
    sorted_elements <- sort(elements)
    return(paste(sorted_elements, collapse = " "))})
  
  # remove any duplicates
  pairs <- unique(pairs)
  
  # separate into columns
  pairs <- data.frame(do.call(rbind, strsplit(pairs$pair, " ")))
  # colnames
  if(any(grepl("geno", pairs$X1))){
    colnames(pairs) <- c(geno_col, "cluster")
  }else{
    colnames(pairs) <- c("cluster", geno_col)
  }
  
  # extract data values
  pairs$cluster <- gsub(".*=(\\d+).*", "\\1", pairs$cluster)
  pairs[[geno_col]] <- gsub(".*=(\\d+).*", "\\1", pairs[[geno_col]])
  
  # note if any genotypes weren't assigned
  missing_genotype <- setdiff(unique(data[[geno_col]]), unique(pairs[[geno_col]]))
  if(length(missing_genotype) > 0){
    message(paste0('No association rules found for genotype: ', missing_genotype, "\n"))
  }

  # merge with hash assignments, rename columns
  HTO_map <- base::merge(pairs, cluster_assignments, by = "cluster", all.x = TRUE)
  colnames(HTO_map) <- c('Kmeans', 'Soup', 'Hash')

  # merge with sample assignments from input csv
  association <- base::merge(hash_table, HTO_map, by='Hash')
  
  # add hash and sample assignment columns to data
  data <- data %>% dplyr::rowwise() %>%
    dplyr::mutate(match_index = base::match(cluster, association$Kmeans),
                  HashAssignment = dplyr::case_when(!is.na(match_index) ~ association$Hash[match_index], TRUE ~ NA_character_),
                  FinalAssignment = dplyr::case_when(!is.na(match_index) ~ association$Sample[match_index], TRUE ~ NA_character_)) %>%
    dplyr::ungroup() 
  
  # add results to seu metadata
  seu_obj[[output_col]] <- data$FinalAssignment
  seu_obj$HashAssignment <- data$HashAssignment
  
  # return additional results tables
  if(res){
    # create final table linking soups to hashes (allowing for multiple soups/hash)
    final <- data.frame()
    for(hto in unique(association$Hash)){
      matched_soups <- association[association$Hash == hto, "Soup"]
      # add unique soups to string
      soups_per_hash <- paste(unique(matched_soups), collapse = ", ")
      final <- rbind(final, data.frame(Hash=hto, Soup=soups_per_hash))
    }
    
    # assoc_output
    association_out <- association[, c("Hash", "Sample", "Soup", "Kmeans")] 
    # if a kmeans cluster wasn't assigned, add NA
    missing_clusters <- setdiff(1:optimal_k, unique(association_out$Kmeans))
    for(clust in missing_clusters){ 
      new_row <- nrow(association_out) + 1
      # add missing kmeans cluster
      association_out[new_row,'Kmeans'] <- clust
      # add hash from kmeans cluster
      association_out[new_row,'Hash'] <- cluster_assignments %>% dplyr::filter(cluster==clust) %>% 
        dplyr::select(representative_HTO)
      # fill sample and soup with NA in quotes for aesthetics, since no association rules were found
      association_out[new_row,] <- replace(association_out[new_row,], is.na(association_out[new_row,]), 'NA')
    }
    # if genotype not assigned, add NA
    for(geno in missing_genotype){
      new_row <- nrow(association_out) + 1
      # add missing genotype 
      association_out[new_row,'Soup'] <- geno
      # fill rest with NA in quotes for aesthetics, since no association rules were found
      association_out[new_row,] <- replace(association_out[new_row,], is.na(association_out[new_row,]), 'NA')
    }
    
    # merge final and input csv
    assigned_df <- base::merge(hash_table, final, by='Hash')
    # contain output from rules
    inspect_rules <- capture.output(inspect(rules))
    # all results
    return(list(seu_obj, association_out, assigned_df, noquote(inspect_rules), sil, graph))
    
    # return seu obj by default
  }else{
    return(seu_obj)
  }
}


