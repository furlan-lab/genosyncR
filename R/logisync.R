#' @title \code{logisync}
#'
#' @description \code{logisync} To demultiplex pooled samples, runs logistic regression on hash data from an input Seurat object.
#' Binary logistic regression is run for each genotype, and hashes are assigned based on which hash has the largest influence on 
#' the genotype (via average marginal effect). Then, genotypes can be linked to samples via the input hash-sample csv. 
#' Note that binary Firth logistic regression will be run for genotypes with small proportions to account for rare cases and reduce bias.
#'
#' @usage logisync(seu_obj, csv, soup_k, output_col='FinalAssignment', res=FALSE)
#'
#'
#' @param seu_obj The input Seurat object. Must contain hash assay named 'HTO' or 'hto' in dgCMatrix format.
#' Must contain Souporcell genotype assay(s) named GENO or geno with the desired k value (ex: 'GENO5'). 
#' Assumes Souporcell was run renaming all multiplet designations to 'multiplet' (see viewmastR
#' add_souporcell_seurat documentation for rename_assignments).
#'
#' @param csv The input hash-sample csv file path. The csv must contain Hash and Sample columns.
#'
#' @param soup_k The desired Souporcell run, indicating the number of genotypes detected. The appropriate 
#' number of genotypes expected for the data should be determined prior to running kmeansync. 
#'
#' @param output_col The name for the sample assignments column that will be added to the output Seurat object.
#' Default is 'FinalAssignment'.
#'
#' @param res Boolean indicating additional results output. Default is FALSE. If TRUE, Seurat object with sample 
#' and hash labels will be returned, along with a dataframe with False Discovery Rate (FDR) and a statistical significance
#' indicator (* when FDR < 0.1). This table shows statistical evidence that the assigned hash has an effect on the genotype. 
#' Also returns a dataframe linking genotypes to samples, and a graph of the average marginal effects for each genotype.
#'
#' @return seu_obj Returns Seurat object with sample and hash labels by default. Note that multiplets have been excluded.
#' If res=TRUE, additional dataframes and graphs are returned as well.
#'
#'
#'
#' @examples
#' \dontrun{
#' 
#'  output_list4 <- logisync(seu, csv='/path/to/hash_sampleABCD.csv', soup_k=4, res=TRUE)
#'  
#'  seu5 <- logisync(seu, csv='/path/to/hash_sampleDEF.csv', soup_k=5, output_col='Sample_Assignment') 
#'   
#' }
#'
#' @importFrom dplyr filter mutate case_when
#' @import PNWColors
#' @import factoextra
#' @import logistf
#' @importFrom stats glm p.adjust as.formula predict coef vcov pnorm
#' @importFrom Seurat NormalizeData ScaleData 
#' @importFrom SeuratObject DefaultAssay
#' @importFrom utils read.csv
#' @importFrom methods as
#' 
#' @export

logisync <- function(seu_obj, csv, soup_k, output_col='FinalAssignment', res=FALSE){
  
  FDR <- Significance <- df.dx <- Soup <- Hash <- NULL
  
  if(!inherits(seu_obj, 'Seurat')){
    stop('Input must be a seurat object.')}
  
  # ensure seu_obj has HTO assay
  if(!any(names(seu_obj@assays) %in% c('HTO', 'hto'))){
    stop('Input seurat object must contain a hash assay named HTO or hto.')}
  
  # check HTO matrix
  if(inherits(seu_obj@assays$HTO$counts) != "dgCMatrix"){
    stop('HTO assay formatted incorrectly in seurat object. HTO assay counts must be a dgCMatrix.')}
  
  # ensure seu_obj has GENO assay
  pattern <- paste0('GENO', soup_k)
  if(!any(grepl(pattern, names(seu_obj@assays), ignore.case = TRUE))){
    stop('Input seurat object must contain a genotype assay named GENO or geno with the desired k value (ex: GENO5).')}
  
  # read in csv
  hash_table <- read.csv(csv) 
  
  # set assay
  SeuratObject::DefaultAssay(seu_obj) <- "HTO"
  # remove multiplets (avoid noise when assigning hashes)
  geno_col <- paste0('geno', soup_k)
  seu_obj <- subset(seu_obj, cells=which(seu_obj@meta.data[[geno_col]] != 'Multiplet'))
  # normalize data 
  seu_obj <- NormalizeData(seu_obj, normalization.method = "CLR", margin = 1)
  # scale data 
  seu_obj <- ScaleData(seu_obj)
  
  
  geno_col <- paste0('geno', soup_k)
  data <- Seurat::FetchData(seu_obj, c(unique(hash_table$Hash), geno_col))
  
  # one-hot encoding
  encoded_genos <- list()
  data_encoded <- data
  for (genotype in data[[geno_col]]){  
    # add column for each genotype
    data_encoded[[paste0('Soup_', genotype)]] <- dplyr::case_when(data[[geno_col]] == genotype ~ 1, 
                                                                  .default = 0)
    # add to column names list
    encoded_genos <- unique(append(encoded_genos, list(paste0('Soup_', genotype))))
  }
  
  # Check for imbalanced classes
  classes <- table(data[[geno_col]])
  proportions <- lapply(classes, function(col){col/sum(classes)})
  names(proportions) <- paste0('Soup_', names(proportions))
  # threshold expected for balanced classes
  balanced_threshold <- 1 / length(encoded_genos)
  
  # test for HTO associated with genotype: logistic regression for each genotype
  regs <- list()
  margins <- list()
  stderr <- list()
  
  # extract hashes
  hash_form <- paste("~ ", paste(unique(hash_table$Hash), collapse = " + "), sep = "")
  # log reg/firth log reg per genotype
  for(col in encoded_genos){
    model_form <- as.formula(paste(factor(col), hash_form))
    
    # see if class proportion is less than expected from balanced classes
    if(proportions[[col]] < balanced_threshold){
      print(paste0('Firth Logistic Regression selected due to class imbalance for: ', col))
      reg <- logistf::logistf(model_form, data = data_encoded, control=logistf.control(maxit=250)) 
    }else{
      # logistic regression for balanced classes
      reg <- stats::glm(model_form, data = data_encoded, family='binomial')
    }
    # add to named list
    regs[[col]] <- reg
    # predicted probability
    predicted_probs <- predict(reg, type = "response")
    coefficients <- coef(reg)
    # calculate influence of each hash on genotype (AME)
    margins[[col]] <- coefficients[-1] * mean(predicted_probs * (1 - predicted_probs))
    # calculate standard error
    std_errors <- sqrt(diag(vcov(reg)))
    stderr[[col]] <- std_errors[-1] * mean(predicted_probs * (1 - predicted_probs))
  }
  
  # initialize output df
  output <- data.frame()
  marginal_effects_list <- data.frame()
  
  # iterate over named list
  for(soup in names(regs)){
    # interpret influence of each hash on genotype
    df_dx <- data.frame(df.dx = margins[[soup]], SE = stderr[[col]])
    # calculate z values, p values (two-tailed)
    df_dx$z <- df_dx$df.dx / df_dx$SE
    df_dx$P_value <- 2 * (1 - pnorm(abs(df_dx$z)))
    # adjust p value
    df_dx$FDR <- stats::p.adjust(df_dx$P_value, 'fdr')
    # identify hash (largest influence on genotype)
    margin_hash <- rownames(df_dx)[which.max(df_dx$df.dx)]
    sig <- ''
    # check statistical significance   
    if(df_dx[margin_hash, 'FDR'] < 0.1){
      sig <- '*'
    }
    
    output <- rbind(output, data.frame(Soup = soup, Hash = margin_hash, FDR = df_dx[margin_hash, 'FDR'], Significance=sig))
    
    # marginal effects for plotting
    df_dx$Soup <- soup
    df_dx$Hash <- rownames(df_dx)
    marginal_effects_list <- rbind(marginal_effects_list, df_dx)
  }
  
  # convert output genos to numeric
  output$Soup <- sapply(strsplit(output$Soup, "_"), "[[", 2)
  # merge csv and output
  pivoted <- merge(hash_table, output %>% dplyr::select(-FDR, -Significance), by='Hash')
  
  # add hash and sample assignment columns to data
  match_index <- match(data[[geno_col]], pivoted$Soup)
  data <- data %>%
    mutate(HashAssignment = ifelse(!is.na(match_index), pivoted$Hash[match_index], NA_character_),
      FinalAssignment = ifelse(!is.na(match_index), pivoted$Sample[match_index], NA_character_))
  
  # add results to seu metadata
  seu_obj[[output_col]] <- data$FinalAssignment
  seu_obj$HashAssignment <- data$HashAssignment
  
  # Plotting marginal effects: order by number
  marginal_effects_list <- marginal_effects_list %>%
    mutate(Soup = factor(Soup, levels = sort(unique(Soup), method = "radix")))
  marginal_effects_list$Hash <- factor(marginal_effects_list$Hash, 
                                       levels = unique(marginal_effects_list$Hash))
  
  # Plot the marginal effects
  pal <- PNWColors::pnw_palette("Sailboat", length(unique(hash_table$Hash)))
  bar <- ggplot(marginal_effects_list, aes(x = Soup, y = df.dx, fill = Hash)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Marginal Effects of Hashes on Genotype",
         x = "Soup", y = "Marginal Effect") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values=pal, name='Hash')  
  
  
  # return additional results tables
  if(res){
    # create final table linking soups to hashes
    final <- data.frame()
    for(hto in unique(output$Hash)){
      matched_soups <- output[output$Hash == hto, "Soup"]
      # add soups to string
      soups_per_hash <- paste(matched_soups, collapse = ", ")
      final <- rbind(final, data.frame(Hash=hto, Soup=soups_per_hash))
    }
    
    # merge final and input csv
    association <- merge(hash_table, final, by=('Hash'))
    
    return(list(seu_obj, output, association, bar))
    
    # return seu obj by default
  }else{
    return(seu_obj)
  }
}