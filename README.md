# genosyncR: Sync Samples to Genotypes via Explainable Machine Learning Techniques

Demultiplex pooled samples from hash data with confidence. Two methods are incorporated in the genosync function: kmeans combined with apriori association analysis, and logistic regression. The outputs of these methods are benchmarked, increasing confidence in sample assignment results when the two independent methods agree.

Multiple [Souporcell](https://github.com/wheaton5/souporcell) runs are iterated through to determine the number of genotypes that best fits the data. However, if the number of genotypes is known, it can be specified with the `max_soup_run` parameter. Souporcell multiplets are excluded from this analysis. 
    
### kmeansync    
The optimal number of clusters for kmeans is selected with the silhouette method. Hashes are assigned to kmeans clusters based on the average hash enrichment in each cluster. Apriori association rules analysis
links hashes from kmeans clusters to Souporcell genotypes. Then, genotypes are linked to samples via the input hash-sample csv. 
    
### logisync
Binary logistic regression is run for each genotype, and hashes are assigned based on which hash has the largest influence on the genotype (via average marginal effect). Genotypes are linked to samples via the input hash-sample csv. Binary Firth logistic regression will be run for genotypes with small proportions to account for rare cases and reduce bias.

## Installation
```
devtools::install_github('furlan-lab/genosyncR')
library(genosyncR)
```

Or clone the repository:
```
git clone https://github.com/furlan-lab/genosyncR.git
```

```
# install in R
setwd('path/to/cloned/repository')
devtools::install('genosyncR')
library(genosyncR)
```

## Example
```
results <- genosync(seu_obj=seu_ABC, hash_csv='hash_ABC.csv', max_soup_run=6)
results
```

**Inputs:**
* Seurat object with Souporcell and HTO assays (see [viewmastR](https://furlan-lab.github.io/viewmastR/reference/add_souporcell_seurat.html)), subset to desired samples
* Hash-Sample csv

![csv](https://github.com/user-attachments/assets/0196f893-3172-4552-8a08-5c728eb2e59a)
* Optional: Number of Souporcell runs to iterate through


**Outputs:**
* Seurat object, labeled with Samples and Hashes
* Consensus results of kmeansync and logisync
     * If no runs agree, all results for all Souporcell runs



![genosync_res](https://github.com/user-attachments/assets/34178109-0cc3-426e-ae19-5e12a78e0c9a)

**A.** The average marginal effects show the influence each hash has on each Souporcell genotype.

**B.** Dataframe of statistical evidence that the assigned hash has an effect on the genotype. The statistical significance indicator * appears when the False Discovery Rate < 0.1.

**C.** The optimal number of clusters for kmeans, selected via the silhouette method.

**D.** A UMAP of kmeans clustering, shaded with average hash enrichment per cluster. Dimensionality of data is reduced via UMAP prior to kmeans running.

**E.** Apriori association rules, linking kmeans clusters to genotypes. High confidence and lift demonstrate stronger associations. Note that association rules may not be found for all genotypes or kmeans clusters, indicating a poor fit of the Souporcell run.

**F.** Dataframe linking kmeans clusters to Souporcell genotypes.

**G.** Dataframe linking samples to genotypes. Consensus between kmeansync and logisync results was found for Souporcell = 3.





