# 1. Calculate eight distinct matrices to measure the tissue specificity
rule task1:
    shell:
        '''        
        tspex --log tissue_median.csv tissue_median_tau.tsv tau
        tspex --log tissue_median.csv tissue_median_counts.tsv counts
        tspex --log tissue_median.csv tissue_median_gini.tsv gini
        tspex --log tissue_median.csv tissue_median_simpson.tsv simpson
        tspex --log tissue_median.csv tissue_median_shannon_specificity.tsv shannon_specificity
        tspex --log tissue_median.csv tissue_median_roku_specificity.tsv roku_specificity
        tspex --log tissue_median.csv tissue_median_spm_dpm.tsv spm_dpm
        tspex --log tissue_median.csv tissue_median_js_specificity_dpm.tsv js_specificity_dpm
        '''

# 2. Detect tissue-specific genes
rule task2:
    script:
        "Detect_tissue-specific_genes.R"


