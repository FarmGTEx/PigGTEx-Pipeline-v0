# wget http://ftp.ensembl.org/pub/release-100/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.100.gtf.gz
sample = ""  # sample id

# 1. Ancestry estimation (plink v1.90, admixture v1.3.0)
rule ancestry_est:
    shell:
        '''        
        plink --bfile  RNA_7008_SNP_5000 --update-ids rna2.txt --out RNA_7008_SNP_5000 --make-bed
        plink --bfile RNA_7008_SNP_5000 --recode --out newRNA_7008_SNP_5000 --make-bed
        for i in 3 4 5 
        do
           admixture -j24 merge-dcdp.bed $i | tee log $i.out
        done
        '''

# 2. Predict breed labels using RandomForest
rule pred_breed:
    script:
        "Breed_prediction_RF.R"
