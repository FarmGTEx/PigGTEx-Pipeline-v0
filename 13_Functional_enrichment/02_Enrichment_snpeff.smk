# 1. snpeff enrichment with TORUS for cis-QTLs
rule run_enrichment:
    input:
        sample = "name",
        tissue = "Adipose"
    output:
        "*.enrichment.est"
    shell:
        '''        
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_3utr.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_3utr.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_5utr.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_5utr.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_downstream.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_downstream.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_intron.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_intron.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_missense.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_missense.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_nc.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_nc.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_splice_acceptor.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_splice_acceptor.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_splice_donor.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_splice_donor.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_splice_region.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_splice_region.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_stop_gained.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_stop_gained.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_synonymous.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_synonymous.enrichment.est
        torus -d {input.sample}_{input.tissue}.cis.nominal.torus.all.txt.gz -annot {input.sample}_{input.tissue}.cis.nominal.snpeff_upstream.snp.txt.gz -est > {input.sample}_{input.tissue}.snpeff_upstream.enrichment.est
        '''
