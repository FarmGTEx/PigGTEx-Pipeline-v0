# 1. Cis-heritability estimation
rule cis_heri:
    input:
        CHR=1,
        TIS="Adipose"
    shell:
        '''        
        julia Cis_heritability.jl {input.CHR} {input.TIS} 
        '''
