# 1. Extend the top 1000 tissue specific gene's TSS to up and down 20K
rule task1:
    input:
        wkdir="~/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_pig_gtex/top1000/",
        tss_bed="~/genome/pig/fold_enrich_region/TSS_esemble100_colin.bed"
    shell:
        '''        
        cd {input.wkdir}
        mkdir -p "../top1000_20K/"
        ls *.txt | while read id;
        do
            # rm "../top1000_20K/"${{id%%.*}}.bed
            cat $id | while read gene
            do
                grep -w $gene {input.tss_bed} | awk '{{if ($2>=20000) print $1, $2-20000, $3+19999, $4, $5, $6, $7}}' | sed 's/ /\t/g'>> "../top1000_20K/"${{id%%.*}}.bed
            done
        done
        '''

# 2. Enrichment for top1000 tissues specific gene in regulator
rule run_enrichment:
    input:
        wkdir="~/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify",
        dir_top1000_20K="~/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_pig_gtex/top1000_20K/"
    output:
        outfile="fold_enrichment_region_Ajust_TSS.csv"
    shell:
        '''
        cd {input.wkdir}
        echo OUTPUTSAMPLE_pig_five_maker_15 > fold_enrichment_region_Ajust_TSS.txt
        ls *segments.bed | while read id
        do
          ls {input.dir_top1000_20K}/*.bed | while read bedfile
          do
             echo $bedfile
             ls -lh $bedfile
             echo state $id ${{bedfile##*/}} $(basename $bedfile .bed) >> fold_enrichment_region_Ajust_TSS.txt
             for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
             do
             A=$(grep -w $i $id | awk '{{print ($3-$2) }}' |awk '{{sum+=$1}} END {{print sum}}')
             B=$(cat $bedfile | awk '{{print ($3-$2) }}' |awk '{{sum+=$1}} END {{print sum}}')
             C=$(bedtools intersect  -a <(grep -w $i $id) -b $bedfile |awk '{{print ($3-$2) }}' |awk '{{sum+=$1}} END {{print sum}}')
             D=2478444698
             if [ "$C" = "" ]; 
             then 
             let C=0
             echo $C
             fi
             echo $i $id ${{bedfile##*/}} $(bc <<< "scale=10;($C/$A)/($B/$D)") >> fold_enrichment_region_Ajust_TSS.txt
             done
          done
        done
        sed 's/ /\t/g' fold_enrichment_region_Ajust_TSS.txt > {output.outfile}
        '''


# 3. Summary the enrichment
rule summaryRes:
    input:
        dir_top1000_20K="~/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_pig_gtex/top1000_20K/"
    shell:
        '''
        touch fold_enrichment.txt
        ls *segments.bed | while read id
        do
            echo $id
            grep -w $id fold_enrichment_region_Ajust_TSS.csv > 1.txt
            ls {input.dir_top1000_20K}/*.bed | while read bedfile
            do
            grep -w ${{bedfile##*/}} 1.txt | cut -f 4 | paste -s >> 2.txt
            done
            awk '{{for(i=1;i<=NF;i++)a[NR,i]=$i}}END{{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}}' 2.txt > $(basename $id ".bed")_fold_enrich1.txt
            rm 2.txt
            paste state.txt $(basename $id ".bed")_fold_enrich1.txt |sed 's/ /\t/g' > "enrichment/"$(basename $id ".bed")_fold_enrich_TSE1000_20K.txt
            rm $(basename $id ".bed")_fold_enrich1.txt
        done
        '''

# 4. Add col name and row name
rule add_names:
    input:
        wkdir="~/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/enrichment", # from the output of the previous step
        gene_gtf="~/data/Sus_scrofa.Sscrofa11.1.100.simple.gtf.gz"
    output:
        "fold_enrich_TSE1000_20K_combine.csv"
    shell:
        '''
        cd {input.wkdir}
        rm 2.csv
        ls *TSE1000_20K.txt| while read id;
        do
            echo $id
            b=$(basename $id "_fold_signature_2.txt")
            c=${{id%%_*}}
            echo $b
            echo $c
            echo "cat "${{id}}"| sed '1d'| awk '\$(NF+1)=\""$c"\"' >> "2.csv"" > ${{i}}.sh; bash ${{i}}.sh
            rm ${{i}}.sh
        done
        sed -n '1p' Adipose_15_segments_fold_enrich_TSE1000_20K.txt   >1.txt
        echo tissue >2.txt
        paste 1.txt 2.txt > 3.txt
        sed -i 's/ /\t/g' 2.csv
        cat 3.txt 2.csv > {output}
        wait
        '''
