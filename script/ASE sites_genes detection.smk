PAIRED_SAMPLES = ["SAMEA103939830", "SAMEA4436905"] # example sample id
SINGLE_SAMPLES = ["aa","bb"]
SAMPLES = PAIRED_SAMPLES + SINGLE_SAMPLES

PHASER = "/path/to/phASER/phaser.py"
PHASER_GENE_AE = "/path/to/phaser_gene_ae.py"
BINOM_R = "binomial.ASE.r"

FEATURES = "Sus_scrofa.Sscrofa11.1.100.gene.phaser.bed"
HAPLO_COUNT_BLACKLIST = ""   # keep "" if not used
# BED file containing genomic intervals to be excluded from haplotypic counts. Reads from any variants which lie within these regions will not be counted for haplotypic counts.

rule all:
    input:
        expand("results/{sample}/{sample}.RNAim.ASE.sites", sample=SAMPLES),
        "results/phaser_expr.7008.done"
        
# 1. ASE sites detection using phaser (python2 is needed)
rule paired_phASER_ASE:
  input:
      vcf = VCF
      bam = BAM
  output:
        ""
    threads: 5
    shell:
      ```
      python STAR/phaser.py \
        --vcf {input.VCF} \
        --bam {input.BAM} \
        --paired_end 1 \
        --mapq 255 \
        --baseq 10 \
        --sample {sample} \
        --threads 5 \
        --temp_dir /tmp \
        --haplo_count_blacklist {wildcards.haplo_count_blacklit} \
        --gw_phase_vcf 1 \
        --o {sample} \  ## the output nam
        --output_read_ids 1
      ```
rule paired_phASER_ASE:
    input:
      vcf = VCF
      bam = BAM
  output:
        ""
    threads: 5
    shell:
      ```
      python STAR/phaser.py \
        --vcf {input.VCF} \
        --bam {input.BAM} \
        --paired_end 0 \
        --mapq 255 \
        --baseq 10 \
        --sample {sample} \
        --threads 5 \
        --temp_dir /tmp \
        --haplo_count_blacklist {wildcards.haplo_count_blacklit} \
        --gw_phase_vcf 1 \
        --o {sample} \  ## the output nam
        --output_read_ids 1
      ```
# 2. Produce gene level haplotype counts for allelic expression studies
    # phASER Gene AE. Uses output from phASER to produce gene level haplotype counts for allelic expression studies. It does this by summing reads from both single variants and phASER haplotype blocks using their phase for each gene. 
    # Developed by [Stephane E. Castel](mailto:scastel@nygenome.org) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.
    rule gene_level_hap_counts:
        input:
            haplotype_count = "results/{sample}/{sample}.haplotypic_counts.txt"
            features= "Sus_scrofa.Sscrofa11.1.100.gene.phaser.bed"
            #features: File in BED format (0 BASED COORDINATES - chr,start,stop,name) containing the features to produce counts for.
            shell:
                ```
                python2 /phaser-master/phaser_gene_ae/phaser_gene_ae.py \
                    --haplotypic_counts {input.haplotype_count} \
                    --features {input.features} \
                    --o {sample}_phaser.gene_ae.txt
                ```
# 3. Aggregates gene-level haplotypic expression measurement files
    rule aggr_gene_level_expr:
        shell:
            '''
            #phASER-pop
            #Aggregates gene-level haplotypic expression measurement files across samples to produce a single haplotypic expression matrix, where each row is a gene and each column is a sample
            python2 /phaser-master/phaser_pop/phaser_expr_matrix.py \
                --gene_ae_dir path/{sample}_phaser.gene_ae.txt\ # the path contains all sample's phaser.gene_ae.txt
                --features Sus_scrofa.Sscrofa11.1.100.gene.phaser.bed \
                --t 5 \
            ```
            --o phaser_expr.7008
        '''
# 4. filter
    rule filter_ASE_sites
        input:
            raw.ASE = "results/{sample}/{sample}.allelic_counts.txt
        output: 
            filter.ASE = "result/filter/{sample}.filter.sites"
        shell:
            ```
            awk 'NR >=2 && $8>=15 && $6>2 && $7 >2 && $6/$8> 0.02 && $6/$8 < 0.98 {print $1,$2,$4,$5,$6,$7,$8,$6/$8}' OFS="\t" results/{sample}/{sample}.allelic_counts.txt > {output.filter.ASE}
            ```
# 5. binom test
    Rscript binomal.ASE.r 
        
