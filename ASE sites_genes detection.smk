paired_sample= ["SAMEA103939830", "SAMEA4436905"] # example sample id
single_sample = ["aa","bb"]
haplo_count_blacklit = " " 
# BED file containing genomic intervals to be excluded from haplotypic counts. Reads from any variants which lie within these regions will not be counted for haplotypic counts.
phASER_dir = "path/to/phASER", # Path to STAR,
VCF = "results/{sample}/selected.{sample}.filtered.imputed.vcf.gz",
BAM = "results/{sample}/dedupped_{sample}.uniq.WASP.bam"
rule all:
    input:
        expand("results/{sample}/dedupped_{sample}.uniq.WASP.bam", sample=PAIRED_SAMPLES + SINGLE_SAMPLES)
        
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
                    --features Sus_scrofa.Sscrofa11.1.100.gene.phaser.bed \
                    --o {sample}_phaser.gene_ae.txt
        
