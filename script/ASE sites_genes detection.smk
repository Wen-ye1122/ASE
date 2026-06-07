</> Python
import pandas as pd

configfile: "config/config.yaml"

sample_info = pd.read_csv("config/sample_info.tsv", sep="\t").set_index("sample")
SAMPLES = sample_info.index.tolist()


def get_paired_end(wc):
    layout = sample_info.loc[wc.sample, "layout"]
    return 1 if layout == "PE" else 0


# BED file containing genomic intervals to be excluded from haplotypic counts. Reads from any variants which lie within these regions will not be counted for haplotypic counts.

rule all:
    input:
        expand("results/{sample}/{sample}.RNAim.ASE.sites", sample=SAMPLES),
        "results/phaser_expr.7008.done"
        
# 1. ASE sites detection using phaser (python2 is needed)
rule phASER_ASE:
  input:
    vcf = "results/{sample}/selected.{sample}.filtered.imputed.vcf.gz",
    bam = "results/{sample}/dedupped_{sample}.uniq.WASP.bam"
  output:
    hap = "results/{sample}/{sample}.haplotypic_counts.txt",
    allelic = "results/{sample}/{sample}.allelic_counts.txt"
    params:
        outprefix = "results/{sample}/{sample}",
        paired = get_paired_end,
        blacklist = lambda wc: (
            f"--haplo_count_blacklist {config['haplo_count_blacklist']}"
            if config["haplo_count_blacklist"] else ""
        )
    threads: 
        config["threads_phaser"]
    shell:
      ```
        mkdir -p results/{wildcards.sample}
        python STAR/phaser.py \
            --vcf {input.VCF} \
            --bam {input.BAM} \
            --paired_end {params.paired} \
            --mapq 255 \
            --baseq 10 \
            --sample {wildcards.sample} \
            --threads {threads} \
            --temp_dir /tmp \
            {params.blacklist} \
            --gw_phase_vcf 1 \
            --o {params.outprefix} \ 
            --output_read_ids 1
      ```

# 2. Produce gene level haplotype counts for allelic expression studies
    # phASER Gene AE. Uses output from phASER to produce gene level haplotype counts for allelic expression studies. It does this by summing reads from both single variants and phASER haplotype blocks using their phase for each gene. 
    # Developed by [Stephane E. Castel](mailto:scastel@nygenome.org) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.
rule gene_level_hap_counts:
    input:
        haplotype_count = "results/{sample}/{sample}.haplotypic_counts.txt"
        features = config["features"]
    output:
        gene_ae = "results/{sample}/{sample}_phaser.gene_ae.txt" 
    shell:
        ```
        python2 {config[phaser_gene_ae]} \
            --haplotypic_counts {input.haplotype_count} \
            --features {input.features} \
            --o {sample}_phaser.gene_ae.txt
        ```
# 3. Aggregates gene-level haplotypic expression measurement files
    rule aggr_gene_level_expr:
        input:
            expand("results/{sample}/{sample}_phaser.gene_ae.txt", sample=SAMPLES)
        output: 
            touch("results/phaser_expr.7008.done")
    shell:
        '''
        mkdir -p results/gene_ae_all
        cp results/*/*_phaser.gene_ae.txt results/gene_ae_all/
        #phASER-pop
        #Aggregates gene-level haplotypic expression measurement files across samples to produce a single haplotypic expression matrix, where each row is a gene and each column is a sample
        python2 {config[phaser_expr_matrix]} \
            --gene_ae_dir results/gene_ae_all\ # the path contains all sample's phaser.gene_ae.txt
            --features {config[features]} \
            --t 5 \
            --o results/phaser_expr.7008
        '''
# 4. filter
    rule filter_ASE_sites
        input:
            raw.ASE = "results/{sample}/{sample}.allelic_counts.txt"
        output: 
            filter.ASE = "results/{sample}/{sample}.RNAim.tested.sites"
        shell:
            ```
            awk 'BEGIN{{OFS="\\t";
                print "chr","pos","col4","col5","refCount","altCount","totalCount","ref_ratio"}}
                'NR >=2 && $8>=15 && $6>2 && $7 >2 && $6/$8> 0.02 && $6/$8 < 0.98 
                {{print $1,$2,$4,$5,$6,$7,$8,$6/$8}}' {input.allelic} > {output.filtered}
            ```
# 5. binom test
    rule binomial_ase_test:
    input:
        filtered = "results/{sample}/{sample}.RNAim.tested.sites",
        script = "scripts/binomial.ASE.r"
    output:
        ase = "results/{sample}/{sample}.RNAim.ASE.sites"
    shell:
        ```
        Rscript {input.script} {input.filtered} {output.ase}
        ```
