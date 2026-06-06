
paired_sample= ["SAMEA103939830", "SAMEA4436905"] # example sample id
single_sample = ["aa","bb"]

POP_VCF = "RNA_SNPs.7008Samples.chr1_18.imputed.vcf.gz"
REF = "Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
GTF = "Sus_scrofa.Sscrofa11.1.104.gtf"
STAR = "/path/to/STAR"
STAR_GENOME = "star-genome"

rule all:
    input:
        expand("results/{sample}/dedupped_{sample}.uniq.WASP.bam", sample=PAIRED_SAMPLES + SINGLE_SAMPLES)

# 1. extract and filter individual imputed genotype from population-level imputed genotype data 
rule extract_individual_imputed_vcf::
    input:
        pop_vcf = POP_VCF,
        reference = REF
    output:
        raw = "results/{sample}/{sample}.raw.imputed.vcf.gz",
        het = "results/{sample}/{sample}.het.imputed.vcf.gz",
        filt = "results/{sample}/{sample}.filtered.imputed.vcf.gz"
        select = "results/{sample}/selected.{sample}.filtered.imputed.vcf.gz"
    threads: 2
    shell:
        '''
        mkdir -p results/{wildcards.sample}
        bcftools view -s {wildcards.sample} -Oz -o {output.raw} {input.pop_vcf}
        tabix -p vcf {output.raw}
        bcftools view -i 'GT!="0/0" && GT!="0|0"' -Oz -o {output.het} {output.raw}
        tabix -p vcf {output.het}
        gatk --java-options "-Djava.io.tmpdir=/tmp" VariantFiltration -R {input.reference} -V {output.het} -O {output.filt} -window 10 -cluster 3 
        gatk --java-options "-Djava.io.tmpdir=/tmp" SelectVariants -R {input.reference} -V {output.filt} --select-type-to-include SNP --restrict-alleles-to BIALLELIC --exclude-filtered -O {output.select}
        tabix -p vcf {output.select}
        ```

# 2.  mapping bias remove (Star --waspoutmode)
rule star_wasp_map_paired:
    input:
        vcf="results/{sample}/selected.{sample}.filtered.imputed.vcf.gz", 
        fq1 = "fastq/{sample}_1.clean.fq.gz",
        fq2 = "fastq/{sample}_2.clean.fq.gz",
        gtf = GTF
    output:
        bam = "results/{sample}/{sample}-STARAligned.sortedByCoord.out.bam"
    threads: 5
    shell:
        '''
        mkdir -p results/{wildcards.sample}
        {STAR} --runThreadN {threads} \
            --genomeDir {STAR_GENOME}  \
            --sjdbGTFfile {input.gtf} \
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMmapqUnique 255  \
            --readFilesCommand zcat \
            --outFilterMismatchNmax 3 \
            --waspOutputMode SAMtag \
            --varVCFfile {input.vcf} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix results/{wildcards.sample}/{wildcards.sample}-STAR
        
rule star_wasp_map_single:
    input:
        vcf = "results/{sample}/selected.{sample}.filtered.imputed.vcf.gz",
        fq = "fastq/{sample}.clean.fq.gz",
        gtf = GTF
    output:
        bam = "results/{sample}/{sample}-STARAligned.sortedByCoord.out.bam"
    threads: 5
    shell:
        ```
        mkdir -p results/{wildcards.sample}
        {STAR} --runThreadN {threads} \
            --genomeDir  {STAR_GENOME}  \
            --sjdbGTFfile {input.gtf} \
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMmapqUnique 255  \
            --readFilesCommand zcat \
            --outFilterMismatchNmax 3 \
            --varVCFfile {input.vcf} \
            --readFilesIn {input.fq} \
            --outFileNamePrefix --readFilesIn {input.fq}
        '''
# 3. select the unbias bam (samtools version > )
rule select_nonbias_unique_bam:
    input:
        bam="results/{sample}/{sample}-STARAligned.sortedByCoord.out.bam", 
        reference=REF
    output:
        outbam = "results/{sample}/dedupped_{sample}.uniq.WASP.bam",
        bai = "results/{sample}/dedupped_{sample}.uniq.WASP.bai"
    shell:
        ```
        samtools index {input.bam}
        samtools view -bq 225 {input.bam} > results/{wildcards.sample}/{wildcards.sample}.STARAligned.uniq.bam
        samtools index results/{wildcards.sample}/{wildcards.sample}.STARAligned.uniq.bam
        samtools view -b -d vW:1 results/{wildcards.sample}/{wildcards.sample}.STARAligned.uniq.bam -o results/{wildcards.sample}/{wildcards.sample}.STARAligned.uniq.wasp.bam
        samtools index results/{wildcards.sample}/{wildcards.sample}.STARAligned.uniq.wasp.bam
        gatk AddOrReplaceReadGroups -I results/{wildcards.sample}/{wildcards.sample}.STARAligned.uniq.wasp.bam -O results/{wildcards.sample}/rg_added_{wildcards.sample}.STARAligned.uniq.wasp.bam -RGID {wildcards.sample} -RGLB lib1 -RGPL illumina -RGPU run -RGSM {wildcards.sample} -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate
        gatk MarkDuplicates -I results/{wildcards.sample}/rg_added_{wildcards.sample}.STARAligned.uniq.wasp.bam -O {output.outbam} -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT --READ_NAME_REGEX null -M results/{wildcards.sample}/dedupped_{wildcards.sample}.marked_dup_metrics.txt
        ```
