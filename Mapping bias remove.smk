
sample="" # sample id

# 1. extract and filter individual imputed genotype from population-level imputed genotype data 
rule extract/bulit individual genotype vcf:
    input:
        pop.vcf="RNA_SNPs.7008Samples.chr1_18.imputed.vcf.gz", # population-level imputed genotype data, from Farm-GTEx.
        reference="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        out_filter="selected.{sample}.filtered.imputed.vcf.gz" # Please provide the output files names from this step.
    threads: 2
    shell:
        '''
        bcftools view -s {sample} -o {sample}.raw.imputed.vcf.gz -O z {input.pop.vcf}
        zcat {sample}.raw.imputed.vcf.gz |grep -v "0|0" >> {sample}.imputed.vcf.gz
        tabix -p vcf {sample}.imputed.vcf.gz
        gatk --java-options "-Djava.io.tmpdir=/tmp" VariantFiltration -R {input.reference} -V {sample}.vcf.gz -O {sample}.filtered.imputed.vcf.gz -window 10 -cluster 3 
        gatk --java-options "-Djava.io.tmpdir=/tmp" SelectVariants -R {input.reference} -V {output.out_filter} --select-type-to-include SNP --restrict-alleles-to BIALLELIC --exclude-filtered -O {output.out_filter}

# 2.  mapping bias remove (Star --waspoutmode)
rule map:
    input:
        STAR_dir="", # Path to STAR
        gtf="Sus_scrofa.Sscrofa11.1.104.gtf",
    output:
        "{sample}-STARAligned.sortedByCoord.out.bam"
    threads: 5
    shell:
        '''
        cd {sample}/
        
        if [ -e *_2.fastq.gz ]  # Paired-end ???
        then
            ${input.STAR_dir}/STAR --runThreadN {threads} \
                --genomeDir star-genome  \
                --sjdbGTFfile {input.gtf} \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMmapqUnique 255  \
                --readFilesCommand zcat \
                --outFilterMismatchNmax 3 \
                --waspOutputMode SAMtag \
                --varVCFfile selected.{sample}.filtered.imputed.vcf.gz \
                --readFilesIn {sample}_1.clean.fq.gz {sample}_2.clean.fq.gz \
                --outFileNamePrefix {sample}-STAR
        
        else    # Single-ended sequencing ???
            ${STAR_dir}/STAR --runThreadN {threads} \
                --genomeDir  star-genome  \
                --sjdbGTFfile {input.gtf} \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMmapqUnique 255  \
                --readFilesCommand zcat \
                --outFilterMismatchNmax 3 \
                --varVCFfile selected.{sample}.filtered.imputed.vcf.gz \
                --readFilesIn {sample}.clean.fq.gz \
                --outFileNamePrefix {sample}-STAR
        fi
        '''
# 3. select the unbias bam (samtools version > )
rule : select Non-bias unique mapping bam
    input:
        bam="{sample}-STARAligned.sortedByCoord.out.bam ", 
        reference="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        outbam="dedupped_{sample}.uniq.WASP.bam" # Please provide the output files names from this step.
    shell:
        ```
        samtools index {input.sample}
        samtools view -bq 225 {input.sample} > {sample}.STARAligned.uniq.bam
        samtools index {sample}.STARAligned.uniq.bam
        samtools view -d vW:1 -o {sample}.STARAligned.uniq.wasp.bam
        samtools index {sample}.STARAligned.uniq.wasp.bam
        gatk AddOrReplaceReadGroups -I {sample}.STARAligned.uniq.wasp.bam -O rg_added_{sample}.STARAligned.uniq.wasp.bam -RGID 4 -RGLB lib1 -RGPL illumina -RGPU run -RGSM 20 -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate
        gatk MarkDuplicates -I rg_added_{sample}.STARAligned.uniq.wasp.bam -O {output.outbam} -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT --READ_NAME_REGEX null -M dedupped_{sample}.marked_dup_metrics.txt
        ```
