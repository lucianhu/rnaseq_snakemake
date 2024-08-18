import pandas as pd

# Load the config file
configfile: "config/config.yaml"

# Load the sample sheet from the config file
samplesheet = pd.read_csv(config["samples"])

# Access the trimming setting
trimming_active = config["trimming"]["activate"]

# Create a new column that combines 'sample', 'unit', and 'condition'
samplesheet['sample_name'] = samplesheet['sample'] + "_" + samplesheet['unit'] + "_" + samplesheet['condition']

# Ensure that 'sample_name' column is properly recognized
sample_list = samplesheet['sample_name'].tolist()

def is_paired_end(sample_name):
    # Select the relevant row(s) for the given sample_name
    sample_units = samplesheet[samplesheet['sample_name'] == sample_name]
    
    # Check if 'fq2' is null (i.e., missing) for all the rows corresponding to this sample
    fq2_null = sample_units["fq2"].isnull()
    
    # Determine if this sample is paired-end
    paired = ~fq2_null
    
    # Check if all associated rows are either paired-end or single-end
    all_paired = paired.all()
    all_single = (~paired).all()
    
    # Assert that the sample is either completely paired-end or single-end
    assert (
        all_single or all_paired
    ), f"Invalid units for sample {sample_name}, must be all paired-end or all single-end."
    
    return all_paired

def get_paired_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.sample].iloc[0]
    if is_paired_end(wildcards.sample):
        return [sample_data['fq1'], sample_data['fq2']]  # Return a list of file paths
    else:
        raise ValueError(f"Sample {wildcards.sample} is not paired-end.")

def get_single_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.sample].iloc[0]
    if not is_paired_end(wildcards.sample):
        return sample_data['fq1']  # Return the file path as a string
    else:
        raise ValueError(f"Sample {wildcards.sample} is paired-end.")

def get_paired_trimmed_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.sample].iloc[0]
    if is_paired_end(wildcards.sample):
        if config["trimming"]["activate"]:
            return [
                f"results/trimmed/{wildcards.sample}_trimmed_R1.fq.gz",
                f"results/trimmed/{wildcards.sample}_trimmed_R2.fq.gz"
            ]
        else:
            return [sample_data['fq1'], sample_data['fq2']]
    else:
        raise ValueError(f"Sample {wildcards.sample} is not paired-end.")

def get_single_trimmed_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.sample].iloc[0]
    if not is_paired_end(wildcards.sample):
        if config["trimming"]["activate"]:
            return [
                f"results/trimmed/{wildcards.sample}_trimmed_R1.fq.gz"
            ]
        else:
            return [sample_data['fq1']]
    else:
        raise ValueError(f"Sample {wildcards.sample} is paired-end.")

# Split the samples into paired-end and single-end
paired_end_samples = [sample for sample in sample_list if is_paired_end(sample)]
single_end_samples = [sample for sample in sample_list if not is_paired_end(sample)]

# Add this rule to define the final target (all) for Snakemake
rule all:
    input:
        # Directories and files related to quality control and trimming results
        expand("results/fastqc/{sample}", sample=paired_end_samples + single_end_samples),
        expand("results/trimmed/{sample}", sample=paired_end_samples + single_end_samples),
        
        # Expanding file patterns for original BAM, sorted BAM and related statistics
        expand("results/align/bam_original/{sample}/{sample}.Aligned.out.bam", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_original/{sample}/{sample}.Aligned.toTranscriptome.out.bam", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam.bai", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam.stats", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam.flagstat", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam.idxstats", sample=paired_end_samples + single_end_samples),

        # Expanding file patterns for BAM markup results
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.MarkDuplicates.metrics.txt", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.bai", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.stats", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.flagstat", sample=paired_end_samples + single_end_samples),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.idxstats", sample=paired_end_samples + single_end_samples),

        # Expanding file patterns for StringTie results and related files
        expand("results/stringtie/{sample}/{sample}.transcripts.gtf", sample=paired_end_samples + single_end_samples),
        expand("results/stringtie/{sample}/{sample}.gene.abundance.txt", sample=paired_end_samples + single_end_samples),
        expand("results/stringtie/{sample}/{sample}.coverage.gtf", sample=paired_end_samples + single_end_samples),
        expand("results/stringtie/{sample}/{sample}.ballgown", sample=paired_end_samples + single_end_samples),

        # Expanding file patterns for feature count results
        expand("results/feature_counts/{sample}/{sample}.featurecounts.txt", sample=paired_end_samples + single_end_samples),
        expand("results/feature_counts/{sample}/{sample}.featurecounts.txt.summary", sample=paired_end_samples + single_end_samples),
        expand("results/feature_counts/{sample}/{sample}.biotype_counts_rrna_mqc.tsv", sample=paired_end_samples + single_end_samples),

        expand("results/qualimap/{sample}", sample=paired_end_samples + single_end_samples),
        expand("results/dupradar/{sample}", sample=paired_end_samples + single_end_samples),

        expand("results/rseqc/bam_stat/{sample}/{sample}.bam_stat.txt", sample=paired_end_samples + single_end_samples),
        #expand("results/rseqc/inner_distance/{sample}/{sample}.inner_distance_plot.pdf", sample=paired_end_samples + single_end_samples),
        expand("results/rseqc/infer_experiment/{sample}/{sample}.infer_experiment.txt", sample=paired_end_samples + single_end_samples),
        expand("results/junction_annotation/{sample}/{sample}.junction.bed", sample=paired_end_samples + single_end_samples),
        #expand("results/junction_saturation/{sample}/{sample}.junction_saturation.pdf", sample=paired_end_samples + single_end_samples),
        #expand("results/read_distribution/{sample}/{sample}.read_distribution.txt", sample=paired_end_samples + single_end_samples),
        #expand("results/read_duplication/{sample}/{sample}.read_duplication.pdf", sample=paired_end_samples + single_end_samples)

        
# Rule for Preparing the Reference Genome
rule import_genome:
    input:
        config["ref"]["genome"]
    output:
        "resources/genome.fa"
    shell:
        "cp {input} {output}"

rule import_annotation:
    input:
        config["ref"]["annotation"]
    output:
        "resources/annotation.gtf"
    shell:
        "cp {input} {output}"

rule filter_gtf:
    input:
        gtf="resources/annotation.gtf",
        fasta="resources/genome.fa"
    output:
        gtf_filtered="resources/annotation.filtered.gtf"
    message: 
        "Refining the input GTF"
    log:
        "logs/filter_gtf.log"
    conda:
        "envs/python.yaml"
    shell:
        """
        python scripts/filter_gtf.py --gtf {input.gtf} --fasta {input.fasta} --prefix resources/annotation 2>> {log}
        """

rule gtf_to_bed:
    input:
        gtf="resources/annotation.gtf"
    output:
        bed="resources/annotation.filtered.bed"
    message:
        "Converting a GTF file to a BED file format."
    log:
        "logs/gtf_to_bed.log"
    conda:
        "envs/perl.yaml"
    shell:
        "scripts/gtf2bed {input.gtf} > {output}"

rule rsem_prepare_reference:
    input:
        gtf="resources/annotation.filtered.gtf",
        fasta="resources/genome.fa"
    output:
        "resources/genome.transcripts.fa"
    message:
        "Preparing a reference transcript sequence for RNA-Seq quantification"
    params:
        prefix="rsem_genome"
    log:
        "logs/rsem_prepare_reference.log"
    conda:
        "envs/rsem.yaml"
    threads: 12
    shell:
        """
        rsem-prepare-reference \\
            --gtf {input.gtf} \\
            --num-threads {threads} \\
            {input.fasta} \\
            {params.prefix} 2>> {log} &&
        cp {params.prefix}.transcripts.fa {output}
        """

rule samtools_faidx_sizes:
    input:
        fasta="resources/genome.fa"
    output:
        sizes="resources/genome.fa.sizes"
    message:
        "Listing the sizes of all sequences in the fasta file"
    log:
        "logs/samtools_faidx_sizes.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools faidx {input.fasta} 2>> {log}
        cut -f 1,2 {input.fasta}.fai > {output}
        """

rule star_genome_generate:
    input:
        gtf="resources/annotation.filtered.gtf",
        fasta="resources/genome.fa"
    output:
        directory("resources/star_index")
    message:
        "Creating a reference genome index"
    log:
        "logs/star_genome_generate.log"
    conda:
        "envs/star.yaml"
    threads: 12
    shell:
        """
        samtools faidx {input.fasta} 2>> {log}
        NUM_BASES=$(gawk '{{sum = sum + $2}}END{{if ((log(sum)/log(2))/2 - 1 > 14) {{printf "%.0f", 14}} else {{printf "%.0f", (log(sum)/log(2))/2 - 1}}}}' {input.fasta}.fai)
        STAR --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} --runThreadN {threads} --genomeSAindexNbases $NUM_BASES --limitGenomeGenerateRAM 77209411328 &>> {log}
        """

# Rule for Read QC
rule fastqc:
    input:
        fq1=lambda wildcards: samplesheet.loc[samplesheet['sample_name'] == wildcards.sample, 'fq1'].values[0],
        fq2=lambda wildcards: [] if pd.isna(samplesheet.loc[samplesheet['sample_name'] == wildcards.sample, 'fq2'].values[0]) else samplesheet.loc[samplesheet['sample_name'] == wildcards.sample, 'fq2'].values[0]
    output:
        fastqc_dir=directory("results/fastqc/{sample}")
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/fastqc/{sample}.fastqc.log"
    message:
        "{wildcards.sample}: Performing quality control on raw reads"
    threads: 6
    shell:
        """
        mkdir -p {output.fastqc_dir}

        fastqc --quiet --threads {threads} {input.fq1} {input.fq2:q} \
        --outdir {output.fastqc_dir} 2>> {log}
        """

# Rule for Adapter and quality trimming
rule trim_galore_paired:
    input:
        get_paired_fq
    output:
        trimmed_dir=directory("results/trimmed/{sample}"),
        fq1="results/trimmed/{sample}_trimmed_R1.fq.gz",
        fq2="results/trimmed/{sample}_trimmed_R2.fq.gz"
    message:
        "{wildcards.sample}: Trimming and performing quality control on paired-end FASTQ files"
    log:
        "logs/trim_galore/{sample}.trimmed.log"
    threads: 8
    conda:
        "envs/trimgalore.yaml"
    params:
        extra="--fastqc_args '-t 12'",
    shell:
        """
        trim_galore {params.extra} --cores {threads} \
            --paired \
            --gzip \
            --output_dir {output.trimmed_dir} \
            {input[0]} {input[1]} \
            2> {log}
        mv {output.trimmed_dir}/*_val_1.fq.gz {output.fq1}
        mv {output.trimmed_dir}/*_val_2.fq.gz {output.fq2}
        """

rule trim_galore_single:
    input:
        fq1=get_single_fq
    output:
        trimmed_dir=directory("results/trimmed/{sample}"),
        fq1="results/trimmed/{sample}_trimmed_R1.fq.gz"
    message:
        "{wildcards.sample}: Trimming and performing quality control on single-end FASTQ files"
    log:
        "logs/trim_galore/{sample}.trimmed.log"
    threads: 8
    conda:
        "envs/trimgalore.yaml"
    params:
        extra="--fastqc_args '-t 12'",
    shell:
        """
        trim_galore {params.extra} --cores {threads} \
            --gzip \
            --output_dir {output.trimmed_dir} \
            {input} \
            2> {log}
        mv {output.trimmed_dir}/*_trimmed.fq.gz {output.fq1}
        """

# Rule for Alignment
rule align_paired:
    input:
        unpack(get_paired_trimmed_fq),
        gtf="resources/annotation.filtered.gtf",
        star_index="resources/star_index"
    output:
        "results/align/bam_original/{sample}/{sample}.Aligned.out.bam",
        "results/align/bam_original/{sample}/{sample}.Aligned.toTranscriptome.out.bam"
    message:
        "{wildcards.sample}: Aligning on paired-end FASTQ files"
    log:
        "logs/star_align/{sample}.star_align.log"
    conda:
        "envs/star.yaml"
    threads: 12
    shell:
        """
        # Run STAR alignment
        STAR --genomeDir {input.star_index} \
             --readFilesIn {input[0]} {input[1]} \
             --runThreadN {threads} \
             --genomeLoad NoSharedMemory \
             --outFileNamePrefix results/align/bam_original/{wildcards.sample}/{wildcards.sample}. \
             --sjdbGTFfile {input.gtf} \
             --outSAMattrRGline 'ID:{wildcards.sample}' 'SM:{wildcards.sample}' \
             --outSAMtype BAM Unsorted \
             --quantMode TranscriptomeSAM \
             --twopassMode Basic \
             --readFilesCommand zcat \
             --runRNGseed 0 \
             --alignSJDBoverhangMin 1 \
             --outSAMattributes NH HI AS NM MD \
             --outSAMstrandField intronMotif \
             --quantTranscriptomeSAMoutput BanSingleEnd \
             --outFilterMultimapNmax 1 \
             --outFilterMismatchNmax 3 \
             --peOverlapNbasesMin 10 \
             --alignSplicedMateMapLminOverLmate 0.5 \
             --chimSegmentMin 10 \
             --chimOutType WithinBAM SoftClip \
             --chimJunctionOverhangMin 10 \
             --chimScoreMin 1 \
             --chimScoreDropMax 30 \
             --chimScoreJunctionNonGTAG 0 \
             --chimScoreSeparation 1 \
             --alignSJstitchMismatchNmax 5 -1 5 5 \
             --chimSegmentReadGapMax 3 2>> {log}
        """

rule align_single:
    input:
        unpack(get_single_trimmed_fq),
        gtf="resources/annotation.filtered.gtf",
        star_index="resources/star_index"
    output:
        "results/align/bam_original/{sample}/{sample}.Aligned.out.bam",
        "results/align/bam_original/{sample}/{sample}.Aligned.toTranscriptome.out.bam"
    message:
        "{wildcards.sample}: Aligning on single-end FASTQ files"
    log:
        "logs/star_align/{sample}.star_align.log"
    conda:
        "envs/star.yaml"
    threads: 12
    shell:
        """
        # Run STAR alignment
        STAR --genomeDir {input.star_index} \
             --readFilesIn {input[0]} \
             --runThreadN {threads} \
             --genomeLoad NoSharedMemory \
             --outFileNamePrefix results/align/bam_original/{wildcards.sample}/{wildcards.sample}. \
             --sjdbGTFfile {input.gtf} \
             --outSAMattrRGline 'ID:{wildcards.sample}' 'SM:{wildcards.sample}' \
             --outSAMtype BAM Unsorted \
             --quantMode TranscriptomeSAM \
             --twopassMode Basic \
             --readFilesCommand zcat \
             --runRNGseed 0 \
             --alignSJDBoverhangMin 1 \
             --outSAMattributes NH HI AS NM MD \
             --outSAMstrandField intronMotif \
             --quantTranscriptomeSAMoutput BanSingleEnd \
             --outFilterMultimapNmax 1 \
             --outFilterMismatchNmax 3 \
             --peOverlapNbasesMin 10 \
             --alignSplicedMateMapLminOverLmate 0.5 \
             --chimSegmentMin 10 \
             --chimOutType WithinBAM SoftClip \
             --chimJunctionOverhangMin 10 \
             --chimScoreMin 1 \
             --chimScoreDropMax 30 \
             --chimScoreJunctionNonGTAG 0 \
             --chimScoreSeparation 1 \
             --alignSJstitchMismatchNmax 5 -1 5 5 \
             --chimSegmentReadGapMax 3 2>> {log}
        """

# Rule for Sort and index alignments 
rule samtools_process:
    input:
        bam="results/align/bam_original/{sample}/{sample}.Aligned.out.bam",
        fasta="resources/genome.fa"
    output:
        sorted_bam="results/align/bam_original/{sample}/{sample}.sorted.bam",
        sorted_bam_bai="results/align/bam_original/{sample}/{sample}.sorted.bam.bai",
        stats="results/align/bam_original/{sample}/{sample}.sorted.bam.stats",
        flagstat="results/align/bam_original/{sample}/{sample}.sorted.bam.flagstat",
        idxstats="results/align/bam_original/{sample}/{sample}.sorted.bam.idxstats"
    message:
        "{wildcards.sample}: Sorting and indexing alignments"
    log:
        "logs/samtools_process/{sample}.samtools_process.log"
    conda:
        "envs/samtools.yaml"
    threads: 8
    shell:
        """
        # Sorts alignments in BAM file
        samtools sort -@ {threads} -m 8G -T tmp -o {output.sorted_bam} {input.bam} 2>> {log}
        
        # Indexes the sorted BAM file
        samtools index -@ 1 {output.sorted_bam} 2>> {log}
        
        # Generates statistics for the sorted BAM file
        samtools stats --threads 1 --reference {input.fasta} {output.sorted_bam} > {output.stats} 2>> {log}
        
        # Counts the number of alignments in the BAM file for each FLAG type
        samtools flagstat --threads 1 {output.sorted_bam} > {output.flagstat} 2>> {log}
        
        # Produces index statistics for the sorted BAM file
        samtools idxstats {output.sorted_bam} > {output.idxstats} 2>> {log}
        """

# Rule for Duplicate read marking
rule picard_mark_duplicates:
    input:
        bam="results/align/bam_original/{sample}/{sample}.sorted.bam",
        fasta="resources/genome.fa"
    output:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        metrics="results/align/bam_markup/{sample}/{sample}.markdup.sorted.MarkDuplicates.metrics.txt"
    message:
        "{wildcards.sample}: Marking duplicate reads"
    log:
        "logs/markdup/{sample}.markduplicates.log"
    conda:
        "envs/picard_markduplicates.yaml"
    resources:
        mem_mb=4096
    threads: 12
    shell:
        """
        picard MarkDuplicates \
        --SORTING_COLLECTION_SIZE_RATIO 0.12 \
        --ASSUME_SORTED true \
        --REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR tmp \
        --INPUT {input.bam} \
        --OUTPUT {output.bam} \
        --REFERENCE_SEQUENCE {input.fasta} \
        --METRICS_FILE {output.metrics} \
        2>> {log}
        """

rule samtools_markdup_processing:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        fasta="resources/genome.fa"
    output:
        bam_bai="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.bai",
        bam_stats="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.stats",
        bam_flagstat="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.flagstat",
        bam_idxstats="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.idxstats"
    message:
        "{wildcards.sample}: Running Samtools to mark duplicate reads"
    log:
        "logs/markdup/{sample}.markduplicates.samtools.log"
    conda:
        "envs/samtools.yaml"
    threads: 1
    shell:
        """
        # Index the BAM file
        samtools index -@ {threads} {input.bam} {output.bam_bai} 2>> {log}
        
        # Get BAM statistics
        samtools stats --threads {threads} -r {input.fasta} {input.bam} > {output.bam_stats} 2>> {log}
        
        # Get BAM flagstat
        samtools flagstat --threads {threads} {input.bam} > {output.bam_flagstat} 2>> {log}
        
        # Get BAM idxstats
        samtools idxstats {input.bam} > {output.bam_idxstats} 2>> {log}
        """

# Rule for Transcript assembly and quantification
rule stringtie:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf"
    output:
        gtf="results/stringtie/{sample}/{sample}.transcripts.gtf",
        abundance="results/stringtie/{sample}/{sample}.gene.abundance.txt",
        coverage="results/stringtie/{sample}/{sample}.coverage.gtf",
        ballgown_dir=directory("results/stringtie/{sample}/{sample}.ballgown")
    message:
        "{wildcards.sample}: Running Stringtie to assemble and quantify transcripts"
    log:
        "logs/stringtie/{sample}.stringtie.log"
    conda:
        "envs/stringtie.yaml"
    threads: 6
    shell:
        """
        stringtie {input.bam} \
                  --rf \
                  -G {input.gtf} \
                  -o {output.gtf} \
                  -A {output.abundance} \
                  -C {output.coverage} \
                  -b {output.ballgown_dir} \
                  -p {threads} \
                  -v \
                  -e 2>> {log}
        """

# Rule for Read counting relative to gene biotype
rule feature_counts:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf",
        header="templates/biotypes_header.txt"
    output:
        counts="results/feature_counts/{sample}/{sample}.featurecounts.txt",
        summary="results/feature_counts/{sample}/{sample}.featurecounts.txt.summary",
        stats="results/feature_counts/{sample}/{sample}.biotype_counts_rrna_mqc.tsv"
    message:
        "{wildcards.sample}: Read counting relative to gene biotype"
    log:
        count_log="logs/feature_counts/{sample}.feature_counts.log",
        stats_log="logs/feature_counts/{sample}.feature_counts.stats.log"
    conda:
        "envs/featurecounts.yaml"
    threads: 6
    shell:
        """
        featureCounts -B -C \\
            -g gene_biotype \\
            -t exon \\
            -p \\
            -T {threads} \\
            -a {input.gtf} \\
            -s 2 \\
            -o {output.counts} \\
            {input.bam} &> {log.count_log} &&
        
        cut -f 1,7 {output.counts} \\
            | tail -n +3 \\
            | cat {input.header} - \\
            > {output.stats} &&
        
        python scripts/mqc_features_stat.py \\
            {output.stats} \\
            -s {wildcards.sample} \\
            -f rRNA \\
            -o {output.stats} &> {log.stats_log}
        """

# Rule for Qualimap
rule qualimap:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf",
    output:
        directory("results/qualimap/{sample}")
    message:
        "{wildcards.sample}: Running Qualimap"
    log:
        "logs/qualimap/{sample}.qualimap.log"
    conda:
        "envs/qualimap.yaml"
    threads: 6
    shell:
        """
        unset DISPLAY &&
        mkdir -p tmp &&
        export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp &&
        qualimap --java-mem-size=29491M rnaseq \\
        -bam {input.bam} \\
        -gtf {input.gtf} \\
        -p strand-specific-reverse \\
        -pe \\
        -outdir {output} 2>> {log}
        """

# Rule for Assessment of technical / biological read duplication
rule dupradar:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf"
    output:
        dupradar_dir=directory("results/dupradar/{sample}")
    log:
        "logs/dupradar/{sample}.dupradar.log"
    params:
        is_paired=lambda wc: "paired" if is_paired_end(wc.sample) else "single"
    message:
        "{wildcards.sample}: Running dupRadar to evaluate technical and biological read duplication"
    conda:
        "envs/dupradar.yaml"
    threads: 4
    shell:
        """
        mkdir -p {output.dupradar_dir}
        scripts/dupradar.r {input.bam} {output.dupradar_dir} {input.gtf} 2 {params.is_paired} {threads} 2>> {log}
        """

rule bam_stat:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam"
    output:
        "results/rseqc/bam_stat/{sample}/{sample}.bam_stat.txt"
    log:
        "logs/rseqc/{sample}/{sample}.bam_stat.log"
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running BAM stat"
    shell:
        """
        bam_stat.py -i {input.bam} > {output} 2> {log}
        """

rule inner_distance:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        pdf="results/rseqc/inner_distance/{sample}/{sample}.inner_distance_plot.pdf"
    log:
        "logs/rseqc/{sample}/{sample}.inner_distance.log"
    conda:
        "envs/rseqc.yaml"
    threads: 1
    shell:
        """
        inner_distance.py -i {input.bam} -r {input.bed} -o results/rseqc/inner_distance/{wildcards.sample}/
        """

rule infer_experiment:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        "results/rseqc/infer_experiment/{sample}/{sample}.infer_experiment.txt"
    log:
        "logs/rseqc/{sample}/{sample}.infer_experiment.log"
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running Infer experiment"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        """

rule junnction_annotation:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        "results/junction_annotation/{sample}/{sample}.junction.bed",
    priority: 1
    log:
        "logs/rseqc/{sample}.junction_annotation.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log[0]} 2>&1"
        """

rule junction_saturation:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        pdf="results/junction_saturation/{sample}/{sample}.junction_saturation.pdf"
    log:
        "logs/junction_saturation/{sample}.junction_saturation.log"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        junction_saturation.py -i {input.bam} -r {input.bed} -o results/junction_saturation/{wildcards.sample}/
        """

rule read_distribution:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        "results/read_distribution/{sample}.read_distribution.txt"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.bed} > {output}
        """

rule read_duplication:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam"
    output:
        pdf="results/read_duplication/{sample}/{sample}.read_duplication.pdf"
    log:
        "logs/read_duplication/{sample}.read_duplication.log"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        read_duplication.py -i {input.bam} -o results/read_duplication/{wildcards.sample}/
        """


