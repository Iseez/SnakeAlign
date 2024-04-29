#Configuration file with defined variables
configfile: "config.yaml"

#Parse arguments
samples = config["samples"]
reference_file = config["reference_file"]
reference_name = config["reference_name"]
quality = config["quality"]
SIZE = config["size"]
STEP = config["step"]

rule clean_fastq:
    input:
        r1 = "data/reads/{sample}_R1.fastq",
        r2 = "data/reads/{sample}_R2.fastq"
    output:
        report_html = "data/{sample}/reports/{sample}_fastp.html",
        report_json = "data/{sample}/reports/{sample}_fastp.json",
        unpaired1 = "data/{sample}/clean/{sample}_R1_unpaired.fastq.gz",
        unpaired2 = "data/{sample}/clean/{sample}_R2_unpaired.fastq.gz",
        r1_out = "data/{sample}/clean/{sample}_R1_clean.fastq.gz",
        r2_out = "data/{sample}/clean/{sample}_R2_clean.fastq.gz"
    params:
        threads = 16,
        compression = 9
    shell:
        """
        module load fastp/0.20.0
        fastp -i {input.r1} -I {input.r2} -o {output.r1_out} -O {output.r2_out} \
            --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} \
            -w {params.threads} -y -x -z {params.compression} \
            -j {output.report_json} -h {output.report_html}
        """

rule align:
    input:
        ref = reference_file,
        r1 = rules.clean_fastq.output.r1_out,
        r2 = rules.clean_fastq.output.r2_out
    params:
        threads = 10
    output:
        bam = f"data/{{sample}}/bam/{{sample}}_{reference_name}.bam"
    shell:
        """
        module load bwa/0.7.4 htslib/1.2.1 gcc/5.1.0 samtools/1.10
        bwa mem -M -t {params.threads} {input.ref} {input.r1} {input.r2} | samtools view -hbS - \
            | samtools sort -o {output.bam} -
        """
rule mark_duplicates:
    input:
        bam = rules.align.output.bam
    output:
        rmdup = temp(f"data/{{sample}}/bam/{{sample}}_{reference_name}.rmdup.bam")
    params:
        metrics = f"data/{{sample}}/bam/{{sample}}_{reference_name}_duplicateMatrix"
    shell:
        """
        module load picard/2.6.0
        picard MarkDuplicates INPUT= {input.bam} OUTPUT= {output.rmdup} \
            METRICS_FILE= {params.metrics} VALIDATION_STRINGENCY=LENIENT
        rm {params.metrics}
        """

rule groups:
    input:
        bam = rules.mark_duplicates.output.rmdup
    output:
        groups = f"data/{{sample}}/bam/{{sample}}_{reference_name}.rmdup.addgp.bam"
    params:
        platform  = config["platform"]
    shell:
        """
        module load picard/2.6.0
        picard AddOrReplaceReadGroups I= {input.bam} O= {output.groups} \
            LB={wildcards.sample} PL=illumina PU={wildcards.sample} \
            SM={wildcards.sample} VALIDATION_STRINGENCY=LENIENT
        """

rule index:
    input:
        bam = rules.groups.output.groups
    output:
        index = f"data/{{sample}}/bam/{{sample}}_{reference_name}.rmdup.addgp.bam.bai"
    shell:
        """
        module load samtools/1.10
        samtools index {input.bam}
        """

rule depth:
    input:
        bam = rules.groups.output.groups,
        bai = rules.index.output.index
    output:
        dp_file = f"data/{{sample}}/reports/{{sample}}_{reference_name}_Q{quality}.depth"
    params:
        qual = quality
    shell:
        """
        module load htslib/1.2.1 gcc/5.1.0 samtools/1.10
        samtools depth -aa -Q {params.qual} {input.bam} > {output.dp_file}
        """

rule plot_coverage:
    input:
        dp = rules.depth.output.dp_file
    output:
        PerContig = f"data/{{sample}}/stats/{{sample}}_{reference_name}_{quality}_{SIZE}by{STEP}_perContig.csv",
        Windows = f"data/{{sample}}/stats/{{sample}}_{reference_name}_{quality}_{SIZE}by{STEP}_windows.csv",
        PerRef = f"data/{{sample}}/stats/{{sample}}_{reference_name}_{quality}_{SIZE}by{STEP}_perReference.csv",
        imgdpth = f"data/{{sample}}/stats/{{sample}}_{reference_name}_{quality}_{SIZE}by{STEP}.pdf",
        imgdpth_mt = f"data/{{sample}}/stats/{{sample}}_{reference_name}_{quality}_{SIZE}by{STEP}_mt.pdf"
    script:
        "scripts/make_depth_plot.R"

rule normalized_coverage:
    input:
        PerRef = rules.plot_coverage.output.PerRef,
        Windows = rules.plot_coverage.output.Windows
    params:
        strain = lambda w: w.get("sample"),
        ref = reference_name,
        folder = lambda w: f"data/{w.sample}/stats/"
    output:
        contigs = f"data/{{sample}}/stats/{{sample}}_{reference_name}_{quality}_{SIZE}by{STEP}_normalized_contigs.csv"
    script:
        "scripts/normalized_coverage.R"

rule final:
    input:
        expand(f"data/{{sample}}/stats/{{sample}}_{reference_name}_{quality}_{SIZE}by{STEP}_normalized_contigs.csv",sample = samples)
    output:
        "data/done.txt"
    shell:
        """
        echo Finished processing samples. > {output}
        """
