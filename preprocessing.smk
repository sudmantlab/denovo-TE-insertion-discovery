
rule uBAMtoFastq:
    # Helper rule that takes an unaligned BAM file and transforms it into a fastq.
    # This transformation retains the read quality (rq) + MM & ML (methylation information) tags for each read.
    input:
        "data/{specimen}/{lane}/{smrtcell}.bam"
    output:
        "data/{specimen}/{lane}/{smrtcell}.fastq.gz"
    conda: "../envs/environment.yml"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 10
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T MM,ML {input} -0 {output}
        """


