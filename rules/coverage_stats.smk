rule genomicCoverage:
    input:
        bam = 'output/alignment/scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam',
    output:
        genomicCoverageTab = "output/alignment/scaffolded/minimap2/standard/coverage_stats/{specimen}.coverage.tab",
    conda:
        "../envs/environment.yml"
    threads: 5
    shell:
        """
        samtools depth --threads 5 -a {input.bam} |  awk '{{sum+=$$3}} END {{ print "Average = ",(sum/NR)*2}}' > {output.genomicCoverageTab}
        """
