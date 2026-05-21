rule minimap2:
    input:
        hifi = "data/{specimen}/{lane}/{smrtcell}.fastq.gz"
    output: 
        temp("output/alignment/reference/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/environment.yml"
    threads: 14
    shell: 
        """
        minimap2 --version && minimap2 {params.refgenome} {input.hifi} -t {threads} -ax map-hifi -Y -y -L --eqx --cs --MD -R '{params.readgroup}' | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}
        """

rule minimap2miniSV:
    input:
        hifi = "data/{specimen}/{lane}/{smrtcell}.fastq.gz",
        fa = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta"
    output:
        refPaf = temp("output/alignment/reference/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.paf"),
        selfPaf = temp("output/alignment/scaffolded/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.paf")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/environment.yml"
    threads: 14
    shell:
        """
        minimap2 --version && minimap2 {params.refgenome} {input.hifi}  -cxmap-hifi -s50 --ds -t {threads} > {output.refPaf}
        minimap2 --version && minimap2 {input.fa} {input.hifi}  -cxmap-hifi -s50 --ds -t {threads} > {output.selfPaf} 
        """

rule minimap2_to_scaffolded:
    input:
        hifi = "data/{specimen}/{lane}/{smrtcell}.fastq.gz",
        fa = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta"
    output:
        temp("output/alignment/scaffolded/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        readgroup = config['minimap2']['readgroup'],
        # minQ = config['samtools']['minQ'] # skip minQ for multimapping
    conda: "../envs/environment.yml"
    threads: 14
    shell:
        """
        minimap2 --version && minimap2 -t {threads} -ax map-hifi -Y -y -L --eqx --cs -I8g --MD -R '{params.readgroup}' {input.fa} {input.hifi} | samtools view -b > {output}
        """
