### Germline SVs ###
'''
rule sniffles_standard:
    input:
        bam = "output/alignment/{reference}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{reference}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/alignment/{reference}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/{reference}/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/{reference}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        reference = "[A-Za-z0-9]+",
    conda:
        '../envs/environment.yml'
    threads:
        10
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/{reference}/{mapper}/standard/variants/sniffles_standard/{specimen}.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {params.refgenome} \
        --tandem-repeats {params.repeats} \
        --threads {threads} \
        --mapq {params.mapq} \
        --output-rnames &> {log}
        """


rule sniffles_standard_scaffolded:
    # The same rule as sniffles_standard, except it doesn't use a tandem repeat annotation file
    # and uses the self assembly fasta as a reference.
    input:
        bam = "output/alignment/scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
    output:
        vcf='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        '../envs/environment.yml'
    threads:
        10
    params:
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --threads {threads} \
        --mapq {params.mapq} \
        --output-rnames &> {log}
        """

### Mosaic SVs (low-frequency SVs) ###

rule sniffles_mosaic:
    input:
        bam = "output/alignment/{reference}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{reference}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
        #repeats = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.all_simple_repeats.bed",
    output:
        vcf='output/alignment/{reference}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/{reference}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/{reference}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        reference = "[A-Za-z0-9]+",
    conda:
        '../envs/environment.yml'
    threads:
        10
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/{reference}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {params.refgenome} \
        --tandem-repeats {params.repeats} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min {params.mosaic_af_min} \
        --mosaic-af-max {params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} &> {log}
        """

rule sniffles_mosaic_scaffolded:
    # The same rule as sniffles_mosaic, except it 
    # uses the self assembly fasta as a reference.
    input:
        bam = "output/alignment/scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
        repeats = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.all_simple_repeats.bed"
    output:
        vcf='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        '../envs/environment.yml'
    threads:
        10
    params:
        # Remove limiter on mapQ due to multimapping
        minsupport = config['sniffles']['minsupport'],
        mapq = 0,
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.log"
    shell:
        """
        # Remove limiter on mapQ due to multimapping
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --tandem-repeats {input.repeats} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min {params.mosaic_af_min} \
        --mosaic-af-max {params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} &> {log}
        """

'''

rule sniffles_mosaic_scaffolded_qc_all:
    # The same rule as sniffles_mosaic, except it doesn't use a tandem repeat annotation file (only compatible with hg38)
    # and uses the self assembly fasta as a reference.
    # Yields all candidates without filtering.
    input:
        bam = "output/alignment/scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
        repeats = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.all_simple_repeats.bed",
    output:
        vcf='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz',
        snf='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.snf',
        tbi='output/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz.tbi'
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        '../envs/environment.yml'
    threads:
        10
    params:
        minsupport = 0,
        mapq = 0,
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.log"
    shell:
        """
        # Yield all candidates, regardless of QC status
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --tandem-repeats {input.repeats} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min={params.mosaic_af_min} \
        --mosaic-af-max={params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} \
        --dev-no-qc &> {log}
        """

rule longcall_scaffolded:
    input:
        bam = "output/alignment/scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/scaffolded/{mapper}/standard/variants/longcall/{specimen}.vcf'
    threads:
        12
    params:
        longcallPath = config['longcall_path'],
        library = config['repeatmasker']['custom_lib']        
    shell:
        """
        {params.longcallPath}/longcallD call -s -t12 {input.fasta} {input.bam} -T {params.library} > {output.vcf}
        """



rule miniSV_ref: 
    input:
        pafSelf = "output/alignment/scaffolded/minimap2/standard/mapped/{specimen}.paf",
        pafRef = "output/alignment/reference/minimap2/standard/mapped/{specimen}.paf",
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta"
    output:
        rsv = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.withref.rsv',
        msv = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.withref.msv',
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.withref.vcf'
    params:
        miniSVPath =  config['miniSV_path'],
        repeats = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.all_simple_repeats.bed",
    resources:
        mem_mb = 80000
    shell:
        """
        export PATH="$PATH:/global/scratch/users/landen_gozashti/tools/minisv/"
        export PATH="$PATH:/global/scratch/users/landen_gozashti/tools/k8-1.2/"
        {params.miniSVPath}/minisv.js e -0b {params.repeats} {input.pafRef} {input.pafSelf} | bash > {output.rsv}
        cat {output.rsv} | sort -k1,1 -k2,2 -S4g | {params.miniSVPath}/minisv.js merge -c1 -s0 -e1 - | awk '$5<2' | grep "SVTYPE=INS" > {output.msv}
        {params.miniSVPath}/minisv.js genvcf {output.msv} > {output.vcf}
        """ 

rule miniSV:
    input:
        pafSelf = "output/alignment/scaffolded/minimap2/standard/mapped/{specimen}.paf",
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
        snifflesVcf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz'
    output:
        rsv = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.rsv',
        msv = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.msv',
        vcf = temp('output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.vcf.pre'),
        header = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.header.txt' 
    params:
        miniSVPath =  config['miniSV_path'],
        repeats = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.all_simple_repeats.bed",
    resources:
        mem_mb = 80000
    shell:
        """
        export PATH="$PATH:/global/scratch/users/landen_gozashti/tools/minisv/"
        export PATH="$PATH:/global/scratch/users/landen_gozashti/tools/k8-1.2/"
        {params.miniSVPath}/minisv.js e -0b {params.repeats} {input.pafSelf} | bash > {output.rsv}
        cat {output.rsv} | sort -k1,1 -k2,2 -S4g | {params.miniSVPath}/minisv.js merge -c1 -s0 -e1 - | awk '$5<2' | grep "SVTYPE=INS" > {output.msv}
        {params.miniSVPath}/minisv.js genvcf {output.msv} > {output.vcf}
        gzip -c -d {input.snifflesVcf} | grep '##contig=' > {output.header}
        """

rule filterMiniSV:
    input:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.vcf.pre',
        header = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.header.txt'
    output:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.vcf'
    run:
        with open(input.header,'r') as f:
            header = f.read()
      
        def replace_ins_with_seq(in_vcf, out_vcf, header):
            
            #Read a VCF, and for any record whose ALT is '<INS>',
            #replace ALT with the sequence from SEQ=... in the INFO field.
            #Write the modified VCF to out_vcf.
            
            with open(in_vcf) as inp, open(out_vcf, 'w') as out:
                for line in inp:
                    if line.startswith('#'):
                        out.write(line)
                        if '##fileformat' in line:
                            out.write(header)
                        continue

                    cols = line.rstrip('\n').split('\t')
                    if len(cols) < 8:
                        out.write(line)
                        continue
                    chrom = cols[0]
                    if 'h2tg' in chrom or 'h1tg' in chrom:
                        continue
                    alt = cols[4]
                    info = cols[7]

                    if alt != '<INS>':
                        out.write(line)
                        continue

                # find SEQ= in INFO
                    seq = None
                    for field in info.split(';'):
                        if field.startswith('SEQ='):
                            seq = field.split('=', 1)[1]
                        if field.startswith('END='):
                            oldEnd = field
                            newEnd = 'END=' + str(int(field.split('=', 1)[1]) + 2)

                    newInfo = info.replace(oldEnd,newEnd)
                    cols[7] = newInfo
                    cols[3] = 'A'
                    if seq is not None and seq != '':
                        cols[4] = seq
 
                    out.write('\t'.join(cols) + '\n')


        replace_ins_with_seq(input.vcf, output.vcf, header)


rule sniffles_mosaic_reference_qc_all:
    # The same rule as sniffles_mosaic, except it doesn't use a tandem repeat annotation file (only compatible with hg38)
    # and uses the self assembly fasta as a reference.
    # Yields all candidates without filtering.
    input:
        bam = "output/alignment/reference/minimap2/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/reference/minimap2/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
    output:
        vcf='output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz',
        snf='output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.snf',
        tbi='output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz.tbi'
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        '../envs/environment.yml'
    threads:
        10
    params:
        minsupport = 0,
        mapq = 0,
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.log"
    shell:
        """
        # Yield all candidates, regardless of QC status
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --tandem-repeats {input.repeats} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min={params.mosaic_af_min} \
        --mosaic-af-max={params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} \
        --dev-no-qc &> {log}
        """

rule longcall_reference:
    input:
        bam = "output/alignment/reference/minimap2/standard/mapped/{specimen}.sorted.merged.bam",
        fasta = config['reference']['fasta']
    output:
        vcf='output/alignment/reference/minimap2/standard/variants/longcall/{specimen}.vcf'
    threads:
        12
    params:
        longcallPath = config['longcall_path'],
        library = config['repeatmasker']['custom_lib']
    shell:
        """
        {params.longcallPath}/longcallD call -s -t12 {input.fasta} {input.bam} -T {params.library} > {output.vcf}
        """
