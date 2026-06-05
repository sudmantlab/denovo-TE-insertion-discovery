#Moved these functions from samtools utils
def get_temp_bams_per_sample(wildcards, sample_table = config['sample_table']) -> List[str]:
    """
    A helper function for creating a list of expected per-smrtcell bam outputs per sample to be merged into a single bam.
    Adjusted to optionally include {hap} in the basepath if wildcards.hap exists.
    """
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)

    # Check if '{hap}' is in the basepath and format accordingly
    if hasattr(wildcards, 'hap'):
        basepath = "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}.{hap}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam"
        input_samples = [basepath.format(refalias=wildcards.refalias, mapper=wildcards.mapper, specimen=s[0], hap=wildcards.hap, lane=s[2], smrtcell=s[3]) for s in samples]
    else:
        basepath = "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam"
        input_samples = [basepath.format(refalias=wildcards.refalias, mapper=wildcards.mapper, specimen=s[0], lane=s[2], smrtcell=s[3]) for s in samples]

    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

def get_fastqs_per_sample(wildcards, sample_table = config['sample_table'], basepath = config['fastq_basepath']) -> List:
    """
    A helper function for creating a list of fastq.gz paths for each sample ahead of hifiasm assembly.
    """
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    input_samples = [basepath.format(specimen=s[0], group = s[1], lane=s[2], smrtcell = s[3]) for s in samples]
    #print(input_samples)
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples




rule hifiasm:
    input: get_fastqs_per_sample
    output:
        p_ctg_hap1 = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.gfa",
        p_ctg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.noseq.gfa",
        p_utg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.gfa",
        p_utg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.lowQ.bed",
        p_utg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.noseq.gfa",
        r_utg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.gfa",
        r_utg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.lowQ.bed",
        r_utg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.noseq.gfa",
        ec = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ec.bin",
        ovlp_reverse = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ovlp.reverse.bin",
        ovlp_source = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ovlp.source.bin"
    params:
        prefix = "output/assembly/hifiasm/{specimen}/{specimen}.asm"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 12
    log: 
        "logs/assembly/hifiasm/{specimen}/{specimen}.asm.log"
    conda: 
        "../envs/environment.yml"
    shell: 
        """
        hifiasm -o {params.prefix} -t {threads} {input} > {log} 2>&1
        """


rule gfaToFa:
    input: "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.{hap}.p_ctg.gfa"
    output: "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    log: "logs/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa.log"
    threads: 1
    conda: 
        "../envs/environment.yml"
    shell: 
        "gfatools gfa2fa {input} > {output} 2> {log}"

rule ragtag_scaffold:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    output:
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.fasta"
        #paf = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.asm.paf"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 10
    params:
        refgenome = config['reference']['fasta'],
        outdir = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}"
    shell:
        """
        mkdir -p {params.outdir}

        ragtag.py scaffold {params.refgenome} \
        {input} \
        -u -w --aligner minimap2 -t {threads} \
        -o {params.outdir}

        # Rename the ragtag output files to include specimen and hap info
        find {params.outdir} -name "ragtag.scaffold.*" \
        -exec sh -c 'for f do dir=$(dirname "$f"); \
        base=$(basename "$f"); suffix=${{base#ragtag.scaffold.}}; \
        mv "$f" "$dir/{wildcards.specimen}.{wildcards.hap}.scaffold.$suffix"; \
        done' sh {{}} +

        # Modify all FASTA headers to include hap information
        sed -i 's/_RagTag$/_RagTag_{wildcards.hap}/' {output.fasta}
        """

rule update_scaffold_agp:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.agp"
    output:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.updated.agp"
    shell:
        """
        awk -F'\t' '{{sub(/_RagTag$/, "_RagTag_{wildcards.hap}", $1);
        print}}' OFS='\t' {input} > {output}
        """

rule cat_scaffolds:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/hap1/{specimen}.hap1.scaffold.fasta",
        "output/assembly/hifiasm/{specimen}/scaffolded/hap2/{specimen}.hap2.scaffold.fasta"
    output:
        fa = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
        fai = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta.fai"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 1
    conda:
        "../envs/environment.yml"
    shell:
        """
        cat {input} > {output.fa}
        samtools faidx {output.fa}
        """
'''
rule cat_contigs:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.hap1.fa"
        "output/assembly/hifiasm/{specimen}/{specimen}.hap2.fa"
    output:
        fasta = "output/assembly/hifiasm/{specimen}/{specimen}.diploid.fasta",
        fai = "output/assembly/hifiasm/{specimen}/{specimen}.diploid.fasta.fai"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 1
    conda:
        "../envs/environment.yml"
    shell:
        """
        cat {input} > {output.fasta}
        samtools faidx {output.fasta}
        """
'''

### QC: QUAST (assembly completeness) ###

rule quast_raw:
    # Optional: Run QUAST on the raw assembly fasta without scaffolding.
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.hap1.fa",
        "output/assembly/hifiasm/{specimen}/{specimen}.hap2.fa"
    output:
        "output/assembly/hifiasm/{specimen}/quast/raw/report.html"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 6
    params:
        outdir = "output/assembly/hifiasm/{specimen}/quast/raw"
    shell:
        """
        quast.py {input} \
        -o {params.outdir} \
        --large --est-ref-size 3100000000 --no-icarus
        """

use rule quast_raw as quast_scaffolded with:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/hap1/{specimen}.hap1.scaffold.fasta",
        "output/assembly/hifiasm/{specimen}/scaffolded/hap2/{specimen}.hap2.scaffold.fasta"
    output:
        "output/assembly/hifiasm/{specimen}/quast/scaffolded/report.html"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 6
    params:
        outdir = "output/assembly/hifiasm/{specimen}/quast/scaffolded"

### QC: Flagger (assembly error estimation) ###

rule prepare_flagger_inputs:
    # Flagger requires a BED file that describes the assembly contigs + coordinates.
    input:
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
        fai = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta.fai"
    output:
        bed = "output/assembly/flagger/{specimen}/whole_genome.bed",
        json = "output/assembly/flagger/{specimen}/annotations_path.json"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.bed})
        
        # Create whole genome bed from fai
        cat {input.fai} | awk '{{print $1"\t0\t"$2}}' > {output.bed}
        
        # Create annotations JSON
        echo '{{"whole_genome": "{output.bed}"}}' > {output.json}
        """

'''
#Commenting out flagger for now...

use rule prepare_flagger_inputs as HG002_unscaffolded_prepare with:
    # Use the unscaffolded version of the HG002 assembly as a reference/control.
    input:
        fasta = "output/assembly/hifiasm/HG002/HG002.unscaffolded.diploid.fa",
        fai = "output/assembly/hifiasm/HG002/HG002.unscaffolded.diploid.fa.fai"
    output:
        bed = "output/assembly/flagger/HG002/unscaffolded/HG002_whole_genome.bed",
        json = "output/assembly/flagger/HG002/unscaffolded/HG002_annotations_path.json"


rule flagger_bam_to_coverage:
    input:
        bam = "output/alignment/scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam",
        bed = "output/assembly/flagger/{specimen}/whole_genome.bed",
        json = "output/assembly/flagger/{specimen}/annotations_path.json"
    output:
        cov = "output/assembly/flagger/{specimen}/coverage_file.cov.gz"
    threads: 16
    resources:
        mem_mb = 32000
    params:
        base_dir = config['workdir']
    shell:
        """
        singularity exec \
            --bind {params.base_dir}:{params.base_dir} \
            docker://mobinasri/flagger:v1.1.0 \
            bam2cov \
                --bam {input.bam} \
                --output {output.cov} \
                --annotationJson {input.json} \
                --threads {threads} \
                --baselineAnnotation whole_genome
        """

rule run_flagger:
    input:
        cov = "output/assembly/flagger/{specimen}/coverage_file.cov.gz",
        alpha_tsv = "config/packages/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv"
    output:
        "output/assembly/flagger/{specimen}/prediction_summary_final.tsv"
    threads: 16
    resources:
        mem_mb = 32000
    params:
        base_dir = config['workdir']
    shell:
        """        
        singularity exec \
            --bind {params.base_dir}:{params.base_dir} \
            docker://mobinasri/flagger:v1.1.0 \
            hmm_flagger \
                --input {input.cov} \
                --outputDir $(dirname {output}) \
                --alphaTsv {input.alpha_tsv} \
                --labelNames Err,Dup,Hap,Col \
                --threads {threads}
        """

rule summarize_flagger:
    input:
        expand("output/assembly/flagger/{specimen}/prediction_summary_final.tsv", specimen = [x for x in specimens if x != '900'])
    output:
        "output/assembly/flagger/all_prediction_summary_final.tsv"
    threads: 1
    shell:
        """
        echo -e "specimen\tErr\tDup\tHap\tCol\tUnk" > {output}
        find output/assembly/flagger -name "prediction_summary_final.tsv" | while read file; do
            specimen=$(echo "$file" | awk -F'/' '{{print $(NF-1)}}')
            awk -v specimen="$specimen" '$1=="PREDICTION" && $2=="base_level" && $3=="percentage" && $4=="annotation" && $5=="whole_genome" {{print specimen "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}}' "$file" >> {output}
        done
        """
'''


### Optional: RepeatMasker annotation of individual assemblies ###

rule filter_canonical_chromosomes:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.fasta"
    output:
        fasta = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.canonical.fasta",
        temp_chroms = temp("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/canonical_chroms.txt")
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 1
    conda:
        "../envs/environment.yml"
    shell:
        """
        samtools faidx {input}
        
        grep "^chr" {input}.fai | grep -v -E "random|chrUn|_h[12]tg" | grep -v "\." | cut -f1 > {output.temp_chroms}
        
        samtools faidx {input} -r {output.temp_chroms} > {output.fasta}
        """

'''
rule filter_hg38_canonical_chromosomes:
    input:
        config['reference']['fasta']
    output:
        fasta = config['reference']['fasta'].replace('.fa.gz', '.canonical.fa'),
        temp_chroms = temp(config['reference']['fasta'].replace('.fa.gz', '.canonical.chroms.txt'))
    threads: 1
    conda:
        "../envs/environment.yml"
    shell:
        """
        samtools faidx {input}
        
        awk '$1 ~ /^chr[1-9]$/ || $1 ~ /^chr[1-2][0-9]$/ || $1 ~ /^chr[XYM]$/' {input}.fai | cut -f1 > {output.temp_chroms}
        
        if [ ! -s {output.temp_chroms} ]; then
            echo "No chromosomes matched the pattern" >&2
            exit 1
        fi
        
        samtools faidx {input} -r {output.temp_chroms} > {output.fasta}
        """


rule split_canonical_fasta_OLD:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.canonical.fasta"
    output:
        expand("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa", allow_missing = True, chr = chrs) # for all chrs
        #directory("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas/")
    conda:
        "../envs/environment.yml"
    params:
        outdir = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas"
    shell:
        """
        mkdir -p {params.outdir}
        awk -v outdir="{params.outdir}" '
        /^>/ {{
            if (file) {{
                close(file)
            }}
            filename = substr($0, 2)
            sub(/_RagTag.*$/, "", filename)
            file = outdir "/" filename ".fa"
            print $0 > file
            next
        }}
        {{ if (file) print > file }}' {input}
        
        """

'''

rule split_canonical_fasta:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/{specimen}.{hap}.scaffold.canonical.fasta"
    output:
        #directory("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas")
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa"
    conda:
        "../envs/environment.yml"
    params:
        outdir = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas"
    shell:
        r"""
        mkdir -p {params.outdir}
        awk -v outdir="{params.outdir}" '
        /^>/ {{
            if (file) {{
                close(file)
            }}
            # drop leading ">"
            header = substr($0, 2)
            # remove trailing _RagTag... if present
            sub(/_RagTag.*$/, "", header)
            # drop everything after first space (description)
            split(header, b, /[ \t]/)
            header = b[1]
            # if header starts with chr..., keep only up to first "_"
            if (header ~ /^chr/) {{
                split(header, a, "_")
                filename = a[1]
            }} else {{
                filename = header
            }}
            file = outdir "/" filename ".fa"
            print $0 > file
            next
        }}
        {{ if (file) print > file }}' {input}
        if [ ! -f {output} ]; then touch {output}; fi
        """


'''

rule repeatmasker_per_chr:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa"
    output:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.fa.out"
    log:
        "logs/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{chr}.log"
    conda:
        "../envs/environment.yml"
    params:
        engine = config['repeatmasker']['engine'],
        species = config['repeatmasker']['species'],
        custom_lib = config['repeatmasker']['custom_lib'],
        outdir = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr"
    threads: 4
    resources:
        mem_mb = 24000
    shell:
        """
        #mkdir -p {params.outdir}
        
        if [ -n "{params.custom_lib}" ]; then
            # Run with -lib if param is not empty
            RepeatMasker -pa {threads} \
                -engine {params.engine} \
                -nocut -gff \
                -lib {params.custom_lib} \
                -dir {params.outdir} \
                {input} &> {log}
        else
            # Run with -species otherwise
            RepeatMasker -pa {threads} \
                -engine {params.engine} \
                -nocut -gff \
                -species {params.species} \
                -dir {params.outdir} \
                {input} &> {log}
        fi
        """

        #"""
        #RepeatMasker -pa {threads} \
        #    -engine {params.engine} \
        #    -nocut -gff \
        #    -species {params.species} \
        #    -dir {params.outdir} \
        #    {input} &> {log}
        #"""
'''

'''
rule trf_per_chr:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa"
    output:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.trf.bed"
    log:
        "logs/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{chr}.trf.log"
    params:
        outdir = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr"
    conda:
        "../envs/other_masking.yaml"
    threads: 12    
    resources:
        mem_mb = 200000
    shell:
        """
        mkdir -p {params.outdir}
        if [[ "{input}" == *"chrY"* ]]; then
            trf {input} 2 7 7 80 10 50 2000 -l 30 -h -ngs | awk 'BEGIN {{OFS="\\t"}} /@/ {{chr=$1; gsub(/@/, "", chr); next}} /^[0-9]/ {{print chr, $1-1, $2, "TRF_repeat", $8, "+"}}' > {output} || touch {output}
        else
            trf {input} 2 7 7 80 10 50 2000 -l 30 -h -ngs | awk 'BEGIN {{OFS="\\t"}} /@/ {{chr=$1; gsub(/@/, "", chr); next}} /^[0-9]/ {{print chr, $1-1, $2, "TRF_repeat", $8, "+"}}' > {output}
        fi
        """
'''


rule ULTRA_per_chr:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa"
    output:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.ultra.bed"
    log:
        "logs/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{chr}.ultra.log"
    params:
        outdir = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr",
        ultraDir = config["ultra_path"]
    conda:
        "../envs/other_masking.yaml"
    threads: 12
    resources:
        mem_mb = 200000
    shell:
        """
        mkdir -p {params.outdir}
        if [[ "{input}" == *"chrY"* ]]; then
            {params.ultraDir}/ultra -t 12 --bed -o {output} {input} || touch {output}
        else
            {params.ultraDir}/ultra -t 12 --bed -o {output} {input}
        fi
        """

rule longdust_per_chr:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa"
    output:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.longdust.bed"
    log:
        "logs/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{chr}.longdust.log"
    conda:
        "../envs/longdust.yaml"
    resources:
        mem_mb = 50000
    shell:
        """
        if [[ "{input}" == *"chrY"* ]]; then
            longdust {input} > {output} || touch {output}
        else
            longdust {input} > {output}
        fi
        """



'''
rule repeatmasker_to_bed:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.fa.out"
    output:
        temp("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.bed")
    threads: 1
    shell:
        """
        tail -n +4 {input} | \
        awk 'BEGIN{{OFS="\t"}} 
        {{
            print $5, $6-1, $7, $10"#"$11, $1, $9
        }}' > {output}
        """
'''

rule combine_trf_outs:
    input:
        expand("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.trf.bed", allow_missing = True, chr = chrs, hap = ['hap1', 'hap2']) #changed from True
    output:
        trfOut = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.trf.bed",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output.trfOut}
        """

rule combine_ultra_outs:
    input:
        expand("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.ultra.bed", allow_missing = True, chr = chrs, hap = ['hap1', 'hap2']) #changed from T>
    output:
        ultraOut = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.ultra.bed",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output.ultraOut}
        """

rule combine_longdust_outs:
    input:
        expand("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.longdust.bed", allow_missing = True, chr = chrs, hap = ['hap1', 'hap2']) #changed from True
    output:
        longDustOut = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.longdust.bed",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output.longDustOut}
        """

'''
rule combine_longdust_and_trf_outs:
    input:
        longDustIn = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.longdust.bed",
        trfIn = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.trf.bed",
    output:
        allSimpleOut = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.all_simple_repeats.bed",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input.longDustIn} {input.trfIn} | sort -k1,1 -k2,2n > {output.allSimpleOut}
        """


rule combine_longdust_and_trf_final:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/hap1/repeatmasker/{specimen}.hap1.all_simple_repeats.bed",
        "output/assembly/hifiasm/{specimen}/scaffolded/hap2/repeatmasker/{specimen}.hap2.all_simple_repeats.bed",
    output:
        final = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.all_simple_repeats.bed",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output.final}
        """
'''

rule combine_longdust_and_ultra_outs:
    input:
        longDustIn = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.longdust.bed",
        ultraIn = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.ultra.bed",
    output:
        allSimpleOut = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.all_simple_repeats.bed",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input.longDustIn} {input.ultraIn} | sort -k1,1 -k2,2n > {output.allSimpleOut}
        """


rule combine_longdust_and_ultra_final:
    input:
        "output/assembly/hifiasm/{specimen}/scaffolded/hap1/repeatmasker/{specimen}.hap1.all_simple_repeats.bed",
        "output/assembly/hifiasm/{specimen}/scaffolded/hap2/repeatmasker/{specimen}.hap2.all_simple_repeats.bed",
    output:
        final = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.all_simple_repeats.bed",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output.final}
        """

'''
rule combine_repeatmasker_outs:
    input:
        expand("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.fa.out", allow_missing = True, chr = chrs, hap = ['hap1', 'hap2'])
    output:
        rmout = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.repeatmasker.out",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input} > {output.rmout}
        """


rule combine_repeatmasker_beds:
    input:
        expand("output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/per_chr/{chr}.bed", allow_missing = True, chr = chrs, hap = ['hap1', 'hap2'])
    output:
        bed = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.repeatmasker.bed.gz",
        tbi = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}/repeatmasker/{specimen}.{hap}.repeatmasker.bed.gz.tbi"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/environment.yml"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bgzip > {output.bed}
        tabix -p bed {output.bed}
        """
'''
