
rule make_insertion_fasta_noqcall:
    input:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.vcf.gz'
    output:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}_insertions.fa'
    shell:
        """
        gzip -c -d {input} | grep -v "^#" | grep "INS" | awk '{{print ">"$1"_"$2"\\n"$5}}' > {output}
        """


rule repeatmasker_insertions_noqcall:
    input:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}_insertions.fa'
    output:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/repeatmasker/{specimen}_insertions.fa.out'
    conda:
        "../envs/environment.yml"
    params:
        engine = config['repeatmasker']['engine'],
        species = config['repeatmasker']['species'],
        custom_lib = config['repeatmasker']['custom_lib'],
        outdir = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/repeatmasker/'
    threads: 8
    resources:
        mem_mb = 24000
    shell:
        """
        mkdir -p {params.outdir}

        if [ -n "{params.custom_lib}" ]; then
            # Run with -lib if param is not empty
            RepeatMasker -pa {threads} \
                -s \
                -xsmall \
                -engine {params.engine} \
                -nocut -gff \
                -lib {params.custom_lib} \
                -dir {params.outdir} \
                {input}
        else
            # Run with -species otherwise
            RepeatMasker -pa {threads} \
                -s \
                -xsmall \
                -engine {params.engine} \
                -nocut -gff \
                -species {params.species} \
                -dir {params.outdir} \
                {input} 
        fi
        """

rule make_insertion_fasta_qcall:
    input:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz'
    output:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}_insertions.qc_all.fa'
    shell:
        """
        gzip -c -d {input} | grep -v "^#" | grep "INS" | awk '{{print ">"$1"_"$2"\\n"$5}}' > {output}
        """


rule repeatmasker_insertions_qcall:
    input:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}_insertions.qc_all.fa'
    output:
        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/repeatmasker/{specimen}_insertions.qc_all.fa.out'
    conda:
        "../envs/environment.yml"
    params:
        engine = config['repeatmasker']['engine'],
        species = config['repeatmasker']['species'],
        custom_lib = config['repeatmasker']['custom_lib'],
        outdir = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/repeatmasker/'
    threads: 8
    resources:
        mem_mb = 24000
    shell:
        """
        mkdir -p {params.outdir}

        if [ -n "{params.custom_lib}" ]; then
            # Run with -lib if param is not empty
            RepeatMasker -pa {threads} \
                -s \
                -xsmall \
                -engine {params.engine} \
                -nocut -gff \
                -lib {params.custom_lib} \
                -dir {params.outdir} \
                {input}
        else
            # Run with -species otherwise
            RepeatMasker -pa {threads} \
                -s \
                -xsmall \
                -engine {params.engine} \
                -nocut -gff \
                -species {params.species} \
                -dir {params.outdir} \
                {input}
        fi
        """

rule filterVcfByCoverage:
    input:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz'
    output:
        filteredVcf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.vcf.gz',
        vcfTemp = temp('output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.temp.vcf.gz')
    params:
        minCov = 9,
        svlen = 16000,
        minlength = 100,
        max_freq = 0.1,
    conda:
        "../envs/environment.yml"
    shell:
        """
        bcftools filter -i 'INFO/SVTYPE="INS" & INFO/SVLEN <= {params.svlen} & INFO/SVLEN >= {params.minlength} & (DR + DV) > {params.minCov} & (DV <= {params.max_freq} * (DR + DV))' {input.vcf} -o {output.vcfTemp} --write-index 
        gzip -c -d {output.vcfTemp} | sed 's/0\/0/1\/1/g' | grep -v "h1tg\|h2tg"  | gzip > {output.filteredVcf}
        """

rule simpleRepeatIntersect:
    input:
        filteredVcf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.vcf.gz',
        repeats = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.all_simple_repeats.bed"
    output:
        SVsInsimpleRepeats = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.simpleRepeats.bed'
    conda:
        "../envs/environment.yml"        
    shell: 
        """
        gzip -c -d {input.filteredVcf} |  bedtools intersect -wo -a - -b {input.repeats} > {output.SVsInsimpleRepeats}
        """

rule filterLongCall:
    input:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/longcall/{specimen}.vcf'
    output:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/longcall/{specimen}_filtered.vcf'
    conda:
        "../envs/environment.yml"
    run:
        with open(input.vcf,'r') as f:
            with open(output.vcf,'w') as outF:
                for line in f:
                    if '#' in line[0]:
                        outF.write(line)
                    else:
                        info = line.strip().split('\t')[-3].split(';')
                        if info[0] == 'SOMATIC' and info[1] == 'MEI':
                            if int(info[4].split('=')[-1]) >= 200:
                                outF.write(line.strip() + '\n')

rule graffiti_sniffles:
    input:
        filteredVcf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.vcf.gz'
    output:
        graffitiOut = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out/3_TSD_search/pangenome.vcf'
    params:
        workDir = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/work',
        outDir = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out',
        fa = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
        custom_lib = config['repeatmasker']['custom_lib'],
        graffitiImage = "/global/scratch/users/landen_gozashti/tools/GraffiTE/graffite_latest.sif"
    handover: False
    threads: 
        4
    shell:
        """
        module load nextflow/24.10.4
        nextflow run -r v1.1dev  cgroza/GraffiTE  -profile standard -w {params.workDir} --cores 4  --vcf {input.filteredVcf} --reference {params.fa}  --TE_library {params.custom_lib} --out {params.outDir} --genotype false  -with-singularity {params.graffitiImage}
        """

rule filterGraffiti_sniffles:
    input:
        graffitiVcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out/3_TSD_search/pangenome.vcf'
    output:
        graffitiFilteredVcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out/3_TSD_search/pangenome_filtered.vcf'
    run:
        def filter_vcf(input_file, output_file):
            """
            Filter VCF file based on:
            - n_hits=1
            - "Satellite" not in the line
            - SVLEN > 250
            """
            
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    # Skip header lines (lines starting with #)
                    if line.startswith('#'):
                        outfile.write(line)
                        continue
                    
                    # Check if line contains "Satellite" - if yes, skip it
                    if "Satellite" in line or "Simple_repeat" in line:
                        continue
                    
                    # Split the line to access the INFO field (column 7, 0-indexed)
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue
        
                    if fields[6].upper() != "PASS":
                        continue
                    
                    info_field = fields[7]
                    
                    # Check for n_hits=1
                    n_hits_found = False
                    svlen_value = 0
                    
                    # Parse INFO field
                    info_parts = info_field.split(';')
                    for part in info_parts:
                        if part.startswith('n_hits='):
                            n_hits_value = part.split('=')[1]
                            if n_hits_value == '1':
                                n_hits_found = True
                        elif part.startswith('SVLEN='):
                            svlen_value = int(part.split('=')[1])
                    
                    # Apply filters
                    if n_hits_found and svlen_value > 250 and "TSD" in info_field:
                        outfile.write(line)
        
        filter_vcf(input.graffitiVcf, output.graffitiFilteredVcf)

rule graffiti_miniSV:
    input:
        filteredVcf = 'output/alignment/scaffolded/minimap2/standard/variants/miniSV/{specimen}.vcf'
    output:
        graffitiOut = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/out/3_TSD_search/pangenome.vcf'
    params:
        workDir = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/work',
        outDir = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/out',
        fa = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta",
        custom_lib = config['repeatmasker']['custom_lib'],
        graffitiImage = "/global/scratch/users/landen_gozashti/tools/GraffiTE/graffite_latest.sif"
    handover: False
    threads:
        4
    shell:
        """
        module load nextflow/24.10.4
        nextflow run -r v1.1dev  cgroza/GraffiTE  -profile standard -w {params.workDir} --cores 4  --vcf {input.filteredVcf} --reference {params.fa}  --TE_library {params.custom_lib} --out {params.outDir} --genotype false  -with-singularity {params.graffitiImage}
        """

rule filterGraffiti_miniSV:
    input:
        graffitiVcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/out/3_TSD_search/pangenome.vcf'
    output:
        graffitiFilteredVcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/out/3_TSD_search/pangenome_filtered.vcf'
    run:
        def filter_vcf(input_file, output_file):
            """
            Filter VCF file based on:
            - n_hits=1
            - "Satellite" not in the line
            - SVLEN > 250
            """

            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    # Skip header lines (lines starting with #)
                    if line.startswith('#'):
                        outfile.write(line)
                        continue

                    # Check if line contains "Satellite" - if yes, skip it
                    if "Satellite" in line or "Simple_repeat" in line:
                        continue

                    # Split the line to access the INFO field (column 7, 0-indexed)
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue

                    if fields[6].upper() != "PASS":
                        continue

                    info_field = fields[7]

                    # Check for n_hits=1
                    n_hits_found = False
                    svlen_value = 0

                    # Parse INFO field
                    info_parts = info_field.split(';')
                    for part in info_parts:
                        if part.startswith('n_hits='):
                            n_hits_value = part.split('=')[1]
                            if n_hits_value == '1':
                                n_hits_found = True
                        elif part.startswith('SVLEN='):
                            svlen_value = int(part.split('=')[1])

                    # Apply filters
                    if n_hits_found and svlen_value > 250 and "TSD" in info_field:
                        outfile.write(line)

        filter_vcf(input.graffitiVcf, output.graffitiFilteredVcf)


rule filterForYoungSniffles:
    input:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out/3_TSD_search/pangenome_filtered.vcf',
    output:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out/3_TSD_search/pangenome_filtered_young.vcf'
    params:
        rmDir = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out/2_Repeat_Filtering/*/repeatmasker_dir/indels.fa.out'
    run:
        import glob
        
        def annotate_and_filter_vcf_with_rm_text(
            vcf_in,
            vcf_out,
            rm_pattern,
            threshold,
            info_tag="RM_PCTDIV",
        ):
            
            #vcf_in:     path to input VCF
            #vcf_out:    path to output (filtered) VCF
            #rm_pattern: glob pattern for RepeatMasker .out files,
            #            e.g. 'PATH/*/repeatmasker_dir/indels.fa.out'
            #threshold:  keep variants with length-weighted percent divergence < threshold
            #info_tag:   INFO key to store the length-weighted percent divergence
            
        
            def parse_repeatmasker_out(path_pattern):
                hits = {}
                for out_fn in glob.glob(path_pattern):
                    with open(out_fn) as fh:
                        for line in fh:
                            line = line.strip()
                            if not line or line.startswith("SW") or line.startswith("score"):
                                continue
                            parts = line.split()
                            if len(parts) < 11:
                                continue
                            try:
                                perc_div = float(parts[1])
                            except ValueError:
                                continue
                            qid = parts[4]
                            q_begin = int(parts[5])
                            q_end = int(parts[6])
                            class_family = parts[10]
                            length = q_end - q_begin + 1
                            # Store percent divergence directly
                            hits.setdefault(qid, []).append((class_family, length, perc_div))
                return hits
        
            def length_weighted_percent_div(hits_for_id, target_classes):
                num = 0.0
                den = 0
                for cls, length, pdiv in hits_for_id:
                    if cls in target_classes:
                        num += pdiv * length
                        den += length
                if den == 0:
                    return None
                return num / den
        
            # Load RepeatMasker hits
            rm_hits = parse_repeatmasker_out(rm_pattern)
        
            header_lines = []
            info_header_present = False
        
            # First pass: read header lines only
            with open(vcf_in) as fin:
                for line in fin:
                    if not line.startswith("#"):
                        break
                    header_lines.append(line.rstrip("\n"))
                    if line.startswith("##INFO=") and f"ID={info_tag}," in line:
                        info_header_present = True
        
            # Prepare output
            with open(vcf_in) as fin, open(vcf_out, "w") as fout:
                # Write header (with added INFO if needed)
                for hl in header_lines:
                    fout.write(hl + "\n")
        
                if not info_header_present:
                    fout.write(
                        f'##INFO=<ID={info_tag},Number=1,Type=Float,'
                        f'Description="Length-weighted percent divergence from RepeatMasker for matching_classes">\n'
                    )
        
                # Rewind to start, skip old header lines, then process records
                for line in fin:
                    if line.startswith("#"):
                        continue
                    break  # line is first data record or EOF
        
                # If the first non-header line was already read, process it
                if not line.startswith("#"):
                    first_record_line = line
                else:
                    first_record_line = None
        
                def process_record(rec_line):
                    rec_line = rec_line.rstrip("\n")
                    if not rec_line or rec_line.startswith("#"):
                        return  # safety
                    cols = rec_line.split("\t")
                    if len(cols) < 8:
                        return
        
                    chrom, pos, vid, ref, alt, qual, flt, info_str = cols[:8]
                    vcf_id = vid
        
                    # Parse INFO into dict
                    info_dict = {}
                    if info_str and info_str != ".":
                        for field in info_str.split(";"):
                            if not field:
                                continue
                            if "=" in field:
                                k, v = field.split("=", 1)
                                info_dict[k] = v
                            else:
                                info_dict[field] = True
        
                    # Get matching_classes
                    mc_val = info_dict.get("matching_classes")
                    if mc_val is None:
                        return
        
                    # matching_classes may have comma-separated values
                    classes = set(mc_val.split(","))
        
                    # No RepeatMasker hits for this ID
                    if vcf_id not in rm_hits:
                        return
        
                    lw_pdiv = length_weighted_percent_div(rm_hits[vcf_id], classes)
                    if lw_pdiv is None:
                        return
        
                    info_dict[info_tag] = f"{lw_pdiv:.3f}"
        
                    # Filter
                    print(vcf_id, lw_pdiv)
                    if lw_pdiv >= threshold:
                        return
        
                    # Rebuild INFO string
                    info_items = []
                    for k, v in info_dict.items():
                        if v is True:
                            info_items.append(k)
                        else:
                            info_items.append(f"{k}={v}")
                    new_info_str = ";".join(info_items) if info_items else "."
        
                    cols[7] = new_info_str
                    fout.write("\t".join(cols) + "\n")
        
                if first_record_line:
                    process_record(first_record_line)
        
                for line in fin:
                    process_record(line)

        annotate_and_filter_vcf_with_rm_text(
        vcf_in=input.vcf,
        vcf_out=output.vcf,
        rm_pattern=params.rmDir,
        threshold=3,           
        info_tag="RM_PCTDIV",
        )


rule filterForYoungMinisv:
    input:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/out/3_TSD_search/pangenome_filtered.vcf',
    output:
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/out/3_TSD_search/pangenome_filtered_young.vcf'
    params:
        rmDir = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/minisv/out/2_Repeat_Filtering/*/repeatmasker_dir/indels.fa.out'
    run:
        import glob
        
        def annotate_and_filter_vcf_with_rm_text(
            vcf_in,
            vcf_out,
            rm_pattern,
            threshold,
            info_tag="RM_PCTDIV",
        ):
            
            #vcf_in:     path to input VCF
            #vcf_out:    path to output (filtered) VCF
            #rm_pattern: glob pattern for RepeatMasker .out files,
            #            e.g. 'PATH/*/repeatmasker_dir/indels.fa.out'
            #threshold:  keep variants with length-weighted percent divergence < threshold
            #info_tag:   INFO key to store the length-weighted percent divergence
            
        
            def parse_repeatmasker_out(path_pattern):
                hits = {}
                for out_fn in glob.glob(path_pattern):
                    with open(out_fn) as fh:
                        for line in fh:
                            line = line.strip()
                            if not line or line.startswith("SW") or line.startswith("score"):
                                continue
                            parts = line.split()
                            if len(parts) < 11:
                                continue
                            try:
                                perc_div = float(parts[1])
                            except ValueError:
                                continue
                            qid = parts[4]
                            q_begin = int(parts[5])
                            q_end = int(parts[6])
                            class_family = parts[10]
                            length = q_end - q_begin + 1
                            # Store percent divergence directly
                            hits.setdefault(qid, []).append((class_family, length, perc_div))
                return hits
        
            def length_weighted_percent_div(hits_for_id, target_classes):
                num = 0.0
                den = 0
                for cls, length, pdiv in hits_for_id:
                    if cls in target_classes:
                        num += pdiv * length
                        den += length
                if den == 0:
                    return None
                return num / den
        
            # Load RepeatMasker hits
            rm_hits = parse_repeatmasker_out(rm_pattern)
        
            header_lines = []
            info_header_present = False
        
            # First pass: read header lines only
            with open(vcf_in) as fin:
                for line in fin:
                    if not line.startswith("#"):
                        break
                    header_lines.append(line.rstrip("\n"))
                    if line.startswith("##INFO=") and f"ID={info_tag}," in line:
                        info_header_present = True
        
            # Prepare output
            with open(vcf_in) as fin, open(vcf_out, "w") as fout:
                # Write header (with added INFO if needed)
                for hl in header_lines:
                    fout.write(hl + "\n")
        
                if not info_header_present:
                    fout.write(
                        f'##INFO=<ID={info_tag},Number=1,Type=Float,'
                        f'Description="Length-weighted percent divergence from RepeatMasker for matching_classes">\n'
                    )
        
                # Rewind to start, skip old header lines, then process records
                for line in fin:
                    if line.startswith("#"):
                        continue
                    break  # line is first data record or EOF
        
                # If the first non-header line was already read, process it
                if not line.startswith("#"):
                    first_record_line = line
                else:
                    first_record_line = None
        
                def process_record(rec_line):
                    rec_line = rec_line.rstrip("\n")
                    if not rec_line or rec_line.startswith("#"):
                        return  # safety
                    cols = rec_line.split("\t")
                    if len(cols) < 8:
                        return
        
                    chrom, pos, vid, ref, alt, qual, flt, info_str = cols[:8]
                    vcf_id = vid
        
                    # Parse INFO into dict
                    info_dict = {}
                    if info_str and info_str != ".":
                        for field in info_str.split(";"):
                            if not field:
                                continue
                            if "=" in field:
                                k, v = field.split("=", 1)
                                info_dict[k] = v
                            else:
                                info_dict[field] = True
        
                    # Get matching_classes
                    mc_val = info_dict.get("matching_classes")
                    if mc_val is None:
                        return
        
                    # matching_classes may have comma-separated values
                    classes = set(mc_val.split(","))
        
                    # No RepeatMasker hits for this ID
                    if vcf_id not in rm_hits:
                        return
        
                    lw_pdiv = length_weighted_percent_div(rm_hits[vcf_id], classes)
                    if lw_pdiv is None:
                        return
        
                    info_dict[info_tag] = f"{lw_pdiv:.3f}"
        
                    # Filter
                    print(vcf_id, lw_pdiv)
                    if lw_pdiv >= threshold:
                        return
        
                    # Rebuild INFO string
                    info_items = []
                    for k, v in info_dict.items():
                        if v is True:
                            info_items.append(k)
                        else:
                            info_items.append(f"{k}={v}")
                    new_info_str = ";".join(info_items) if info_items else "."
        
                    cols[7] = new_info_str
                    fout.write("\t".join(cols) + "\n")
        
                if first_record_line:
                    process_record(first_record_line)
        
                for line in fin:
                    process_record(line)

        annotate_and_filter_vcf_with_rm_text(
        vcf_in=input.vcf,
        vcf_out=output.vcf,
        rm_pattern=params.rmDir,
        threshold=3,           
        info_tag="RM_PCTDIV",
        )



rule liftover:
    input:
        #vcf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.vcf.gz',
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/sniffles/out/3_TSD_search/pangenome_filtered.vcf',
        fa = "output/assembly/hifiasm/{specimen}/scaffolded/{specimen}.diploid.fasta"
    output:
        liftedBed = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.qc_all.covfiltered.lifted.bed',
        bed = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.qc_all.covfiltered.bed',
        header = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}_header.txt'
    params:
        refgenome = config['reference']['fasta'],
        paf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.paf',
    conda:
        "../envs/environment.yml"
    threads:
        14
    shell:
        """
        minimap2 -t14 -c  {params.refgenome} {input.fa} > {params.paf}
        if grep -v -q "#" {input.vcf} ; then
        # Action if match is found
            grep -v "^#"  {input.vcf} | awk '{{print $1"\\t"$2-1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}}' > {output.bed} 
        else
        # Action if NO match is found
            touch {output.bed} 
        fi
        #grep -v "^#"  {input.vcf} | awk '{{print $1"\\t"$2-1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}}' > {output.bed}
        grep "^#" {input.vcf}  > {output.header}
        paftools.js liftover {params.paf} {output.bed} > {output.liftedBed} || true
        """
'''

rule liftover:
    input:
        pafHap1 = "output/assembly/hifiasm/{specimen}/scaffolded/hap1/{specimen}.hap1.scaffold.asm.paf",
        pafHap2 = "output/assembly/hifiasm/{specimen}/scaffolded/hap2/{specimen}.hap2.scaffold.asm.paf",
        vcf = 'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{specimen}/out/3_TSD_search/pangenome_filtered.vcf',
    output:
        liftedBed = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.qc_all.covfiltered.lifted.bed',
        bed = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.qc_all.covfiltered.bed',
        header = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}_header.txt'
    params:
        refgenome = config['reference']['fasta'],
        paf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.paf',
    conda:
        "../envs/environment.yml"
    threads:
        8
    shell:
        """
        cat {input.pafHap1} {input.pafHap2} > {params.paf}
        grep -v "^#"  {input.vcf} | awk '{{print $1"\\t"$2-1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}}' > {output.bed}
        grep "^#" {input.vcf}  > {output.header}
        paftools.js liftover {params.paf} {output.bed} > {output.liftedBed}
        """
'''


rule bedToVcf:
    input:
        liftedBed = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.qc_all.covfiltered.lifted.bed',
        bed = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.qc_all.covfiltered.bed',
        header = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}_header.txt'
    output:
        liftedVcf = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.qc_all.covfiltered.lifted.vcf'
    conda:
        "../envs/environment.yml"
    run:
        import sys
        locToVcfEntry = {}

        with open(input.bed, 'r') as f:
            for line in f:
                if '#' not in line:
                    sp = line.strip().split('\t')
                    entry = '_'.join(sp[:3])
                    vcfLine = sp[3:-1] + ['1/1:' + ':'.join(sp[-1].split(':')[1:])]
                    locToVcfEntry[entry] = vcfLine

        liftoverDic = {}
        with open(input.liftedBed, 'r') as f:
            for line in f:
                sp = line.strip().split('\t')
                entry = sp[3]
                liftoverDic[entry] = sp[:2]

        with open(input.header,'r') as f:
            header = f.read()

        with open(output.liftedVcf,'w') as outF:
            outF.write(header)
            for entry in liftoverDic:
                outF.write('\t'.join(liftoverDic[entry]) + '\t' + '\t'.join(locToVcfEntry[entry]) + '\n')
                




rule survivorMerge:
    input:
        vcfs = expand(f'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{{specimen}}.qc_all.covfiltered.lifted.vcf', specimen=specimens)
    output:
        mergedVCF = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/all_merged.vcf'
    params:
        vcfList = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/lifted_vcfs.txt',
        SURVIVOR_PATH = config["SURVIVOR_path"]
    conda:
        "../envs/environment.yml"
    shell:
        """
        ls {input.vcfs} > {params.vcfList}
        {params.SURVIVOR_PATH}/SURVIVOR merge {params.vcfList} 1000 0 1 0 0 80  {output.mergedVCF}
        """

rule newLiftoverHeader:
    input:
        refgenome = config['reference']['fasta']
    output:
        liftedHeader = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/newHeader.txt' 
    params:
        refgenome = config['reference']['fasta']
    conda:
        "../envs/environment.yml"
    shell:
        """
        samtools faidx {params.refgenome}
        awk 'BEGIN{{OFS=""}} {{print "##contig=<ID=" $1 ",length=" $2 ">"}}' {params.refgenome}.fai > {output.liftedHeader}
        """


rule removeMultiSampleSitesCustom:
    input:
        liftedHeader = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/newHeader.txt',
        mergedVCF = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/all_merged.vcf'
    output:
        mergedVCFFiltered = 'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/all_merged.final.vcf'
    conda:
        "../envs/environment.yml"
    run:
        first = True
        with open(input.liftedHeader,'r') as f:
            header = f.read()

        with open(output.mergedVCFFiltered,'w') as outF:
            with open(input.mergedVCF,'r') as f:
                for line in f:
                    if '#' in line:
                        if "##contig=" in line:
                            if first == True:
                                outF.write(header)
                                first = False
                        else:
                            outF.write(line)
                    else:
                        genotypes = [i.split(':')[0] for i in line.strip().split('\t')[9:]]
                        alleleCount = max([int(i.split('/')[0]) for i in genotypes if i.split('/')[0] !='.'])
                        if genotypes.count('1/1') < 2 and alleleCount < 2:
                            outF.write(line)



