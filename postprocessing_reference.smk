refalias = config['reference']['alias']

rule filterVcfByCoverage_reference:
    input:
        vcf = 'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz'
    output:
        filteredVcf = 'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.vcf.gz',
        vcfTemp = temp('output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.covfiltered.temp.vcf.gz')
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


rule filterLongCall_reference:
    input:
        vcf = 'output/alignment/reference/minimap2/standard/variants/longcall/{specimen}.vcf'
    output:
        vcf = 'output/alignment/reference/minimap2/standard/variants/longcall/{specimen}_filtered.vcf'
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

rule graffiti_reference:
    input:
        filteredVcf = f'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{{specimen}}.qc_all.covfiltered.vcf.gz'
    output:
        graffitiOut = f'output/alignment/reference/minimap2/standard/variants/graffiti/{{specimen}}/out/3_TSD_search/pangenome.vcf'
    params:
        workDir = f'output/alignment/reference/minimap2/standard/variants/graffiti/{{specimen}}/work',
        outDir = f'output/alignment/reference/minimap2/standard/variants/graffiti/{{specimen}}/out',
        fa = config['reference']['fasta'],
        custom_lib = config['repeatmasker']['custom_lib'],
        graffitiImage = "/global/scratch/users/landen_gozashti/tools/GraffiTE/graffite_latest.sif"
    handover: False
    shell:
        """
        module load nextflow/24.10.4
        nextflow run -r v1.1dev  cgroza/GraffiTE  -profile standard -w {params.workDir} --cores 4  --vcf {input.filteredVcf} --reference {params.fa}  --TE_library {params.custom_lib} --out {params.outDir} --genotype false  -with-singularity {params.graffitiImage}
        """

rule filterGraffiti_reference:
    input:
        graffitiVcf = 'output/alignment/reference/minimap2/standard/variants/graffiti/{specimen}/out/3_TSD_search/pangenome.vcf'
    output:
        graffitiFilteredVcf = 'output/alignment/reference/minimap2/standard/variants/graffiti/{specimen}/out/3_TSD_search/pangenome_filtered.vcf'
    run:
        def filter_vcf(input_file, output_file):
            
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


                



rule survivorMerge_reference:
    input:
        vcfs = expand(f'output/alignment/reference/minimap2/standard/variants/graffiti/{{specimen}}/out/3_TSD_search/pangenome_filtered.vcf',specimen=specimens)
    output:
        mergedVCF = 'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/all_merged.vcf'
    params:
        vcfList = 'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/lifted_vcfs.txt',
        SURVIVOR_PATH = config["SURVIVOR_path"]
    conda:
        "../envs/environment.yml"
    shell:
        """
        ls {input.vcfs} > {params.vcfList}
        {params.SURVIVOR_PATH}/SURVIVOR merge {params.vcfList} 1000 0 1 0 0 80  {output.mergedVCF}
        """



rule removeMultiSampleSitesCustom_reference:
    input:
        mergedVCF = 'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/all_merged.vcf'
    output:
        mergedVCFFiltered = 'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/all_merged.final.vcf'
    conda:
        "../envs/environment.yml"
    run:
        first = True

        with open(output.mergedVCFFiltered,'w') as outF:
            with open(input.mergedVCF,'r') as f:
                for line in f:
                    if '#' in line:
                        outF.write(line)
                    else:
                        genotypes = [i.split(':')[0] for i in line.strip().split('\t')[9:]]
                        alleleCount = max([int(i.split('/')[0]) for i in genotypes if i.split('/')[0] !='.'])
                        if genotypes.count('1/1') < 2 and alleleCount < 2:
                            outF.write(line)


