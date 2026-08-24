from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
from pathlib import Path
import glob
import os
import numpy as np
import pandas as pd

configfile: "/global/scratch/users/landen_gozashti/projects/Sperm_diversity/stacy_pipeline/analysis/macaque_sperm/config/pipelineConfigs/template_config.yml"
workdir: config['workdir'] 

refalias= config['reference']['alias'] 

### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)

specimens = sample_table['specimen']
groups =sample_table["group"]
ref_fasta = config['reference']['fasta']



bad_groups = [g for g in groups if not re.match(r'^4X[0-9]+$', g)] 
print("Bad groups:", bad_groups)

chrs = []
with open(ref_fasta) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        # Take the first token after '>'
        name = line[1:].split()[0]  # e.g. "chr1" or "chrCAXLPS010000009.1"
        # Keep only "main" chromosomes: chr + number, X, Y (case-insensitive)
        if not name.lower().startswith("chr"):
            continue
        core = name[3:]  # part after 'chr'
        if core.isdigit() or core.upper() in {"X", "Y"}:
            chrs.append(name)

# Optionally deduplicate while preserving order
chrs = list(dict.fromkeys(chrs))

# preprocessing steps
include: "rules/preprocessing.smk"

# Assembly and QC
include: "rules/assembly.smk"

# Alignment (and realignment)
include: "rules/minimap2.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"

# postprocessing
include: "rules/postprocessing.smk"


rule all:
    input:
        # self-alignment: assembly + QC, variant calls through qc_all stage
        expand(f"output/assembly/hifiasm/{{group}}/{{group}}.asm.bp.p_ctg.gfa", group=groups),
        expand("output/alignment/scaffolded/minimap2/standard/coverage_stats/{group}_{specimen}.coverage.tab",zip, group = groups, specimen = specimens),
        expand(f"output/assembly/hifiasm/{{group}}/scaffolded/{{hap}}/{{group}}.{{hap}}.scaffold.canonical.fasta",group= groups, hap = ['hap1','hap2']),


        
        expand(f"output/assembly/hifiasm/{{group}}/scaffolded/{{hap}}/repeatmasker/split_fastas/{{chr}}.fa", group= groups, hap = ['hap1','hap2'], chr = chrs),
        expand(f"output/assembly/hifiasm/{{group}}/scaffolded/{{hap}}/repeatmasker/per_chr/{{chr}}.longdust.bed", group= groups, hap = ['hap1','hap2'], chr = chrs,allow_missing = True),
        expand(f"output/assembly/hifiasm/{{group}}/scaffolded/{{hap}}/repeatmasker/per_chr/{{chr}}.ultra.bed", group= groups, hap = ['hap1','hap2'], chr = chrs,allow_missing = True),
        expand('output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{group}_{specimen}.qc_all.vcf.gz', zip, group = groups, specimen = specimens),
        expand('output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{group}_{specimen}.qc_all.covfiltered.vcf.gz',zip, group = groups, specimen = specimens),
        expand('output/alignment/scaffolded/minimap2/standard/variants/sniffles_standard/{group}_{specimen}.vcf.gz',zip, group = groups, specimen = specimens),
  
        expand('output/alignment/scaffolded/minimap2/standard/variants/graffiti/{group}_{specimen}/sniffles/out/3_TSD_search/pangenome_filtered.vcf',zip, group = groups, specimen = specimens),
  
        expand('output/alignment/scaffolded/minimap2/standard/variants/graffiti/{group}_{specimen}/sniffles/out/3_TSD_search/pangenome_filtered_young_final.vcf',zip, group = groups,  specimen = specimens),

        'output/alignment/scaffolded/minimap2/standard/variants/sniffles_standard/intersection/mergedISEC.vcf',
        "output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/all_merged.final.vcf",
        
   
