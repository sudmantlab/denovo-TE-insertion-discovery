from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
from pathlib import Path
import glob
import os
import numpy as np
import pandas as pd

configfile: "/global/scratch/users/landen_gozashti/projects/Sperm_diversity/stacy_pipeline/analysis/20251205_sperm_trial/config/pipelineConfigs/template_config.yml"
workdir: config['workdir'] 

refalias= config['reference']['alias'] #correct syntax?

### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen'].unique()
ref_fasta = config['reference']['fasta']

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
include: "rules/postprocessing_reference.smk"



rule all:
    input:
        # self-alignment: assembly + QC, variant calls through qc_all stage
        expand(f"output/alignment/scaffolded/minimap2/standard/coverage_stats/{{specimen}}.coverage.tab", specimen = specimens),
        expand(f"output/assembly/hifiasm/{{specimen}}/scaffolded/{{hap}}/{{specimen}}.{{hap}}.scaffold.canonical.fasta",specimen= specimens, hap = ['hap1','hap2']),
        expand(f"output/assembly/hifiasm/{{specimen}}/scaffolded/{{hap}}/repeatmasker/split_fastas/{{chr}}.fa", specimen= specimens, hap = ['hap1','hap2'], chr = chrs),
        expand(f"output/assembly/hifiasm/{{specimen}}/scaffolded/{{hap}}/repeatmasker/per_chr/{{chr}}.longdust.bed", specimen= specimens, hap = ['hap1','hap2'], chr = chrs,allow_missing = True),
        expand(f"output/assembly/hifiasm/{{specimen}}/scaffolded/{{hap}}/repeatmasker/per_chr/{{chr}}.ultra.bed", specimen= specimens, hap = ['hap1','hap2'], chr = chrs,allow_mising= True),
        expand(f'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{{specimen}}.qc_all.vcf.gz', specimen = specimens),
        expand(f'output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/{{specimen}}.qc_all.covfiltered.vcf.gz', specimen = specimens),
        #expand(f"output/alignment/scaffolded/minimap2/standard/variants/miniSV/{{specimen}}.vcf", specimen = specimens),
        expand(f'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{{specimen}}/sniffles/out/3_TSD_search/pangenome_filtered.vcf', specimen = specimens),
        #expand(f'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{{specimen}}/minisv/out/3_TSD_search/pangenome_filtered.vcf', specimen = specimens),
        expand(f'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{{specimen}}/sniffles/out/3_TSD_search/pangenome_filtered_young_final.vcf', specimen = specimens),
        #expand(f'output/alignment/scaffolded/minimap2/standard/variants/graffiti/{{specimen}}/minisv/out/3_TSD_search/pangenome_filtered_young_final.vcf', specimen = specimens),
        #expand(f'output/alignment/scaffolded/minimap2/standard/variants/longcall/{{specimen}}_filtered.vcf', specimen = specimens),
        "output/alignment/scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/all_merged.final.vcf",
        
        #for reference alignments
        #expand(f'output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/{{specimen}}.qc_all.covfiltered.vcf.gz', specimen = specimens),
        #expand(f'output/alignment/reference/minimap2/standard/variants/graffiti/{{specimen}}/out/3_TSD_search/pangenome_filtered.vcf', specimen = specimens),
        ##expand(f'output/alignment/reference/minimap2/standard/variants/longcall/{{specimen}}_filtered.vcf', specimen = specimens),
        #"output/alignment/reference/minimap2/standard/variants/sniffles_mosaic/all_merged.final.vcf"
