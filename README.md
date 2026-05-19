# denovo-TE-insertion-discovery
Pipeline to identify de novo TE insertions from PacBio long-read data

# Getting started

To download the pipeline: 
``` git clone ```

Then,
``` cd denovo-TE-insertion-discovery ```

This directory has the following structure:
├── env/
│   ├── environment.yml
│   ├── longdust.yaml
│   └──other_masking.yaml
├── rules/
│   ├── sniffles.smk
│   ├── minimap2.smk
│   ├── samtools_utils.smk
│   ├── postprocessing.smk
│   ├── postprocessing_reference.smk
│   ├── coverage_stats.smk
│   ├── preprocessing.smk
│   └── assembly.smk
└── config/
    ├── snakemake
        └── config.yaml
    ├── pipelineConfigs/
        ├── template_config.yml
        └── template_samples.tsv
