# denovo-TE-insertion-discovery
Pipeline to identify de novo TE insertions from PacBio long-read data

# Getting started

To download the pipeline: 
``` git clone ```

Then,
``` cd denovo-TE-insertion-discovery ```

This directory has the following structure:
```
.
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
```

## Downloading necessary packages that cannot be installed with snakemake

### SURVIVOR
```
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug
make
```

### GraffiTE dev version

First pull the singlularity image:
```
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
apptainer pull --arch amd64 graffite_latest.sif library://cgroza/collection/graffite:latest
```

Then, 
``` git clone https://github.com/cgroza/GraffiTE.git ```
```nextflow pull -r v1.1dev https://github.com/cgroza/GraffiTE ```

## Preparing configs and inputs  

You can leave ```snakemake/config.yaml``` as is, but you'll need to generate dataset-specific template_config.yml and template_samples.tsv before running the pipeline.

For each dataset you plan on analyzing with the pipeline, you'll need to generate a new working directory (this directory can be anywhere). Then specify the full path to this directory in the your pipeline config file (e.g. template_config.yml). In this working directory, you also need to make a directory named ```data```. This directory will how the input data for the pipeline with following structure: ```data/{specimen}/{lane}/{smrtcell}.fastq.gz```. These wildcard (speciment, lane, smrtcell) need to be specified in the sample tsv file (e.g. template_samples.tsv). The columns in this file provide each of these wildcards:

| specimen | group | lane | smrtcell | 
| :--- | :--- | :--- | :--- |

In the config file you'll need to edit ```workdir```, ```sample_table```, ```my_tmp_path```, ```SURVIVOR_path```, ```reference```, and ```custom_lib```.

* workdir: full path to working directory
* sample_table: full path to the sample table
* my_tmp_path: full path to a directory where there is a ton of space that is accessible to compute nodes
* SURVIVOR_path: full path to directory with the survivor executable
* reference: full path to the reference assembly for your species
* custom_lib: full path to a lineage-specific curated TE library for your species

You'll also have to edit one path in rules/postprocessing.smk. (I know this is annoying).

```graffitiImage = "/global/scratch/users/landen_gozashti/tools/GraffiTE/graffite_latest.sif"``` needs to be edited to your local graffite singlularity image.
