#!/bin/bash                                                                                                               
#SBATCH --job-name=cellRanger --output=cellRanger.out --error=cellRanger.err --time=16:00:00 --qos=normal --nodes=1 --mem-per-cpu=32G --ntasks-per-node=4

#id = output folder name
#path is the path to fastq files
#sample is the prefix of the sample name

/home/groups/oro/software/cellranger-4.0.0/cellranger count --id=D3_WT_cellRanger \
                   --transcriptome=/home/groups/oro/software/refdata-gex-GRCh38-2020-A/ \
                   --fastqs=/oak/stanford/groups/oro/anncoll/scRNA/D3_WT/ \
                   --sample=WTD3 \
