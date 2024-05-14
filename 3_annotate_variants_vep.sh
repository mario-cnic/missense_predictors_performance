#########################################
#    Anota utilizando vep local las     # 
#   variantes en el archivo con los     #
#     predictores de interés            # 
#                                       #
#########################################

#!/bin/bash

in_file='/home/mruizp/Proyectos/alphamissense/VEP_format/sarcomeric_pre_vep_hg38.csv'
out_file='/home/mruizp/Proyectos/alphamissense/VEP_format/vep_results/sarcomeric_vep_results_hg38_local.tsv'
data3='/home/mruizp/NetVolumes'
ams_path='/home/mruizp/Proyectos/alphamissense/_in/raw/AlphaMissense_hg38.tsv.gz'
# conda deactivate
source activate ensembl-vep

vep --fork 6 --offline --cache -species homo_sapiens --cache_version 109 --dir ${data3}/data3/genomes/homo_sapiens/annotations/VEP \
    --fasta ${data3}/data3/genomes/Homo_sapiens/GATK_bundle/BROAD_hg38/v0/Homo_sapiens_assembly38.fasta --assembly GRCh38 \
    -i ${in_file} \
    -o ${out_file} \
    --tab \
    --force_overwrite \
    --check_existing \
    --everything \
    --no_stats \
    --gene_phenotype \
    --pick \
    --custom file=${data3}/data3/genomes/homo_sapiens/annotations/gnomAD/v4/gnomAD.v4.concat.vcf.gz,short_name=gnomADv4,format=vcf,type=exact,coords=0,fields=AF%AF_AFR%AF_AMR%AF_ASJ%AF_EAS%AF_FIN%AF_NFE%AF_OTH \
    --plugin CADD,"${data3}/data3/genomes/homo_sapiens/annotations/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz" \
    --plugin REVEL,${data3}/data3/genomes/Homo_sapiens/annotations/REVEL/GRCh38/new_tabbed_revel_grch38.tsv.gz \
    --plugin Mastermind,${data3}/data3/genomes/Homo_sapiens/annotations/Mastermind/mastermind_cited_variants_reference-2019.09.27-grch38.vcf.gz \
    --plugin LoFtool,${data3}/data3/genomes/Homo_sapiens/annotations/LoFtool/LoFtool_scores.txt \
    --plugin dbscSNV,${data3}/data3/genomes/Homo_sapiens/annotations/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
    --plugin dbNSFP,"${data3}/data3/genomes/homo_sapiens/annotations/dbNSFP/dbNSFP/dbNSFP4.3a_grch38.gz,ALL" \
    # --plugin AlphaMissense,${ams_path}
    # --plugin FATHMM,"python ${data3}/data3/genomes/Homo_sapiens/annotations/FATHMM/fathmm.py" \ da errores y no lo vamos a usar

    
