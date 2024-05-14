
# ruta='/home/mruizp/NetVolumes/data3/genomes/homo_sapiens/annotations/alphamissense'
ruta2='/home/mruizp/Proyectos/alphamissense/_in/raw/'

gsutil -m cp \
  "gs://dm_alphamissense/AlphaMissense_hg19.tsv.gz" \
  "gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz" \
  ${ruta2}

# crea archivo de index para VEP solo si se va a anotar con el plugin de VEP
tabix -s 1 -b 2 -e 2 -f -S 1 ${ruta2}AlphaMissense_hg19.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 ${ruta2}AlphaMissense_hg38.tsv.gz
