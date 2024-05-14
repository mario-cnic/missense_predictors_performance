#########################################
#   Da formato a las variantes          #
# a analizar para poder pasarlas por    # 
#           VEP web                     #
#########################################

import pandas as pd

in_file = '_in/sarcomeric_filtered_hg38.csv'
out_file = 'VEP_format/sarcomeric_pre_vep_hg38.csv'

data = pd.read_csv(in_file,header=0,sep=',')

variants = data.iloc[:,:4]

variants.insert(loc = 2,
            column = 'ID',
            value = '.')

variants.insert(loc = 5,
            column = 'FILTER',
            value = '.')

variants.insert(loc = 6,
            column = 'INFO',
            value = '.')

variants.insert(loc = 7,
            column = 'FORMAT',
            value = '.')

# Para que pueda ser aceptado por VEP
variants.rename(columns={'CHROM':'#CHROM'},inplace=True)


variants.to_csv(out_file,sep=' ',index=None)
print(f'Se ha creado el archivo {out_file}')


