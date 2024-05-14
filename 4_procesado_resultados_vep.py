#########################################
# Procesa el archivo txt obtenido de    # 
# analizar las variantes mediante VEP   # 
#                                       #
#########################################


import pandas as pd
import numpy as np

def read_header(file,only_header = False):
	values = None
	with open(file, 'r') as f:
		for line in f:
			listline = line.strip().split('\t')
			# if len(listline) != 740 and line[:2] != '##':
			# 	print(line)
			if line.startswith('#Uploaded_variation'):
				listline[0] = listline[0][1:]
				if only_header:
					return listline
				# creamos diccionario
				values = {x:[] for x in listline}
			elif not line.startswith('##') and values is not None:
				for k, v in zip(values.keys(), listline):
					values[k].append(v)
	if values is not None:
        # Convert the dictionary to a DataFrame
		df = pd.DataFrame(values)
		df.replace('-','',inplace=True)
		return df
	else:
		return None 
     
directorio = 'VEP_format/vep_results/'
in_file = directorio+'sarcomeric_vep_results_hg38_local.tsv'
out_sufix = 'vep_prediction_scores.csv'

# importamos 
#vep_variants = read_header(in_file)

# vep_variants = pd.read_csv(in_file, delimiter='\t', comment='#', header=None)
vep_variants = pd.read_csv('VEP_format/vep_results/sarcomeric_vep_results_hg38_local.tsv',delimiter='\t',comment = '#',low_memory=False,header=None)
vep_variants.columns = read_header(in_file,only_header=True)
vep_variants.to_csv('_other_files/proper_format_vep.csv',index=None)
# obtenemos el cromosoma y la posición
vep_variants[['CHROM','POS']] = vep_variants['Location'].str.extract(r'^(\d+):(\d+)$')
# obtenemos la frecuencia alélica máxima de gnomAD exomas
vep_variants['MAX_AF'] = vep_variants.loc[:,['AF', 'gnomADe_AF', 'gnomADe_AFR_AF', 'gnomADe_AMR_AF',
       'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF',
       'gnomADe_OTH_AF', 'gnomADe_SAS_AF']].max(axis='columns')


# seleccionamos las columnas de interés
columns = ['CHROM','POS','Allele','PolyPhen','MAX_AF','gnomAD_exomes_AF','gnomADe_AF','CADD_phred','CADD_RAW','FATHMM_score','MutPred_score',
       'MutationTaster_score','Polyphen2_HDIV_score', 'Polyphen2_HVAR_score','REVEL_score','SIFT4G_score','SIFT_score','DANN_score',
	   'MetaLR_score','MetaRNN_score','PrimateAI_score','gnomADv4_AF']

# # escogemos el valor máximo de las columnas en forma de string (intuyo que devuelven 1 resultado por transcrito aunque todas arrojan lo mismo)
parse_cols = ['REVEL_score','MutationTaster_score','FATHMM_score','SIFT4G_score','SIFT_score','Polyphen2_HVAR_score','Polyphen2_HDIV_score','MetaRNN_score']
for col in parse_cols:
	vep_variants[col] = vep_variants[col].str.strip()
	vep_variants[col] = vep_variants[col].str.split(',')
	vep_variants[col] = vep_variants[col].apply(lambda x: max([float(i) if i not in ['.','-'] else np.nan for i in x ]))

df = vep_variants.loc[:,columns]
df.rename(columns={'Allele':'ALT'},inplace=True)

# print(df.head())

df.to_csv(directorio+out_sufix,index=None)
print(f'Anotaciones guardadas en {directorio}{out_sufix}')