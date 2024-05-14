import pandas as pd
import numpy as np
import gzip


def annotate_variants(ams,sarc):
    """Recorre el archivo tsv de alphamissense y hace un inner join con mis datos de evaluación"""
    sarco = pd.read_csv(sarc,header=0)
    chunksize = 100000
    result = []
    for chunk in pd.read_csv(ams,chunksize=chunksize,skiprows=3,header=0,compression='gzip',sep='\t'):
        chunk.rename(columns={'#CHROM':'CHROM'},inplace=True)
        chunk['CHROM'] = chunk['CHROM'].str.replace('chr','')
        chromvals = list(chunk['CHROM'].unique())
        # print(chromvals,any(x in chromvals for x in ['X','Y','M']))
        if any(x in chromvals for x in ['X','Y','M']): 
            chunk = chunk[(chunk['CHROM'] != 'X') & (chunk['CHROM'] != 'Y') & (chunk['CHROM'] != 'M')] # M = mitocondrial
            print('ELiminando cromosomas sexuales que no son números, como estos no los usamos se ignora (no es lo ideal)')
        chunk['CHROM'] = np.int32(chunk['CHROM']) # debería trabajar con string en lugar de quitarme los chr sexuales
        df = sarco.merge(chunk,on=['CHROM','POS','REF','ALT'],how='inner')
        if (len(df)>0):
            print(f'Anotando variantes en cromosoma {chromvals}')
            result.append(df)
    result = pd.concat(result)
    data = check_missing_variants(sarco,result)
    return data

def check_missing_variants(sarc: pd.DataFrame, res: pd.DataFrame):
    """Comprueba cuántas variantes no se encuentran anotadas por alphamissense y las añade con valores Na al dataframe"""
    data = sarc.merge(res,how='left',on=['CHROM','POS','REF','ALT','Class','Freq','Gene'])
    missing = data[data['am_class'].isna()]
    new_df = pd.concat([res,missing])
    print(f'Hay {len(missing)} variantes que no están en alphamissense, se incluyen también')
    return new_df
if __name__ == '__main__':
    assembly = 'hg38'
    ruta_ams = f'_in/raw/AlphaMissense_{assembly}.tsv.gz'
    ruta_sarcomericos = '_in/sarcomeric_filtered_hg38.csv'
    out = '_in/sarcomeric_annotated_ams.csv' 
    data = annotate_variants(ruta_ams,ruta_sarcomericos)
    data.to_csv(out,header=True,index=False)
    print(f'Archivo escrito en {out}')