#########################################
# Procesa los archivos con las variantes# 
#                                       #
#########################################
import pandas as pd
import numpy as np
import os
import re
from liftOver.convert_coordinates import liftO

def read_csv(ruta):
    """Write into pd.dataframe the file passed or csv files in folder.
    Returns dataframe"""
    if os.path.isfile(ruta):
        raw_data = pd.read_csv(ruta,header=0,sep=';')
    else:
        dfs = []
        for file in os.listdir(ruta):
            df = pd.read_csv(ruta+'/'+file,header=0,sep=';')
            dfs.append(df)
            print(file)
        raw_data = pd.concat(dfs,ignore_index=True)
        raw_data = raw_data.drop_duplicates()

    return raw_data
    

def data_preprocessing(raw_data, filter_cols = True, inclusion_criteria = 'S',spanish = False):
    """inclusion criteria (Strong/Weak) if strong only include pathogenic benign, if weak include likeky pathogenic/likely benign"""

    raw_data.rename(columns = 
                    {'Chr':'CHROM',
                    'Chromosomic name':'POS-REF-ALT',
                    'Pat.':'Class',
                    'Freq.':'Freq',
                    'Gene name': 'Gene'}, inplace=True)
    if spanish:
        raw_data.rename(columns = 
                    {'Chr':'CHROM',
                    'Nombre cromosómico':'POS-REF-ALT',
                    'Pat.':'Class',
                    'Frec.':'Freq',
                    'Nombre del gen': 'Gene'}, inplace=True)

    
    raw_data['Class'] = raw_data['Class'].str.strip()
    # eliminamos variantes sin clasificar
    nullv = raw_data['Class'].isna()
    # print(raw_data.loc[raw_data['Class'].isna(),['Class']].head())
    raw_data = raw_data.loc[~raw_data['Class'].isna(),:]
    print(f'Se han eliminado {sum(nullv)} variantes que no tenían clase')
    # selccionamos variantes 
    if inclusion_criteria == 'S':
        # P & B
        raw_data['Class'] = np.where(raw_data['Class'].str.startswith('--'),
                            0,
                            np.where(raw_data['Class'].str.startswith('+++'),
                                    1,
                                    0.5)
                            )
    elif inclusion_criteria == 'W':
        # LP/P & LB/B
        raw_data.loc[:,'Class'] = np.where(raw_data['Class'].str.contains(r'-[-\|?]'),
                            0,
                            np.where(raw_data['Class'].str.startswith('++'),
                                    1,
                                    0.5)
                            )    
    else:
        raise ValueError(f'El valor {inclusion_criteria} no está definido')
    # si la variante tiene nombre diferente: NC_000014.8:g.23898525_23898527delCAGinsAAC la eliminamos (de momento son solo 4, TODO: corregir)
    rare_data = raw_data['POS-REF-ALT'].str.startswith('NC')
    raw_data = raw_data.loc[~rare_data,:]
    print(f'Se han eliminado {sum(rare_data)} variantes por tener nombres raros, esto no debería ocurrir')

    # separamos Posición alelo ref y alelo alt en columnas independientes
    raw_data[['POS-REF','ALT']] = raw_data['POS-REF-ALT'].str[2::].str.split('>',expand = True)
    raw_data['REF']= raw_data['POS-REF'].str[-1]
    raw_data['POS'] = raw_data['POS-REF'].str[:-1]

    # simplificamos la columna de frecuencias
    raw_data['Freq'] = raw_data.Freq.str.strip()
    raw_data['Freq'] = np.where(raw_data['Freq'].str.startswith('?'),None,
                                np.where(raw_data['Freq'].str.startswith('>'),'>0.01',
                                np.where(raw_data['Freq'].str.startswith('<'),'<0.01', '0')
                                )
                            )
    # eliminamos columnas que ya han sido parseadas
    raw_data = raw_data.drop(columns=['POS-REF-ALT','POS-REF'])
    # obtenemos el nombre del gen
    raw_data['Gene'] = raw_data['Gene'].str.strip()
    raw_data['Gene'] = raw_data['Gene'].str.split(' ',expand=True)[0]
    if filter_cols:
        # seleccionamos columnas de interés
        raw_data = raw_data.loc[:,['CHROM','POS','REF','ALT','Class','Freq','Gene']]
        
    # ordenamos por cromosoma y posición
    raw_data.sort_values(by=['CHROM','POS'],inplace=True)
    # eliminamos indice de filas
    return raw_data

def save_dataframe(df,name = '_in/file.csv'):
    print(f'Archivo guardado en {os.getcwd()+'/_in'}')
    df.to_csv(name,index=False)


if __name__ == '__main__':

    path_sarc = '_in/raw/variants_sarcomeric'
    out_name = '_in/sarcomeric_filtered.csv'
    inclusion_criteria = 'W' # incluir solo patogénicas y benignas o también las likely
    spanish = False # si el archivo original tiene headers en español o inglés
    filter_cols = True # incluye todas las columnas o solo las de interés
    df = read_csv(path_sarc)
    clean_df = data_preprocessing(df,filter_cols,inclusion_criteria,spanish)
    # df_hg38 = liftO(clean_df) # convierte a coordenadas hg38
    save_dataframe(clean_df,name = out_name)
