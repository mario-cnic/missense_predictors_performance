import pandas as pd 
import os

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

path_myh7 = '_in/raw/variants_sarcomeric'
sarc_genes = read_csv(path_myh7)
# get carriers for each variant
carriers = sarc_genes['Phe'].str.extract(r'Total carriers: \s?(\d+) carriers?').astype(float)

print(f'The total number of carriers (non distinct) is {carriers.sum()[0]}')

