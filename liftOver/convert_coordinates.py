import pandas as pd
from pyliftover import LiftOver

def liftO(data: str | pd.DataFrame, old: str= 'hg19', new: str = 'hg38', pos: int = 0, change_pos = True) -> (pd.DataFrame | None):
    """Convierte entre ensamblajes de genomas para SNV

    @Params:
        path: ruta del archivo a convertir o dataframe
        old: assembly del que proviene las variantes
        new: assembly destino
        pos: 1-based o 0-based start
        change_pos: Si `True` crea el mismo archivo, cambiando las coordenadas. 
        Si `False` añade al final del archivo las nuevas coordenadas
    @Returns:
        pd.DataFrame: dataframe con la posición en el nuevo assembly
        """

    lo = LiftOver(old,new)
    if isinstance(data,str):
        hg19_data = pd.read_csv(data,header=0,sep=',')
    else:
        hg19_data = data
    print(hg19_data.head())
    result = hg19_data.apply(lambda row: lo.convert_coordinate('chr'+str(row.CHROM),row.POS)[0],axis=1,result_type='expand')
    result.columns = ['CHROM38','POS38','STRAND','CONV_SCORE']
    result['CHROM38'] = pd.to_numeric(result['CHROM38'].str.replace('chr',''))
    result_merged = pd.merge(hg19_data,result,how='left',left_index=True,right_index=True)

    if change_pos:
        result_merged['POS'] = result_merged['POS38']
        result_merged = result_merged.drop(['STRAND','CONV_SCORE','CHROM38','POS38'],axis=1)
    else:
        result_merged = result_merged.drop(['STRAND','CONV_SCORE','CHROM38'],axis=1)

    return result_merged

if __name__ == '__main__':
    dire = '_in/'
    path = "sarcomeric_filtered.csv"
    new_pos = liftO(dire+path)
    new_pos.to_csv(f'{dire}/{path.replace('.csv','')}_hg38.csv',index=False)