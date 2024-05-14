import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import plotly.express as px
import re
import numpy as np

def sensitivity(validator: pd.Series,predictor: pd.Series, binary_class: list = ['likely_benign','likely_pathogenic']) -> int:
    """Calcula la sensibilidad de un predictor.

    Parameters:
        validator: pd.Series
        predictor: pd.Series
        binary_class: list = lista de tamaño 2 donde el primer elemento indica la clase 0 y el segundo la clase 1

    Returns:
        Float
    """
    true_pos = (validator == 0) & (predictor == binary_class[0])
    true_pos_sum = sum(true_pos)
    false_neg = (validator == 0) & (predictor == binary_class[1])
    print(true_pos_sum)
    false_neg_sum = sum(false_neg)
    sens = true_pos_sum / (true_pos_sum + false_neg_sum)

    # print(f"La sensibilidad es:\n{round(sens,2)}")
    return sens

def specificity(validator: pd.Series, predictor: pd.Series, binary_class: list = ['likely_benign','likely_pathogenic']) -> int:
    """Calcula la especificidad de un predictor.
    
    Params:
        validator: pd.Series
        predictor: pd.Series
        binary_class: list = lista de tamaño 2 donde el primer elemento indica la clase 0 y el segundo la clase 1
        
    Returns:
        Float
    """
    true_neg = (validator == 1) & (predictor == binary_class[1])
    true_neg_sum = sum(true_neg)
    false_pos = (validator == 1) & (predictor == binary_class[0])
    false_pos_sum = sum(false_pos)
    spec = true_neg_sum / (true_neg_sum + false_pos_sum)
    # print(f"La especificidad es:\n{round(spec,2)}")
    return spec


def plot_curves(data,curve_type: str = 'roc'):
    i = 0
    # Dynamically choose colors
    colors = plt.cm.rainbow(np.linspace(0, 1, len(data)))
    plt.figure(figsize=(10, 10))
    for k,v in data.items():
        if curve_type == 'pr':
            y,x,threshold = precision_recall_curve(v.iloc[:,0],v.iloc[:,1])
            title = 'Precision-Recall (PR) Curve for classifying variants in sarcomeric genes'
            x_lab = 'Precision'
            y_lab = 'Recall'
            line = [[0, 1], [0.5, 0.5]]
        elif curve_type == 'roc':
            x,y,threshold = roc_curve(v.iloc[:,0],v.iloc[:,1])
            title = 'ROC curve for classifying variants in sarcomeric genes'
            x_lab = 'False Positive Rate (1 - specificity)'
            y_lab = 'True Positive Rate (sensibility)'
            line = [[0, 1], [0, 1]]
        pr_auc = auc(x,y)
        label_name = re.sub("_.+$","",k) if k != 'am_pathogenicity' else 'AlphaMissense'
        plt.plot(x,y,color=colors[i], lw=2, label='{}(area = {:.2f})'.format(label_name,pr_auc))
        i += 1
    else:
        plt.plot(line[0],line[1], color='navy', lw=2, linestyle='--')
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.title(title)
        plt.legend(loc='lower right')
        plt.savefig(f'analysis/out/{curve_type}.png',bbox_inches='tight')
        plt.show()

def performance_scores(data):
    perf_scores = dict()
    auc_score = dict()
    for k,v in data.items():
        fpr,tpr,t = roc_curve(v.iloc[:,0],v.iloc[:,1])
        rc,pr,t2 = precision_recall_curve(v.iloc[:,0],v.iloc[:,1])
        # print(k,len(t),len(t))
        t2 =np.append(t2,1)
        perf_scores[k] = {'FPR': fpr,'TPR':tpr,'Thresholds': t,'Precision':pr,'Recall':rc,'Thresholds_pr':t2}
        auc_score[k] = {'roc':auc(fpr,tpr),'pr':auc(pr,rc)}
        print(k)
        # print(len(pr),len(t2))
    auc_df = pd.DataFrame.from_dict(auc_score)
    
    roc_df = pd.DataFrame()
    pr_df = pd.DataFrame()

    for predictor, output in perf_scores.items():
        # print(predictor,
        #       len(output['FPR']),len(output['TPR']),len(output['Thresholds']),
        #       len(output['Precision']),len(output['Recall']),len(output['Thresholds_pr']))
        # print(output['Thresholds'][0],output['Thresholds'][-1],type(output['Thresholds_pr']))
        temp_df = pd.DataFrame({
            'Predictor': predictor,
            'FPR': output['FPR'],
            'TPR': output['TPR'],
            'Thresholds': output['Thresholds'],
            
        })
        temp2_df = pd.DataFrame({
            'Predictor': predictor,
            'Precision': output['Precision'],
            'Recall': output['Recall'],
            'Thresholds': output['Thresholds_pr']
        })

        roc_df = pd.concat([roc_df, temp_df])
        pr_df = pd.concat([pr_df,temp2_df])

    return roc_df,pr_df,auc_df

def plotly_curves(data,curve_type: str = 'roc'):
    curves = []  # List to store curve traces
    colors = px.colors.qualitative.Plotly
    i = 0
    for k,v in data.items():
        if curve_type == 'pr':
            y,x,threshold = precision_recall_curve(v.iloc[:,0],v.iloc[:,1])
            title = 'Precision-Recall (PR) Curve for classifying variants in sarcomeric genes'
            x_lab = 'Precision'
            y_lab = 'Recall'
            line = [[0, 1], [0.5, 0.5]]
        elif curve_type == 'roc':
            x,y,threshold = roc_curve(v.iloc[:,0],v.iloc[:,1])
            title = 'ROC curve for classifying variants in sarcomeric genes'
            x_lab = 'False Positive Rate (1 - specificity)'
            y_lab = 'True Positive Rate (sensitivity)'
            line = [[0, 1], [0, 1]]
        pr_auc = auc(x,y)
        label_name = re.sub("_.+$","",k) if k != 'am_pathogenicity' else 'AlphaMissense'
        # Create ROC curve trace
        curve_trace = go.Scatter(x=x, y=y,
                                    mode='lines',
                                    line=dict(color=colors[i % len(colors)]),
                                    name='{} ({:.3f})'.format(label_name, pr_auc)
                                    )
        curves.append(curve_trace)  # Append current curve trace to the list
        i+=1

    # Create layout
    layout = go.Layout(title=title,
                    xaxis=dict(title=x_lab),
                    yaxis=dict(title=y_lab),
                    showlegend=True)
    # Create figure
    fig = go.Figure(data=curves, layout=layout)

    fig.update_layout(title_text=title,
                    xaxis_title=x_lab,
                    yaxis_title=y_lab,
                    legend_title="AUC",
                    width = 800,
                    height = 700,
                    # margin=dict(l=20, r=20, t=20, b=20),
                    template='simple_white')
    fig.show()
    # fig.write_html('out/')