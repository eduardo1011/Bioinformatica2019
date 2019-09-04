import os
import requests
import sys
import urllib.request
from urllib.request import urlopen
import re
from pandas import DataFrame
import pandas as pd
from pandas.compat import StringIO
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
import json
from PIL import Image
import webbrowser
#from wordcloud import WordCloud, STOPWORDS
import subprocess
from datetime import datetime
import networkx as nx
import nxviz as nxv

#####
from matplotlib import cm
words = ['Biological Process', 'Molecular Function', 'Cellular Component']
pie_colors = {'Set3':cm.Set3(np.arange(12)/12.),
              'Set2':cm.Set2(np.arange(8)/8.),
              'Set1':cm.Set1(np.arange(9)/9.),
              'Pastel2':cm.Pastel2(np.arange(8)/8.),
              'Pastel1':cm.Pastel1(np.arange(9)/9.),
              'Dark2':cm.Dark2(np.arange(8)/8.),
              'Paired':cm.Paired(np.arange(12)/12.),
              'Accent':cm.Accent(np.arange(8)/8.),
              'Spectral':cm.Spectral(np.arange(11)/11.),
              'tab20':cm.tab20(np.arange(20)/20.),
              'tab20b':cm.tab20b(np.arange(20)/20.),
              'tab20c':cm.tab20c(np.arange(20)/20.)}
colors = {'#8DD3C7':pie_colors['Set3'][0:1], '#FFFFB3':pie_colors['Set3'][1:2],
         '#BEBADA':pie_colors['Set3'][2:3], '#FB8072':pie_colors['Set3'][3:4],
         '#80B1D3':pie_colors['Set3'][4:5], '#FDB462':pie_colors['Set3'][5:6],
         '#B3DE69':pie_colors['Set3'][6:7], '#FCCDE5':pie_colors['Set3'][7:8],
         '#D9D9D9':pie_colors['Set3'][8:9], '#BC80BD':pie_colors['Set3'][9:10],
         '#CCEBC5':pie_colors['Set3'][10:11], '#FFED6F':pie_colors['Set3'][11:12]}
sizes = ['3x3','4x4','5x5','6x6','7x7','8x8','9x9','10x10']
escala_uno = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]


#from __future__ import division
#import matplotlib.pyplot as plt
#from matplotlib import colors as mcolors
#colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
# lista de colores python
common_colors = {'blue': 'blue',
    'green': 'green',
    'red': 'red',
    'cyan': 'cyan',
    'magenta': 'magenta',
    'yellow': 'yellow',
    'black': 'black',
    'white': 'white',
    'aliceblue': '#F0F8FF',
    'antiquewhite': '#FAEBD7',
    'aqua': '#00FFFF',
    'aquamarine': '#7FFFD4',
    'azure': '#F0FFFF',
    'beige': '#F5F5DC',
    'bisque': '#FFE4C4',
    'black': '#000000',
    'blanchedalmond': '#FFEBCD',
    'blue': '#0000FF',
    'blueviolet': '#8A2BE2',
    'brown': '#A52A2A',
    'burlywood': '#DEB887',
    'cadetblue': '#5F9EA0',
    'chartreuse': '#7FFF00',
    'chocolate': '#D2691E',
    'coral': '#FF7F50',
    'cornflowerblue': '#6495ED',
    'cornsilk': '#FFF8DC',
    'crimson': '#DC143C',
    'cyan': '#00FFFF',
    'darkblue': '#00008B',
    'darkcyan': '#008B8B',
    'darkgoldenrod': '#B8860B',
    'darkgray': '#A9A9A9',
    'darkgreen': '#006400',
    'darkgrey': '#A9A9A9',
    'darkkhaki': '#BDB76B',
    'darkmagenta': '#8B008B',
    'darkolivegreen': '#556B2F',
    'darkorange': '#FF8C00',
    'darkorchid': '#9932CC',
    'darkred': '#8B0000',
    'darksalmon': '#E9967A',
    'darkseagreen': '#8FBC8F',
    'darkslateblue': '#483D8B',
    'darkslategray': '#2F4F4F',
    'darkslategrey': '#2F4F4F',
    'darkturquoise': '#00CED1',
    'darkviolet': '#9400D3',
    'deeppink': '#FF1493',
    'deepskyblue': '#00BFFF',
    'dimgray': '#696969',
    'dimgrey': '#696969',
    'dodgerblue': '#1E90FF',
    'firebrick': '#B22222',
    'floralwhite': '#FFFAF0',
    'forestgreen': '#228B22',
    'fuchsia': '#FF00FF',
    'gainsboro': '#DCDCDC',
    'ghostwhite': '#F8F8FF',
    'gold': '#FFD700',
    'goldenrod': '#DAA520',
    'gray': '#808080',
    'green': '#008000',
    'greenyellow': '#ADFF2F',
    'grey': '#808080',
    'honeydew': '#F0FFF0',
    'hotpink': '#FF69B4',
    'indianred': '#CD5C5C',
    'indigo': '#4B0082',
    'ivory': '#FFFFF0',
    'khaki': '#F0E68C',
    'lavender': '#E6E6FA',
    'lavenderblush': '#FFF0F5',
    'lawngreen': '#7CFC00',
    'lemonchiffon': '#FFFACD',
    'lightblue': '#ADD8E6',
    'lightcoral': '#F08080',
    'lightcyan': '#E0FFFF',
    'lightgoldenrodyellow': '#FAFAD2',
    'lightgray': '#D3D3D3',
    'lightgreen': '#90EE90',
    'lightgrey': '#D3D3D3',
    'lightpink': '#FFB6C1',
    'lightsalmon': '#FFA07A',
    'lightseagreen': '#20B2AA',
    'lightskyblue': '#87CEFA',
    'lightslategray': '#778899',
    'lightslategrey': '#778899',
    'lightsteelblue': '#B0C4DE',
    'lightyellow': '#FFFFE0',
    'lime': '#00FF00',
    'limegreen': '#32CD32',
    'linen': '#FAF0E6',
    'magenta': '#FF00FF',
    'maroon': '#800000',
    'mediumaquamarine': '#66CDAA',
    'mediumblue': '#0000CD',
    'mediumorchid': '#BA55D3',
    'mediumpurple': '#9370DB',
    'mediumseagreen': '#3CB371',
    'mediumslateblue': '#7B68EE',
    'mediumspringgreen': '#00FA9A',
    'mediumturquoise': '#48D1CC',
    'mediumvioletred': '#C71585',
    'midnightblue': '#191970',
    'mintcream': '#F5FFFA',
    'mistyrose': '#FFE4E1',
    'moccasin': '#FFE4B5',
    'navajowhite': '#FFDEAD',
    'navy': '#000080',
    'oldlace': '#FDF5E6',
    'olive': '#808000',
    'olivedrab': '#6B8E23',
    'orange': '#FFA500',
    'orangered': '#FF4500',
    'orchid': '#DA70D6',
    'palegoldenrod': '#EEE8AA',
    'palegreen': '#98FB98',
    'paleturquoise': '#AFEEEE',
    'palevioletred': '#DB7093',
    'papayawhip': '#FFEFD5',
    'peachpuff': '#FFDAB9',
    'peru': '#CD853F',
    'pink': '#FFC0CB',
    'plum': '#DDA0DD',
    'powderblue': '#B0E0E6',
    'purple': '#800080',
    'rebeccapurple': '#663399',
    'red': '#FF0000',
    'rosybrown': '#BC8F8F',
    'royalblue': '#4169E1',
    'saddlebrown': '#8B4513',
    'salmon': '#FA8072',
    'sandybrown': '#F4A460',
    'seagreen': '#2E8B57',
    'seashell': '#FFF5EE',
    'sienna': '#A0522D',
    'silver': '#C0C0C0',
    'skyblue': '#87CEEB',
    'slateblue': '#6A5ACD',
    'slategray': '#708090',
    'slategrey': '#708090',
    'snow': '#FFFAFA',
    'springgreen': '#00FF7F',
    'steelblue': '#4682B4',
    'tan': '#D2B48C',
    'teal': '#008080',
    'thistle': '#D8BFD8',
    'tomato': '#FF6347',
    'turquoise': '#40E0D0',
    'violet': '#EE82EE',
    'wheat': '#F5DEB3',
    'white': '#FFFFFF',
    'whitesmoke': '#F5F5F5',
    'yellow': '#FFFF00',
    'yellowgreen': '#9ACD32'}
#####


###
def barras(df = DataFrame([]), column = 0, dim = 111, title = '', row_num = 10, color = 'salmon',
           size_x = 9, size_y = 15, xlabel = 10, size_title = 25, size_bartxt = 12, sep = 3):
    if len(df) == 0:
        print('Data frame sin datos')
    else:
        plt.subplot(dim, facecolor= 'white')
        barWidth = 0.9
        if row_num == len(df):
            ejey = list(df.iloc[0:len(df),3])
            val = max(ejey)
            ejex = list(df.iloc[0:len(df),column])
            borde = list(np.repeat('white', len(ejey)))
            colores = list(np.repeat(color, len(ejey)))
            linea = list(np.repeat(0, len(ejey)))
        if row_num < len(df):
            ejey = list(df.iloc[0:row_num,3]) + [df.iloc[row_num:len(df),3].sum()]
            val = max(ejey)
            ejex = list(df.iloc[0:row_num,column]) + ['Others']
            borde = list(np.repeat('white', len(ejey)-1)) + ['black']
            colores = list(np.repeat(color, len(ejey)-1)) + ['linen']
            linea = list(np.repeat(0, len(ejey)-1)) + [1]
        if row_num > len(df):
            ejey = list(df.iloc[0:len(df),3])
            val = max(ejey)
            ejex = list(df.iloc[0:len(df),column])
            borde = list(np.repeat('white', len(ejey)))
            colores = list(np.repeat(color, len(ejey)))
            linea = list(np.repeat(0, len(ejey)))
    
        for i, j, k, l, m in zip(ejex , ejey, borde, colores, linea):
            plt.barh(i, j,
                     color= l,
                     align='center',
                     height= 0.7,linewidth = m,
                     alpha = 1, edgecolor= k)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_position(('data',-0.6))
        plt.gca().spines['left'].set_visible(False)
        plt.title(title,size= size_title,loc='left')
        plt.tick_params(axis="y", color="gray")
        plt.yticks(size = size_y )

        #if max(list(ejey)) >= 151:
        #    escala = 50
        #if 101 < max(list(ejey)) <= 150:
        #    escala = 20
        #if 51 <= max(list(ejey)) <= 100:
        #    escala = 10
        #if max(list(ejey)) <= 50:
        #    escala = 5
            
        v1 = -50
        v2 = 0
        v3 = 0
        for i in range(10000):
            v1 += 50
            v2 += 50
            v3 += 10
            if v1 <= max(list(ejey)) < v2:
                #print(v3, v1, val, v2)
                escala = v3
        

        plt.xticks(range(0,val,escala),size=size_x) #fontweight='bold'
        plt.xlabel("Number of Proteins",size= xlabel)

        for j, k in zip(ejey,range(0,len(ejey))):
            plt.text(j+sep,k-0.2,j,
                     size=size_bartxt,ha='left',color= 'black')
###

def pastel(df = DataFrame([]), column = 0, dim = 111, title = '',
           row_num = 7, angle = 0, legend = 'Box', # lagend options = 'Box' y 'Labels'
           size_text = 10, size_title = 20, color = pie_colors['Set2'],
           explode = 0, open_center = 0.4, box_col = 1, rot_lab = True, alpha = 1,
              lab_dist = 0.1):
    if len(df) < 1:
        print('Data frame sin datos')
    else:
        if row_num == len(df):
            ejey = list(df.iloc[0:len(df),3])
            val = max(ejey)
            ejex = list(df.iloc[0:len(df),column])
            colores = list(color[0:len(ejey)])
        if row_num < len(df):
            cc1 = list(color) * 20
            ejey = list(df.iloc[0:row_num,3]) + [df.iloc[row_num:len(df),3].sum()]
            val = max(ejey)
            ejex = list(df.iloc[0:row_num,column]) + ['Others']
            colores =  list(cc1[0:row_num]) + ['linen']
        if row_num > len(df):
            ejey = list(df.iloc[0:len(df),3])
            val = max(ejey)
            ejex = list(df.iloc[0:len(df),column])
            colores = list(color[0:len(ejey)])
        
        
        explo = explode
        porcentaje = np.round(list((np.array(ejey) / sum(ejey)) * 100), 1)
        labb = ejex
        labels_box = [i+' ('+str(j)+'%)' for i, j in zip(labb, porcentaje)]
        labels_lab = ejex
        opcion_leyenda = {'Box':[labels_box, None, None, (0.95 + (explo / 2), 0.95 + (explo / 2)), size_text],
                          'Labels':['', labels_lab, '%1.1f%%', (0, 0), 0],
                          'Rotate':['', labels_box, None, (0, 0), 0]}
        #plt.subplot(dim)
        explode0 = np.repeat(explode, len(ejey))
        plt.pie(ejey, labels = opcion_leyenda[legend][1],
                labeldistance= 1,
                wedgeprops=dict(edgecolor=None, alpha = alpha),
                autopct = opcion_leyenda[legend][2], startangle=angle, radius = 1,
                colors = colores, rotatelabels = rot_lab,
                explode = explode0, textprops=dict(size = size_text))
            
        plt.legend(opcion_leyenda[legend][0], bbox_to_anchor= opcion_leyenda[legend][3],
                   loc=2, labelspacing = 0.1,
                   borderaxespad=None, fontsize = opcion_leyenda[legend][4],
                   ncol= box_col, frameon=False)
        centre_circle = plt.Circle((0,0),open_center,fc='white')
        #plt.title('', y=0.95 + (explo / 2),size=20,fontweight='bold')
        plt.gca().add_artist(centre_circle)

# Network
from networkx import path_graph, random_layout
import matplotlib.pyplot as plt
import networkx as nx

def net_plot(df = DataFrame([]), layout = 'Spring', label = 'none', column = 0, label_size = 5,diam_nodos = 10, espe_edges = 0.1,
             inter = 10, color_inter_min = 'k',color_inter_max = 'blue',
             edge_alpha_min = 0.3, edge_alpha_max = 0.3, k_num = 3, color_nodo = 'red', node_alpha = 0.7, backg = 'white',
             label_color = 'black'):
    #
    df1 = df[['GO', 'Entry', 'Term', 'Short_Term']].merge(df, on = 'Entry', how = 'left').drop_duplicates()
    df2 = DataFrame(df[['GO', 'Term', 'Short_Term', 'Entry']].drop_duplicates().groupby(['GO','Term', 'Short_Term']).Entry.count()).reset_index()
    #### >>>>>>>>>>>>>>>>
    #### A partir de una matriz de datos extrae valores no redundantes
    matrix = df1.pivot_table(values='Entry',index=['GO_x', 'Term_x', 'Short_Term_x'],aggfunc=len,columns=['GO_y', 'Term_y', 'Short_Term_y'])
    ###
    df_mat = []
    n = -1
    for i in list(matrix.columns.values):
        n += 1
        new = DataFrame(matrix.iloc[n:len(matrix)][i])
        nn = -1
        for index, row in new.iterrows():
            nn += 1
            df_mat.append([index, i, new.iloc[nn][i]])
        nn = 0
    ###
    df_mat = DataFrame(df_mat, columns = ['go0', 'go1', 'val']).dropna()
    ###
    nodos = []
    for index, row in df_mat.iterrows():
        if row.go0 == row.go1:
            #print(row.go0, row.go1)
            continue
        else:
            #print(row.go0, row.go1)
            nodos.append([row.go0, row.go1, row.val])
    nodos = DataFrame(nodos)
    columnas = {0:'GO', 1:'Term', 2:'Short_Term'}
    nodos = DataFrame([[i[column] for i in nodos[0]], [i[column] for i in nodos[1]], nodos[2]]).T
    #### >>>>>>>>>>>>>>>>
    # si interacciona con mas uno, eliminar la redundancia, y si no interacciona con ninguno, dejar el nodo
    # y su valor, este se verá en la red como un nodo aislado
    aislado = [i for i in matrix.columns if len(matrix[[i]].dropna()) == 1]
    aislado = [df_mat[df_mat.go0 == i] for i in aislado]
    if len(aislado) > 0:
        aislado = pd.concat(aislado)
        aislado.columns = [0, 1, 2]
        aislado = DataFrame([[i[column] for i in aislado[0]], [i[column] for i in aislado[1]], aislado[2]]).T
        nodos = pd.concat([nodos, aislado])
    else:
        pass
    ####################
    # https://networkx.github.io/documentation/networkx-2.3/auto_examples/drawing/plot_weighted_graph.html#sphx-glr-auto-examples-drawing-plot-weighted-graph-py
    G=nx.Graph()
    for index, row in nodos.iterrows():
        G.add_edge(row[0], row[1],weight = row[2])
    elarge=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] >= inter]
    esmall=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] < inter]
    ###
    layouts = {'Circular':nx.circular_layout,
          'Random':nx.random_layout,
          'Shell':nx.shell_layout,
          'Spectral':nx.spectral_layout,
          'Spring':nx.spring_layout,
          'KK':nx.kamada_kawai_layout}
    #circular_layout
    #random_layout
    #shell_layout
    #spring_layout
    #spectral_layout
    #pos=nx.spring_layout(G, k = k_num) # positions for all nodes
    if layouts[layout] == nx.spring_layout:
        pos=layouts[layout](G, k = k_num)
    else:
        pos=layouts[layout](G)
    #pos=layout(G)
    # nodes
    #------------------------------------------------------------------------------
    # ordenar los valores para representarlos en el tama;o del nodo
    order = []
    for index, row in nodos.iterrows():
        order.append(row[0])
        order.append(row[1])
    orden3 = DataFrame(order).drop_duplicates(keep = 'first').reset_index(drop = True)
    orden3.columns = [columnas[column]]
    orden4 = pd.merge(orden3, df2, on = [columnas[column]], how = 'left')
    #------------------------------------------------------------------------------
    # https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.drawing.nx_pylab.draw_networkx_edges.html
    nx.draw_networkx_nodes(G,pos,node_size= np.array(orden4.Entry) * diam_nodos,
                           node_color= color_nodo,alpha= node_alpha)
    # edges
    nx.draw_networkx_edges(G,pos,edgelist=esmall,
                           width = np.array([i[2] for i in esmall]) * espe_edges,
                           alpha= edge_alpha_min,edge_color= color_inter_min,style='-')
    nx.draw_networkx_edges(G,pos, edgelist=elarge,
                           width = np.array([i[2] for i in elarge]) * espe_edges,
                           alpha= edge_alpha_max,edge_color= color_inter_max,style= '-')

    # labels
    posicion = {} ## posicion de las etiquetas, ligeramente arriba
    for key, value in pos.items():
        posicion[key] = value + 0.05
    # arreglo de las posiciones de los nodos en el plano cartesiano
    arr = np.array([[i for i in value] for key, value in pos.items()])
    
    # labels
    if label == 'label':
        nx.draw_networkx_labels(G,posicion,font_size=label_size, font_color=label_color) # ,font_weight='bold'
        if label == 'label':
            plt.axis([arr[:,0].min() - 0.3, arr[:,0].max() + 0.3,
                      arr[:,1].min() - 0.3, arr[:,1].max() + 0.3])
        #plt.axis('off')
        #plt.show() # display
    if label == 'none':
        plt.axis([arr[:,0].min() - 0.2, arr[:,0].max() + 0.2,
                  arr[:,1].min() - 0.2, arr[:,1].max() + 0.2])
        #plt.axis('off')
        #plt.show() # display
    plt.gca().set_facecolor(backg)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)   
    #return print(len(nodos),'Connections')
##################
import seaborn as sns
#https://xkcd.com/color/rgb/
def heatmap_plot(df = DataFrame([]), colors = 'hsv', label_x = 'GO', label_y = 'Term', xticks_size = 12, yticks_size = 12, ylabel_size = 18):
    if len(df) == 0:
        print('Data frame sin valores')
    else:
        df1 = df.groupby(['Aspect']).get_group(df.Aspect.iloc[0])[['GO', 'Term', 'Entry']].drop_duplicates().merge(df, on = 'Entry', how = 'left').drop_duplicates()
        new_df1 = []
        for index, row in df1.iterrows():
            if row.GO_x == row.GO_y:
                continue
            if row.Term_x == row.Term_y:
                continue
            else:
                new_df1.append([row.GO_x, row.Term_x, row.Entry, row.GO_y, row.Term_y, row.Aspect])
        df2 = DataFrame(new_df1, columns =['GO_x', 'Term_x', 'Entry', 'GO_y', 'Term_y', 'Aspect'])
        matrix = df2.pivot_table(values='Entry',index=['GO_x', 'Term_x'],aggfunc=len,columns=['GO_y', 'Term_y'])
        ###
        #plt.subplots(figsize=(size_plot,size_plot))
        plt.gca().set_facecolor('whitesmoke')
        plt.xticks(size = xticks_size) # fontweight='bold'
        plt.yticks(size = yticks_size) # fontweight='bold'

        if label_x == 'GO':
            xlab = [i[0] for i in matrix.index.values]
        if label_x == 'Term':
            xlab = [i[1] for i in matrix.index.values]
        if label_y == 'GO':
            ylab = [i[0] for i in matrix.columns.values]
        if label_y == 'Term':
            ylab = [i[1] for i in matrix.columns.values]
        
        sns.heatmap(np.array(matrix),  square=True, annot = False, cmap = colors,
                mask = np.zeros_like(np.array(matrix)),
                #annot_kws={"size": 13, 'fontweight': 'bold'},fmt= '.0f',
                linewidths=0.005, linecolor='w',
                cbar=True, cbar_kws={"shrink": 0.5}, # 'label': '
                xticklabels= xlab,
                yticklabels= ylab)
        plt.gca().figure.axes[-1].set_ylabel('Interaction degrees\n(Proteins)', size= ylabel_size)
        plt.gca().set_facecolor('xkcd:light grey')
##################
# ##################
# ##################
# ##################
import networkx as nx
import nxviz as nxv

def circos(df = DataFrame([]), label = True, node_color = 'None', column = 1, inter = 25, size = 5,
           fontsize = 10, type_edges = 'Normal', increase = 2):
    df1 = df[['GO', 'Entry', 'Term', 'Short_Term']].merge(df, on = 'Entry', how = 'left').drop_duplicates()
    df2 = DataFrame(df[['GO', 'Term', 'Short_Term', 'Entry']].drop_duplicates().groupby(['GO','Term', 'Short_Term']).Entry.count()).reset_index()
    #### >>>>>>>>>>>>>>>>
    #### A partir de una matriz de datos extrae valores no redundantes
    matrix = df1.pivot_table(values='Entry',index=['GO_x', 'Term_x', 'Short_Term_x'],aggfunc=len,columns=['GO_y', 'Term_y', 'Short_Term_y'])
    ###
    df_mat = []
    n = -1
    for i in list(matrix.columns.values):
        n += 1
        new = DataFrame(matrix.iloc[n:len(matrix)][i])
        nn = -1
        for index, row in new.iterrows():
            nn += 1
            df_mat.append([index, i, new.iloc[nn][i]])
        nn = 0
    ###
    df_mat = DataFrame(df_mat, columns = ['go0', 'go1', 'val']).dropna()
    ###
    nodos = []
    for index, row in df_mat.iterrows():
        if row.go0 == row.go1:
            #print(row.go0, row.go1)
            continue
        else:
            #print(row.go0, row.go1)
            nodos.append([row.go0, row.go1, row.val])
    nodos = DataFrame(nodos)
    columnas = {0:'GO', 1:'Term', 2:'Short_Term'}
    nodos = DataFrame([[i[column] for i in nodos[0]], [i[column] for i in nodos[1]], nodos[2]]).T
    #### >>>>>>>>>>>>>>>>
    # si interacciona con mas uno, eliminar la redundancia, y si no interacciona con ninguno, dejar el nodo
    # y su valor, este se verá en la red como un nodo aislado
    aislado = [i for i in matrix.columns if len(matrix[[i]].dropna()) == 1]
    aislado = [df_mat[df_mat.go0 == i] for i in aislado]
    if len(aislado) > 0:
        aislado = pd.concat(aislado)
        aislado.columns = [0, 1, 2]
        aislado = DataFrame([[i[column] for i in aislado[0]], [i[column] for i in aislado[1]], aislado[2]]).T
        nodos = pd.concat([nodos, aislado])
        aislado.to_csv('datos/nodes_without_interactions.tsv', sep = '\t', index = None)
    else:
        pass
    nodos.columns = ['Source','Target','Weight']
    edges = nodos
    order = []
    for index, row in nodos.iterrows():
        order.append(row[0])
        order.append(row[1])
    orden3 = DataFrame(order).drop_duplicates(keep = 'first').reset_index(drop = True)
    orden3.columns = [columnas[column]]
    nodes = pd.merge(orden3, df2, on = [columnas[column]], how = 'left')
    
    def make_graph(nodes, edges):
        g = nx.Graph()

        for i,row in nodes.iterrows():
            keys = row.index.tolist()
            values = row.values
            # The dict contains all attributes
            g.add_node(row[nodes.columns[0]], **dict(zip(keys,values)))

        for i,row in edges.iterrows():
            keys = row.index.tolist()
            values = row.values
            g.add_edge(row['Source'], row['Target'], **dict(zip(keys,values)))
        return g
    
    g = make_graph(nodes, edges)
    for i,row in nodes.iterrows():
        if row['Entry'] >= inter:
            g.add_node(row[nodes.columns[0]], umbral='up')
        if row['Entry'] < inter:
            g.add_node(row[nodes.columns[0]], umbral='down')
    color_nodo = {'Uniques':nodes.columns[0],
              'Umbral':'umbral',
              'None':False}
    edges_values_modificados = {'Normal':'Weight',
                                'log10':list(np.log10(edges['Weight'].tolist()) * increase)}
    c = nxv.CircosPlot(g,
                   node_color= color_nodo[node_color], # nodes.columns[0],
                   node_grouping= color_nodo[node_color],
                   node_labels=label,
                   node_label_layout='rotation',
                   edge_width= edges_values_modificados[type_edges],
                   #edge_color = 'umbral',
                   figsize=(size,size),
                   fontsize = fontsize)
    return c.draw()
    ####################
####
def venn2_plot(set1 = set(),
               set2 = set(),
               lab_set1 = 'Set1',
               lab_set2 = 'Set2',
               linewidth = 1,
               color_line = 'black',
               alpha_sets = 0.3,
               font_sets = False, # False o 'bold'
               size_vals_sets = 12,
               alpha_inter = 0.3,
               font_inter = False, # False o 'bold'
               size_vals_inter = 12,
               size_label = 12,
               font_label = False): # False o 'bold'
    v = venn2_unweighted(subsets = (set1, set2), set_labels = (lab_set1, lab_set2))
    c = venn2_circles(subsets = (1, 1, 1),
                    linestyle='--', linewidth= linewidth, color=color_line)
    v.get_patch_by_id('10').set_alpha(0)
    partes = ['10', '01', '11']
    partes2 = ['10', '01']
    venn_info = [[i, j] for i, j in zip(v.subset_labels, partes)]
    for i in venn_info:
        if i[0] != None:
            if i[1] in partes2:
                v.get_patch_by_id(i[1]).set_alpha(alpha_sets) # i[1] = el conjunto creado,  0 = alpha del conjunto
                v.get_label_by_id(i[1]).set_fontweight(font_sets)
                v.get_label_by_id(i[1]).set_fontsize(size_vals_sets)
            if i[1] == '11': # configurar la intersección independientemente '111'
                v.get_patch_by_id('11').set_alpha(alpha_inter) # i[1] = el conjunto creado,  0 = alpha del conjunto
                v.get_label_by_id('11').set_fontweight(font_inter)
                v.get_label_by_id('11').set_fontsize(size_vals_inter)    
    for text in v.set_labels:
        text.set_fontsize(size_label)
        text.set_fontweight(font_label)
def venn3_plot(set1 = set(),
               set2 = set(),
               set3 = set(),
               lab_set1 = 'Set1',
               lab_set2 = 'Set2',
               lab_set3 = 'Set3',
               linewidth = 1,
               color_line = 'black',
               alpha_sets = 0.3,
               font_sets = False, # False o 'bold'
               size_vals_sets = 12,
               alpha_inter = 0.3,
               font_inter = False, # False o 'bold'
               size_vals_inter = 12,
               size_label = 12,
               font_label = False): # False o 'bold'
    v = venn3_unweighted(subsets = (set1, set2, set3), set_labels = (lab_set1, lab_set2, lab_set3))
    c = venn3_circles(subsets = (1, 1, 1, 1, 1, 1, 1),
                      linestyle='--', linewidth = linewidth, color = color_line)
    partes = ['100', '010', '110', '001', '101', '011', '111']
    partes2 = ['100', '010', '110', '001', '101', '011']
    venn_info = [[i, j] for i, j in zip(v.subset_labels, partes)]
    for i in venn_info:
        if i[0] != None:
            if i[1] in partes2:
                v.get_patch_by_id(i[1]).set_alpha(alpha_sets) # i[1] = el conjunto creado,  0 = alpha del conjunto
                v.get_label_by_id(i[1]).set_fontweight(font_sets)
                v.get_label_by_id(i[1]).set_fontsize(size_vals_sets)
            if i[1] == '111': # configurar la intersección independientemente '111'
                v.get_patch_by_id('111').set_alpha(alpha_inter) # i[1] = el conjunto creado,  0 = alpha del conjunto
                v.get_label_by_id('111').set_fontweight(font_inter)
                v.get_label_by_id('111').set_fontsize(size_vals_inter)    
    for text in v.set_labels:
        text.set_fontsize(size_label)
        text.set_fontweight(font_label)
