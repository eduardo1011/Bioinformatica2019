"""
=========================

Funciones, curso bioinformática

1. barras(df = DataFrame([]), column = 1, dim = 111, title = 'Title',
            row_num = 10, color = '#ff7f0e', size_x = 15, size_y = 15,
            xlabel = 20, size_title = 25, size_bartxt = 12)

2. pastel(df = DataFrame([]), column = 0, dim = 111, title = 'Title',
            row_num = 3, angle = 0, size_text = 15, size_title = 20)

3. word_cloud(df = DataFrame([]), size = 10)

4. GO_DAG(go_terms = [], term = '')

5. venn3_plot(set1 = set(), set2 = set(), set3 = set(), lab_set1 = 'Set1',
            lab_set2 = 'Set2', lab_set3 = 'Set3', size_label = 20,
            size_vals = 20)

6. net_plot(df = DataFrame([]),label = 'none',diam_nodos = 10, espe_edges = 0.1,
            inter = 10, color_inter_min = 'k',color_inter_max = 'blue',
            edge_alpha_min = 0.3, edge_alpha_max = 0.3, k_num = 3,
            color_nodo = 'red', node_alpha = 0.7)

7. get_UniProtKB_info1(ids = [])

8. get_UniProtKB_info2(id_organism = 0)

9. get_UniProtKB_info3(ids = [])

10. get_UniProtKB_info4(ids = [])

11. heatmap_plot(df = DataFrame([]), colors = 'Spectral', label_x = 'GO',
            label_y = 'Term', size_plot = 10, xticks_size = 12,
            yticks_size = 12, ylabel_size = 18)

12. go_file()

13. run_Blast_Windows(x_p = '', file = '', db = '')

14. run_Blast_Linux(x_p = '', file = '', db = '')

15. run_blast_jupyter(x_p = '', file = '', db = '')

=========================

"""
##
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
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
import json
from PIL import Image
import webbrowser
import subprocess
from datetime import datetime
# nuevo
import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, interact_manual,\
Button, HBox, VBox, IntSlider, Label, IntRangeSlider
from ipywidgets import Checkbox, RadioButtons
import random
from matplotlib import cm
###
files = [i for i in os.listdir('datos/') if re.search('.tsv', i)]
fs = widgets.Dropdown(options =  files)
def box(fs):
    print('\n        ',fs)
    df = pd.read_csv('datos/'+fs, sep = '\t', header = None).dropna()
    return display(df.head()), print(len(df),'rows')
output = widgets.interactive_output(box, {'fs':fs})
file = widgets.HBox([VBox([fs]), output])
def frame():
    df = pd.read_csv('datos/'+fs.value, sep = '\t', header = None).dropna()
    return df
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
              'Spectral':cm.Spectral(np.arange(11)/11.)}
colors = {'#8DD3C7':pie_colors['Set3'][0:1], '#FFFFB3':pie_colors['Set3'][1:2],
         '#BEBADA':pie_colors['Set3'][2:3], '#FB8072':pie_colors['Set3'][3:4],
         '#80B1D3':pie_colors['Set3'][4:5], '#FDB462':pie_colors['Set3'][5:6],
         '#B3DE69':pie_colors['Set3'][6:7], '#FCCDE5':pie_colors['Set3'][7:8],
         '#D9D9D9':pie_colors['Set3'][8:9], '#BC80BD':pie_colors['Set3'][9:10],
         '#CCEBC5':pie_colors['Set3'][10:11], '#FFED6F':pie_colors['Set3'][11:12]}
sizes = ['3x3','4x4','5x5','6x6','7x7','8x8','9x9','10x10']
escala_uno = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
common_colors = ['red','orange','yellow','green','blue','purple','brown','magenta','tan',
                 'cyan','olive','maroon','navy','aquamarine','turquoise','silver',
                 'lime','teal','indigo','violet','pink','black','white']
#####

code = {200:'The request was processed successfully.',
        400:'Bad request. There is a problem with your input.',
        404:'Not found. The resource you requested doesn’t exist.',
        410:'Gone. The resource you requested was removed.',
        500:'Internal server error. Most likely a temporary problem, but if the problem persists please contact us.',
        503:'Service not available. The server is being updated, try again later.'}
###
asp = {'P':'Biological Process',
       'F':'Molecular Function',
       'C':'Cellular Component'}
###
if os.path.exists('datos/go-basic.obo'):
    print('\nYa está descargado el archivo go-basic.obo\n')
    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    #print(urllib.request.urlopen(url).headers)
else:
    # Método 1: urllib.request.urlretrieve('http://purl.obolibrary.org/obo/go.obo', 'datos/go.obo')
    # Método 2:
    print('\nDescargando la Ontología: go-basic.obo\n')
    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    go_file = 'datos/go-basic.obo'
    with open(go_file, 'wb') as f:
        #print ("Downloading %s" % file_name)
        response = requests.get(url, stream=True)
        print(code[response.status_code],'\n')
        total_length = response.headers.get('content-length')
        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=8192):
                dl += len(data)
                f.write(data)
                done = int(40 * dl / total_length)
                sys.stdout.write("\rDescargando archivo go-basic.obo [%s%s] %s MB" % ('■' * done, ' ' * (40-done), round(dl/1000000,2)), )    
                sys.stdout.flush()
    # información de la base de datos
    #print('')
    #print('')
    #print(urllib.request.urlopen(url).headers)
###
def go_file():
    if os.path.exists('datos/go-basic.obo'):
        url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        #print(urllib.request.urlopen(url).headers)
        with open('datos/go-basic.obo', 'r') as g:
            go_obo = g.read()
        go1 = go_obo.split('[Term]')
    else:
        # Método 1: urllib.request.urlretrieve('http://purl.obolibrary.org/obo/go.obo', 'datos/go.obo')
        # Método 2:
        url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        go_file = 'datos/go-basic.obo'
        with open(go_file, 'wb') as f:
            #print ("Downloading %s" % file_name)
            response = requests.get(url, stream=True)
            print(code[response.status_code])
            total_length = response.headers.get('content-length')
            if total_length is None: # no content length header
                f.write(response.content)
            else:
                dl = 0
                total_length = int(total_length)
                for data in response.iter_content(chunk_size=8192):
                    dl += len(data)
                    f.write(data)
                    done = int(40 * dl / total_length)
                    sys.stdout.write("\rDescargando archivo go-basic.obo [%s%s] %s MB" % ('■' * done, ' ' * (40-done), round(dl/1000000,2)), )    
                    sys.stdout.flush()
        with open('datos/go-basic.obo', 'r') as g:
            go_obo = g.read()
            go1 = go_obo.split('[Term]')
        # información de la base de datos
        #print(urllib.request.urlopen(url).headers)
    return go1
###
ontology_file = go_file()
aspect = {'biological_process':'P', 'molecular_function':'F', 'cellular_component':'C'}
items = []
for i in ontology_file[1:len(ontology_file)]:
    items.append([i.split('\n')[1].split(': ')[1],
                 i.split('\n')[2].split(': ')[1],
                 aspect[i.split('\n')[3].split(': ')[1]]])
ontologia = DataFrame(items, columns = ['GO', 'Term', 'Aspect'])
def ontology():
    return ontologia
###
from ipywidgets import Button, Layout
goslim = [re.search('go\w+', i).group() for i in ontology_file[0].split('\n') if re.search(' go\w+', i)]
x = widgets.Dropdown(options=goslim,description='GO Slim:',disabled=False,
                    layout=Layout(width='40%', height='25px'))
def go_slim(x):
    slim = []
    for i in ontology_file[1:len(ontology_file)]:
        if re.findall('subset: '+x, i):
            slim.append(i.split('\n')[1].split(': ')[1])
        else:
            continue
    slim = DataFrame(slim, columns = ['GO'])
    slim = slim.merge(ontologia, on = 'GO', how = 'left').dropna()
    return display(slim.head()), print(len(slim))
out = widgets.interactive_output(go_slim, {'x':x})
GOSLIM = widgets.HBox([x, out])
def frame_goslim():
    slim = []
    for i in ontology_file[1:len(ontology_file)]:
        if re.findall('subset: '+x.value, i):
            slim.append(i.split('\n')[1].split(': ')[1])
        else:
            continue
    slim = DataFrame(slim, columns = ['GO'])
    slim = slim.merge(ontologia, on = 'GO', how = 'left').dropna()
    return slim



###
# un dict
iiss = {}
for i in ontology_file[1:len(ontology_file)]:
    cero = i.split('\n')[1].split(': ')[1]
    uno = re.findall('\nis_a: GO:.*!',i)
    is_a = []
    for j in uno:
        if re.search('GO:\d+', j):
            is_a.append(re.search('GO:\d+', j).group())
        else:
            continue
    iiss[cero] = is_a
    #if len(is_a) >= 1: # eliminamos términos sin "is_a"
    #    iiss[cero] = is_a
    #else:
    #    continue
def IS_A():
    is_a = []
    for i in iiss:
        for j in list(iiss[i]):
            is_a.append([i,j])
    is_a = DataFrame(is_a, columns = ['go', 'GO'])
    is_a = pd.merge(is_a, ontologia, on = 'GO', how = 'left')
    is_a = is_a[is_a['Term'].str.contains('biological_process|molecular_function|cellular_component') == False]
    is_a.columns = ['GO', 'is_a', 'Term','Aspect'] # el Term y Aspect pertenecen a "is_a"
    return is_a
###
# boocle para buscar super términos, esto sustituye la instalación de Orange
def super_terms(terms):
    terms = [terms] if isinstance(terms, str) else terms
    uno = set()
    dos = set(terms)
    while dos:
        term = dos.pop()
        uno.add(term)
        terminos = set(iiss[term]) # el haber hecho un "dict" potenció drásticamente el tiempo en buscar supertérminos 
        dos.update(terminos - uno)
    return uno
def mapping_super_terms(list_gos = []):
    sup1 = {}
    for i in list_gos:
        sup1[i] = list(super_terms(i))
    df = []
    for i in sup1:
        for j in list(sup1[i]):
            df.append([i,j])
    return DataFrame(df, columns = ['GO', 'is_a'])
#--------------------------------------------------------------------------------------
###
###
def barras(df = DataFrame([]), column = 1, dim = 111, title = 'Title', row_num = 10, color = '#ff7f0e',
           size_x = 15, size_y = 15, xlabel = 20, size_title = 25, size_bartxt = 12, sep = 3):
    if len(df) == 0:
        print('Data frame sin datos')
    else:
        plt.subplot(dim, facecolor= 'white')
        barWidth = 0.9
        val = max(list(df.iloc[0:row_num,3]))
        plt.barh(list(df.iloc[0:row_num,column]) ,list(df.iloc[0:row_num,3]),
                 color= color,
                 align='center',
                 height= 0.7,
                 edgecolor='white')
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_position(('data',-1))
        plt.gca().spines['left'].set_visible(False)
        plt.title(title,size= size_title,loc='left', fontweight='bold')
        plt.tick_params(axis="y", color="gray")
        plt.yticks(size = size_y )
        #plt.xticks(size = size_x ,fontweight='bold')
        if max(list(df.iloc[0:row_num,3])) > 200:
            escala = 50
        if max(list(df.iloc[0:row_num,3])) < 200:
            escala = 20
        plt.xticks(range(0,val,escala),size=size_x) #fontweight='bold'
        plt.xlabel("Proteins",size= xlabel)

        for j, k in zip(list(df.iloc[0:row_num,3]),range(0,len(list(df.iloc[0:row_num,1])))):
            plt.text(j+sep,k-0.2,j,
                     size=size_bartxt,ha='left',color= 'black')
###


def pastel(df = DataFrame([]), column = 0, dim = 111, title = 'Title', row_num = 5, angle = 0,
           size_text = 15, size_title = 20, color = pie_colors['Set1']):
    if len(df) < 1:
        print('Data frame sin datos')
    else:
        plt.subplot(dim)
        labels = list(df.iloc[0:row_num,column])
        sizes = list(df.iloc[0:row_num,3])
        explode = np.repeat(0.03,row_num)
        plt.pie(sizes,labels=labels, autopct='%1.1f%%', startangle=angle, radius = 1,
                colors = color,
                explode = explode,textprops=dict(size = size_text))

        centre_circle = plt.Circle((0,0),0,fc='white')
        #centre_circle = plt.Circle((0,0),0.50,fc='white')
        plt.title(title,size=size_title,fontweight='bold')
        plt.gca().add_artist(centre_circle)

###
def GO_DAG(go_terms = [], term = ''):
    if len(go_terms) == 0:
        print('Lista vacía')
    else:
        if term == 'none':
            formato1 = []
            for i in go_terms:
                formato = {i:{"title":i, "body":DataFrame(data = {'GO': [i]}).merge(ontologia, on = 'GO', how = 'left').Term[0], "fill":"gold", "font":"black", "border":"black"}}
                formato1.append(['"'+i+'"', json.dumps(formato[i])])
            dag = re.sub('}:','},', '{'+':'.join([':'.join(i) for i in formato1])+'}')
        
        else:
            formato1 = []
            for i in go_terms:
                if i == term:
                    formato = {i:{"title":i, "body":DataFrame(data = {'GO': [i]}).merge(ontologia, on = 'GO', how = 'left').Term[0], "fill":"cyan", "font":"black", "border":"black"}}
                    formato1.append(['"'+i+'"', json.dumps(formato[i])])
                else:
                    formato = {i:{"title":i, "body":DataFrame(data = {'GO': [i]}).merge(ontologia, on = 'GO', how = 'left').Term[0], "fill":"gold", "font":"black", "border":"black"}}
                    formato1.append(['"'+i+'"', json.dumps(formato[i])])
            dag = re.sub('}:','},', '{'+':'.join([':'.join(i) for i in formato1])+'}')
        ## codifica una url a partir de formato utf8
        #import urllib.parse
        query = dag
        convert_utf8_to_url = urllib.parse.quote(query)
        ##
        formato_png = 'http://amigo.geneontology.org/visualize?term_data='+convert_utf8_to_url+'&term_data_type=json&mode=amigo&inline=false&format=png'
        formato_svg = 'http://amigo.geneontology.org/visualize?inline=false&mode=amigo&term_data='+convert_utf8_to_url+'&format=svg&term_data_type=json'
        ## decodifica la url
        #from urllib.parse import unquote
        #url = 'http://amigo.geneontology.org/visualize?term_data='+convert_utf8_to_url+'&term_data_type=json&mode=amigo&inline=false&format=png'
        #unquote(url)
        #####
        #import webbrowser
        #chrome_path = 'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'
        #webbrowser.get(chrome_path).open(formato_png)
        return [dag, formato_png, formato_svg]
###
def venn2_plot(set1, set2, label1, label2,size_label,size_num):
    v = venn2([set1, set2],
              ([label1, label2]))
    c=venn2_circles(subsets = [set1, set2],
                    linestyle='dashed', linewidth=2, color="k")
    v.get_patch_by_id('11').set_color('w')
    v.get_patch_by_id('11').set_alpha(1)
    v.get_label_by_id('11').set_color('k')
    v.get_label_by_id('11').set_fontweight('bold')
    for text in v.set_labels:
        text.set_fontsize(size_label)
        text.set_fontweight('bold')
    for text in v.subset_labels:
        text.set_fontsize(size_num)
    #plt.show()
###
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
################
###
# Network
from networkx import path_graph, random_layout
import matplotlib.pyplot as plt
import networkx as nx

def net_plot(df = DataFrame([]), layout = 'Spring', label = 'none', column = 0, label_size = 5,diam_nodos = 10, espe_edges = 0.1,
             inter = 10, color_inter_min = 'k',color_inter_max = 'blue',
             edge_alpha_min = 0.3, edge_alpha_max = 0.3, k_num = 3, color_nodo = 'red', node_alpha = 0.7):
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
          'Spring':nx.spring_layout}
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
        nx.draw_networkx_labels(G,posicion,font_size=label_size) # ,font_weight='bold'
        if label == 'label':
            plt.axis([arr[:,0].min() - 0.2, arr[:,0].max() + 0.2,
                      arr[:,1].min() - 0.2, arr[:,1].max() + 0.2])
        #plt.axis('off')
        #plt.show() # display
    if label == 'none':
        plt.axis([arr[:,0].min() - 0.1, arr[:,0].max() + 0.1,
                  arr[:,1].min() - 0.1, arr[:,1].max() + 0.1])
        #plt.axis('off')
        #plt.show() # display
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    #return print(len(nodos),'Connections')
###
def get_UniProtKB_info0(ids = []):
    sets_ids = []
    for k in range(0, len(ids),250):
        sets_ids.append(ids[k:k+250])
    print('Grupos:', len(sets_ids))
    if len(sets_ids) == 0:
        print('Lista vacía')
    else:
        uniprot = []
        n = 0
        for i in sets_ids:
            n += 1
            params = {'from':'ACC','to':'ACC','format':'tab','query': ' '.join(i)}
            data = urllib.parse.urlencode(params)
            url = 'https://www.uniprot.org/uploadlists/'
            response = requests.get(url, data)
            print('Set '+str(n)+' ('+str(len(i))+' IDs):', code[response.status_code])
            respuesta = response.content.decode()
            names = ['ACC', 'Entry']
            df = pd.read_csv(StringIO(respuesta),sep='\t',header=None).drop(index = [0])
            df.columns = names
            uniprot.append(df)
        uniprotkb0 = pd.concat(uniprot).fillna('NA')
        return uniprotkb0
###
def get_UniProtKB_info1(ids = []):
    code = {200:'The request was processed successfully.',
        400:'Bad request. There is a problem with your input.',
        404:'Not found. The resource you requested doesn’t exist.',
        410:'Gone. The resource you requested was removed.',
        500:'Internal server error. Most likely a temporary problem, but if the problem persists please contact us.',
        503:'Service not available. The server is being updated, try again later.'}
    sets_ids = []
    for k in range(0, len(ids),200):
        sets_ids.append(ids[k:k+200])
    print('Grupos:', len(sets_ids))
    if len(sets_ids) == 0:
        print('Lista vacía')
    else:
        uniprot = []
        n = 0
        for i in sets_ids:
            n += 1
            url = 'https://www.uniprot.org/uniprot/'
            uu = '%20OR%20id:'.join(i)
            uu = re.sub('^', '?query=id:',uu)
            columnas = '&format=tab&columns=id,entry%20name,protein%20names,organism,genes,length,go-id,go,organism-id'
            new_url = url + uu + columnas
            response = requests.get(new_url)
            print('Set '+str(n)+' ('+str(len(i))+' IDs):', code[response.status_code])
            respuesta = response.content.decode()
            names = ['Entry', 'Entry_name', 'Protein_name', 'Organism', 'Gene', 'Length', 'GO', 'Terms', 'Tax_ID']
            df = pd.read_csv(StringIO(respuesta),sep='\t')
            df.columns = names
            uniprot.append(df)
        uniprotkb = pd.concat(uniprot).fillna('NA')
        return uniprotkb
###
def get_UniProtKB_info2(id_organism = 0):
    if id_organism == 0:
        print('ID de organismo incorrecto')
    else:
        link = 'https://uniprot.org/uniprot/?query=organism:'+str(id_organism)+'&format=tab&columns=id,entry name,protein names,organism,genes,length,go-id,go,organism-id'
        file_name = 'annotation_'+str(id_organism)
        with open(file_name, 'wb') as f:
            #print ("Downloading %s" % file_name)
            response = requests.get(link, stream=True)
            total_length = response.headers.get('X-Total-Results')
            if total_length is None: # no content length header
                f.write(response.content)
            else:
                dl = 0
                total_length = int(total_length)
                for data in response.iter_content(chunk_size=8192):
                    dl += len(data)
                    f.write(data)
                    done = int(0.1 * dl / total_length)
                    sys.stdout.write("\rLoading [%s%s] %s MB" % ('■' * done, ' ' * (5-done), round(dl/1000000,2)), ) 
                    sys.stdout.flush()
        ###
        acc_uniprot_GO_id=pd.read_csv('annotation_'+str(id_organism),sep='\t')
        names = ['Entry', 'Entry_name', 'Protein_name', 'Organism', 'Gene', 'Length', 'GO', 'Terms', 'Tax_ID']
        acc_uniprot_GO_id.columns = names
        return acc_uniprot_GO_id
###
def get_UniProtKB_info3(ids = []):
    sets_ids = []
    for k in range(0, len(ids),250):
        sets_ids.append(ids[k:k+250])
    print('Grupos:', len(sets_ids))
    if len(sets_ids) == 0:
        print('Lista vacía')
    else:
        uniprot = []
        n = 0
        for i in sets_ids:
            n += 1
            params = {'from':'ACC','to':'ACC','format':'tab','query': ' '.join(i),
                     'columns':'id,entry name,protein names,organism,genes,length,go-id,go,organism-id'}
            data = urllib.parse.urlencode(params)
            url = 'https://www.uniprot.org/uploadlists/'
            response = requests.get(url, data)
            print('Set '+str(n)+' ('+str(len(i))+' IDs):', code[response.status_code])
            respuesta = response.content.decode()
            names = ['Entry', 'Entry_name', 'Protein_name', 'Organism', 'Gene', 'Length', 'GO', 'Terms', 'Tax_ID']
            df = pd.read_csv(StringIO(respuesta),sep='\t',header=None).drop(columns = [9]).drop(index = [0])
            df.columns = names
            uniprot.append(df)
        uniprotkb = pd.concat(uniprot).fillna('NA')
        return uniprotkb
###
import urllib.parse
import urllib.request
def get_UniProtKB_info4(ids = []):
    sets_ids = []
    for k in range(0, len(ids),250):
        sets_ids.append(ids[k:k+250])
    print('Grupos:', len(sets_ids))
    if len(sets_ids) == 0:
        print('Lista vacía')
    else:
        uniprot = []
        n = 0
        for i in sets_ids:
            n += 1
            params = {'from':'ACC','to':'ACC','format':'tab','query': ' '.join(i),
                     'columns':'id,entry name,protein names,organism,genes,length,go-id,go,organism-id'}
            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            url = 'https://www.uniprot.org/uploadlists/'
            req = urllib.request.Request(url, data)
            with urllib.request.urlopen(req) as f:
                response = f.read().decode('utf-8')
            print('Set '+str(n)+' ('+str(len(i))+' IDs):', code[f.getcode()])
            names = ['Entry', 'Entry_name', 'Protein_name', 'Organism', 'Gene', 'Length', 'GO', 'Terms', 'Tax_ID']
            df = pd.read_csv(StringIO(response),sep='\t',header=None).drop(columns = [9]).drop(index = [0])
            df.columns = names
            uniprot.append(df)
        uniprotkb = pd.concat(uniprot).fillna('NA')
        return uniprotkb
###
import seaborn as sns
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
        #plt.xticks(rotation=70)

        #plt.show()
###

###
def run_Blast_Windows(x_p = '', file = '', db = '', evalue = ''):
    if file == '':
        print('Intentar nuevamente')
    else:
        urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/runBlast.py', 'runBlast.py')
        blast = open('runBlast.py','r')
        blast = blast.read()
        blast = re.sub('xxxxx',"'"+file+"'", blast)
        blast = re.sub('yyyyy',"'"+x_p+"'", blast)
        blast = re.sub('zzzzz',"'"+db+"'", blast)
        blast = re.sub('uuuuu',"'"+str(evalue)+"'", blast)
        save= open('run.py','w')
        save.write(blast)
        save.close()
        runblast = os.system("start cmd /k python run.py")
###
def run_Blast_Linux(x_p = '', file = '', db = '', ):
    if file == '':
        print('Intentar nuevamente')
    else:
        urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/runBlast.py', 'runBlast.py')
        blast = open('runBlast.py','r')
        blast = blast.read()
        blast = re.sub('xxxxx',"'"+file+"'", blast)
        blast = re.sub('yyyyy',"'"+x_p+"'", blast)
        blast = re.sub('zzzzz',"'"+db+"'", blast)
        blast = re.sub('#wwwww', 'print("")\nprint("Finalizado")\ntime.sleep(60)', blast)
        save= open('run.py','w')
        save.write(blast)
        save.close()
        os.system('gnome-terminal -e "python3 run.py"')
###
def run_blast_jupyter(x_p = '', file = '', db = '', evalue = 0.1):
    if file == '':
        print('Intentar nuevamente')
    else:
        print('\nCorriendo '+x_p+'\n')
        version = subprocess.check_output([x_p, '-version'])
        print(version.decode())
        # 
        fasta1 = open(file,'r')
        filename = file.split('.')[0]
        name = x_p+'_'+filename.split('/')[-1]+'.tab'
        #
        print(x_p,'-db', db,'-query', file,'-evalue',str(evalue),'-outfmt\n',
        "6 qacc sacc qlen slen length score bitscore evalue pident nident\n",
        "mismatch positive gaps gapopen stitle",'-max_target_seqs','10\n',
        '-max_hsps','1','-out', name+'\n')
        #
        fasta1 = fasta1.read()
        # 
        len(re.findall('>', fasta1))
        ###
        out_text = []
        n = 0
        dld = 0
        tim = datetime.now()
        xx = datetime.now()
        secs = len(fasta1.split('>')) - 1
        sin_resultados = []
        for i in fasta1.split('>')[1:int(float(secs)) + 1]:
            i = re.sub('\n$', '', i)
            if re.search('[|]', i):
                identifier = i.split('|')[1]
            else:
                identifier = re.search('\w+', i).group()
            #identifier = re.search('\w+', i).group()
            ###
            save1= open(identifier,'w')
            save1.write('>'+i)
            save1.close()
            ###
            blast = subprocess.call([x_p,'-db',db,'-query', identifier,'-evalue',str(evalue),'-outfmt',
                                "6 qacc sacc qlen slen length score bitscore evalue pident nident mismatch positive gaps gapopen stitle",
                                '-max_target_seqs','10','-max_hsps','1','-out', identifier+'.txt'])
            ###
            out = open(identifier+'.txt','r')
            out = out.read()
            n += 1
            if len(out) == 0:
                out_bytes = len(out)
                sin_resultados.append('Secuencia '+str(n)+' | '+str(identifier))
                os.remove(identifier)
                os.remove(identifier+'.txt')
                continue
            else:
                out_text.append(pd.read_csv(StringIO(out),sep='\t',header=None))
                out_bytes = len(out)
                ###
                dif = max([i for i in range(1, out_bytes+100,100)]) - out_bytes
                total_length = int(out_bytes)
                dld += out_bytes
                dl = 0
                for dat in [i for i in range(1, out_bytes+100,100)]:
                    tim = datetime.now() - xx
                    dl = dat - dif
                    done = int(30 * dl / total_length)
                    sys.stdout.write('\rSecuencia '+str(n)+' | %s | %s' % ('{}'.format(tim).split('.')[0], identifier)) 
                    sys.stdout.flush()
            os.remove(identifier)
            os.remove(identifier+'.txt')
        ###
        if len(out_text) == 0:
            print('Secuencias sin resultados: '+str(len(sin_resultados)))
        else:
            resultados_blast = pd.concat(out_text)
            header = ('qacc','Entry','qlen','slen','length','score','bitscore','evalue','pident','nident',
                              'mismatch','positive','gaps','gapopen','stitle')
            resultados_blast.columns = header
            resultados_blast.to_csv(name, sep = '\t',index = None)
            print('\n\nTiempo: {}'.format(tim).split('.')[0], '(h:m:s)')
            print('\nResultado: ', name, '\n')
            print('Secuencias sin resultados: '+str(len(sin_resultados)))
        #for i in sin_resultados:
        #    print(i)
###



def explorar_terms(Term = random.choices(population=list(ontologia.GO.drop_duplicates()), k=50)):
    for i in ontology_file[1:len(ontology_file)]:
        if re.findall('id: '+Term, i):
            print('\n[Term]', i)
    return Term