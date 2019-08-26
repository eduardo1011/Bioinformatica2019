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
from wordcloud import WordCloud, STOPWORDS
import subprocess
from datetime import datetime
import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, interact_manual, Button, HBox, VBox, IntSlider, Label, IntRangeSlider
from ipywidgets import Checkbox, RadioButtons
from ipywidgets import Button, Layout
import random
from matplotlib import cm
from IPython.display import display
from networkx import path_graph, random_layout
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import nxviz as nxv


#modulo = urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/bar_pas_net_heat_cir.py', 'bar_pas_net_heat_cir.py')
#import bar_pas_net_heat_cir
#from bar_pas_net_heat_cir import *

modulo = urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/venn_bar_pas_net_heat_cir.py', 'venn_bar_pas_net_heat_cir.py')
import venn_bar_pas_net_heat_cir
from venn_bar_pas_net_heat_cir import *

#-------------------------------------------------------------
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
    print('\nThe go-basic.obo file is already downloaded')
    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    #print(urllib.request.urlopen(url).headers)
else:
    # Método 1: urllib.request.urlretrieve('http://purl.obolibrary.org/obo/go.obo', 'datos/go.obo')
    # Método 2:
    print('\nDownloading the Ontology: go-basic.obo')
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
                sys.stdout.write("\rDownloading Ontology (go-basic.obo) [%s%s] %s MB" % ('■' * done, ' ' * (40-done), round(dl/1000000,2)), )    
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

#>>>>>>>>>>>>>>>>



#---------------------------------------------------------------
# uniprot
###
code = {200:'The request was processed successfully.',
        400:'Bad request. There is a problem with your input.',
        404:'Not found. The resource you requested doesn’t exist.',
        410:'Gone. The resource you requested was removed.',
        500:'Internal server error. Most likely a temporary problem, but if the problem persists please contact us.',
        503:'Service not available. The server is being updated, try again later.'}
def get_UniProtKB_info0(ids = []):
    sets_ids = []
    for k in range(0, len(ids),250):
        sets_ids.append(ids[k:k+250])
    #print('Grupos:', len(sets_ids))
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
            #print('Set '+str(n)+' ('+str(len(i))+' IDs):', code[response.status_code])
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
    #print('Grupos:', len(sets_ids))
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
            #print('Set '+str(n)+' ('+str(len(i))+' IDs):', code[response.status_code])
            respuesta = response.content.decode()
            names = ['Entry', 'Entry_name', 'Protein_name', 'Organism', 'Gene', 'Length', 'GO', 'Terms', 'Tax_ID']
            df = pd.read_csv(StringIO(respuesta),sep='\t')
            df.columns = names
            uniprot.append(df)
        uniprotkb = pd.concat(uniprot).fillna('NA')
        return uniprotkb
#---------------------------------------------------------------
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
#---------------------------------------------------------------
def inputs():
    lista_seleccionada = pd.read_csv('lista_seleccionada.tsv', sep = '\t') 
    goslim_seleccionado = pd.read_csv('goslim_seleccionado.tsv', sep = '\t')
    return lista_seleccionada, goslim_seleccionado



def RUN(lista = [], GosliM = DataFrame([])):
    print('■ 4. Runing...')
    uniprotkb = get_UniProtKB_info1(ids = get_UniProtKB_info0(ids = lista).Entry)
    #uniprotkb.to_csv('uniprotkb.tsv', sep = '\t', index = None)
    anotation = []
    for index, row in uniprotkb.iterrows():
        for j in row.GO.split('; '):
            if re.search('GO:\d+', j):
                anotation.append([row.Entry, re.search('GO:\d+', j).group()])
            else:
                anotation.append([row.Entry, j])
    anotation = DataFrame(anotation, columns = ['Entry', 'GO'])
    anotation = anotation[anotation['GO'].str.contains('NA') == False]
    ontologia = ontology()
    ###
    slim = GosliM
    is_a = IS_A()
    dfsup = mapping_super_terms(list_gos = anotation.GO.drop_duplicates())
    asp = {'P':'Biological Process', 'F':'Molecular Function', 'C':'Cellular Component'}
    df2 = anotation.merge(dfsup, on = 'GO', how = 'inner')
    df3 = df2.merge(ontologia, on = 'GO', how = 'left')
    df3.columns = ['Entry', 'Original GO', 'is_a','Original Term', 'Original Aspect']
    df4 = df3.merge(is_a[['is_a','Term','Aspect']], on = 'is_a', how = 'left').drop_duplicates().reset_index(drop = True).dropna()
    df4 = df4.rename(columns={'is_a':'GO'}) # is_a se convierte en GO, para mapearlos con el GO Slim seleccionado
    slim['ASPECT'] = [asp[i] for i in slim.Aspect]
    GO_Slim = slim.merge(df4, on = ['GO', 'Term', 'Aspect'],
                         how = 'left').drop_duplicates().reset_index(drop = True).dropna().reset_index(drop = True)
    etiquetas = []    
    for i in GO_Slim.Term:
        if len(i.split(' ')) <= 3:
            etiquetas.append(re.sub(' ', '\n', i))
        else:
            etiquetas.append(''.join(i.split(' ')[0])+'\n'+''.join(i.split(' ')[1])+''.join(i.split(' ')[2])+'...') # +'\n'+i.split(' ')[2]
    GO_Slim['Short_Term'] = etiquetas 
    # exportar este data frame, será usado por las funciones interactivas
    #GO_Slim.to_csv('GO_Slim.tsv', sep = '\t', index = None)
    #print('\n\nColumnas del GO Slim: [GO, Term, Aspect, Short_Term]\n')
    print('     Finished')
    print('■ 5. REPORT')
    print('     1. Unique Original GO Terms:',GO_Slim['Original GO'].drop_duplicates().count())
    print('       • P:', GO_Slim[GO_Slim['Original Aspect'] == 'P']['Original Term'].drop_duplicates().count())
    print('       • F:', GO_Slim[GO_Slim['Original Aspect'] == 'F']['Original Term'].drop_duplicates().count())
    print('       • C:', GO_Slim[GO_Slim['Original Aspect'] == 'C']['Original Term'].drop_duplicates().count())
    print('     2. Unique GO Slim Terms:',GO_Slim.GO.drop_duplicates().count())
    print('       • P:', GO_Slim[GO_Slim['Aspect'] == 'P']['Term'].drop_duplicates().count())
    print('       • F:', GO_Slim[GO_Slim['Aspect'] == 'F']['Term'].drop_duplicates().count())
    print('       • C:', GO_Slim[GO_Slim['Aspect'] == 'C']['Term'].drop_duplicates().count())
    print('    3. Unique mapped identifiers (Proteins):',GO_Slim.Entry.drop_duplicates().count())
    print('       • P:', GO_Slim[GO_Slim['Aspect'] == 'P'].Entry.drop_duplicates().count())
    print('       • F:', GO_Slim[GO_Slim['Aspect'] == 'F'].Entry.drop_duplicates().count())
    print('       • C:', GO_Slim[GO_Slim['Aspect'] == 'C'].Entry.drop_duplicates().count())
    print('■ 6. Results Visualization')
    uniprotkb.to_csv('uniprotkb_mapping.tsv', sep = '\t', index = None)
    GO_Slim.to_csv('GO_Slim_summarize.tsv', sep = '\t', index = None)
    return uniprotkb, GO_Slim
#--------------------------------------------------------------

def outputs():
    uniprotkb_mapping = pd.read_csv('uniprotkb_mapping.tsv', sep = '\t') 
    GO_Slim_summarize = pd.read_csv('GO_Slim_summarize.tsv', sep = '\t')
    return uniprotkb_mapping, GO_Slim_summarize

#---------------------------------------------------------------
import networkx as nx
def dfUniprot_dfGOSlim(dfUniprot = DataFrame([]), dfGOSlim = DataFrame([])):    
    from ipywidgets import Button, Layout
    from matplotlib import cm
    columnas = ['GO', 'Entry', 'Aspect', 'Term', 'ASPECT', 'Short_Term']
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
    siz = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    escala_uno = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    common_colors = ['red','orange','yellow','green','blue','purple','brown','magenta','tan',
                     'cyan','olive','maroon','navy','aquamarine','turquoise','silver',
                     'lime','teal','indigo','violet','pink','black','white']
    common_colors.sort()
    label_options = {'GO':0, 'Term':1, 'Short_Term':2}
    import networkx as nx
    layouts = {'Circular':nx.circular_layout,
              'Random':nx.random_layout,
              'Shell':nx.shell_layout,
              'Spectral':nx.spectral_layout,
              'Spring':nx.spring_layout,
              'KK':nx.kamada_kawai_layout}
    design = list(layouts.keys())

    #..........................................................................................
    #..........................................................................................
    #..........................................................................................

    a = widgets.ToggleButtons(options=words,description='',disabled=False,button_style='')
    b = widgets.Dropdown(options=list(range(1,101)),value=10,description='rows_num:',disabled=False,
                        layout=Layout(width='95%', height='25px'))
    c = widgets.Dropdown(options=list(range(1,21)),value=10,description='size_x:',disabled=False,
                        layout=Layout(width='95%', height='25px'))
    d = widgets.Dropdown(options=list(range(1,21)),value=10,description='size_y:',disabled=False,
                        layout=Layout(width='95%', height='25px'))
    e = widgets.Dropdown(options=list(range(1,21)),value=10,description='xlabel:',disabled=False,
                        layout=Layout(width='95%', height='25px'))
    f = widgets.Dropdown(options=list(range(1,21)),value=10,description='size_title:',disabled=False,
                        layout=Layout(width='95%', height='25px'))
    g = widgets.Dropdown(options=list(range(1,21)),value=10,description='size_bartxt:',disabled=False,
                        layout=Layout(width='95%', height='25px'))
    #h = widgets.Dropdown(options=list(colors.keys()),description='color:',disabled=False) # colors
    h = widgets.Dropdown(options=common_colors,value='teal',description='color:',disabled=False,
                        layout=Layout(width='95%', height='25px')) # colors
    n = widgets.Dropdown(options=['GO', 'Term', 'Short_Term'],value='GO',
                          description='label_option:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    o = widgets.Dropdown(options=list(range(1,21)),value=3,description='sep:',disabled=False,
                        layout=Layout(width='95%', height='25px'))
    v = widgets.SelectionSlider(options=sizes,value='5x5',description='size_plot',disabled=False,
                              continuous_update=False,orientation='horizontal',readout=True)
    #z = widgets.Checkbox(value=False,description='Save',disabled=False)
    z = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                             tooltip='Description')
    w = widgets.VBox([v, n, b, c, d, e, f, g, h, o, z])

    def box(a, b, c, d, e, f, g, h, n, o, z, v):
        #print((a, b, c, d, e, f, g, h))
        frame1 = dfGOSlim[columnas].groupby(['ASPECT']).get_group(a).drop_duplicates()
        frame2 = DataFrame(frame1[['GO', 'Term', 'Short_Term', 'Entry']].drop_duplicates().groupby(['GO','Term', 'Short_Term']).Entry.count()).reset_index()
        frame3 = frame2.sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        aa = widgets.ToggleButtons(options=['# GO Terms  =  '+str(frame3.GO.count())],
                                      description='',
                                      disabled=False,
                                      button_style='danger')# 'success', 'info', 'warning', 'danger' or ''
        bb = widgets.ToggleButtons(options=['# Proteins  =  '+str(frame1.Entry.drop_duplicates().count())],
                                      description='',
                                      disabled=False,
                                      button_style='info')# 'success', 'info', 'warning', 'danger' or ''
        ###
        frame_uniprot = frame1[['Entry', n]].merge(dfUniprot[['Entry', 'Protein_name', 'Gene']],
                                                              on = 'Entry', how = 'left')
        frame4 = frame3.iloc[0:b].sort_values(by ='Entry',ascending=True).reset_index(drop=True)
        prot_por_term = {}
        for ii in frame4[n].drop_duplicates():
            prot_name = []
            for x, y in zip(list(frame_uniprot[frame_uniprot[n] == ii].Entry.drop_duplicates()),
                            list(frame_uniprot[frame_uniprot[n] == ii].Protein_name.drop_duplicates())):
                prot_name.append(x+' = '+y)                
            prot_por_term[ii] = prot_name
        qq = widgets.Dropdown(options=list(prot_por_term.keys()),description='Term:',disabled=False)
        xx = widgets.VBox([qq])
        def fb(qq):
            #print(qq)
            ww = widgets.Dropdown(options=prot_por_term[qq],
                              description='',disabled=False)
            display(ww)
        out = widgets.interactive_output(fb, {'qq':qq})
        display(widgets.HBox([aa, bb, VBox([xx, out])]))
        ###
    
        plt.subplots(figsize=(int(float(v.split('x')[0])),int(float(v.split('x')[0]))))
        barras(df = frame3,
               column = label_options[n],
               dim = 111,
               title = a,
               row_num = b,
               color = h,  # colors[h]
               size_x = c,
               size_y = d,
               xlabel = e,
               size_title = f,
               size_bartxt = g,
               sep = o)
        if z == True:
            plt.savefig('barras_'+a+'.png', dpi = 600, bbox_inches='tight')
        else:
            pass
    out = widgets.interactive_output(box, {'a':a, 'b': b, 'c': c, 'd':d, 'e':e, 'f':f, 'g':g,
                                           'h':h, 'n':n, 'o':o, 'z':z, 'v':v})
    def barras_interactive_plot():
        display(HBox([a]),HBox([w, out]))
    
    #..........................................................................................
    #..........................................................................................
    #..........................................................................................

    escala_uno1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    a1 = widgets.ToggleButtons(options=words,description='',disabled=False,button_style='')
    b1 = widgets.Dropdown(options=list(range(0,101)),value=5,description='rows_num:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    c1 = widgets.Dropdown(options=list(range(0,21)),value=10,description='size_text:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    d1 = widgets.Dropdown(options=list(range(0,361)),value=0,description='angle:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    e1 = widgets.Dropdown(options=list(pie_colors.keys()),value='Pastel1',description='colors:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    f1 = widgets.Dropdown(options=list(range(0,21)),value=10,description='size_title:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    g1 = widgets.SelectionSlider(options=sizes,value='5x5',description='size_plot',disabled=False,
                              continuous_update=False,orientation='horizontal',readout=True)
    n1 = widgets.Dropdown(options=['GO', 'Term', 'Short_Term'],value='GO',description='label_option:',
                          disabled=False,layout=Layout(width='95%', height='25px'))
    o1 = widgets.SelectionSlider(options=list([0.0])+list(np.round(np.linspace(0.01,0.2,10), 2)),value=0.0,
                        description='explode', disabled=False, continuous_update=False,
                        orientation='horizontal', readout=True)
    p1 = widgets.Dropdown(options=['Labels', 'Box'],value='Labels',description='legend:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    q1 = widgets.SelectionSlider(options=escala_uno1,value=0.5,description='center:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    z1 = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                             tooltip='Description')
    w1 = widgets.VBox([g1, p1, q1, n1, b1, c1, d1, f1, e1, o1, z1])

    def box1(a1, b1, c1, d1, e1, f1, g1, n1, o1, p1, q1, z1):
        #print((a, b, c, d, e, f, g, h))
        frame1 = dfGOSlim[columnas].groupby(['ASPECT']).get_group(a1).drop_duplicates()
        frame2 = DataFrame(frame1[['GO', 'Term','Short_Term', 'Entry']].drop_duplicates().groupby(['GO','Term', 'Short_Term']).Entry.count()).reset_index()
        frame3 = frame2.sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        aa1 = widgets.ToggleButtons(options=['# GO Terms = '+str(frame3.GO.count())],
                                      description='',
                                      disabled=False,
                                      button_style='danger')# 'success', 'info', 'warning', 'danger' or ''
        bb1 = widgets.ToggleButtons(options=['# Proteins  =  '+str(frame1.Entry.drop_duplicates().count())],
                                      description='',
                                      disabled=False,
                                      button_style='info')# 'success', 'info', 'warning', 'danger' or ''
        ###
        frame_uniprot = frame1[['Entry', n1]].merge(dfUniprot[['Entry', 'Protein_name', 'Gene']],
                                                              on = 'Entry', how = 'left')
        frame4 = frame3.iloc[0:b1].sort_values(by ='Entry',ascending=True).reset_index(drop=True)
        prot_por_term = {}
        for jj in frame4[n1].drop_duplicates():
            prot_name = []
            for x, y in zip(list(frame_uniprot[frame_uniprot[n1] == jj].Entry.drop_duplicates()),
                            list(frame_uniprot[frame_uniprot[n1] == jj].Protein_name.drop_duplicates())):
                prot_name.append(x+' = '+y)                
            prot_por_term[jj] = prot_name
        qq1 = widgets.Dropdown(options=list(prot_por_term.keys()),description='Term:',disabled=False)
        xx1 = widgets.VBox([qq1])
        def fb(qq1):
            #print(qq)
            ww1 = widgets.Dropdown(options=prot_por_term[qq1],
                              description='',disabled=False)
            display(ww1)
        out1 = widgets.interactive_output(fb, {'qq1':qq1})
        display(widgets.HBox([aa1, bb1, VBox([xx1, out1])]))
        ###
        plt.subplots(figsize=(int(float(g1.split('x')[0])),int(float(g1.split('x')[0]))))
        pastel(df = frame3,
               column = label_options[n1],
               dim = 111,
               title = a1,
               row_num = b1,
               angle = d1,
               legend = p1,
               size_text = c1,
               size_title = f1,
               color = pie_colors[e1],
               explode = o1,
               open_center = q1)
        if z1 == True:
            plt.savefig('pastel_'+a1+'.png', dpi = 600, bbox_inches='tight')
        else:
            pass
    out1 = widgets.interactive_output(box1, {'a1':a1, 'b1': b1, 'c1': c1, 'd1':d1, 'e1':e1,
                                             'f1':f1, 'g1':g1, 'n1':n1, 'o1':o1, 'p1':p1, 'q1':q1, 'z1':z1})

    def pastel_interactive_plot():
        display(HBox([a1]), HBox([w1, out1]))

    #..........................................................................................
    #..........................................................................................
    #..........................................................................................

    a2 = widgets.ToggleButtons(options=words,description='',disabled=False,button_style='')
    b2 = widgets.Dropdown(options=list(range(1,51)),value=10,description='diam_nodos:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    c2 = widgets.Dropdown(options=list(range(1,101)),value=50,description='inter:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    d2 = widgets.Dropdown(options=escala_uno,value=0.2,description='wid_edges:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    e2 = widgets.Dropdown(options=list(range(0,21)),value=3,description='k_num:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    f2 = widgets.Dropdown(options=['none', 'label'],value='none',description='node_label:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    g2 = widgets.Dropdown(options=common_colors,value='silver',description='col_in_min:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    h2 = widgets.Dropdown(options=common_colors,value='red',description='col_in_max:',disabled=False,
                         layout=Layout(width='95%', height='25px')) # colors
    i2 = widgets.Dropdown(options=common_colors,value='teal',description='col_node:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    j2 = widgets.Dropdown(options=escala_uno,value=0.7,description='alpha_min:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    k2 = widgets.Dropdown(options=escala_uno,value=0.5,description='alpha_max:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    l2 = widgets.Dropdown(options=escala_uno,value=0.7,description='node_alpha:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    m2 = widgets.Dropdown(options=list(range(1,21)),value=7,description='label_size:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    n2 = widgets.Dropdown(options=['GO', 'Term', 'Short_Term'],value='GO',
                          description='label_option:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    o2 = widgets.Dropdown(options=design,value='Spring',
                          description='design:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    p2 = widgets.Dropdown(options=common_colors,value='white',description='back_col:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    q2 = widgets.Dropdown(options=common_colors,value='black',description='lab_col:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    v2 = widgets.SelectionSlider(options=sizes,value='5x5',description='size_plot',disabled=False,
                              continuous_update=False,orientation='horizontal',readout=True)

    z2 = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                             tooltip='Description')
    w2 = widgets.VBox([v2, o2, f2, n2, m2, b2, c2, d2, e2, g2, h2, i2, j2, k2, l2, p2, q2, z2])

    def box2(a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2, l2, m2, n2, o2, p2, q2, z2, v2):
        #print((a, b, c, d, e, f, g, h))
        frame1 = dfGOSlim[columnas].groupby(['ASPECT']).get_group(a2).drop_duplicates()
        frame2 = DataFrame(frame1[['GO', 'Term','Short_Term', 'Entry']].drop_duplicates().groupby(['GO','Term', 'Short_Term']).Entry.count()).reset_index()
        frame3 = frame2.sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        aa2 = widgets.ToggleButtons(options=['# GO Terms  =  '+str(frame1.GO.drop_duplicates().count())],
                                      description='',
                                      disabled=False,
                                      button_style='danger')# 'success', 'info', 'warning', 'danger' or ''
        bb2 = widgets.ToggleButtons(options=['# Proteins  =  '+str(frame1.Entry.drop_duplicates().count())],
                                      description='',
                                      disabled=False,
                                      button_style='info')# 'success', 'info', 'warning', 'danger' or ''
        ###
        frame_uniprot = frame1[['Entry', n2]].merge(dfUniprot[['Entry', 'Protein_name', 'Gene']],
                                                              on = 'Entry', how = 'left')
        #frame4 = frame3.iloc[0:b1].sort_values(by ='Entry',ascending=True).reset_index(drop=True)
        prot_por_term = {}
        for jj in frame3[n2].drop_duplicates():
            prot_name = []
            for x, y in zip(list(frame_uniprot[frame_uniprot[n2] == jj].Entry.drop_duplicates()),
                            list(frame_uniprot[frame_uniprot[n2] == jj].Protein_name.drop_duplicates())):
                prot_name.append(x+' = '+y)               
            prot_por_term[jj] = prot_name
        qq2 = widgets.Dropdown(options=list(prot_por_term.keys()),description='Term:',disabled=False)
        xx2 = widgets.VBox([qq2])
        def fb(qq2):
            #print(qq)
            ww2 = widgets.Dropdown(options=prot_por_term[qq2],
                          description='',disabled=False)
            display(ww2)
        out2 = widgets.interactive_output(fb, {'qq2':qq2})
        display(widgets.HBox([aa2, bb2, VBox([xx2, out2])]))
        ###
        plt.subplots(figsize=(int(float(v2.split('x')[0])),int(float(v2.split('x')[0]))))
        net_plot(df = frame1,
                 label = f2,
                 layout = o2,
                 column = label_options[n2],
                 label_size = m2,
                 diam_nodos = b2,
                 espe_edges = d2,
                 inter = c2,
                 color_inter_min = g2, #
                 color_inter_max = h2, #
                 edge_alpha_min = j2,
                 edge_alpha_max = k2,
                 k_num = e2,
                 color_nodo = i2, #
                 node_alpha = l2,
                 backg = p2,
                 label_color = q2)
        if z2 == True:
            plt.savefig('net_'+a2+'.png', dpi = 600, bbox_inches='tight')
        else:
            pass
    out2 = widgets.interactive_output(box2, {'a2':a2, 'b2': b2, 'c2': c2, 'd2':d2, 'e2':e2, 'f2':f2,
                                            'g2':g2, 'h2':h2, 'i2':i2, 'j2':j2, 'k2':k2, 'l2':l2, 'm2':m2,
                                            'n2':n2, 'o2':o2, 'p2':p2, 'q2':q2, 'z2':z2, 'v2':v2, 'z2':z2})
    def network_interactive_plot():
        display(HBox([a2]),HBox([w2, out2]))
    
    #..........................................................................................
    #..........................................................................................
    #..........................................................................................
    heatmap_colors = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues',
                  'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                  'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg',
                  'gist_gray', 'gray', 'bone', 'pink', 'spring', 'summer', 'autumn', 'winter',
                  'cool', 'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper', 'PiYG', 'PRGn',
                  'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr',
                  'seismic', 'twilight', 'twilight_shifted', 'hsv']
    heatmap_colors.sort()
    a3 = widgets.ToggleButtons(options=words,description='',disabled=False,button_style='')
    d3 = widgets.Dropdown(options=list(range(1,21)),value=8,description='xticks_size:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    e3 = widgets.Dropdown(options=list(range(1,21)),value=8,description='yticks_size:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    f3 = widgets.Dropdown(options=list(range(1,21)),value=12,description='ylabel_size:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    z3 = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                                tooltip='Description')
    b3 = widgets.Dropdown(options=['GO', 'Term'],value='GO',
                            description='label_x:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    c3 = widgets.Dropdown(options=['GO', 'Term'],value='Term',
                            description='label_y:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    i3 = widgets.Dropdown(options=heatmap_colors,value='viridis',description='color:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    v3 = widgets.SelectionSlider(options=sizes,value='5x5',description='size_plot',disabled=False,
                              continuous_update=False,orientation='horizontal',readout=True)
    w3 = widgets.VBox([v3, b3, c3, d3, e3, f3, i3, z3])

    def box3(a3, d3, e3, f3, z3, b3, c3, i3, v3):
        frame1 = dfGOSlim[columnas].groupby(['ASPECT']).get_group(a3).drop_duplicates()
        frame2 = DataFrame(frame1[['GO', 'Term','Short_Term', 'Entry']].drop_duplicates().groupby(['GO','Term', 'Short_Term']).Entry.count()).reset_index()
        frame3 = frame2.sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        aa3 = widgets.ToggleButtons(options=['# GO Terms  =  '+str(frame1.GO.drop_duplicates().count())],
                                        description='',
                                        disabled=False,
                                        button_style='danger')# 'success', 'info', 'warning', 'danger' or ''
        bb3 = widgets.ToggleButtons(options=['# Proteins  =  '+str(frame1.Entry.drop_duplicates().count())],
                                        description='',
                                        disabled=False,
                                        button_style='info')# 'success', 'info', 'warning', 'danger' or ''
        frame_uniprot = frame1[['Entry', c3]].merge(dfUniprot[['Entry', 'Protein_name', 'Gene']],
                                                                on = 'Entry', how = 'left')
        #frame4 = frame3.iloc[0:b1].sort_values(by ='Entry',ascending=True).reset_index(drop=True)
        prot_por_term = {}
        for jj in frame3[c3].drop_duplicates():
            prot_name = []
            for x, y in zip(list(frame_uniprot[frame_uniprot[c3] == jj].Entry.drop_duplicates()),
                            list(frame_uniprot[frame_uniprot[c3] == jj].Protein_name.drop_duplicates())):
                prot_name.append(x+' = '+y)               
            prot_por_term[jj] = prot_name
        qq3 = widgets.Dropdown(options=list(prot_por_term.keys()),description='Term:',disabled=False)
        xx3 = widgets.VBox([qq3])
        def fb(qq3):
            #print(qq)
            ww3 = widgets.Dropdown(options=prot_por_term[qq3],
                            description='',disabled=False)
            display(ww3)
        out3 = widgets.interactive_output(fb, {'qq3':qq3})
        display(widgets.HBox([aa3, bb3, VBox([xx3, out3])]))
        ###
        plt.subplots(figsize=(int(float(v3.split('x')[0])),int(float(v3.split('x')[0]))))
        heatmap_plot(df = frame1,
                    colors = i3,
                    label_x = b3,
                    label_y = c3,
                    xticks_size = d3,
                    yticks_size = e3,
                    ylabel_size = f3)
        if z3 == True:
                plt.savefig('heatmap_'+a3+'.png', dpi = 600, bbox_inches='tight')
        else:
            pass
    out3 = widgets.interactive_output(box3, {'a3':a3, 'd3':d3, 'e3':e3, 'f3':f3, 'z3':z3, 'b3':b3, 'c3':c3,
                                            'i3':i3, 'v3':v3})
    def heatmap_interactive_plot():
        display(HBox([a3]),HBox([w3, out3]))

    #..........................................................................................
    #..........................................................................................
    #..........................................................................................
    color_nodo = ['Uniques', 'Umbral', 'None']

    a4 = widgets.ToggleButtons(options=words,description='',disabled=False,button_style='')
    b4 = widgets.Dropdown(options=list(range(1,21)),value=10,description='fontsize:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    c4 = widgets.Dropdown(options=list(range(5,101)),value=20,description='node _Prot >=:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    d4 = widgets.Dropdown(options=[True, False],value=False,description='label:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    e4 = widgets.Dropdown(options=color_nodo,value='None',description='node_color:',disabled=False,
                         layout=Layout(width='95%', height='25px'))
    z4 = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                                tooltip='Description')
    v4 = widgets.SelectionSlider(options=sizes,value='5x5',description='size_plot',disabled=False,
                              continuous_update=False,orientation='horizontal',readout=True)
    
    w4 = widgets.VBox([v4, b4, d4, c4, e4, z4])

    def box4(a4, b4, d4, e4, c4, v4, z4):
        frame1 = dfGOSlim[columnas].groupby(['ASPECT']).get_group(a4).drop_duplicates()        
        ###
        circos(df = frame1,
                node_color = e4,
                size = int(float(v4.split('x')[0])),
                label = d4,
                inter = c4,
                fontsize = b4)
        if z4 == True:
                plt.savefig('circos_'+a4+'.png', dpi = 600, bbox_inches='tight')
        else:
            pass
    out4 = widgets.interactive_output(box4, {'a4':a4,'b4':b4,'d4':d4, 'e4':e4,'c4':c4,'v4':v4,'z4':z4})
    def circos_interactive_plot():
        display(HBox([a4]),HBox([w4, out4]))
    
    #..........................................................................................
    #..........................................................................................
    #..........................................................................................

    words2 = {'Biological Process':'P', 'Molecular Function':'F', 'Cellular Component':'C'} # nuevo
    escala_uno1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1] # nuevo
    a5 = widgets.ToggleButton(value=True,description='Biological Process',disabled=False,button_style='',
                                tooltip='Description')
    a51 = widgets.ToggleButton(value=True,description='Molecular Function',disabled=False,button_style='',
                                tooltip='Description')
    a52 = widgets.ToggleButton(value=True,description='Cellular Component',disabled=False,button_style='',
                                tooltip='Description')
    v5 = widgets.SelectionSlider(options=sizes,value='5x5',description='size_plot:',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True)
    z5 = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                                tooltip='Description')
    f5 = widgets.Dropdown(options=['none', 'Aspects', 'Abbreviation', 'Name_with jump'],value='none',description='node_label:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    g5 = widgets.SelectionSlider(options=[0, 1, 2, 3, 4, 5, 6, 7],value=1,description='linewidth:',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True)
    h5 = widgets.Dropdown(options=escala_uno1,value=0.5,description='alpha_sets:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    i5 = widgets.Dropdown(options=[False, 'bold'],value=False,description='font_sets:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    j5 = widgets.Dropdown(options=list(range(1,21)),value=15,description='size_sets:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    k5 = widgets.Dropdown(options=escala_uno1,value=0.5,description='alpha_inter:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    l5 = widgets.Dropdown(options=[False, 'bold'],value=False,description='font_inter:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    m5 = widgets.Dropdown(options=list(range(1,31)),value=15,description='size_inter:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    n5 = widgets.Dropdown(options=list(range(1,51)),value=17,description='size_label:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    o5 = widgets.Dropdown(options=[False, 'bold'],value=False,description='font_label:',disabled=False,
                            layout=Layout(width='95%', height='25px'))
    p5 = widgets.Dropdown(options=common_colors,value='black',description='color_border:',disabled=False,
                        layout=Layout(width='95%', height='25px')) # colors
    uxu = Button(layout=Layout(width='95%', height='5px'), disabled=True)
    uxu.style.button_color = 'purple'

    w5 = widgets.HBox([a5, a51, a52])
    w51 = widgets.VBox([v5, f5, n5, o5, uxu, h5, i5, j5, uxu, k5, l5, m5, uxu, g5, p5, z5])

    def box5(a5, a51, a52, f5,  g5, h5, i5, j5, k5, l5, m5, n5, o5, p5, v5, z5):
        seleccion = [a5, a51, a52]
        if seleccion.count(True) < 2:
            print('!!!Select at least two groups!!!')
        else:
            sel = [] # guarda las categorias seleccionadas
            if a5:
                sel.append(words2['Biological Process'])
            if a51:
                sel.append(words2['Molecular Function'])
            if a52:
                sel.append(words2['Cellular Component'])
            #print(sel)
            intersecciones = {}
            for i in sel: # usa la lista de categorias seleccionadas para extraer info del df final
                intersecciones[i] = set(dfGOSlim[dfGOSlim.Aspect == i].Entry)
        
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if len(sel) == 2:
                #
                if f5 == 'none':
                    lab_venn = {sel[0]:'Set1', sel[1]:'Set2'}
                if f5 == 'Aspects':
                    lab_venn = {'P':'Biological Process', 'F':'Molecular Function', 'C':'Cellular Component'} # nuevo
                if f5 == 'Abbreviation':
                    lab_venn = {'P':'P', 'F':'F', 'C':'C'} # nuevo
                if f5 == 'Name_with jump':
                    lab_venn = {'P':'Biological\nProcess', 'F':'Molecular\nFunction', 'C':'Cellular\nComponent'} # nuevo
                #
                plt.subplots(figsize=(int(float(v5.split('x')[0])),int(float(v5.split('x')[0]))))
                venn2_plot(intersecciones[sel[0]], intersecciones[sel[1]],
                        lab_set1 = lab_venn[sel[0]],
                        lab_set2 = lab_venn[sel[1]],
                        linewidth = g5,
                        color_line = p5,
                        alpha_sets = h5,
                        font_sets = i5,
                        size_vals_sets = j5,
                        alpha_inter = k5,
                        font_inter = l5,
                        size_vals_inter = m5,
                        size_label = n5,
                        font_label = o5)
                if z5 == True:
                    plt.savefig('Venn_'+'_'.join(sel)+'.png', dpi = 600, bbox_inches='tight')
                else:
                    pass
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if len(sel) == 3:
                #
                if f5 == 'none':
                    lab_venn = {sel[0]:'Set1', sel[1]:'Set2', sel[2]:'Set3'}
                if f5 == 'Aspects':
                    lab_venn = {'P':'Biological Process', 'F':'Molecular Function', 'C':'Cellular Component'} # nuevo
                if f5 == 'Abbreviation':
                    lab_venn = {'P':'P', 'F':'F', 'C':'C'} # nuevo
                if f5 == 'Name_with jump':
                    lab_venn = {'P':'Biological\nProcess', 'F':'Molecular\nFunction', 'C':'Cellular\nComponent'} # nuevo
                #
                plt.subplots(figsize=(int(float(v5.split('x')[0])),int(float(v5.split('x')[0]))))
                venn3_plot(intersecciones[sel[0]], intersecciones[sel[1]], intersecciones[sel[2]],
                        lab_set1 = lab_venn[sel[0]],
                        lab_set2 = lab_venn[sel[1]],
                        lab_set3 = lab_venn[sel[2]],
                        linewidth = g5,
                        color_line = p5,
                        alpha_sets = h5,
                        font_sets = i5,
                        size_vals_sets = j5,
                        alpha_inter = k5,
                        font_inter = l5,
                        size_vals_inter = m5,
                        size_label = n5,
                        font_label = o5)
                #>>>>>>>>>>>>>>
                if z5 == True:
                    plt.savefig('Venn_'+'_'.join(sel)+'.png', dpi = 600, bbox_inches='tight')
                else:
                    pass
    out5 = widgets.interactive_output(box5, {'a5':a5, 'a51':a51, 'a52':a52, 'f5':f5, 'g5':g5,
                                            'h5':h5, 'i5':i5, 'j5':j5, 'k5':k5, 'l5':l5, 'm5':m5,
                                            'n5':n5, 'o5':o5, 'p5':p5, 'v5':v5, 'z5':z5})

    #..........................................................................................
    #..........................................................................................
    #..........................................................................................

    www = Button(layout=Layout(width='94%', height='15px'), disabled=True)
    www.style.button_color = 'white'
    uuu = Button(layout=Layout(width='94%', height='5px'), disabled=True)
    uuu.style.button_color = 'black'
    graficas = {'VENN CHART': (w5, HBox([out5, w51])),
                'BAR CHART': (HBox([a]), HBox([w, out])),
                'PIE CHART': (HBox([a1]), HBox([w1, out1])),
                'NETWORK CHART': (HBox([a2]),HBox([w2, out2])),
                'HEATMAP':(HBox([a3]),HBox([w3, out3])),
                'CIRCOS CHART': (HBox([a4]), HBox([out4, w4]))}
    exe = widgets.ToggleButtons(options=list(graficas.keys()),disabled=False,button_style='')
    def ss(exe):
        display(graficas[exe][0], graficas[exe][1])
    ou = widgets.interactive_output(ss, {'exe':exe})

    return display(VBox([www, uuu, exe, uuu, www, ou]))

#------------------------------------------------------


#----------------------------------------------------------------


import tkinter as tk
from tkinter import filedialog

up = widgets.Checkbox(description='Upload file',value=False,disabled=False)
goslim = [re.search('go\w+', i).group() for i in ontology_file[0].split('\n') if re.search(' go\w+', i)]
gs = widgets.Dropdown(options=goslim,description='GO Slim:',disabled=False,
                    layout=Layout(width='40%', height='25px'))
yyy = Button(layout=Layout(width='94%', height='5px'), disabled=True)
yyy.style.button_color = 'red'
xxx = Button(layout=Layout(width='94%', height='5px'), disabled=True)
xxx.style.button_color = 'white'
run = widgets.Checkbox(description='Run GO Slim',value=False,disabled=False)
def file_goslim(up):
    if up == True:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        file_path = filedialog.askopenfilename()
        root.destroy()
        if file_path == '':
            print('\n!!!!!!! File not found !!!!!!!')
            sys.exit()
        else:
            print('■ 1.',file_path)
            df = pd.read_csv(file_path, sep = '\t', header = None).dropna()
            df.to_csv('lista_seleccionada.tsv', sep = '\t', index = None)
            ###
            new_url = 'https://www.uniprot.org/uniprot/?query=id:'+df.iloc[0][0]+'&format=tab&columns=id,organism'
            respuesta = requests.get(new_url).content.decode()
            dff = pd.read_csv(StringIO(respuesta),sep='\t')
            print('■ 2. Organism identified: ',dff.Organism[0])
salida = widgets.interactive_output(file_goslim, {'up':up})
##
def file_goslim2(gs):
    slim = []
    for i in ontology_file[1:len(ontology_file)]:
        if re.findall('subset: '+gs, i):
            slim.append(i.split('\n')[1].split(': ')[1])
        else:
            continue
    slim = DataFrame(slim, columns = ['GO'])
    slim = slim.merge(ontologia, on = 'GO', how = 'left').dropna()
    slim = slim[slim['Term'].str.contains('biological_process|molecular_function|cellular_component') == False]
    slim.to_csv('goslim_seleccionado.tsv', sep = '\t', index = None)
    print('■ 3.',gs)
salida2 = widgets.interactive_output(file_goslim2, {'gs':gs})
##
def file_goslim3(run):
    if run == True:
        result = RUN(lista = inputs()[0]['0'], GosliM = inputs()[1])
        dfUniprot_dfGOSlim(dfUniprot = result[0], dfGOSlim = result[1])
salida3 = widgets.interactive_output(file_goslim3, {'run':run})

in1 = VBox([yyy, xxx, up, salida, xxx])
in2 = VBox([yyy, xxx, gs, salida2, xxx])
in3 = VBox([yyy, xxx, run, salida3, xxx, yyy])

GO_Slim_Annotation = VBox([in1, in2, in3])
#--------------------------------------------------------------
