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

9. heatmap_plot(df = DataFrame([]), colors = 'Spectral', label_x = 'GO',
            label_y = 'Term', size_plot = 10, xticks_size = 12,
            yticks_size = 12, ylabel_size = 18)

10. go_file()

11. run_Blast_Windows(x_p = '', file = '', db = '')

12. run_Blast_Linux(x_p = '', file = '', db = '')

13. run_blast_jupyter(x_p = '', file = '', db = '')

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
from wordcloud import WordCloud, STOPWORDS
import subprocess
from datetime import datetime

import Orange
from orangecontrib.bio import go

code = {200:'The request was processed successfully.',
        400:'Bad request. There is a problem with your input.',
        404:'Not found. The resource you requested doesn’t exist.',
        410:'Gone. The resource you requested was removed.',
        500:'Internal server error. Most likely a temporary problem, but if the problem persists please contact us.',
        503:'Service not available. The server is being updated, try again later.'}
###
if os.path.exists('datos/go-basic.obo'):
    print('\nYa está descargado el archivo go-basic.obo\n')
    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    print(urllib.request.urlopen(url).headers)
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
    print(urllib.request.urlopen(url).headers)
###
# definir la ontologia, ubicación del archivo go-basic.obo
ontology = go.Ontology('datos/go-basic.obo')
###
def barras(df = DataFrame([]), column = 1, dim = 111, title = 'Title', row_num = 10, color = '#ff7f0e',
           size_x = 15, size_y = 15, xlabel = 20, size_title = 25, size_bartxt = 12):
    if len(df) == 0:
        print('Data frame sin datos')
    else:
        plt.subplot(dim, facecolor= 'white')
        barWidth = 0.9
        val = max(list(df.iloc[0:row_num,2]))
        plt.barh(list(df.iloc[0:row_num,column]) ,list(df.iloc[0:row_num,2]),
                 color= color,
                 align='center',
                 height= 0.95,
                 edgecolor='white')
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_position(('data',-1))
        plt.gca().spines['left'].set_visible(False)
        plt.title(title,size= size_title,fontweight='bold',loc='left')
        plt.tick_params(axis="y", color="gray")
        plt.yticks(size = size_y ,fontweight='bold')
        #plt.xticks(size = size_x ,fontweight='bold')
        if max(list(df.iloc[0:row_num,2])) > 200:
            escala = 50
        if max(list(df.iloc[0:row_num,2])) < 200:
            escala = 20
        plt.xticks(range(0,val,escala),size=size_x,fontweight='bold')
        plt.xlabel("Proteins",size= xlabel,fontweight='bold')

        for j, k in zip(list(df.iloc[0:row_num,2]),range(0,len(list(df.iloc[0:row_num,1])))):
            plt.text(j+3,k-0.2,j,
                     size=size_bartxt,ha='left',color= 'black',fontweight='bold')
###
def pastel(df = DataFrame([]), column = 0, dim = 111, title = 'Title', row_num = 3, angle = 0,
           size_text = 15, size_title = 20):
    if len(df) < 1:
        print('Data frame sin datos')
    else:
        plt.subplot(dim)
        labels = list(df.iloc[0:row_num,column])
        sizes = list(df.iloc[0:row_num,2])
        explode = np.repeat(0.03,row_num)
        plt.pie(sizes,labels=labels, autopct='%1.1f%%', startangle=angle, radius = 1,
                explode = explode,textprops=dict(size = size_text,fontweight='bold'))

        centre_circle = plt.Circle((0,0),0,fc='white')
        #centre_circle = plt.Circle((0,0),0.50,fc='white')
        plt.title(title,size=size_title,fontweight='bold')
        plt.gca().add_artist(centre_circle)
###
def word_cloud(df = DataFrame([]), size = 10):
    if len(df) < 3:
        print('Datos insificuentes, mínimo 3 patrones diferentes')
    else:
        frecuencia = []
        for index, row in df.iterrows():
            frec = list(np.repeat(re.sub(' |, |-','_', row.Term), row.Entry))
            frecuencia.append(' '.join(frec))
        text = ' '.join(frecuencia)
        ovalo = urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/ovalo.jpg', 'ovalo.jpg')
        wave_mask = np.array(Image.open( "ovalo.jpg"))
        wordcloud = WordCloud(mask=wave_mask,collocations=False,max_font_size=1000,
                              #contour_color='firebrick',contour_width=10, stopwords = set(STOPWORDS),
                              max_words=2000, background_color="white").generate(text)
        plt.subplots(figsize=(size,20))
        plt.imshow(wordcloud, interpolation="bilinear")
        plt.axis("off")
        #plt.show()
###
def GO_DAG(go_terms = [], term = ''):
    if len(go_terms) == 0:
        print('Lista vacía')
    else:
        if term == 'none':
            formato1 = []
            for i in go_terms:
                formato = {ontology[i].id:{"title":ontology[i].id, "body":ontology[i].name, "fill":"gold", "font":"black", "border":"black"}}
                formato1.append(['"'+ontology[i].id+'"', json.dumps(formato[ontology[i].id])])
            dag = re.sub('}:','},', '{'+':'.join([':'.join(i) for i in formato1])+'}')
        
        else:
            formato1 = []
            for i in go_terms:
                if ontology[i].id == term:
                    formato = {ontology[i].id:{"title":ontology[i].id, "body":ontology[i].name, "fill":"cyan", "font":"black", "border":"black"}}
                    formato1.append(['"'+ontology[i].id+'"', json.dumps(formato[ontology[i].id])])
                else:
                    formato = {ontology[i].id:{"title":ontology[i].id, "body":ontology[i].name, "fill":"gold", "font":"black", "border":"black"}}
                    formato1.append(['"'+ontology[i].id+'"', json.dumps(formato[ontology[i].id])])
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
def venn3_plot(set1 = set(), set2 = set(), set3 = set(), lab_set1 = 'Set1', lab_set2 = 'Set2',
               lab_set3 = 'Set3', size_label = 20, size_vals = 20):
    if len(set1 & set2 & set3) <= 0:
        print('No hay intersecciones entre los datos')
    else:
        v = venn3([set1, set2, set3],
                  ([lab_set1, lab_set2, lab_set3]))
        c = venn3_circles(subsets = [set1, set2, set3],
                        linestyle='dashed', linewidth=2, color="k")
        v.get_patch_by_id('111').set_color('w')
        v.get_patch_by_id('111').set_alpha(1)
        v.get_label_by_id('111').set_color('k')
        v.get_label_by_id('111').set_fontweight('bold')
        for text in v.set_labels:
            text.set_fontsize(size_label)
            text.set_fontweight('bold')
        for text in v.subset_labels:
            text.set_fontsize(size_vals)
        #plt.show()
###
# Network
from networkx import path_graph, random_layout
import matplotlib.pyplot as plt
import networkx as nx

def net_plot(df = DataFrame([]),label = 'none',diam_nodos = 10, espe_edges = 0.1, inter = 10, color_inter_min = 'k',color_inter_max = 'blue',
             edge_alpha_min = 0.3, edge_alpha_max = 0.3, k_num = 3, color_nodo = 'red', node_alpha = 0.7):
    #df = df[df.Aspect == df.Aspect.iloc[0]].drop_duplicates()
    df1 = df.groupby(['Aspect']).get_group(df.Aspect.iloc[0])[['GO', 'Entry']].drop_duplicates().merge(df, on = 'Entry', how = 'left').drop_duplicates()
    df2 = DataFrame(df[['GO', 'Term', 'Entry']].drop_duplicates().groupby(['GO','Term']).Entry.count()).reset_index()
    #### >>>>>>>>>>>>>>>>
    #### A partir de una matriz de datos extrae valores no redundantes
    matrix = df1.pivot_table(values='Entry',index=['GO_x'],aggfunc=len,columns='GO_y')
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
    df_mat = DataFrame(df_mat, columns = ['go0', 'go1', 'val'])
    ###
    nodos = []
    for index, row in df_mat.iterrows():
        if row.go0 == row.go1:
            #print(row.go0, row.go1)
            continue
        else:
            #print(row.go0, row.go1)
            nodos.append([row.go0, row.go1, row.val])
    #### >>>>>>>>>>>>>>>>
    ###
    nodos = DataFrame(nodos).dropna()
    ####################
    # https://networkx.github.io/documentation/networkx-2.3/auto_examples/drawing/plot_weighted_graph.html#sphx-glr-auto-examples-drawing-plot-weighted-graph-py
    G=nx.Graph()
    for index, row in nodos.iterrows():
        G.add_edge(row[0], row[1],weight = row[2])
    elarge=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] >= inter]
    esmall=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] < inter]
    ###
    #circular_layout
    #random_layout
    #shell_layout
    #spring_layout
    #spectral_layout
    pos=nx.spring_layout(G, k = k_num) # positions for all nodes
    # nodes
    #------------------------------------------------------------------------------
    # ordenar los valores para representarlos en el tama;o del nodo
    order = []
    for index, row in nodos.iterrows():
        order.append(row[0])
        order.append(row[1])
    orden3 = DataFrame(order).drop_duplicates(keep = 'first').reset_index(drop = True)
    orden3.columns = ['GO']
    orden4 = pd.merge(orden3, df2, on = 'GO', how = 'left')
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
    if label == 'label':
        nx.draw_networkx_labels(G,pos,font_size=10,font_weight='bold') # ,font_weight='bold'
        plt.axis('off')
        #plt.show() # display
    if label == 'none':
        plt.axis('off')
        #plt.show() # display
###

###
def get_UniProtKB_info1(ids = []):
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
            params = {'from':'ACC','to':'ACC','format':'tab','query': ' '.join(i),'columns':'id,entry name,protein names,length,go-id,organism'}
            data = urllib.parse.urlencode(params)
            url = 'https://www.uniprot.org/uploadlists/'
            response = requests.get(url, data)
            print('Set '+str(n)+' ('+str(len(i))+' IDs):', code[response.status_code])
            respuesta = response.content.decode()
            names = ['Entry', 'Entry_name', 'Protein_name', 'Length', 'GO', 'Organism']
            df = pd.read_csv(StringIO(respuesta),sep='\t',header=None).drop(columns = [6]).drop(index = [0])
            df.columns = names
            uniprot.append(df)
        uniprotkb = pd.concat(uniprot).fillna('NA')
        return uniprotkb
###
def get_UniProtKB_info2(id_organism = 0):
    if id_organism == 0:
        print('ID de organismo incorrecto')
    else:
        link = 'https://uniprot.org/uniprot/?query=organism:'+str(id_organism)+'&format=tab&columns=id,entry name,protein names,length,go-id,organism'
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
                    done = int(0.25 * dl / total_length)
                    sys.stdout.write("\rLoading [%s%s] %s MB" % ('■' * done, ' ' * (5-done), round(dl/1000000,2)), ) 
                    sys.stdout.flush()
        ###
        acc_uniprot_GO_id=pd.read_csv('annotation_'+str(id_organism),sep='\t')
        names = ['Entry', 'Entry_name', 'Protein_name', 'Length', 'GO', 'Organism']
        acc_uniprot_GO_id.columns = names
        return acc_uniprot_GO_id
###
import seaborn as sns
def heatmap_plot(df = DataFrame([]), colors = 'Spectral', label_x = 'GO', label_y = 'Term', size_plot = 10, xticks_size = 12, yticks_size = 12, ylabel_size = 18):
    if len(df) == 0:
        print('Data frame sin valores')
    else:
        df1 = df.groupby(['Aspect']).get_group(df.Aspect.iloc[0])[['GO', 'Entry']].drop_duplicates().merge(df, on = 'Entry', how = 'left').drop_duplicates()
        new_df1 = []
        for index, row in df1.iterrows():
            if row.GO_x == row.GO_y:
                continue
            else:
                new_df1.append([row.GO_x, row.Entry, row.GO_y, row.Aspect, row.Term])
        df1 = DataFrame(new_df1, columns =['GO_x','Entry','GO_y','Aspect','Term'])
        matrix = df1.pivot_table(values='Entry',index=['GO_x'],aggfunc=len,columns='GO_y')
        matrix2 = df1.pivot_table(values='Entry',index=['Term'],aggfunc=len,columns='GO_y')
        ###
        plt.subplots(figsize=(size_plot,size_plot))
        plt.gca().set_facecolor('whitesmoke')
        plt.xticks(size = xticks_size) # fontweight='bold'
        plt.yticks(size = yticks_size) # fontweight='bold'

        if label_x == 'GO':
            xlab = list(matrix.columns)
        if label_x == 'Term':
            xlab = list(matrix2.index)
        if label_y == 'GO':
            ylab = list(matrix.columns)
        if label_y == 'Term':
            ylab = list(matrix2.index)
        
        sns.heatmap(np.array(matrix),  square=True, annot = False, cmap = colors,
                #annot_kws={"size": 13, 'fontweight': 'bold'},fmt= '.0f',
                linewidths=0.005, linecolor='w',
                cbar=True, cbar_kws={"shrink": 0.5}, # 'label': '
                xticklabels= xlab,
                yticklabels= ylab)
        plt.gca().figure.axes[-1].set_ylabel('Interaction degrees\n(Proteins', size= ylabel_size)
        #plt.xticks(rotation=70)

        #plt.show()
###
def go_file():
    if os.path.exists('datos/go-basic.obo'):
        print('\nYa está descargado el archivo go-basic.obo\n')
        url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        print(urllib.request.urlopen(url).headers)
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
        print(urllib.request.urlopen(url).headers)
    return go1

###
def run_Blast_Windows(x_p = '', file = '', db = ''):
    if file == '':
        print('Intentar nuevamente')
    else:
        urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/runBlast.py', 'runBlast.py')
        blast = open('runBlast.py','r')
        blast = blast.read()
        blast = re.sub('xxxxx',"'"+file+"'", blast)
        blast = re.sub('yyyyy',"'"+x_p+"'", blast)
        blast = re.sub('zzzzz',"'"+db+"'", blast)
        save= open('run.py','w')
        save.write(blast)
        save.close()
        runblast = os.system("start cmd /k python run.py")
###
def run_Blast_Linux(x_p = '', file = '', db = ''):
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
def run_blast_jupyter(x_p = '', file = '', db = ''):
    if file == '':
        print('Intentar nuevamente')
    else:
        print('\nCorriendo '+x_p+'\n')
        # 
        fasta1 = open(file,'r')
        filename = file.split('.')[0]
        name = x_p+'_'+filename.split('/')[-1]+'.tab'
        #
        print(x_p,'-db', db,'-query', file,'-evalue','1E-6','-outfmt\n',
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
            blast = subprocess.call([x_p,'-db',db,'-query', identifier,'-evalue','1E-6','-outfmt',
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
                    sys.stdout.write('\rSecuencia '+str(n)+' | %s | %s [%s%s] Bytes %s' % ('{}'.format(tim).split('.')[0], identifier,'■' * done, ' ' * (30 - done), dl)) 
                    sys.stdout.flush()
            os.remove(identifier)
            os.remove(identifier+'.txt')
        ###
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
