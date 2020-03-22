import datetime
import re
from pandas import Series, DataFrame 
from pandas import DataFrame
import pandas
version = pandas.__version__
if float(''.join(re.findall('[0-9]{1,3}[.][0-9]{1,3}', version))) >= 0.25:
    from io import StringIO
elif float(''.join(re.findall('[0-9]{1,3}[.][0-9]{1,3}', version))) < 0.25:
    from pandas.compat import StringIO
import pandas as pd
import csv
import pathlib
from urllib.request import urlopen
import urllib.request
import requests

import webbrowser
import shutil, os
import numpy as np
from time import sleep
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import datetime
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

if os.path.exists('tax'): shutil.rmtree('tax')
os.makedirs('img',exist_ok=True)
os.makedirs('tax',exist_ok=True)



def df_input(dfUniprot = DataFrame([])):

    num_phylum =  len(dfUniprot.PHYLUM.drop_duplicates())
    from matplotlib import cm
    set3 = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Set3(np.arange(12)/12.)]
    set2 = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Set2(np.arange(8)/8.)]
    set1 = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Set1(np.arange(9)/9.)]
    pastel2 = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Pastel2(np.arange(8)/8.)]
    pastel1 = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Pastel1(np.arange(9)/9.)]
    dark2 = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Dark2(np.arange(8)/8.)]
    paired = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Paired(np.arange(12)/12.)]
    accent = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Accent(np.arange(8)/8.)]
    spectral = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.Spectral(np.arange(11)/11.)]
    tab20 = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.tab20(np.arange(20)/20.)]
    tab20b = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.tab20b(np.arange(20)/20.)]
    tab20c = [matplotlib.colors.rgb2hex(tuple(i)) for i in cm.tab20c(np.arange(20)/20.)]

    Colors1 = set2 + set1 + dark2 + paired + accent + spectral + tab20 + tab20b  + tab20c
    Colors2 = accent + spectral + tab20 + tab20b  + tab20c + set1 + set2 + dark2 + paired
    Colors3 = dark2 + paired + accent + spectral + tab20 + tab20b  + tab20c + set1 + set2
    Colors4 = tab20b  + tab20c + set1 + set2 + dark2 + paired + accent + spectral + tab20
    Colors5 = spectral + tab20 + tab20b  + tab20c + set1 + set2 + dark2 + paired + accent



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
    circle_colors = {'Colors1':Colors1[0:num_phylum],
                     'Colors2':Colors2[0:num_phylum],
                     'Colors3':Colors3[0:num_phylum],
                     'Colors4':Colors4[0:num_phylum],
                     'Colors5':Colors5[0:num_phylum]}

    def tax_colors(color_list = circle_colors['Colors1'], taxx = dfUniprot):
        tax_cols = ['Entry', 'Tax_ID','KINGDOM','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES', 'Organism']
        new2 = taxx[tax_cols].drop_duplicates()
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>
        phylum0 = new2.groupby(['PHYLUM']).Entry.count().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        asign_color = {}
        for i, j in zip(phylum0.PHYLUM, color_list):
            if i == 'NA':
                asign_color[i] = 'black'
            else:
                asign_color[i] = j
        phylum0['phy_col'] = list(asign_color.values())
        # distribución de Class
        phylum1 = new2.groupby(['PHYLUM', 'CLASS']).Entry.count().reset_index()
        class0 = []
        class0_colors = []
        for i in phylum0.PHYLUM:
            for j in phylum1.PHYLUM:
                if i == j:
                    class0_colors.append(asign_color[j])
                    class0.append(phylum1[phylum1.PHYLUM == i].sort_values(by ='Entry',ascending=False).reset_index(drop=True))
                else:
                    pass
        class1 = pd.concat(class0).drop_duplicates()
        class1['class_col'] = class0_colors
        class0_colors_corregido = []
        for index, row in class1.iterrows():
            if row.PHYLUM == 'NA':
                if row.CLASS == 'NA':
                    class0_colors_corregido.append(row.class_col)
                else:
                    class0_colors_corregido.append('grey')
            else:
                if row.CLASS == 'NA':
                    class0_colors_corregido.append('black')
                else:
                    class0_colors_corregido.append(row.class_col)
        class1['class_col'] = class0_colors_corregido
        class11 = class1.groupby(['CLASS']).Entry.sum().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        class11 = class11.merge(class1[['CLASS','class_col']].drop_duplicates(), on = 'CLASS', how = 'left')
        # distribución de Order
        phylum2 = new2.groupby(['PHYLUM', 'CLASS', 'ORDER']).Entry.count().reset_index()
        order0 = []
        order0_colors = []
        for i in phylum0.PHYLUM:
            for j in phylum2.PHYLUM:
                if i == j:
                    order0_colors.append(asign_color[j])
                    order0.append(phylum2[phylum2.PHYLUM == i].sort_values(by ='Entry',ascending=False).reset_index(drop=True))
                else:
                    pass
        order1 = pd.concat(order0).drop_duplicates()
        order1['order_col'] = order0_colors
        order0_colors_corregido = []
        for index, row in order1.iterrows():
            if row.PHYLUM == 'NA':
                if row.ORDER == 'NA':
                    order0_colors_corregido.append(row.order_col)
                else:
                    order0_colors_corregido.append('grey')
            else:
                if row.ORDER == 'NA':
                    order0_colors_corregido.append('black')
                else:
                    order0_colors_corregido.append(row.order_col)
        order1['order_col'] = order0_colors_corregido
        order11 = order1.groupby(['ORDER']).Entry.sum().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        order11 = order11.merge(order1[['ORDER','order_col']].drop_duplicates(), on = 'ORDER', how = 'left')
        # distribución de Genus
        phylum3 = new2.groupby(['PHYLUM', 'CLASS', 'ORDER', 'GENUS']).Entry.count().reset_index()
        genus0 = []
        genus0_colors = []
        for i in phylum0.PHYLUM:
            for j in phylum3.PHYLUM:
                if i == j:
                    genus0_colors.append(asign_color[j])
                    genus0.append(phylum3[phylum3.PHYLUM == i].sort_values(by ='Entry',ascending=False).reset_index(drop=True))
                else:
                    pass
        genus1 = pd.concat(genus0).drop_duplicates()
        genus1['genus_col'] = genus0_colors
        genus0_colors_corregido = []
        for index, row in genus1.iterrows():
            if row.PHYLUM == 'NA':
                if row.GENUS == 'NA':
                    genus0_colors_corregido.append(row.genus_col)
                else:
                    genus0_colors_corregido.append('grey')
            else:
                if row.GENUS == 'NA':
                    genus0_colors_corregido.append('black')
                else:
                    genus0_colors_corregido.append(row.genus_col)
        genus1['genus_col'] = genus0_colors_corregido
        genus11 = genus1.groupby(['GENUS']).Entry.sum().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        genus11 = genus11.merge(genus1[['GENUS','genus_col']].drop_duplicates(), on = 'GENUS', how = 'left')
        # distribución de Organism
        phylum4 = new2.groupby(['PHYLUM', 'CLASS', 'ORDER', 'GENUS', 'Organism']).Entry.count().reset_index()
        org0 = []
        org0_colors = []
        for i in phylum0.PHYLUM:
            for j in phylum4.PHYLUM:
                if i == j:
                    org0_colors.append(asign_color[j])
                    org0.append(phylum4[phylum4.PHYLUM == i].sort_values(by ='Entry',ascending=False).reset_index(drop=True))
                else:
                    pass
        org1 = pd.concat(org0).drop_duplicates()
        org1['org_col'] = org0_colors
        org0_colors_corregido = []
        for index, row in org1.iterrows():
            if row.PHYLUM == 'NA':
                if row.Organism == 'NA':
                    org0_colors_corregido.append(row.org_col)
                else:
                    org0_colors_corregido.append('grey')
            else:
                if row.Organism == 'NA':
                    org0_colors_corregido.append('black')
                else:
                    org0_colors_corregido.append(row.org_col)
        org1['org_col'] = org0_colors_corregido
        org11 = org1.groupby(['Organism']).Entry.sum().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True)
        org11 = org11.merge(org1[['Organism','org_col']].drop_duplicates(), on = 'Organism', how = 'left')
        os.makedirs('tax',exist_ok=True)
        return phylum0.to_csv('tax/phylum0.tsv', sep = '\t', index = None),\
            class1.to_csv('tax/class1.tsv', sep = '\t', index = None),\
            class11.to_csv('tax/class11.tsv', sep = '\t', index = None),\
            order1.to_csv('tax/order1.tsv', sep = '\t', index = None),\
            order11.to_csv('tax/order11.tsv', sep = '\t', index = None),\
            genus1.to_csv('tax/genus1.tsv', sep = '\t', index = None),\
            genus11.to_csv('tax/genus11.tsv', sep = '\t', index = None),\
            org1.to_csv('tax/org1.tsv', sep = '\t', index = None),\
            org11.to_csv('tax/org11.tsv', sep = '\t', index = None)

    alfas = {'Lineage*':[1, 1, 1, 1, 1],
            'Phylum':[1, 0.3, 0.3, 0.3, 0.3],
            'Class':[0.3, 1, 0.3, 0.3, 0.3],
            'Order':[0.3, 0.3, 1, 0.3, 0.3],
            'Genus':[0.3, 0.3, 0.3, 1, 0.3],
            'Species':[0.3, 0.3, 0.3, 0.3, 1],
            'Gradient1*':[1, 0.85, 0.7, 0.55, 0.4],
            'Gradient2*':[0.4, 0.55, 0.7, 0.85, 1],
            'Attenuate*':[0.3, 0.3, 0.3, 0.3, 0.3],
            'None*':[0, 0, 0, 0, 0]}

    
    def circle_lineage(alphas = alfas['Phylum']):
        #fig, ax = plt.subplots(111, facecolor= 'white')
        #fig, ax = plt.subplot(111)
        phylum0 = pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA')
        class1 = pd.read_csv('tax/class1.tsv', sep = '\t').fillna('NA')
        order1 = pd.read_csv('tax/order1.tsv', sep = '\t').fillna('NA')
        genus1 = pd.read_csv('tax/genus1.tsv', sep = '\t').fillna('NA')
        org1 = pd.read_csv('tax/org1.tsv', sep = '\t').fillna('NA')   

        radio = 0.5

        linaje = [phylum0, class1, order1, genus1, org1]
        #colores = [list(asign_color.values()), class0_colors, order0_colors, genus0_colors, org0_colors] 
        colores = ['phy_col', 'class_col', 'order_col', 'genus_col', 'org_col']
        pat = []
        size=-.205
        for i, j, k in zip(linaje, colores, alphas):
            size+=.205
            patches, texts = plt.pie(i.Entry, radius=radio+size, labels = None,
                labeldistance=0.8, rotatelabels =True,
                colors = i[j], # new_colors(valor = len(i.Entry), col = 'nipy_spectral'),
                wedgeprops=dict(width= 0.2, edgecolor='white', alpha = k),
                textprops=dict(size = 10))
            pat.append(patches)

        #plt.legend(pat[0], df_phylum.PHYLUM, loc=2,fontsize=13,labelspacing = 0.4,
        #          bbox_to_anchor=(1.05, 1),frameon=False)

        plt.gca().set(aspect='equal')
        plt.title('Root',fontsize=10,x = 0.5, y = 0.465)
        plt.text(-1.8, 1.35, 'Lineage', fontsize = 15, ha='left', va='center',
                        color='black')
        #plt.title('Lineage',fontsize=20, fontweight='bold', x = -0.17, y = 1)
        #plt.text(1.1, 1.35, linaje_seleccionado, fontsize = 15, ha='left', va='center',
        #                    color='black')
        #>>>>>>>>>>>>>>>>>>>>>>>
        #### insetplot
        #ax2 = plt.axes([0.1, 0.66, 0.13, 0.14])
        ax2 = plt.axes([-0.07, 1.71, 0.17, 0.18])
        logo = [20, 20, 20, 20, 20, 20, 20, 20]
        logo_col = ['white','white','black','white','white','white','white','white']
        logo_col1 = ['white','white','black','black','black','black','black','black']
        radio = 0.5
        linaje = [logo,logo,logo,logo,logo]
        colores = [logo_col1,logo_col,logo_col,logo_col,logo_col]
        name_linaje = ['Phylum', 'Class', 'Order', 'Genus', 'Species']

        pat = []
        size=-.44
        pos = -.18
        for i, j, k, l in zip(linaje, colores, name_linaje, alphas):
            pos += .47
            size+=.44
            ax2.pie(i, radius=radio+size, labels = None,
                colors = j,wedgeprops=dict(width= 0.35, edgecolor='white', alpha = l),
                textprops=dict(size = 10))
            ax2.text(0.1, pos, k, fontsize = 9, ha='left', va='center', fontweight='bold',
                            alpha = l) #color='black'

    def barras_tax(df = DataFrame([]), column = 0, dim = 111, title = '', row_num = 10, color = ['#ff7f0e'],
            size_x = 8, size_y = 10, ylabel_text = '', xlabel = 10, ylabel = 10, size_title = 15, size_bartxt = 10, sep = 1.2):
        if len(df) == 0:
            print('Data frame sin datos')
        else:
            #plt.subplot(dim, facecolor= 'white')
            barWidth = 0.9
            if row_num == len(df):
                ejey = list(df.iloc[0:len(df),1])
                val = max(ejey)
                ejex = list(df.iloc[0:len(df),column])
                colores = list(df.iloc[0:len(df),2])
                borde = list(np.repeat('white', len(df.iloc[0:row_num,column])))
                linea = list(np.repeat(0, len(df.iloc[0:row_num,column])))
            if row_num < len(df):
                ejey = list(df.iloc[0:row_num,1]) + [df.iloc[row_num:len(df),1].sum()]
                val = max(ejey)
                ejex = list(df.iloc[0:row_num,column]) + ['Others']
                borde = list(np.repeat('white', len(df.iloc[0:row_num,column]))) + ['black']
                colores =  list(df.iloc[0:row_num,2]) + ['linen']
                linea = list(np.repeat(0, len(df.iloc[0:row_num,column]))) + [1]
            if row_num > len(df):
                ejey = list(df.iloc[0:len(df),1])
                val = max(ejey)
                ejex = list(df.iloc[0:len(df),column])
                borde = list(np.repeat('white', len(df.iloc[0:row_num,column])))
                colores = list(df.iloc[0:len(df),2])
                linea = list(np.repeat(0, len(df.iloc[0:row_num,column])))

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
            plt.ylabel(ylabel_text,size= ylabel)
            plt.xlabel("Number of Proteins",size= xlabel)
            #plt.tick_params(top = 'on', bottom = 'on', right = 'on', left = 'on')
            #plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

            for j, k in zip(ejey,range(0,len(ejey))):
                plt.text(j+sep,k-0.2,j,
                        size=size_bartxt,ha='left',color= 'black')

    import ipywidgets as widgets
    from ipywidgets import interact, interactive, fixed, interact_manual, Button, HBox, VBox, IntSlider, Label, IntRangeSlider
    from ipywidgets import Checkbox, RadioButtons
    from ipywidgets import Button, Layout
    alfas = {'Lineage*':[1, 1, 1, 1, 1],
            'Phylum':[1, 0.3, 0.3, 0.3, 0.3],
            'Class':[0.3, 1, 0.3, 0.3, 0.3],
            'Order':[0.3, 0.3, 1, 0.3, 0.3],
            'Genus':[0.3, 0.3, 0.3, 1, 0.3],
            'Species':[0.3, 0.3, 0.3, 0.3, 1],
            'Gradient1*':[1, 0.85, 0.7, 0.55, 0.4],
            'Gradient2*':[0.4, 0.55, 0.7, 0.85, 1],
            'Attenuate*':[0.3, 0.3, 0.3, 0.3, 0.3],
            'None*':[0, 0, 0, 0, 0]}
    plotss = ['Phylum', 'Class', 'Order', 'Genus', 'Species']
    posicion_subplots = []
    n = 0.9
    while n < 2:
        n+= 0.1
        posicion_subplots.append(np.around(n,1))

    color_a6 = widgets.Dropdown(options=list(circle_colors.keys()),value='Colors1',description='Colors:',disabled=False,button_style='',
                                layout=Layout(width='20%', height='25px'))
    a6 = widgets.Dropdown(options=list(alfas.keys()),description='Chart 1:', value = 'Phylum',
                        disabled=False,button_style='',
                        layout=Layout(width='20%', height='25px'))
    a61 = widgets.Dropdown(options=plotss,description='Chart 2:',disabled=False,button_style='',
                        layout=Layout(width='20%', height='25px'))
    pos_sub1 = widgets.Dropdown(options=posicion_subplots,value=1.3,description='xloc1:',
                            disabled=False,layout=Layout(width='15%', height='25px'))
    pos_sub2 = widgets.Dropdown(options=posicion_subplots,value=1.3,description='xloc2:',
                            disabled=False,layout=Layout(width='15%', height='25px'))
    b6 = widgets.Dropdown(options=list(range(0,101)),value=10,description='rows1:',disabled=False,
                            layout=Layout(width='15%', height='25px'))
    c6 = widgets.Dropdown(options=list(range(0,101)),value=10,description='rows2:',disabled=False,
                            layout=Layout(width='15%', height='25px'))
    z6 = widgets.ToggleButton(value=False,description='Save Chart',disabled=False,button_style='',
                                tooltip='Description')
    o6 = widgets.Dropdown(options=[0, 0.25, 0.5, 0.75] + list(range(0,201)),value=3,description='sep1:',disabled=False,
                            layout=Layout(width='15%', height='25px'))
    o61 = widgets.Dropdown(options=[0, 0.25, 0.5, 0.75] + list(range(0,201)),value=3,description='sep2:',disabled=False,
                            layout=Layout(width='15%', height='25px'))
    
    d6 = widgets.Dropdown(options=list(range(0,51)),value=8,description='size_y1:',disabled=False,
                        layout=Layout(width='15%', height='25px'))
    d61 = widgets.Dropdown(options=list(range(0,51)),value=8,description='size_y2:',disabled=False,
                        layout=Layout(width='15%', height='25px'))
    g6 = widgets.Dropdown(options=list(range(0,51)),value=8,description='bartxt1:',disabled=False,
                        layout=Layout(width='15%', height='25px'))
    g61 = widgets.Dropdown(options=list(range(0,51)),value=8,description='bartxt2:',disabled=False,
                        layout=Layout(width='15%', height='25px'))

    xxx = Button(layout=Layout(width='5%', height='25px'), disabled=True)
    xxx.style.button_color = 'white'
    yyy = Button(layout=Layout(width='94%', height='5px'), disabled=True)
    yyy.style.button_color = 'red'

    ww = widgets.HBox([color_a6, xxx, z6])
    w6 = widgets.HBox([a6, b6,  o6, d6, g6, pos_sub1,])
    w7 = widgets.HBox([a61, c6,  o61, d61, g61, pos_sub2,])
    w8 = widgets.VBox([w6, w7, yyy])



    ######

    def col(color_a6):
        tax_colors(color_list = circle_colors[color_a6], taxx = dfUniprot)
    out7 = widgets.interactive_output(col, {'color_a6':color_a6})


    def box1(a6, a61, pos_sub1, pos_sub2, b6, c6, z6, o6, o61, d6, d61, g6, g61):
        yetiquetas_plot1 = {'Lineage*':'Phylum',
                        'Phylum':'Phylum',
                        'Class':'Class',
                        'Order':'Order',
                        'Genus':'Genus',
                        'Species':'Species',
                        'Gradient1*':'Phylum',
                        'Gradient2*':'Phylum',
                        'Attenuate*':'Phylum',
                        'None*':'Phylum'}
        plots1 = {'Lineage*':pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA'),
            'Phylum':pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA'),
            'Class':pd.read_csv('tax/class11.tsv', sep = '\t').fillna('NA'),
            'Order':pd.read_csv('tax/order11.tsv', sep = '\t').fillna('NA'),
            'Genus':pd.read_csv('tax/genus11.tsv', sep = '\t').fillna('NA'),
            'Species':pd.read_csv('tax/org11.tsv', sep = '\t').fillna('NA'),
            'Gradient1*':pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA'),
            'Gradient2*':pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA'),
            'Attenuate*':pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA'),
            'None*':pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA')}
        plots2 = {'Phylum':pd.read_csv('tax/phylum0.tsv', sep = '\t').fillna('NA'),
            'Class':pd.read_csv('tax/class11.tsv', sep = '\t').fillna('NA'),
            'Order':pd.read_csv('tax/order11.tsv', sep = '\t').fillna('NA'),
            'Genus':pd.read_csv('tax/genus11.tsv', sep = '\t').fillna('NA'),
            'Species':pd.read_csv('tax/org11.tsv', sep = '\t').fillna('NA')}
        ax3 = plt.axes([pos_sub2, .97, .3, 0.55])
        ##>>>>>>>>>>> grafico circular
        ax = plt.axes([0, 1, 0.9, 1])
        circle_lineage(alphas = alfas[a6])
        ##>>>>>>>>>>> grafico 1
        #ax2 = plt.axes([pos_sub1, 1.51, .3, 0.55])
        ax2 = plt.axes([pos_sub1, 1.63, .3, 0.4])  #>>>>>>>>>>

        barras_tax(plots1[a6],
        #barras_tax(tax_colors(color_list = circle_colors['Spectral'])[0],
                row_num = b6,
                color = plots1[a6].iloc[0:b6,2],
                sep = o6,
                size_y = d6,
                size_bartxt = g6,
                ylabel_text = yetiquetas_plot1[a6],
                ylabel = 10)

        ##>>>>>>>>>>> grafico 2
        ax3 = plt.axes([pos_sub2, .97, .3, 0.55])

        barras_tax(plots2[a61],
        #barras_tax(tax_colors(color_list = circle_colors['Spectral'])[0],
                row_num = c6,
                color = plots2[a61].iloc[0:b6,2],
                sep = o61,
                size_y = d61,
                size_bartxt = g61,
                ylabel_text = yetiquetas_plot1[a61],
                ylabel = 10)

        ##>>>>>>>>>>>> save
        if z6 == True:
            import datetime
            plt.savefig('img/Lineage'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches='tight')
        else:
            pass
    out6 = widgets.interactive_output(box1, {'a6':a6, 'a61':a61, 'pos_sub1':pos_sub1, 'pos_sub2':pos_sub2,
                                            'b6':b6, 'c6':c6, 'z6':z6, 'o6':o6, 'o61':o61,
                                            'd6':d6, 'd61':d61, 'g6':g6, 'g61':g61})
    import warnings
    warnings.filterwarnings("ignore")
    return display(VBox([yyy, ww, w8, out6]))


def inputs():
    taxonomia_cargada = pd.read_csv('tax_file.tsv', sep = '\t').fillna('NA')
    return taxonomia_cargada


import tkinter as tk
from tkinter import filedialog
yyy = Button(layout=Layout(width='94%', height='5px'), disabled=True)
yyy.style.button_color = 'red'
xxx = Button(layout=Layout(width='94%', height='5px'), disabled=True)
xxx.style.button_color = 'white'
upload_tax = widgets.Checkbox(description='Upload file',value=False,disabled=False)
def file_tax(upload_tax):
    if upload_tax == True:
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
            df = pd.read_csv(file_path, sep = '\t')
            df.to_csv('tax_file.tsv', sep = '\t', index = None)
            #df_input(dfUniprot = inputs())
            
salida = widgets.interactive_output(file_tax, {'upload_tax':upload_tax})
in1 = VBox([yyy, xxx, upload_tax, salida, xxx])

run = widgets.Checkbox(description='Show result',value=False,disabled=False)
def file_goslim3(run):
    if run == True:
        df_input(dfUniprot = inputs())
salida3 = widgets.interactive_output(file_goslim3, {'run':run})

in2 = VBox([yyy, xxx, run, salida3, xxx, yyy])
import warnings
warnings.filterwarnings("ignore")
TAX = VBox([in1, in2])


#--------------------------------------------------------------

# extracción de taxonomia desde NCBI


import tkinter as tk
from tkinter import filedialog

# funcion para extraer informacion taxonomica desde NCBI
# columns = ['Tax_ID','kingdom','phylum','class','order','genus']
def Get_Lineage(): # acepta un id
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
        from datetime import datetime 
        inicio_total = datetime.now()
        frame = pd.read_csv(file_path, sep = '\t')
        frame['Tax_ID'] = [str(i) for i in frame.Tax_ID]
        tax_id = frame.iloc[0:len(frame), 2].drop_duplicates().tolist()
        taxa = []
        for i in tax_id:
            ox = requests.get('https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id='+i+'&lvl=3&lin=f&keep=1&srchmode=1&unlock').content.decode()
            txid = ''.join(re.findall('NCBI:txid[0-9]{1,100}',ox))
            txid = re.sub('NCBI:txid','',txid)
            kingdom = ''.join(re.findall('TITLE="kingdom">[\w+ ]{0,100}',ox))
            kingdom = re.sub('TITLE="kingdom">','',kingdom)
            phylum = ''.join(re.findall('TITLE="phylum">[\w+ ]{0,100}',ox))
            phylum = re.sub('TITLE="phylum">','',phylum)
            clase = ''.join(re.findall('TITLE="class">[\w+ ]{0,100}',ox))
            clase = re.sub('TITLE="class">','',clase)
            order = ''.join(re.findall('TITLE="order">[\w+ ]{0,100}',ox))
            order = re.sub('TITLE="order">','',order)
            family = ''.join(re.findall('TITLE="order">[\w+ ]{0,100}',ox))
            family = re.sub('TITLE="order">','',order)
            genus = ''.join(re.findall('TITLE="genus">[\w+ ]{0,100}',ox))
            genus = re.sub('TITLE="genus">','',genus)
            species = ''.join(re.findall('TITLE="species">[\w+ ]{0,100}',ox))
            species = re.sub('TITLE="species">','',species)
            tax = []
            for j in [txid,kingdom,phylum,clase,order,family,genus,species,species]:
                if j:
                    tax.append(j)
                else:
                    j = 'NA'
                    tax.append(j)
            taxa.append(tax)
        lapso_total = datetime.now() - inicio_total
        print('■ 2. Time (hh:mm:ss.ms) {}'.format(lapso_total))
        taxonomy = DataFrame(taxa,columns = ['Tax_ID','KINGDOM','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES','Organism'])
        frame1 = frame.merge(taxonomy, on = 'Tax_ID', how = 'left')
        frame1.to_csv('taxonomy_classification_'+file_path.split('/')[-1], sep = '\t', index = None)
        print('■ 3. Output file: ', 'taxonomy_classification_'+file_path.split('/')[-1])
        return #frame1.head()
#

# funcion para extraer informacion taxonomica desde NCBI
# columns = ['Tax_ID','kingdom','phylum','class','order','genus']
def Get_Taxonomic_Lineage_one(tax_id = ['']): # acepta un id
    if ''.join(tax_id) == '':
        print('No TaxID submit')
    else:
        taxa = []
        for i in tax_id:
            ox = requests.get('https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id='+str(i)+'&lvl=3&lin=f&keep=1&srchmode=1&unlock').content.decode()
            txid = ''.join(re.findall('NCBI:txid[0-9]{1,100}',ox))
            txid = re.sub('NCBI:txid','',txid)
            kingdom = ''.join(re.findall('TITLE="kingdom">[\w+ ]{0,100}',ox))
            kingdom = re.sub('TITLE="kingdom">','',kingdom)
            phylum = ''.join(re.findall('TITLE="phylum">[\w+ ]{0,100}',ox))
            phylum = re.sub('TITLE="phylum">','',phylum)
            clase = ''.join(re.findall('TITLE="class">[\w+ ]{0,100}',ox))
            clase = re.sub('TITLE="class">','',clase)
            order = ''.join(re.findall('TITLE="order">[\w+ ]{0,100}',ox))
            order = re.sub('TITLE="order">','',order)
            family = ''.join(re.findall('TITLE="family">[\w+ ]{0,100}',ox))
            family = re.sub('TITLE="family">','',family)
            genus = ''.join(re.findall('TITLE="genus">[\w+ ]{0,100}',ox))
            genus = re.sub('TITLE="genus">','',genus)
            species = ''.join(re.findall('TITLE="species">[\w+ ]{0,100}',ox))
            species = re.sub('TITLE="species">','',species)
            tax = []
            for j in [txid,kingdom,phylum,clase,order,family,genus,species]:
                if j:
                    tax.append(j)
                else:
                    j = 'NA'
                    tax.append(j)
            taxa.append(tax)
        taxonomy = DataFrame(taxa,columns = ['Tax_ID','KINGDOM','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES'])
        return taxonomy
