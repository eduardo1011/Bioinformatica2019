#!/usr/bin/env python
# coding: utf-8

# # GLOBAL




import os, re, pandas
import json
from functools import reduce
from pandas import DataFrame
import pandas as pd
version = pandas.__version__
if float(re.sub('[.]$', '', version[0:4])) >= 0.25:
    from io import StringIO
elif float(re.sub('[.]$', '', version[0:4])) < 0.25:
    from pandas.compat import StringIO
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, _plot_dendrogram, average


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import matplotlib
import warnings
warnings.filterwarnings("ignore")

import seaborn as sns
from ipywidgets import Button, GridBox, Layout, ButtonStyle, Box, HBox, VBox
 
import ipywidgets as widgets 
from IPython.display import clear_output, display 

import urllib.request

import itertools

os.makedirs('Plots16S', exist_ok = True)



# https://seaborn.pydata.org/tutorial/color_palettes.html
# https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
qualitative_colors = {'Pastel1':9, 'Pastel1_r':9,
                      'Pastel2':8, 'Pastel2_r':8,
                      'Paired':12, 'Paired_r':12,
                      'Accent':8, 'Accent_r':8,
                      'Dark2':8, 'Dark2_r':8,
                      'Set1':9, 'Set1_r':9,
                      'Set2':8, 'Set2_r':8,
                      'Set3':12, 'Set3_r':12,
                      'tab10':10, 'tab10_r':10,
                      'tab20':20, 'tab20_r':20,
                      'tab20b':20, 'tab20b_r':20,
                      'tab20c':20, 'tab20c_r':20}

va = []
n = 0.1
for i in range(100):
    va.append(round(n, 3))
    n += 0.1
va = [0, 0.01, 0.05] + va

mark = ['s',  'o', 'p', '*', '+', 'X', '^', 'v', '<', '>', '8', 'h', 'H', 'D', 'd', 'P']

def barcolor(lista = []):
    mpl.rcParams.update(mpl.rcParamsDefault)
    customPalette = lista
    sns.set_style('white')
    sns.set_palette(customPalette)
    sns.palplot(customPalette)
    plt.gcf().set_size_inches(2, 0.25)
    plt.gca().axis('off')
    plt.show()
def barcolor2(color = ''):
    mpl.rcParams.update(mpl.rcParamsDefault)
    if color in list(qualitative_colors.keys()):
        palette = sns.color_palette(color, qualitative_colors[color])
        sns.palplot(palette)
        plt.gcf().set_size_inches(2, 0.25)
        plt.gca().axis('off')
        plt.show()
    else:
        palette = sns.color_palette(color, 25)
        sns.palplot(palette)
        plt.gcf().set_size_inches(2, 0.25)
        plt.gca().axis('off')
        plt.show()


QUALITATIVE_colors = {}
for i in qualitative_colors:
    QUALITATIVE_colors[i] = [matplotlib.colors.to_hex(i) for i in plt.get_cmap(i)(np.arange(qualitative_colors[i]))]


import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import *


def RampaS(list_color = list()):
    mpl.rcParams.update(mpl.rcParamsDefault)
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))


    root= tk.Tk()
    root.title("Colors")
    root.geometry("430x460")
    root.configure(background='white')
    root.attributes("-topmost", True)

    lab1 = Label(root, text="Python colors", bg = 'white',
                 font=("Arial", 8, "bold"), fg = 'black')


    lencol = len(list_color)

    figure, axes = plt.subplots(figsize=(3.7,4.2), nrows=lencol)
    figure.subplots_adjust(top=0.95, bottom=0.01, left=0.25, right=0.99)


    bar1 = FigureCanvasTkAgg(figure, root)
    bar1.get_tk_widget().grid(row=1, column = 0, columnspan = 4)

    for ax, name in zip(axes, list_color):
        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3]/2.
        figure.text(x_text, y_text, name, va='center', ha='right', fontsize=9)
        ax.axis('off')
        plt.close()

    root.mainloop()



metricas = ['euclidean','braycurtis','canberra','chebyshev','cityblock','correlation','cosine','dice',
            'hamming','jaccard','jensenshannon','kulsinski','mahalanobis',
            'matching','minkowski','rogerstanimoto','russellrao','seuclidean','sokalmichener',
            'sokalsneath','sqeuclidean','yule']
metodos = ['complete','single','average','weighted','centroid','ward']



sampledata = """
Sample	Coffee_Variety	Processing	Cultivation	Time_Dry	Kit	OTA
1R0	Robusta	Dry	Conventional	0	DPS	1.25
2R3	Robusta	Dry	Conventional	3	DPS	0.70
3R8	Robusta	Dry	Conventional	8	DPS	0.74
4R10	Robusta	Dry	Conventional	10	DPS	3.11
5A0	Arabica	Dry	Conventional	0	DPS	0.86
6A6	Arabica	Dry	Conventional	6	DPS	0.70
7A15	Arabica	Dry	Conventional	15	DPS	15.83
8A0G	Arabica	Wet	Organic	0	DPS	3.19
9A3G	Arabica	Wet	Organic	3	DPS	0.64
10A7G	Arabica	Wet	Organic	7	DPS	0.00
11C0	Arabica	Wet	Conventional	0	DPS	0.73
12C3	Arabica	Wet	Conventional	3	DPS	0.00
13C7	Arabica	Wet	Conventional	7	DPS	0.59
140G	Arabica	Wet	Organic	0	DUCM	3.19
153G	Arabica	Wet	Organic	3	DUCM	0.64
167G	Arabica	Wet	Organic	7	DUCM	0.00
17R0	Robusta	Dry	Conventional	0	DUCM	1.25
18R3	Robusta	Dry	Conventional	3	DUCM	0.70
19R8	Robusta	Dry	Conventional	8	DUCM	0.74
20R10	Robusta	Dry	Conventional	10	DUCM	3.11
21A0	Arabica	Dry	Conventional	0	DUCM	0.86
22A6	Arabica	Dry	Conventional	6	DUCM	0.70
23A15	Arabica	Dry	Conventional	15	DUCM	15.83
240C	Arabica	Wet	Conventional	0	DUCM	0.73
253C	Arabica	Wet	Conventional	3	DUCM	0.00
267C	Arabica	Wet	Conventional	7	DUCM	0.59
"""

Sampledata = pd.read_csv(StringIO(sampledata), sep='\t')
Sampledata['Time_Dry'] = [str(i) for i in Sampledata.Time_Dry]
Sampledata['OTA'] = [str(i) for i in Sampledata.OTA]

variables = list(Sampledata.iloc[:, 1:].columns)

def merge_tabla(DF = DataFrame([]),
                varcol=QUALITATIVE_colors['Pastel1'],
                procol=QUALITATIVE_colors['Pastel1_r'],
                culcol=QUALITATIVE_colors['Pastel2'],
                timcol=QUALITATIVE_colors['Pastel2_r'],
                kitcol=QUALITATIVE_colors['Paired'],
                otacol=QUALITATIVE_colors['tab20']):
    
    names1 = DF.merge(Sampledata, on = 'Sample', how = 'left')

    # variedad
    var = names1.Coffee_Variety.unique()
    variedad = dict(zip(var, varcol[0:len(var)]))
    #variedad = {'Arabica':"salmon",'Robusta':"limegreen"}
    
    # Processing
    pro = names1.Processing.unique()
    procesado = dict(zip(pro, procol[0:len(pro)]))
    #procesado = {'Wet':'olive','Dry':'violet'}

    # Cultivation
    cul = names1.Cultivation.unique()
    cultivo = dict(zip(cul, culcol[0:len(cul)]))
    #cultivo = {'Organic':"purple",'Conventional':"orange"}

    # Time_Dry
    tim = names1.Time_Dry.unique()
    tiempo_secado = dict(zip(tim, timcol[0:len(tim)]))
    #tiempo_secado = {'0':Pastel1[0],'3':Pastel1[1],'6':Pastel1[2],'7':Pastel1[3],'8':Pastel1[4],'10':Pastel1[5],'15':Pastel1[6]}
    
    # kit
    kitt = names1.Kit.unique()
    kit = dict(zip(kitt, kitcol[0:len(kitt)]))
    #kit = {'DPS':"cyan",'DUCM':"darkred"}

    otaa = names1.OTA.unique()
    ota = dict(zip(otaa, otacol[0:len(otaa)]))

    return variedad, procesado, cultivo, tiempo_secado, kit, ota, names1


NCBI_RDP_SILVA_summary_ASVs = pd.read_csv('Anexos16S/ASVs_Taxonomy_Counts.tab', sep = '\t')
muestras_ASVs = {}
for i in list(NCBI_RDP_SILVA_summary_ASVs.columns[15:]):
    co = i
    df = NCBI_RDP_SILVA_summary_ASVs[['#OTU ID', co, 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    df = df[df[co] > 0]
    df = df.drop_duplicates()
    df['Per'] = ((np.array(df[co].tolist()) / np.sum(df[co].tolist())) * 100)
    muestras_ASVs[co] = df
NCBI_RDP_SILVA_summary_OTUs = pd.read_csv('Anexos16S/OTUs_Taxonomy_Counts.tab', sep = '\t')
muestras_OTUs = {}
for i in list(NCBI_RDP_SILVA_summary_OTUs.columns[15:]):
    co = i
    df = NCBI_RDP_SILVA_summary_OTUs[['#OTU ID', co, 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    df = df[df[co] > 0]
    df = df.drop_duplicates()
    df['Per'] = ((np.array(df[co].tolist()) / np.sum(df[co].tolist())) * 100)
    muestras_OTUs[co] = df




ordenado4 = list(set([int(re.search('[0-9]{1,2}', i).group()) for i in list(NCBI_RDP_SILVA_summary_ASVs.columns[15:])]))
ordenado3 = {}
for q in ordenado4:
    for w in list(NCBI_RDP_SILVA_summary_ASVs.columns[15:]):
        if re.search('^\d+', w[0:2]): 
            if q == int(re.search('^\d+', w[0:2]).group()):
                ordenado3[q] = w



#ASVs
TAX_DFS_ASVs = {}
for TAX in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
    tax = TAX
    dat0 = []
    for i in list(NCBI_RDP_SILVA_summary_ASVs.columns[15:]):
        co = i
        ww = NCBI_RDP_SILVA_summary_ASVs[[tax, co]][NCBI_RDP_SILVA_summary_ASVs[co] >= 2] # con mas de 2 cuentas
        ee = pd.pivot_table(ww[[tax, co]].drop_duplicates(), values=co, index=[tax], aggfunc=sum).reset_index()
        ff = ee.sort_values(by =co,ascending=True).reset_index(drop=True)
        dat0.append(ff)
    # algunos ASVs se encuentran en varias especies, con esto todos los reads son sumados a una sola especie
    APROBADOS_ASVs = reduce(lambda  left,right: pd.merge(left,right,on=[tax], how='outer'), dat0).fillna(0)
    TAX_DFS_ASVs[TAX] = APROBADOS_ASVs
#OTUs
TAX_DFS_OTUs = {}
for TAX in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
    tax = TAX
    dat0 = []
    for i in list(NCBI_RDP_SILVA_summary_OTUs.columns[15:]):
        co = i
        ww = NCBI_RDP_SILVA_summary_OTUs[[tax, co]][NCBI_RDP_SILVA_summary_OTUs[co] >= 2] # con mas de 2 cuentas
        ee = pd.pivot_table(ww[[tax, co]].drop_duplicates(), values=co, index=[tax], aggfunc=sum).reset_index()
        ff = ee.sort_values(by =co,ascending=True).reset_index(drop=True)
        dat0.append(ff)
    # algunos OTUs se encuentran en varias especies, con esto todos los reads son sumados a una sola especie
    APROBADOS_OTUs = reduce(lambda  left,right: pd.merge(left,right,on=[tax], how='outer'), dat0).fillna(0)
    TAX_DFS_OTUs[TAX] = APROBADOS_OTUs




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



def update_plot(mostrar_dendrograma = True, mostrar_circulos = True, mostrar_metadata = True,
               rampa = 'tab20', etiqueta = '', umbral = '', orientacion = 'HBar',
               png = False, svg = False):
    
    mpl.rcParams.update(mpl.rcParamsDefault)
    ########################## archivos
    sumary111 = []
    uno = open('Anexos16S/SUMARY2.txt', 'r')
    for enu, line in enumerate(uno):
        line = line.rstrip()
        if enu == 0:
            header = line.split('\t')
        else:
            sumary111.append(line.split('\t'))
    uno.close()
    sumary2 = DataFrame(sumary111, columns = header)
    sumary2 = sumary2.set_index('Sample')
    sumary2 = sumary2.astype('float64')
    
    names111 = []
    uno = open('Anexos16S/names1.txt', 'r')
    for enu, line in enumerate(uno):
        line = line.rstrip()
        if enu == 0:
            header = line.split('\t')
        else:
            names111.append(line.split('\t'))
    uno.close()
    names1 = DataFrame(names111, columns = header)
    
    with open('Anexos16S/dend.json', 'r') as fp:
        dendro = json.load(fp)
    with open('Anexos16S/ssss.json', 'r') as fp:
        ssss = json.load(fp)
    with open('Anexos16S/dictionaries.json', 'r') as fp:
        oooo = json.load(fp) 
    
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    unicos = [str(x) for x in sorted([int(i) for i in tiempo_secado])]
    tiempo_secado = dict(zip(unicos, [tiempo_secado[j] for j in unicos]))
    
    kit = oooo[4]
    ota = oooo[5]
    unicos2 = [str(x) for x in sorted([float(i) for i in ota])]
    ota = dict(zip(unicos2, [ota[j] for j in unicos2]))
    
    ancho_barra = 0.8
    
    ##########################


    if orientacion == 'VBar':
        if [mostrar_dendrograma, mostrar_circulos, mostrar_metadata] == [False, False, False]:
            mpl.rcParams.update(mpl.rcParamsDefault)

            fig = plt.figure()

            CUADRO = 0.8
            denx2pos = 0

            ax1 = fig.add_axes([denx2pos+0.005, 0, CUADRO, CUADRO])

            sumary2.plot(kind='bar', stacked=True, figsize=(4+0.5, 4), width= ancho_barra,
                                     color = QUALITATIVE_colors[rampa], alpha = 0.7, ax = ax1)

            plt.xticks(size=8, rotation=90, ha = 'right')
            plt.yticks([25, 50, 75, 100], ['25', '50', '75', '100'], size=7)

            for d in [25, 50, 75]:
                plt.plot(np.array(range(-1, len(sumary2)+1)),np.repeat(d, len(sumary2)+2), 
                         marker='o', markeredgewidth=0, zorder=0, linestyle='-',
                                     markersize=0, color='black', linewidth=0.3, alpha = 1)

            plt.setp(ax1.get_legend().get_title(), fontweight='bold')
            plt.setp(ax1.get_legend().get_texts(), fontsize='8')
            ax1.tick_params(bottom=True, right=False, top=False, left=True, width = 0.3, length=2, color='black')
            ax1.spines['left'].set_linewidth(0.3)
            ax1.spines['bottom'].set_linewidth(0.2)
            ax1.spines['left'].set_color('black')
            ax1.spines['bottom'].set_color('black')
            ax1.spines['right'].set_color(None)
            ax1.spines['top'].set_color(None)
            ax1.spines['top'].set_linewidth(0.2)
            ax1.spines['top'].set_linewidth(0.3)

            ax1.set_ylabel('Relative abundance (%)', fontsize=8)
            ax1.set_xlabel('', fontsize=0)

            ax1.legend(title='Level', title_fontsize = 8,  bbox_to_anchor=(1, 0.98), loc=2, ncol = 1,
                       handletextpad=0.5,
                       fancybox=True, framealpha=0.5, shadow=False,
                       handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                       borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                       prop={'size':7})

            if png == True:
                plt.savefig('Plots16S/'+etiqueta+'_Global_Taxonomic_Level_Stacked_VBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
            if svg == True:
                plt.savefig('Plots16S/'+etiqueta+'_Global_Taxonomic_Level_Stacked_VBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
            plt.show()
        else:
            print('To show the vertical graph you have to disable:\nDendrogram, Pie Plots and Metadata')
            
    if orientacion == 'HBar':
    
        mpl.rcParams.update(mpl.rcParamsDefault)

        fig = plt.figure()

        CUADRO = 0.8
        denx2pos = 0
        #########################################################################
        ################   dendrograma       ########################
        #########################################################################
        if mostrar_dendrograma == True:
            denx2pos = 0.25
            ax0 = fig.add_axes([0, 0.005, denx2pos, CUADRO-0.0105])

            ax0.set_facecolor('none')

            matplotlib.rcParams['lines.linewidth'] = 1
            no_labels=True
            color_list = dendro['color_list']
            Z = np.asarray(dendro['dcoord'], order='c')
            Zs = Z.shape
            N = Zs[0] + 1
            mh = max(Z[:, 2])

            _plot_dendrogram(dendro['icoord'], dendro['dcoord'], dendro['ivl'], 0, N, mh, 'left',
                                 no_labels, color_list, leaf_font_size=None,
                                 leaf_rotation=None, contraction_marks=None,
                                 ax=ax0)

            ax0.axis('off')

        #########################################################################
        ##################     stacked      ###########################
        #########################################################################

        ax1 = fig.add_axes([denx2pos+0.005, 0, CUADRO, CUADRO])

        sumary2.plot(kind='barh', stacked=True, figsize=(4+0.5, 4), width=0.8,
                                 color = QUALITATIVE_colors[rampa], alpha = 0.7, ax = ax1)

        for d in [25, 50, 75]:
            plt.plot(np.repeat(d, len(sumary2)+2),np.array(range(-1, len(sumary2)+1)), 
                     marker='o', markeredgewidth=0, zorder=0, linestyle='-',
                                 markersize=0, color='black', linewidth=0.3, alpha = 1)

        plt.yticks(size=8, rotation=0, ha = 'right', fontweight='bold')
        plt.xticks([25, 50, 75, 100], ['25', '50', '75', '100'], size=7)

        for enu, lab in enumerate(ax1.get_yticklabels()):
            ax1.text(101, enu, lab.get_text(), va = 'center', ha = 'left', fontsize=7)
            #print(enu, lab.get_text())

        plt.setp(ax1.get_legend().get_title(), fontweight='bold')
        plt.setp(ax1.get_legend().get_texts(), fontsize='8')
        plt.gca().tick_params(bottom=True, right=False, top=False, left=False, width = 0.3, length=2, color='black')
        plt.gca().spines['left'].set_linewidth(0.3)
        plt.gca().spines['bottom'].set_linewidth(0.2)
        plt.gca().spines['left'].set_color(None)
        plt.gca().spines['bottom'].set_color('black')
        plt.gca().spines['right'].set_color(None)
        #plt.gca().spines['top'].set_color(None)
        plt.gca().spines['top'].set_linewidth(0.3)
        plt.gca().spines['top'].set_linewidth(0.2)

        ax1.set_ylabel('', fontsize=0)
        ax1.set_xlabel('Relative abundance (%)', fontsize=8)
        ax1.set_yticklabels('', fontsize=0)
        ax1.get_legend().set_visible(False)

        #########################################################################
        ####################     circles       ######################
        #########################################################################
        if mostrar_circulos == True:
            ax2 = fig.add_axes([denx2pos+0.005, CUADRO+0.01, CUADRO, 0.2])
            ax2.set_facecolor('none')

            sample_ordenado0 = list(sumary2.index)
            sample_ordenado0 = list(reversed(sample_ordenado0))
            cuadros = CUADRO/len(sample_ordenado0)

            inicio = denx2pos+0.005
            inicioy = CUADRO+0.01 + (cuadros * 5)
            for u in sample_ordenado0:
                ax222 = plt.axes([inicio, inicioy, cuadros, cuadros])

                ax222.text(0, 0, '   '+u, fontsize=7, color = 'black', ha='center',
                                            va = 'bottom', rotation=90)
                for i, k in zip(ssss[u], [QUALITATIVE_colors[rampa][i] for i, j in enumerate(list(sumary2.columns))]):

                    ax22 = plt.axes([inicio, inicioy, cuadros, cuadros])

                    ax22.pie([ssss[u][i][0], ssss[u][i][1]], #autopct ='%1.1f%%',
                            pctdistance = 0.5, labeldistance= 1,
                            startangle = 0, radius = 1.1, rotatelabels = False,
                            colors = [k,'silver'],
                            wedgeprops={'alpha':1, 'linewidth': 0, 'edgecolor':'black'},
                            explode = np.array([0.05, 0.05]), textprops=dict(size = 0))
                    centre_circle = plt.Circle((0,0),0.2,fc = 'white')

                    plt.gca().add_artist(centre_circle)
                    inicioy -= cuadros

                inicio += cuadros
                inicioy = CUADRO+0.01 + (cuadros * 5)

            inicioy2 = CUADRO-0.05 + (cuadros * 5)
            for ta in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
                ax2.text(0, inicioy2, ta+' ', fontsize=7, color = 'black', ha='right',
                         va = 'center', rotation=0)
                inicioy2 -= 0.17

            ax2.axis('off')
        #########################################################################
        #####################       metadata       ########################
        #########################################################################
        if mostrar_metadata == True:
            ax3 = fig.add_axes([(CUADRO)+(denx2pos)+0.045, -0.015, 0.5, CUADRO+0.0415])
            ax3.set_facecolor('none')
            ax3.set_xlim(0, 65)

            level = dict(zip(list(sumary2.columns), QUALITATIVE_colors[rampa]))

            PUNTO = 50
            inicio = 3
            y_inicio = 0
            agregado = 3.7
            posiciones = []
            for index, row in names1[['Sample', 'Coffee_Variety', 'Processing', 'Cultivation', 'Time_Dry', 'Kit', 'OTA']].drop_duplicates().iterrows():
                ax3.scatter(inicio, y_inicio, s = PUNTO * ancho_barra, c = variedad[row.Coffee_Variety])
                ax3.scatter(inicio + (agregado),y_inicio, s = PUNTO, c = procesado[row.Processing])
                ax3.scatter(inicio + (agregado*2),y_inicio, s = PUNTO, c = cultivo[row.Cultivation])
                ax3.scatter(inicio + (agregado*3),y_inicio, s = PUNTO, c = tiempo_secado[row.Time_Dry])
                ax3.scatter(inicio + (agregado*4),y_inicio, s = PUNTO, c = kit[row.Kit])
                ax3.scatter(inicio + (agregado*5),y_inicio, s = PUNTO, c = ota[row.OTA])
                posiciones.append(y_inicio)
                y_inicio += 1

            agre = [0, 1, 2, 3, 4, 5]
            agregado = 3.7
            for ag, VaR in zip(agre, ['Coffee_Variety', 'Processing', 'Cultivation', 'Time_Dry', 'Kit', 'OTA']):
                ax3.text(inicio + (agregado*ag), y_inicio-0.3, VaR, fontsize=6, color = 'black', ha='center',
                                            va = 'bottom', rotation=90)
            xpos = inicio + (agregado*6.2)
            pos = max(list(reversed(posiciones)))
            for i, u in zip([level, variedad, procesado, cultivo, tiempo_secado],
                            ['Level', 'Coffee_Variety', 'Processing', 'Cultivation', 'Time_Dry']):
                ax3.text(xpos, pos, u, fontsize=6, color = 'black', ha='left',
                                va = 'center', weight='bold')
                pos -= 1
                for j in i:
                    if j in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
                        ax3.scatter(xpos+1, pos, s = 13, c = i[j], marker = 's')
                    else:
                        ax3.scatter(xpos+1, pos, s = 20, c = i[j], marker = 'o')
                    ax3.text(xpos+4, pos, j, fontsize=6, color = 'black', ha='left', va = 'center')
                    pos -= 1

            xpos2 = inicio + (agregado*12.5)
            pos2 = max(list(reversed(posiciones)))
            for i, u in zip([kit, ota], ['Kit', 'OTA']):
                ax3.text(xpos2, pos2, u, fontsize=6, color = 'black', ha='left',
                        va = 'center', weight='bold')
                pos2 -= 1
                for j in i:
                    ax3.scatter(xpos2+1, pos2, s = 20, c = i[j], marker = 'o')
                    ax3.text(xpos2+4, pos2, j, fontsize=6, color = 'black', ha='left', va = 'center')
                    pos2 -= 1

            ax3.axis('off')
        else:
            ax1.legend(title='Level', title_fontsize = 8,  bbox_to_anchor=(1.07, 1.01), loc=2, ncol = 1,
                           handletextpad=0.5,
                           fancybox=True, framealpha=0.5, shadow=False,
                           handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                           borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                           prop={'size':7})
            plt.gca().set_ylabel('', fontsize=0)

        if png == True:
            plt.savefig('Plots16S/'+etiqueta+'_Global_Taxonomic_Level_Stacked_HBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
        if svg == True:
            plt.savefig('Plots16S/'+etiqueta+'_Global_Taxonomic_Level_Stacked_HBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')

        plt.show()



def update_plot2(mostrar_dendrograma = True, mostrar_circulos = True, mostrar_metadata = True,
                        orientacion = 'HBar', rampa = 'tab20', etiqueta = '', umbral = '', NUCOL = 1, tax = '',
                png = False, svg = False):
    
    mpl.rcParams.update(mpl.rcParamsDefault)
    ########################## archivos
    sumary111 = []
    uno = open('Anexos16S/AAbundancEE.txt', 'r')
    for enu, line in enumerate(uno):
        line = line.rstrip()
        if enu == 0:
            header = line.split('\t')
        else:
            sumary111.append(line.split('\t'))
    uno.close()
    sumary2 = DataFrame(sumary111, columns = header)
    sumary2 = sumary2.set_index('Sample')
    sumary2 = sumary2.astype('float64')
    
    
    names111 = []
    uno = open('Anexos16S/names111.txt', 'r')
    for enu, line in enumerate(uno):
        line = line.rstrip()
        if enu == 0:
            header = line.split('\t')
        else:
            names111.append(line.split('\t'))
    uno.close()
    names111 = DataFrame(names111, columns = header)
    
    with open('Anexos16S/dendo.json', 'r') as fp:
        dendro = json.load(fp)
    with open('Anexos16S/numeros2.json', 'r') as fp:
        numeros2 = json.load(fp)
    with open('Anexos16S/dictionaries2.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    unicos = [str(x) for x in sorted([int(i) for i in tiempo_secado])]
    tiempo_secado = dict(zip(unicos, [tiempo_secado[j] for j in unicos]))
    
    kit = oooo[4]
    ota = oooo[5]
    unicos2 = [str(x) for x in sorted([float(i) for i in ota])]
    ota = dict(zip(unicos2, [ota[j] for j in unicos2]))
    
    ancho_barra = 0.8
    
    ##########################
    if orientacion == 'VBar':
        if [mostrar_dendrograma, mostrar_circulos, mostrar_metadata] == [False, False, False]:
            mpl.rcParams.update(mpl.rcParamsDefault)

            fig = plt.figure()

            CUADRO = 0.8
            denx2pos = 0

            ax1 = fig.add_axes([denx2pos+0.005, 0, CUADRO, CUADRO])

            sumary2.plot(kind='bar', stacked=True, figsize=(4+0.5, 4), width= ancho_barra,
                                     colormap = plt.get_cmap(rampa), alpha = 0.7, ax = ax1)

            plt.xticks(size=8, rotation=90, ha = 'right')
            plt.yticks([25, 50, 75, 100], ['25', '50', '75', '100'], size=7)

            for d in [25, 50, 75]:
                plt.plot(np.array(range(-1, len(sumary2)+1)),np.repeat(d, len(sumary2)+2), 
                         marker='o', markeredgewidth=0, zorder=0, linestyle='-',
                                     markersize=0, color='black', linewidth=0.3, alpha = 1)

            plt.setp(ax1.get_legend().get_title(), fontweight='bold')
            plt.setp(ax1.get_legend().get_texts(), fontsize='8')
            ax1.tick_params(bottom=True, right=False, top=False, left=True, width = 0.3, length=2, color='black')
            ax1.spines['left'].set_linewidth(0.3)
            ax1.spines['bottom'].set_linewidth(0.2)
            ax1.spines['left'].set_color('black')
            ax1.spines['bottom'].set_color('black')
            ax1.spines['right'].set_color(None)
            ax1.spines['top'].set_color(None)
            ax1.spines['top'].set_linewidth(0.2)
            ax1.spines['top'].set_linewidth(0.3)

            ax1.set_ylabel('Relative abundance (%)', fontsize=8)
            ax1.set_xlabel('', fontsize=0)

            ax1.legend(title=tax, title_fontsize = 6,  bbox_to_anchor=(1, 0.97), loc=2, ncol = NUCOL,
                       handletextpad=0.5,
                       fancybox=True, framealpha=0.5, shadow=False,
                       handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                       borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                       prop={'style':'italic', 'size':5})
            
            if png == True:
                plt.savefig('Plots16S/'+etiqueta+'_Taxonomic_Level_'+tax+'_Stacked_VBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
            if svg == True:
                plt.savefig('Plots16S/'+etiqueta+'_Taxonomic_Level_'+tax+'_Stacked_VBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
            plt.show()
        else:
            print('To show the vertical graph you have to disable:\nDendrogram, Pie Plots and Metadata')
            
            
    if orientacion == 'HBar':
        
        mpl.rcParams.update(mpl.rcParamsDefault)

        fig = plt.figure()

        CUADRO = 0.8
        denx2pos = 0
        #########################################################################
        ################   dendrograma       ########################
        #########################################################################
        if mostrar_dendrograma == True:
            denx2pos = 0.25
            ax0 = fig.add_axes([0, 0.005, denx2pos, CUADRO-0.0105])

            ax0.set_facecolor('none')

            matplotlib.rcParams['lines.linewidth'] = 1
            no_labels=True
            color_list = dendro['color_list']
            Z = np.asarray(dendro['dcoord'], order='c')
            Zs = Z.shape
            N = Zs[0] + 1
            mh = max(Z[:, 2])

            _plot_dendrogram(dendro['icoord'], dendro['dcoord'], dendro['ivl'], 0, N, mh, 'left',
                                 no_labels, color_list, leaf_font_size=None,
                                 leaf_rotation=None, contraction_marks=None,
                                 ax=ax0)

            ax0.axis('off')

        #########################################################################
        ##################     stacked      ###########################
        #########################################################################

        ax1 = fig.add_axes([denx2pos+0.005, 0, CUADRO, CUADRO])

        sumary2.plot(kind='barh', stacked=True, figsize=(4+0.5, 4), width=0.8,
                                 colormap = plt.get_cmap(rampa), alpha = 0.7, ax = ax1)

        for d in [25, 50, 75]:
            plt.plot(np.repeat(d, len(sumary2)+2),np.array(range(-1, len(sumary2)+1)), 
                     marker='o', markeredgewidth=0, zorder=0, linestyle='-',
                                 markersize=0, color='black', linewidth=0.3, alpha = 1)

        plt.yticks(size=8, rotation=0, ha = 'right', fontweight='bold')
        plt.xticks([25, 50, 75, 100], ['25', '50', '75', '100'], size=7)

        for enu, lab in enumerate(ax1.get_yticklabels()):
            ax1.text(101, enu, lab.get_text(), va = 'center', ha = 'left', fontsize=7)
            #print(enu, lab.get_text())

        plt.setp(ax1.get_legend().get_title(), fontweight='bold')
        plt.setp(ax1.get_legend().get_texts(), fontsize='8')
        plt.gca().tick_params(bottom=True, right=False, top=False, left=False, width = 0.3, length=2, color='black')
        plt.gca().spines['left'].set_linewidth(0.3)
        plt.gca().spines['bottom'].set_linewidth(0.2)
        plt.gca().spines['left'].set_color(None)
        plt.gca().spines['bottom'].set_color('black')
        plt.gca().spines['right'].set_color(None)
        #plt.gca().spines['top'].set_color(None)
        plt.gca().spines['top'].set_linewidth(0.3)
        plt.gca().spines['top'].set_linewidth(0.2)

        ax1.set_ylabel('', fontsize=0)
        ax1.set_xlabel('Relative abundance (%)', fontsize=8)
        ax1.set_yticklabels('', fontsize=0)
        ax1.legend(title= tax, title_fontsize = 6,  bbox_to_anchor=(1.66, 1.02), loc = 2, ncol = NUCOL,
                           handletextpad=0.5,
                           fancybox=True, framealpha=0.5, shadow=False,
                           handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                           borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                           prop={'style':'italic', 'size':5})
        
        #########################################################################
        ####################     circles       ######################
        #########################################################################
        if mostrar_circulos == True:
            ax2 = fig.add_axes([denx2pos+0.005, CUADRO+0.01, CUADRO, 0.2])
            ax2.set_facecolor('none')

            sample_ordenado0 = list(sumary2.index)
            sample_ordenado0 = list(reversed(sample_ordenado0))
            cuadros = CUADRO/len(sample_ordenado0)
            
            inicio = denx2pos+0.005
            inicioy = CUADRO+0.01
            for u in sample_ordenado0:
                #ax222 = plt.axes([inicio, inicioy, cuadros, cuadros])

                #ax222.text(0, 0, '   '+u, fontsize=7, color = 'black', ha='center', va = 'bottom', rotation=90)
                
                ax22 = plt.axes([inicio, inicioy, cuadros, cuadros])
                
                ax22.pie([numeros2[u][0], numeros2[u][1]], #autopct ='%1.1f%%',
                            pctdistance = 0.5, labeldistance= 1,
                            startangle = 0, radius = 1.1, rotatelabels = False,
                            colors = ['black','silver'],
                            wedgeprops={'alpha':1, 'linewidth': 0, 'edgecolor':'black'},
                            explode = np.array([0.05, 0.05]), textprops=dict(size = 0))
                centre_circle = plt.Circle((0,0),0.2,fc = 'white')

                plt.gca().add_artist(centre_circle)
                
                ax22.text(0.2, 0, '   '+u, fontsize=7, color = 'black', ha='center', va = 'bottom', rotation=90)
                
                inicio += cuadros
                inicioy = CUADRO+0.01


            ax2.axis('off')
            
        #########################################################################
        #####################       metadata       ########################
        #########################################################################
        if mostrar_metadata == True:
            ax3 = fig.add_axes([(CUADRO)+(denx2pos)+0.045, -0.015, 0.5, CUADRO+0.0415])
            ax3.set_facecolor('none')
            ax3.set_xlim(0, 65)

            PUNTO = 50
            inicio = 3
            y_inicio = 0
            agregado = 3.7
            posiciones = []
            for index, row in names111[['Sample', 'Coffee_Variety', 'Processing', 'Cultivation', 'Time_Dry', 'Kit', 'OTA']].drop_duplicates().iterrows():
                ax3.scatter(inicio, y_inicio, s = PUNTO * ancho_barra, c = variedad[row.Coffee_Variety])
                ax3.scatter(inicio + (agregado),y_inicio, s = PUNTO, c = procesado[row.Processing])
                ax3.scatter(inicio + (agregado*2),y_inicio, s = PUNTO, c = cultivo[row.Cultivation])
                ax3.scatter(inicio + (agregado*3),y_inicio, s = PUNTO, c = tiempo_secado[row.Time_Dry])
                ax3.scatter(inicio + (agregado*4),y_inicio, s = PUNTO, c = kit[row.Kit])
                ax3.scatter(inicio + (agregado*5),y_inicio, s = PUNTO, c = ota[row.OTA])
                posiciones.append(y_inicio)
                y_inicio += 1

            agre = [0, 1, 2, 3, 4, 5]
            agregado = 3.7
            for ag, VaR in zip(agre, ['Coffee_Variety', 'Processing', 'Cultivation', 'Time_Dry', 'Kit', 'OTA']):
                ax3.text(inicio + (agregado*ag), y_inicio-0.3, VaR, fontsize=6, color = 'black', ha='center',
                                            va = 'bottom', rotation=90)
            xpos = inicio + (agregado*6.2)
            pos = max(list(reversed(posiciones)))
            for i, u in zip([variedad, procesado, cultivo, tiempo_secado],
                            ['Coffee_Variety', 'Processing', 'Cultivation', 'Time_Dry']):
                ax3.text(xpos, pos, u, fontsize=6, color = 'black', ha='left',
                                va = 'center', weight='bold')
                pos -= 1
                for j in i:
                    if j in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
                        ax3.scatter(xpos+1, pos, s = 13, c = i[j], marker = 's')
                    else:
                        ax3.scatter(xpos+1, pos, s = 20, c = i[j], marker = 'o')
                    ax3.text(xpos+4, pos, j, fontsize=6, color = 'black', ha='left', va = 'center')
                    pos -= 1

            xpos2 = inicio + (agregado*12.5)
            pos2 = max(list(reversed(posiciones)))
            for i, u in zip([kit, ota], ['Kit', 'OTA']):
                ax3.text(xpos2, pos2, u, fontsize=6, color = 'black', ha='left',
                        va = 'center', weight='bold')
                pos2 -= 1
                for j in i:
                    ax3.scatter(xpos2+1, pos2, s = 20, c = i[j], marker = 'o')
                    ax3.text(xpos2+4, pos2, j, fontsize=6, color = 'black', ha='left', va = 'center')
                    pos2 -= 1

            ax3.axis('off')
        else:
            ax1.legend(title= tax, title_fontsize = 6,  bbox_to_anchor=(1.05, 1.02), loc = 2, ncol = NUCOL,
                           handletextpad=0.5,
                           fancybox=True, framealpha=0.5, shadow=False,
                           handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                           borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                           prop={'style':'italic', 'size':5})
            plt.gca().set_ylabel('', fontsize=0)

        if png == True:
            plt.savefig('Plots16S/'+etiqueta+'_Taxonomic_Level_'+tax+'_Stacked_HBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
        if svg == True:
            plt.savefig('Plots16S/'+etiqueta+'_Taxonomic_Level_'+tax+'_Stacked_HBar_'+umbral+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')

        plt.show()



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cor_size = dict(zip(([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), ([4, 4.2, 4.4, 4.6, 4.8, 5, 5.2, 5.4, 5.6, 5.8])))
cor_size


data = widgets.ToggleButtons(options=['ASV', 'OTU'])
data.layout.width = '10%'
#data.style.font_weight = 'bold'
data.style.button_width = '60px'
boton_data = Box(children=[data], layout= Layout(border='1px solid pink', width='69px', height='63px'))

show_dendrograma = widgets.ToggleButtons(options=[True, False])
show_dendrograma.layout.width = '10%'
show_dendrograma.style.button_width = '60px'
boton_dendrogram = Box(children=[show_dendrograma], layout= Layout(border='1px solid pink', width='69px', height='63px'))

show_pies = widgets.ToggleButtons(options=[True, False])
show_pies.layout.width = '10%'
show_pies.style.button_width = '60px'
boton_pies = Box(children=[show_pies], layout= Layout(border='1px solid pink', width='69px', height='63px'))

show_anotaciones = widgets.ToggleButtons(options=[True, False])
show_anotaciones.layout.width = '10%'
show_anotaciones.style.button_width = '60px'
boton_metadata = Box(children=[show_anotaciones], layout= Layout(border='1px solid pink', width='69px', height='63px'))

orientation = widgets.ToggleButtons(options=['HBar', 'VBar'])
orientation.layout.width = '10%'
orientation.style.button_width = '60px'
hor_ver = Box(children=[orientation], layout= Layout(border='1px solid pink', width='69px', height='63px'))


Percentage =widgets.SelectionSlider(options=va,value=0.01,description='Limit (%):',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True,
                                   layout=Layout(width='95%', height='25px'))
Method = widgets.Dropdown(options=metodos,value='complete',disabled=False,
                         layout=Layout(width='88px', height='25px'))
Metric = widgets.Dropdown(options=metricas,value='euclidean',disabled=False,
                         layout=Layout(width='88px', height='25px'))
############


Rampas = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Paired',
                          disabled=False,
                         layout=Layout(width='98px', height='25px'))

tamano_plot = widgets.SelectionSlider(options=list(cor_size.keys()),value=4,disabled=False,
                                      description = 'Size Plot:',
                                continuous_update=False,orientation='horizontal',readout=True)

################
negra = widgets.Button(layout=Layout(width='977px', height='3px'), disabled=True)
negra.style.button_color = 'black'
blanca = widgets.Button(layout=Layout(width='977px', height='3px'), disabled=True)
blanca.style.button_color = 'white'


#----------------------------------------------
button = widgets.Button(description=" UPDATE PLOT ", icon = 'fa-refresh')
button.style.button_color = 'lime' #'deepskyblue'
button.style.font_weight = 'bold'
output = widgets.Output()

def button_clicked(b):
    with output:
        clear_output(True)
        update_plot(mostrar_dendrograma = show_dendrograma.value, mostrar_circulos = show_pies.value,
                   mostrar_metadata = show_anotaciones.value, rampa = Rampas.value,
                   png = False, svg = False, etiqueta = '', umbral = '', orientacion = orientation.value)

button.on_click(button_clicked)
#----------------------------------------------
boton4 = widgets.Button(description=" SEE COLORS ", icon = 'fa-eye')
boton4.style.button_color = 'silver'
output4 = widgets.Output()

def button_clicked4(b):
    RampaS(list_color = list(QUALITATIVE_colors.keys()))
    
boton4.on_click(button_clicked4)
#----------------------------------------------
boton3 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
boton3.style.button_color = 'gold'
output3 = widgets.Output()

def button_clicked3(b):
    with output3:
        clear_output(True)
        update_plot(mostrar_dendrograma = show_dendrograma.value, mostrar_circulos = show_pies.value,
                   mostrar_metadata = show_anotaciones.value, rampa = Rampas.value,
                   png = True, svg = False, etiqueta = data.value, umbral = str(Percentage.value),
                    orientacion = orientation.value)


        
boton3.on_click(button_clicked3)
#----------------------------------------------
boton5 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
boton5.style.button_color = 'gold'
output5 = widgets.Output()

def button_clicked5(b):
    with output5:
        clear_output(True)
        update_plot(mostrar_dendrograma = show_dendrograma.value, mostrar_circulos = show_pies.value,
                   mostrar_metadata = show_anotaciones.value, rampa = Rampas.value,
                   png = False, svg = True, etiqueta = data.value, umbral = str(Percentage.value),
                    orientacion = orientation.value)
        

        

boton5.on_click(button_clicked5)

altos = []
for i in QUALITATIVE_colors:
    if len(QUALITATIVE_colors[i]) > 11:
        altos.append(i)

CCOOFF = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Pastel1',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
PPRROO = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Pastel2',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
CCUULL = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Paired',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
TTIIMM = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Accent',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
KKIITT = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Dark2',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
OOTTAA = widgets.Dropdown(options=altos,value='tab20',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))



def box1(data, Percentage, Method, Metric, Rampas,
         CCOOFF,PPRROO,CCUULL,TTIIMM,KKIITT,OOTTAA):
    barcolor(lista = QUALITATIVE_colors[Rampas])
    if data == 'ASV':
        muestras = muestras_ASVs
    elif data == 'OTU':
        muestras = muestras_OTUs

    umbral_porcentaje = float(Percentage)   
    ssss = {}
    sumary = []
    for co in muestras:
        df = muestras[co]
        df_up = df[df.Per > umbral_porcentaje]
        df_down = df[df.Per <= umbral_porcentaje]
        sm = []
        numeros0 = {}
        for tax in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
            df2 = df_up[df_up.Per > umbral_porcentaje][tax].drop_duplicates()
            ll = [len(df_up[tax].drop_duplicates()), len(df_down[tax].drop_duplicates())]
            numeros0[tax] = ll
            sm.append(len(df2))
        sm2 = list(np.array(sm) / np.sum(sm) * 100)
        sm2.insert(0, co)
        sumary.append(sm2)
        ssss[co] = numeros0

    SUMARY = DataFrame(sumary, columns = ['Sample','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])

    # metodo y metrica
    dend = dendrogram(linkage(SUMARY.iloc[:, 1:].values, Method, metric=Metric),
                orientation='left',  no_plot = True,
                labels=SUMARY.Sample.tolist(),
                distance_sort='descending',
                show_leaf_counts=True)
    SUMARY1 = DataFrame(dend['ivl'], columns = ['Sample']).merge(SUMARY, on = 'Sample', how = 'left')

    SUMARY2 = SUMARY1.set_index('Sample')

    names1 = DataFrame(dend['ivl'], columns = ['Sample'])
    variedad, procesado, cultivo, tiempo_secado, kit, ota, names1 = merge_tabla(DF = names1,
                                                                           varcol=QUALITATIVE_colors[CCOOFF],
                                                                            procol=QUALITATIVE_colors[PPRROO],
                                                                            culcol=QUALITATIVE_colors[CCUULL],
                                                                            timcol=QUALITATIVE_colors[TTIIMM],
                                                                            kitcol=QUALITATIVE_colors[KKIITT],
                                                                            otacol=QUALITATIVE_colors[OOTTAA])
    SUMARY2.to_csv('Anexos16S/SUMARY2.txt', sep = '\t')
    names1.to_csv('Anexos16S/names1.txt', sep = '\t', index = None)
    with open('Anexos16S/dictionaries.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
    with open('Anexos16S/ssss.json', 'w') as fp:
        json.dump(ssss, fp)
    with open('Anexos16S/dend.json', 'w') as fp:
        json.dump(dend, fp)

OUT1 = widgets.interactive_output(box1, {'data':data, 'Percentage':Percentage, 'Method':Method,
                                         'Metric':Metric, 'Rampas':Rampas, 'CCOOFF':CCOOFF,
                                        'PPRROO':PPRROO,'CCUULL':CCUULL,'TTIIMM':TTIIMM,'KKIITT':KKIITT,
                                        'OOTTAA':OOTTAA})


gris = widgets.Button(layout=Layout(width='0%', height='20px'), disabled=True)
gris.style.button_color = 'white'


def box2(CCOOFF):
    barcolor(lista = QUALITATIVE_colors[CCOOFF])
    with open('Anexos16S/dictionaries.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    variedad = dict(zip(list(variedad.keys()), QUALITATIVE_colors[CCOOFF][0:len(variedad)]))
    with open('Anexos16S/dictionaries.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT2 = widgets.interactive_output(box2, {'CCOOFF':CCOOFF})

def box3(PPRROO):
    barcolor(lista = QUALITATIVE_colors[PPRROO])
    with open('Anexos16S/dictionaries.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    procesado = dict(zip(list(procesado.keys()), QUALITATIVE_colors[PPRROO][0:len(procesado)]))
    with open('Anexos16S/dictionaries.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT3 = widgets.interactive_output(box3, {'PPRROO':PPRROO})

def box4(CCUULL):
    barcolor(lista = QUALITATIVE_colors[CCUULL])
    with open('Anexos16S/dictionaries.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    cultivo = dict(zip(list(cultivo.keys()), QUALITATIVE_colors[CCUULL][0:len(cultivo)]))
    with open('Anexos16S/dictionaries.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT4 = widgets.interactive_output(box4, {'CCUULL':CCUULL})

def box5(TTIIMM):
    barcolor(lista = QUALITATIVE_colors[TTIIMM])
    with open('Anexos16S/dictionaries.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    tiempo_secado = dict(zip(list(tiempo_secado.keys()), QUALITATIVE_colors[TTIIMM][0:len(tiempo_secado)]))
    with open('Anexos16S/dictionaries.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT5 = widgets.interactive_output(box5, {'TTIIMM':TTIIMM})

def box6(KKIITT):
    barcolor(lista = QUALITATIVE_colors[KKIITT])
    with open('Anexos16S/dictionaries.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    kit = dict(zip(list(kit.keys()), QUALITATIVE_colors[KKIITT][0:len(kit)]))
    with open('Anexos16S/dictionaries.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT6 = widgets.interactive_output(box6, {'KKIITT':KKIITT})

def box7(OOTTAA):
    barcolor(lista = QUALITATIVE_colors[OOTTAA])
    with open('Anexos16S/dictionaries.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    ota = dict(zip(list(ota.keys()), QUALITATIVE_colors[OOTTAA][0:len(ota)]))
    with open('Anexos16S/dictionaries.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT7 = widgets.interactive_output(box7, {'OOTTAA':OOTTAA})


box_layout = Layout(flex_flow='column',
                    align_items='stretch',
                    border='1px solid silver',
                    width='170px',
                    height='467px',
                    display='flex')
carousel = Box(children=[HBox([widgets.Label('Variable colors')]),
             HBox([widgets.Label('Coffee_Variety:'), CCOOFF]),
             OUT2,
             HBox([widgets.Label('Processing:'), PPRROO]),
             OUT3,
             HBox([widgets.Label('Cultivation:'), CCUULL]),
             OUT4,
             HBox([widgets.Label('Time_Dry:'), TTIIMM]),
             OUT5,
             HBox([widgets.Label('Kit:'), KKIITT]),
             OUT6,
             HBox([widgets.Label('OTA:'), OOTTAA]),
             OUT7
             ], layout=box_layout)


box_layout1 = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='1px solid lime',
                    width='780px',
                   height='500px')


GLOBAL_TAX = VBox([
    #negra,
      Box([HBox([VBox([widgets.Label('Clustering:'), boton_data]),
                 VBox([widgets.Label('Method/Metric:'), Method, Metric]),
                 VBox([widgets.Label('Dendrog.:'), boton_dendrogram]),
                 VBox([widgets.Label('Circles:'), boton_pies]), 
                 VBox([widgets.Label('Metadata:'), boton_metadata]),
                 VBox([widgets.Label('Orientarion:'), hor_ver]),
                 VBox([widgets.Label('Stacked Color:'), VBox([Rampas, OUT1])])                
                ])],
          
      layout = Layout(border='1px solid silver', width='950px')),
      Box(children = [Percentage], layout = Layout(border='1px solid silver', width='950px')), #negra,
     Box(children = [HBox([button, gris, gris, gris, gris, boton3, boton5])],  layout = Layout(border='1px solid silver', width='950px')),
    HBox([Box(children=[output], layout=box_layout1), VBox([carousel, boton4])])
     ])



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# # individual

INDIVIDUAL_COLORS = ['BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Blues',
                     'BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',
                     'Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Pastel1','Pastel2',
                     'Paired','Accent','Dark2','Set1','Set2','Set3','tab10','tab20','tab20b','tab20c']
INDIVIDUAL_COLORS = list(set(INDIVIDUAL_COLORS + list(QUALITATIVE_colors.keys())))
INDIVIDUAL_COLORS = dict(zip(INDIVIDUAL_COLORS, INDIVIDUAL_COLORS))


Rampas2 = widgets.Dropdown(options=list(INDIVIDUAL_COLORS.keys()),value='tab20',
                          disabled=False,
                         layout=Layout(width='98px', height='25px'))


NumcoL =widgets.SelectionSlider(options=[1,2,3,4,5,6,7,8,9,10],value=1,description='Ncol:',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True)


Percentage2 =widgets.SelectionSlider(options=va,value=1,description='Limit (%):',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True,
                                   layout=Layout(width='95%', height='25px'))


data_i = widgets.ToggleButtons(options=['ASV', 'OTU'])
data_i.layout.width = '10%'
#data.style.font_weight = 'bold'
data_i.style.button_width = '60px'
boton_data_i = Box(children=[data_i], layout= Layout(border='1px solid pink', width='69px', height='63px'))

show_dendrograma_i = widgets.ToggleButtons(options=[True, False])
show_dendrograma_i.layout.width = '10%'
show_dendrograma_i.style.button_width = '60px'
boton_dendrogram_i = Box(children=[show_dendrograma_i], layout= Layout(border='1px solid pink', width='69px', height='63px'))

show_pies_i = widgets.ToggleButtons(options=[True, False])
show_pies_i.layout.width = '10%'
show_pies_i.style.button_width = '60px'
boton_pies_i = Box(children=[show_pies_i], layout= Layout(border='1px solid pink', width='69px', height='63px'))

show_anotaciones_i = widgets.ToggleButtons(options=[True, False])
show_anotaciones_i.layout.width = '10%'
show_anotaciones_i.style.button_width = '60px'
boton_metadata_i = Box(children=[show_anotaciones_i], layout= Layout(border='1px solid pink', width='69px', height='63px'))

orientation_i = widgets.ToggleButtons(options=['HBar', 'VBar'])
orientation_i.layout.width = '10%'
orientation_i.style.button_width = '60px'
hor_ver_i = Box(children=[orientation_i], layout= Layout(border='1px solid pink', width='69px', height='63px'))


Method_i = widgets.Dropdown(options=metodos,value='complete',disabled=False,
                         layout=Layout(width='88px', height='25px'))
Metric_i = widgets.Dropdown(options=metricas,value='euclidean',disabled=False,
                         layout=Layout(width='88px', height='25px'))
############

tamano_plot_i = widgets.SelectionSlider(options=list(cor_size.keys()),value=4,disabled=False,
                                      description = 'Size Plot:',
                                continuous_update=False,orientation='horizontal',readout=True)




CCOOFF_i = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Pastel1',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
PPRROO_i = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Pastel2',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
CCUULL_i = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Paired',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
TTIIMM_i = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Accent',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
KKIITT_i = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Dark2',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))
OOTTAA_i = widgets.Dropdown(options=altos,value='tab20',
                          disabled=False,
                         layout=Layout(width='78px', height='25px'))



Linaje = widgets.Dropdown(options=['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'],value='Order',
                          disabled=False,
                         layout=Layout(width='98px', height='25px'))



def box11(data_i, Linaje, Percentage2, Method_i, Metric_i,
         CCOOFF_i,PPRROO_i,CCUULL_i,TTIIMM_i,KKIITT_i,OOTTAA_i):
    tax = Linaje
    
    if data_i == 'ASV':
        TAX_DFS = TAX_DFS_ASVs
    elif data_i == 'OTU':
        TAX_DFS = TAX_DFS_OTUs

    umbral_porcentaje = float(Percentage2)
    
    numeros2 = {}
    abundancias = []
    for i in list(TAX_DFS[Linaje].iloc[:, 1:].columns):
        co = i
        df = TAX_DFS[Linaje][[tax, co]]
        df = df.sort_values(by =co,ascending=False).reset_index(drop=True)
        df = df[df[co] > 0]
        suma = df[co].sum()
        df[co] = (df[co] / suma) * 100
        ff = df[df[co] > umbral_porcentaje] # umbral %
        kf = df[df[co] <= umbral_porcentaje]
        numeros2[i] = [len(ff), len(kf)]
        kf = DataFrame({tax:['Others'], co:[kf[co].sum()]})
        hf = pd.concat([ff, kf])
        abundancias.append(hf)
    
    Abundance = reduce(lambda  left,right: pd.merge(left,right,on=[tax], how='outer'), abundancias).fillna(0)
    Abundance = Abundance.rename(columns={tax: 'Sample'})
    Abundance = Abundance.set_index('Sample')
    Abundance = Abundance.T
    Abundance.insert(loc = 0, column='Sample', value=Abundance.index)
    
    # metodo y metrica
    dendo = dendrogram(linkage(Abundance.iloc[:, 1:].values, Method_i, metric=Metric_i),
                orientation='left',  no_plot = True,
                labels=Abundance.Sample.tolist(),
                distance_sort='descending',
                show_leaf_counts=True)
    AAbundancEE = DataFrame(dendo['ivl'], columns = ['Sample']).merge(Abundance, on = 'Sample', how = 'left')
    AAbundancEE = AAbundancEE.set_index('Sample')
    
    names11 = DataFrame(dendo['ivl'], columns = ['Sample'])
    variedad_i, procesado_i, cultivo_i, tiempo_secado_i, kit_i, ota_i, names111 = merge_tabla(DF = names11,
                                                                           varcol=QUALITATIVE_colors[CCOOFF_i],
                                                                            procol=QUALITATIVE_colors[PPRROO_i],
                                                                            culcol=QUALITATIVE_colors[CCUULL_i],
                                                                            timcol=QUALITATIVE_colors[TTIIMM_i],
                                                                            kitcol=QUALITATIVE_colors[KKIITT_i],
                                                                            otacol=QUALITATIVE_colors[OOTTAA_i])
    
    print('\33[1m\33[7m'+tax+': '+str(len(AAbundancEE.columns))+'\33[0m')
    
    AAbundancEE.to_csv('Anexos16S/AAbundancEE.txt', sep = '\t')
    names111.to_csv('Anexos16S/names111.txt', sep = '\t', index = None)
    with open('Anexos16S/dictionaries2.json', 'w') as fp:
        json.dump([variedad_i, procesado_i, cultivo_i, tiempo_secado_i, kit_i, ota_i], fp)
    with open('Anexos16S/numeros2.json', 'w') as fp:
        json.dump(numeros2, fp)
    with open('Anexos16S/dendo.json', 'w') as fp:
        json.dump(dendo, fp)
    

    
OUT11 = widgets.interactive_output(box11, {'data_i':data_i, 'Linaje':Linaje, 'Percentage2':Percentage2, 'Method_i':Method_i,
                                         'Metric_i':Metric_i, 'CCOOFF_i':CCOOFF_i,
                                        'PPRROO_i':PPRROO_i,'CCUULL_i':CCUULL_i,'TTIIMM_i':TTIIMM_i,
                                        'KKIITT_i':KKIITT_i,'OOTTAA_i':OOTTAA_i})




#----------------------------------------------
buttonnn = widgets.Button(description=" UPDATE PLOT ", icon = 'fa-refresh')
buttonnn.style.button_color = 'lime' #'deepskyblue'
buttonnn.style.font_weight = 'bold'
outputtt = widgets.Output()

def button_clickeddd(b):
    with outputtt:
        clear_output(True)
        update_plot2(mostrar_dendrograma = show_dendrograma_i.value, mostrar_circulos = show_pies_i.value,
                     mostrar_metadata = show_anotaciones_i.value, orientacion = orientation_i.value,
               rampa = Rampas2.value, etiqueta = '', umbral = '', NUCOL = NumcoL.value, tax = Linaje.value,
            png = False, svg = False)

buttonnn.on_click(button_clickeddd)
#----------------------------------------------
boton333 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
boton333.style.button_color = 'gold'
output333 = widgets.Output()

def button_clicked333(b):
    with output333:
        clear_output(True)
        update_plot2(mostrar_dendrograma = show_dendrograma_i.value, mostrar_circulos = show_pies_i.value,
                     mostrar_metadata = show_anotaciones_i.value, orientacion = orientation_i.value,
               rampa = Rampas2.value, etiqueta = data_i.value, umbral = str(Percentage2.value), NUCOL = NumcoL.value, tax = Linaje.value,
            png = True, svg = False)

        
boton333.on_click(button_clicked333)
#----------------------------------------------
boton444 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
boton444.style.button_color = 'gold'
output444 = widgets.Output()

def button_clicked444(b):
    with output444:
        clear_output(True)
        update_plot2(mostrar_dendrograma = show_dendrograma_i.value, mostrar_circulos = show_pies_i.value,
                     mostrar_metadata = show_anotaciones_i.value, orientacion = orientation_i.value,
               rampa = Rampas2.value, etiqueta = data_i.value, umbral = str(Percentage2.value), NUCOL = NumcoL.value, tax = Linaje.value,
            png = False, svg = True)

        
boton444.on_click(button_clicked444)



def boxrampas2(Rampas2):
    barcolor2(color = INDIVIDUAL_COLORS[Rampas2])
OUTboxrampas2 = widgets.interactive_output(boxrampas2, {'Rampas2':Rampas2})



def box22(CCOOFF_i):
    barcolor(lista = QUALITATIVE_colors[CCOOFF_i])
    with open('Anexos16S/dictionaries2.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    variedad = dict(zip(list(variedad.keys()), QUALITATIVE_colors[CCOOFF_i][0:len(variedad)]))
    with open('Anexos16S/dictionaries2.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT22 = widgets.interactive_output(box22, {'CCOOFF_i':CCOOFF_i})

def box33(PPRROO_i):
    barcolor(lista = QUALITATIVE_colors[PPRROO_i])
    with open('Anexos16S/dictionaries2.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    procesado = dict(zip(list(procesado.keys()), QUALITATIVE_colors[PPRROO_i][0:len(procesado)]))
    with open('Anexos16S/dictionaries2.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT33 = widgets.interactive_output(box33, {'PPRROO_i':PPRROO_i})

def box44(CCUULL_i):
    barcolor(lista = QUALITATIVE_colors[CCUULL_i])
    with open('Anexos16S/dictionaries2.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    cultivo = dict(zip(list(cultivo.keys()), QUALITATIVE_colors[CCUULL_i][0:len(cultivo)]))
    with open('Anexos16S/dictionaries2.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT44 = widgets.interactive_output(box44, {'CCUULL_i':CCUULL_i})

def box55(TTIIMM_i):
    barcolor(lista = QUALITATIVE_colors[TTIIMM_i])
    with open('Anexos16S/dictionaries2.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    tiempo_secado = dict(zip(list(tiempo_secado.keys()), QUALITATIVE_colors[TTIIMM_i][0:len(tiempo_secado)]))
    with open('Anexos16S/dictionaries2.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT55 = widgets.interactive_output(box55, {'TTIIMM_i':TTIIMM_i})

def box66(KKIITT_i):
    barcolor(lista = QUALITATIVE_colors[KKIITT_i])
    with open('Anexos16S/dictionaries2.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    kit = dict(zip(list(kit.keys()), QUALITATIVE_colors[KKIITT_i][0:len(kit)]))
    with open('Anexos16S/dictionaries2.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT66 = widgets.interactive_output(box66, {'KKIITT_i':KKIITT_i})

def box77(OOTTAA_i):
    barcolor(lista = QUALITATIVE_colors[OOTTAA_i])
    with open('Anexos16S/dictionaries2.json', 'r') as fp:
        oooo = json.load(fp) 
    variedad = oooo[0]
    procesado = oooo[1]
    cultivo = oooo[2]
    tiempo_secado = oooo[3]
    kit = oooo[4]
    ota = oooo[5]
    ota = dict(zip(list(ota.keys()), QUALITATIVE_colors[OOTTAA_i][0:len(ota)]))
    with open('Anexos16S/dictionaries2.json', 'w') as fp:
        json.dump([variedad, procesado, cultivo, tiempo_secado, kit, ota], fp)
OUT77 = widgets.interactive_output(box77, {'OOTTAA_i':OOTTAA_i})



box_layout = Layout(flex_flow='column',
                    align_items='stretch',
                    border='1px solid silver',
                    width='170px',
                    height='467px',
                    display='flex')
carousel_i = Box(children=[HBox([widgets.Label('Variable colors')]),
             HBox([widgets.Label('Coffee_Variety:'), CCOOFF_i]),
             OUT22,
             HBox([widgets.Label('Processing:'), PPRROO_i]),
             OUT33,
             HBox([widgets.Label('Cultivation:'), CCUULL_i]),
             OUT44,
             HBox([widgets.Label('Time_Dry:'), TTIIMM_i]),
             OUT55,
             HBox([widgets.Label('Kit:'), KKIITT_i]),
             OUT66,
             HBox([widgets.Label('OTA:'), OOTTAA_i]),
             OUT77
             ], layout=box_layout)

box_layout111 = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='1px solid lime',
                    width='780px',
                   height='467px')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def SELECTSAM(SAM_SELECT = ''):
    import datetime
    sumary111 = []
    uno = open('AnexosITS/AAbundancEE.txt', 'r')
    for enu, line in enumerate(uno):
        line = line.rstrip()
        if enu == 0:
            header = line.split('\t')
        else:
            sumary111.append(line.split('\t'))
    uno.close()
    sumary2 = DataFrame(sumary111, columns = header)
    sumary2 = sumary2.set_index('Sample')
    sumary2 = sumary2.astype('float64')


    fig = plt.figure()
    ax1 = fig.add_axes([0, 0, 1, 1])
    sumary2.plot(kind='barh', stacked=True, colormap = plt.get_cmap(Rampas2.value), alpha = 1, ax = ax1)
    atributos = plt.legend(title='', title_fontsize = 8,  bbox_to_anchor=(1, 0.96), loc=2, ncol = 1,
                           prop={'size':7})
    plt.close()

    atributos2 = atributos.legendHandles
    posibles = [matplotlib.colors.to_hex(i) for i in np.array([i.get_facecolor() for i in atributos2])]
    handles, labels = ax1.get_legend_handles_labels()
    # correspondencia de colores de acuerdo al stacked plot
    CorrepondenciA = dict(zip(labels, posibles))


    df_one = sumary2.T[[SAM_SELECT]][sumary2.T[[SAM_SELECT]][SAM_SELECT] > 0].sort_values(by =SAM_SELECT,ascending=False)

    df_two = Sampledata[Sampledata.Sample == SAM_SELECT]
    df_two = df_two.set_index('Sample').T

    namesmax = max([len(i) for i in df_one.index])+2
    nummax = max([len(str(round(i, 4))) for i in df_one[SAM_SELECT]])
    namesmax2 = max([len(i) for i in df_two.index])+1
    nummax2 = max([len(i) for i in df_two[SAM_SELECT]])

    #**************************************************************************

    root= tk.Tk()
    root.title("Sample Data")
    root.geometry("950x450")
    root.configure(background='white')
    root.attributes("-topmost", True)


    cero = Label(root, text=SAM_SELECT, font=("Arial", 12,  "bold"), fg = 'black', bg = 'lime')
    cero.grid(column=0, row=0, sticky = W+E+S)
    
    
    cero1 = Label(root, text='16S Analysis', font=("Arial", 12,  "bold"), fg = 'silver', bg = 'white')
    cero1.grid(column=4, row=0, sticky = W+S)

    labebl = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 0, row = 1)

    primera = Label(root, text= Linaje.value, font=("Arial", 8,  "bold"), fg = 'black', bg = 'gainsboro')
    primera.grid(column=0, row=2, sticky = W+E+S)

    segunda = Label(root, text='Relative abundance', font=("Arial", 8,  "bold"), fg = 'black', bg = 'khaki')
    segunda.grid(column=1, row=2, sticky = W+E+S)


    #--------------------------------------------------------
    table = tk.Text(root, font=("Arial italic", 8), height=len(df_one), width=namesmax, fg = 'black', bg = 'antiquewhite')
    table.grid(column=0, row=3, sticky =  W+E+N)

    for h, i in enumerate(df_one.index):
        table.insert(tk.INSERT, '  '+i+'\n')
    table.config(state=DISABLED)

    table2 = tk.Text(root, font=("Arial", 8), height=len(df_one), width=nummax, fg = 'black', bg ='lightcyan' )
    table2.grid(column=1, row=3, sticky =  W+E+N)
    for h, i in enumerate(df_one[SAM_SELECT]):
        table2.insert(tk.INSERT, '  '+str(round(i, 4))+'\n')
    table2.config(state=DISABLED)


    ##########################################################
    labebl = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 0, row = 4)


    primera = Label(root, text='Variable', font=("Arial", 8,  "bold"), fg = 'black', bg = 'gainsboro')
    primera.grid(column=0, row=5, sticky = W+E+S)

    segunda = Label(root, text='Type', font=("Arial", 8,  "bold"), fg = 'black', bg = 'khaki')
    segunda.grid(column=1, row=5, sticky = W+E+S)

    #--------------------------------------------------------
    table3 = tk.Text(root, font=("Arial", 8,  "bold"), height=len(df_two), width=namesmax2, fg = 'black', bg = 'antiquewhite')
    table3.grid(column=0, row=6, sticky =  W+E+N)

    for h, i in enumerate(df_two.index):
        table3.insert(tk.INSERT, '  '+i+'\n')
    table3.config(state=DISABLED)

    table4 = tk.Text(root, font=("Arial", 8), height=len(df_two), width=nummax2, fg = 'black', bg ='lightcyan' )
    table4.grid(column=1, row=6, sticky =  W+E+N)
    for h, i in enumerate(df_two[SAM_SELECT]):
        table4.insert(tk.INSERT, '  '+i+'\n')
    table4.config(state=DISABLED)



    #--------------------------------------------------------
    labebl = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 3, row = 4)
    
    def salvarpng():
        mpl.rcParams.update(mpl.rcParamsDefault)
        figure = plt.figure(figsize=(6,3))

        ax = figure.add_axes([0, 0, 1, 1])

        ax.pie(df_one[SAM_SELECT].tolist(), #autopct ='%1.1f%%',
               labels = None,
                                    pctdistance = 1, labeldistance= 1,
                                    startangle = 0, radius = 0.4, rotatelabels = True,frame=True,
               center=(0.5, 0.5),
                                    colors = [CorrepondenciA[i] for i in df_one.index],
                                    wedgeprops={'alpha':1, 'linewidth': 0, 'edgecolor':'black'},
                                    explode = np.array([0.0]*len(df_one)), textprops=dict(size = 8))
        if len(df_one) <= 17:
            plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                       handletextpad=0.5,ncol= 1,title=Linaje.value, title_fontsize = 7,
                                       fancybox=True, framealpha=0.5, shadow=False,
                                       handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                       borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                       prop={'style':'italic', 'size':6.5})
        else:
            plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                       handletextpad=0.5,ncol= 2,title=Linaje.value, title_fontsize = 7,
                                       fancybox=True, framealpha=0.5, shadow=False,
                                       handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                       borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                       prop={'style':'italic', 'size':6.5})
    
        centre_circle = plt.Circle((0.5,0.5),0.2,fc = 'white')
        plt.gca().add_artist(centre_circle)

        ax.set_xlim(0,2)

        ax.axis('off')

        plt.savefig('Plots16S/'+SAM_SELECT+'_'+data_i.value+'_'+Linaje.value+'_'+str(Percentage2.value)+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
    
    labebl2 = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl2.grid(column = 0, row = 7)
    
    boton = Button(root, text=" PNG ", cursor="hand2",
                #activebackground= 'black',activeforeground= 'black',
                bg="gold", fg="black",font=("Arial", 8),
                    command = salvarpng) # newwin.destroy
    boton.grid(column = 0, row = 8)
    
    def salvarsvg():
        mpl.rcParams.update(mpl.rcParamsDefault)
        figure = plt.figure(figsize=(6,3))

        ax = figure.add_axes([0, 0, 1, 1])

        ax.pie(df_one[SAM_SELECT].tolist(), #autopct ='%1.1f%%',
               labels = None,
                                    pctdistance = 1, labeldistance= 1,
                                    startangle = 0, radius = 0.4, rotatelabels = True,frame=True,
               center=(0.5, 0.5),
                                    colors = [CorrepondenciA[i] for i in df_one.index],
                                    wedgeprops={'alpha':1, 'linewidth': 0, 'edgecolor':'black'},
                                    explode = np.array([0.0]*len(df_one)), textprops=dict(size = 8))
        if len(df_one) <= 17:
            plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                       handletextpad=0.5,ncol= 1,title=Linaje.value, title_fontsize = 7,
                                       fancybox=True, framealpha=0.5, shadow=False,
                                       handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                       borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                       prop={'style':'italic', 'size':6.5})
        else:
            plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                       handletextpad=0.5,ncol= 2,title=Linaje.value, title_fontsize = 7,
                                       fancybox=True, framealpha=0.5, shadow=False,
                                       handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                       borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                       prop={'style':'italic', 'size':6.5})
    
        centre_circle = plt.Circle((0.5,0.5),0.2,fc = 'white')
        plt.gca().add_artist(centre_circle)

        ax.set_xlim(0,2)

        ax.axis('off')

        plt.savefig('Plots16S/'+SAM_SELECT+'_'+data_i.value+'_'+Linaje.value+'_'+str(Percentage2.value)+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
    
    
    boton = Button(root, text=" SVG ", cursor="hand2",
                #activebackground= 'black',activeforeground= 'black',
                bg="gold", fg="black",font=("Arial", 8),
                    command = salvarsvg) # newwin.destroy
    boton.grid(column = 1, row = 8)
    
    #-------------------------------------


    group_aspect = LabelFrame(root, text = "", font=("Arial", 8,  "bold"), height = 1)
    group_aspect.grid(column=4, row=2, rowspan = 100, sticky= E+N)
    group_aspect.configure(background='white')


    mpl.rcParams.update(mpl.rcParamsDefault)
    figure = plt.figure(figsize=(6,3))

    bar1 = FigureCanvasTkAgg(figure, group_aspect)
    bar1.draw()
    bar1.get_tk_widget().grid(row=0, column = 0)#, rowspan=7, columnspan = 7


    ax = figure.add_axes([0, 0, 1, 1])

    ax.pie(df_one[SAM_SELECT].tolist(), #autopct ='%1.1f%%',
           labels = None,
                                pctdistance = 1, labeldistance= 1,
                                startangle = 0, radius = 0.4, rotatelabels = True,frame=True,
           center=(0.5, 0.5),
                                colors = [CorrepondenciA[i] for i in df_one.index],
                                wedgeprops={'alpha':1, 'linewidth': 0, 'edgecolor':'black'},
                                explode = np.array([0.0]*len(df_one)), textprops=dict(size = 8))
    if len(df_one) <= 17:
        plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                   handletextpad=0.5,ncol= 1, title=Linaje.value, title_fontsize = 7,
                                   fancybox=True, framealpha=0.5, shadow=False,
                                   handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                   borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                   prop={'style':'italic', 'size':6.5})
    else:
        plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                   title=Linaje.value, title_fontsize = 7,
                   handletextpad=0.5,ncol= 2,
                                   fancybox=True, framealpha=0.5, shadow=False,
                                   handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                   borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                   prop={'style':'italic', 'size':6.5})

    centre_circle = plt.Circle((0.5,0.5),0.2,fc = 'white')
    plt.gca().add_artist(centre_circle)

    ax.set_xlim(0,2)

    ax.axis('off')
    

    root.mainloop()

Sample_Select = widgets.Dropdown(options=list(ordenado3.values()),value=list(ordenado3.values())[0],disabled=False,
                         layout=Layout(width='88px', height='40px'))


def samsel(Sample_Select):

    bot0 = widgets.Button(description=Sample_Select, layout=Layout(width='80px'))
    bot0.style.button_color = 'cyan'
    bot0.style.font_weight = 'bold'
    outbot0 = widgets.Output()
    
    display(bot0, outbot0)
    def button_bot0(b):
        with outbot0:
            clear_output(True)
            SELECTSAM(SAM_SELECT = Sample_Select)

    bot0.on_click(button_bot0)
    
samselOUT = widgets.interactive_output(samsel, {'Sample_Select':Sample_Select})

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TAXONOMI_IND = VBox([
    #negra,
      Box([HBox([VBox([widgets.Label('Tax Level:'), VBox([Linaje, OUT11])]),
                 VBox([widgets.Label('Clustering:'), boton_data_i]),
                 VBox([widgets.Label('Method/Metric:'), Method_i, Metric_i]),
                 VBox([widgets.Label('Dendrog.:'), boton_dendrogram_i]),
                 VBox([widgets.Label('Circles:'), boton_pies_i]), 
                 VBox([widgets.Label('Metadata:'), boton_metadata_i]),
                 VBox([widgets.Label('Orientarion:'), hor_ver_i]),
                 VBox([widgets.Label('Stacked Color:'), VBox([Rampas2, OUTboxrampas2])]),
                 VBox([widgets.Label('Explore Samples:'),
     Box(children=[HBox([Sample_Select, samselOUT])], layout= Layout(border='1px solid pink', width='179px', height='46px'))])             
                ])],
          
      layout = Layout(border='1px solid silver', width='950px')),
    Box(children = [Percentage2], layout = Layout(border='1px solid silver', width='950px')),
    Box(children = [HBox([buttonnn, NumcoL, boton333, boton444])],  layout = Layout(border='1px solid silver', width='950px')),
    HBox([Box(children=[outputtt], layout=box_layout111), VBox([carousel_i, boton4])])
])




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


indices_ASVs = pd.read_csv('Anexos16S/ASVs_Indices.tsv', sep = '\t')
indices_OTUs = pd.read_csv('Anexos16S/OTUs_Indices.tsv', sep = '\t')



idividuosASVs = indices_ASVs[indices_ASVs['Unnamed: 0'] == 'Individuals'].iloc[:, 1:]


ASVs_counts = dict(zip(list(idividuosASVs.columns), list(idividuosASVs.values[0])))


idividuosOTUs = indices_OTUs[indices_OTUs['Unnamed: 0'] == 'Individuals'].iloc[:, 1:]


OTUs_counts = dict(zip(list(idividuosOTUs.columns), list(idividuosOTUs.values[0])))


def indice_plot(data_ind = 'ASV', tipo_indice = 'Taxa_S', show_barras = True,
                show_puntos = True, show_texto = True, ColormaP = 'RdYlGn', tamano = 7,
                png = False, svg = False, linea = True, ancholinea = 1):
    if data_ind == 'ASV':
        Indices = indices_ASVs
        dict_counts = ASVs_counts
    elif data_ind == 'OTU':
        Indices = indices_OTUs
        dict_counts = OTUs_counts
    
    INDICE = tipo_indice
    SIZe = tamano
    LETRAS = SIZe

    EJE_X = np.array(range(len(Indices.iloc[:, 1:].columns)))
    EJE_Y = np.array(list(Indices[Indices['Unnamed: 0'] == INDICE][list(ordenado3.values())].iloc[:, :].values[0]))
    PESOS = np.array([dict_counts[ordenado3[i]] for i in ordenado3])/500

    mpl.rcParams.update(mpl.rcParamsDefault)

    fig, ax = plt.subplots(figsize=(SIZe, ((SIZe * 2.5) /  7)))
    
    if show_puntos == True:
        ax.scatter(EJE_X, EJE_Y,
                   marker='o', c=[matplotlib.colors.to_hex(i) for i in plt.get_cmap(ColormaP)(np.arange(len(Indices.iloc[:, 1:].columns))/len(Indices.iloc[:, 1:].columns))],
                   s=PESOS, alpha=1,linewidth=0)
        ax.scatter(np.array([(len(ordenado3)+1), (len(ordenado3)+1), (len(ordenado3)+1)]), np.array([ax.get_ylim()[1]*0.8,
                                                     ax.get_ylim()[1]*0.69,
                                                     ax.get_ylim()[1]*0.55]),
                   zorder=10, s = np.array([50000/500, 100000/500, 150000/500]),
                       marker = 'o', c = 'black', linewidth=0,
                       edgecolors = 'white',alpha=1)

        plt.text((len(ordenado3)+2), ax.get_ylim()[1]*0.8, ' 50000', color = 'black', fontsize = LETRAS, ha = 'left', va = 'center')
        plt.text((len(ordenado3)+2), ax.get_ylim()[1]*0.69, '100000', color = 'black', fontsize = LETRAS, ha = 'left', va = 'center')
        plt.text((len(ordenado3)+2), ax.get_ylim()[1]*0.55, '150000', color = 'black', fontsize = LETRAS, ha = 'left', va = 'center')
        plt.text((len(ordenado3)+1), ax.get_ylim()[1]*0.92,
                     'Individuals', color = 'black', fontsize = LETRAS+1, ha = 'left', va = 'center', weight = 'bold')
        
    else:
        ax.scatter(EJE_X, EJE_Y,
                   marker='o', c='grey',
                   s=20, alpha=1,linewidth=0)
        
    if linea == True:
        ax.plot(EJE_X, EJE_Y, marker='o', linestyle='-', label = '',
                            markeredgewidth=0,markersize=2, zorder = 2, color='black', linewidth=1, alpha = 0.3)
    if show_barras == True:
        ax.plot([EJE_X, EJE_X],
                 [EJE_Y, np.repeat(ax.get_ylim()[0], len(EJE_Y))],
                 marker='o', linestyle='-', label = '',
                            markeredgewidth=0,markersize=0, zorder = 0,
                color='silver',
                linewidth=ancholinea, alpha = 0.3)
    
    if show_texto == True:
        for x, y in zip(EJE_X, EJE_Y):
            ax.text(x+0.1, y, '     '+str(round(y, 2)), size = LETRAS, ha = 'center', va = 'bottom', rotation = 90)
    #

    plt.yticks(size=7)
    plt.gca().tick_params(which='major', width = 0.3, length=2, color='black')
    plt.gca().spines['left'].set_position(('data',-1))
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['right'].set_color(None)
    plt.gca().spines['top'].set_color(None)

    plt.xticks(list(range(len(ordenado4))), [ordenado3[i] for i in ordenado4], size=LETRAS+1, rotation=90, weight = 'bold')
    plt.yticks(size=LETRAS)

    plt.gca().set_ylabel(INDICE, fontsize=LETRAS+1, weight = 'bold')

    plt.ylim(ax.get_ylim()[0], ax.get_ylim()[1] + (ax.get_ylim()[1]*0.15))
    plt.gca().set_xlabel('', fontsize=0)

    if png == True:
        plt.savefig('Plots16S/'+data_ind+'_Diversity_Index_'+tipo_indice+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
    if svg == True:
        plt.savefig('Plots16S/'+data_ind+'_Diversity_Index_'+tipo_indice+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')

    plt.show()





data_ind = widgets.ToggleButtons(options=['ASV', 'OTU'])
data_ind.layout.width = '10%'
#data.style.font_weight = 'bold'
data_ind.style.button_width = '60px'
boton_data_ind = Box(children=[data_ind], layout= Layout(border='1px solid pink', width='69px', height='63px'))





show_barras = widgets.ToggleButtons(options=[True, False])
show_barras.layout.width = '10%'
show_barras.style.button_width = '60px'
boton_barras = Box(children=[show_barras], layout= Layout(border='1px solid pink', width='69px', height='63px'))





show_puntos = widgets.ToggleButtons(options=[True, False])
show_puntos.layout.width = '10%'
show_puntos.style.button_width = '60px'
boton_show_puntos = Box(children=[show_puntos], layout= Layout(border='1px solid pink', width='69px', height='63px'))





show_linea = widgets.ToggleButtons(options=[True, False])
show_linea.layout.width = '10%'
show_linea.style.button_width = '60px'
boton_show_linea = Box(children=[show_linea], layout= Layout(border='1px solid pink', width='69px', height='63px'))





show_texto = widgets.ToggleButtons(options=[True, False])
show_texto.layout.width = '10%'
show_texto.style.button_width = '60px'
boton_show_texto = Box(children=[show_texto], layout= Layout(border='1px solid pink', width='69px', height='63px'))





tipo_indice = widgets.Dropdown(options=indices_ASVs['Unnamed: 0'].tolist(),value='Taxa_S',disabled=False,
                         layout=Layout(width='88px', height='25px'))





but_ind_plot = widgets.Button(description=" UPDATE PLOT ", icon = 'fa-refresh')
but_ind_plot.style.button_color = 'lime' #'deepskyblue'
but_ind_plot.style.font_weight = 'bold'
out_ind_plot = widgets.Output()

def button_ind_plot(b):
    with out_ind_plot:
        clear_output(True)
        indice_plot(data_ind = data_ind.value, tipo_indice = tipo_indice.value, show_barras = show_barras.value,
                show_puntos = show_puntos.value, show_texto = show_texto.value, ColormaP = Rampas3.value,
                   tamano = cor_size2[tamano_plot_ind.value], png = False, svg = False, linea = show_linea.value,
                   ancholinea = lineancho.value)

but_ind_plot.on_click(button_ind_plot)
#----------------------------------------------
but_ind_plot_save = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
but_ind_plot_save.style.button_color = 'gold'
out_ind_plot_save = widgets.Output()

def button_ind_plot_save(b):
    with out_ind_plot_save:
        clear_output(True)
        indice_plot(data_ind = data_ind.value, tipo_indice = tipo_indice.value, show_barras = show_barras.value,
                show_puntos = show_puntos.value, show_texto = show_texto.value, ColormaP = Rampas3.value,
                   tamano = cor_size2[tamano_plot_ind.value], png = True, svg = False, linea = show_linea.value,
                   ancholinea = lineancho.value)

but_ind_plot_save.on_click(button_ind_plot_save)
#----------------------------------------------
but_ind_plot_save2 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
but_ind_plot_save2.style.button_color = 'gold'
out_ind_plot_save2 = widgets.Output()

def button_ind_plot_save2(b):
    with out_ind_plot_save2:
        clear_output(True)
        indice_plot(data_ind = data_ind.value, tipo_indice = tipo_indice.value, show_barras = show_barras.value,
                show_puntos = show_puntos.value, show_texto = show_texto.value, ColormaP = Rampas3.value,
                   tamano = cor_size2[tamano_plot_ind.value], png = False, svg = True, linea = show_linea.value,
                   ancholinea = lineancho.value)

but_ind_plot_save2.on_click(button_ind_plot_save2)

cor_size2 = dict(zip(([1, 2, 3, 4, 5, 6, 7, 8]), ([7, 7.5, 8, 8.5, 9, 9.5, 10, 11])))


tamano_plot_ind = widgets.SelectionSlider(options=list(cor_size2.keys()),value=4,disabled=False,
                                      description = 'Size Plot:',
                                continuous_update=False,orientation='horizontal',readout=True)

lineancho = widgets.SelectionSlider(options=[2,3,4,5,6,7,8,9,10],value=5,disabled=False,
                                    layout=Layout(width='100px', height='25px'),
                                      #description = 'Size Plot:',
                                continuous_update=False,orientation='horizontal',readout=False)
def lin(lineancho):
    print('      '+str(lineancho))
OUTlin = widgets.interactive_output(lin, {'lineancho':lineancho})



Rampas3 = widgets.Dropdown(options=list(INDIVIDUAL_COLORS.keys()),value='RdYlGn',
                          disabled=False,
                         layout=Layout(width='98px', height='25px'))




def boxrampas3(Rampas3):
    barcolor2(color = INDIVIDUAL_COLORS[Rampas3])
OUTboxrampas3 = widgets.interactive_output(boxrampas3, {'Rampas3':Rampas3})



box_layout1111 = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='1px solid lime',
                    width='950px',
                   height='410px')


INDICES_16S = VBox([Box([HBox([VBox([widgets.Label('Clustering:'), boton_data_ind]),
                               VBox([widgets.Label('Index:'), tipo_indice]),
                               VBox([widgets.Label('Bars:'), boton_barras]),
                               VBox([widgets.Label('Bar width:'), lineancho, OUTlin]),
                               VBox([widgets.Label('Points:'), boton_show_puntos]),
                               VBox([widgets.Label('Text:'), boton_show_texto]),
                               VBox([widgets.Label('Line:'), boton_show_linea]),
                               VBox([widgets.Label('Points Color:'), Rampas3, OUTboxrampas3])
                               ])], layout = Layout(border='1px solid silver', width='950px')),
                    Box([HBox([but_ind_plot, tamano_plot_ind, but_ind_plot_save, but_ind_plot_save2])], 
                        layout = Layout(border='1px solid silver', width='950px')),
                    HBox([Box([out_ind_plot], layout=box_layout1111)])
                    ])











#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





def index_plots(data_ind2 = 'ASV', tipo_indice2 = 'Taxa_S', ColormaP2 = 'RdYlGn', var = 'Processing',
               SIG_val = 0.05, SIG = True, TAMANO = 3.5, showmedia = True, showstd = True, png = False, svg = False):
    if data_ind2 == 'ASV':
        Indices = indices_ASVs
    elif data_ind2 == 'OTU':
        Indices = indices_OTUs
    
    LETRA = (10 * TAMANO) / 4
    LETRA2 = (9 * TAMANO) / 4
    
    Indices_transpose = Indices.T.iloc[1:,:]
    Indices_transpose.columns = list(indices_ASVs['Unnamed: 0'])
    index_remove = []
    for i in Indices_transpose.index:
        if re.search('Lower.*|Upper.*', i):
            index_remove.append(re.search('Lower.*|Upper.*', i).group())
    Indices_transpose = Indices_transpose.drop(index = index_remove)
    Indices_transpose.insert(loc = 0, column='Sample', value=Indices_transpose.index)
    Indices_transpose = Indices_transpose.reset_index(drop = True)
    Indices_transpose = Indices_transpose.merge(Sampledata, on = 'Sample', how = 'left')
    
    indice = tipo_indice2
    t = var
    
    It2 = Indices_transpose[['Sample', indice, t]]
    It2[indice] = It2[indice].astype(float).tolist()
    
    if len(set(It2[t].tolist())) == 2:
        datos = {}
        for E in set(It2[t].tolist()):
            datos[E] = It2[It2[t] == E][indice].tolist()
        T_test = stats.ttest_ind(datos[list(datos.keys())[0]], datos[list(datos.keys())[1]], equal_var=False)
        if T_test.pvalue < SIG_val:
            if T_test.pvalue < 0.009999999999999999999999:
                #print('Samples:', i, 'P-value:', format(T_test.pvalue, '.3e'))
                significancia = format(T_test.pvalue, '.3e')
                ColoRlin = 'black'
            else:
                #print('Samples:', i, 'P-value:', round(T_test.pvalue, 3))
                significancia = round(T_test.pvalue, 3)
                ColoRlin = 'black'
        
        else:
            significancia = ''
            ColoRlin = 'white'
        
        mpl.rcParams.update(mpl.rcParamsDefault)

        sns.set(style="whitegrid")
        fig, ax = plt.subplots(figsize=(TAMANO+1, TAMANO))
        ax = sns.boxplot(x=t, y=indice, data=It2, palette=ColormaP2, hue=None, saturation=1, width=0.4,
                         order = list(set(It2[t].tolist())),
                        fliersize=0, linewidth=1) #showmeans=True
        labels_Y = ax.get_yticks()
        ax.tick_params(which='major', width = 3, length=2, color='black')

        ax = sns.swarmplot(x=t, y=indice, data=It2, color="black", order = list(set(It2[t].tolist())))
        
        for e, d in enumerate(set(It2[t].tolist())):
            ed = (e, d)
            if showmedia == True:
                ax.scatter(e, np.mean(datos[d]), c='red', s=1900, alpha=1,linewidth=1, marker='_')
            if showstd == True:
                ax.scatter(e, (np.mean(datos[d])) - np.std(datos[d], ddof = 1), marker='_',
                           c='blue', s=400, alpha=1,linewidth=1)
                ax.scatter(e, (np.mean(datos[d])) + np.std(datos[d], ddof = 1), marker='_',
                           c='blue', s=400, alpha=1,linewidth=1)
                ax.plot([e, e],[(np.mean(datos[d])) - np.std(datos[d], ddof = 1),
                                (np.mean(datos[d])) + np.std(datos[d], ddof = 1)], '-',
                        linewidth=1, color='blue')
            
        plt.yticks(np.round(np.array(labels_Y)[0:-1].astype(float), 3), np.round(np.array(labels_Y)[0:-1].astype(float), 3), size=LETRA2)
        plt.xticks(size=LETRA)
        
        if SIG == True:
            limy = ax.get_ylim()[1] + (ax.get_ylim()[1] * 0.05) # mas el 10%
            ax.plot([0, 1], [limy, limy], '-', linewidth=0.5, color=ColoRlin)
            ax.plot([0, 0], [limy, limy-(limy * 0.02)], '-', linewidth=0.5, color=ColoRlin)
            ax.plot([1, 1], [limy, limy-(limy * 0.02)], '-', linewidth=0.5, color=ColoRlin)
            ax.text(1+0.1, limy, significancia, ha = 'left', va = 'center', fontsize = LETRA2-1)
        
        
        plt.gca().spines['right'].set_color(None)
        plt.gca().spines['top'].set_color(None)
        plt.gca().set_xlabel(var, fontsize=LETRA+1)
        plt.gca().set_ylabel(tipo_indice2, fontsize=LETRA+1)
        
        if png == True:
            plt.savefig('Plots16S/'+data_ind2+'_Diversity_Index_'+var+'_'+tipo_indice2+'_'+str(SIG_val)+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
        if svg == True:
            plt.savefig('Plots16S/'+data_ind2+'_Diversity_Index_'+var+'_'+tipo_indice2+'_'+str(SIG_val)+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
        plt.show()
        
    else:
        if var == 'Time_Dry':
            unicos = Sampledata.Time_Dry.unique().tolist()
            orden_var = sorted([int(i) for i in unicos])
            orden_var = [str(i) for i in orden_var]
        if var == 'OTA':
            unicos = Sampledata.OTA.unique().tolist()
            orden_var = sorted([float(i) for i in unicos])
            orden_var = [str(i) for i in orden_var]


        significativos = []
        datos1 = {}
        for E in orden_var:
            datos1[E] = It2[It2[t] == E][indice].tolist()
        for i in list(itertools.combinations(list(datos1.keys()), 2)):
            T_test_m = stats.ttest_ind(datos1[i[0]], datos1[i[1]], equal_var=False)
            
            if T_test_m.pvalue < SIG_val:
                if T_test_m.pvalue < 0.009999999999999999999999:
                    #print('Samples:', i, 'P-value:', format(T_test_m.pvalue, '.3e'))
                    significativos.append([i, format(T_test_m.pvalue, '.3e')])
                else:
                    #print('Samples:', i, 'P-value:', round(T_test_m.pvalue, 3))
                    significativos.append([i, round(T_test_m.pvalue, 3)])
        
        mpl.rcParams.update(mpl.rcParamsDefault)

        sns.set(style="whitegrid")
        fig, ax = plt.subplots(figsize=((TAMANO*2)-1.5, TAMANO))

        ax = sns.boxplot(x=t, y=indice, data=It2, palette=ColormaP2, hue=None, saturation=1, width=0.4,
                         order = orden_var,
                        fliersize=0, linewidth=0.5) #showmeans=True
        labels_Y = ax.get_yticks()
        
        ax.tick_params(which='major', width = 3, length=2, color='black')

        ax = sns.swarmplot(x=t, y=indice, data=It2, color="black",
                           order = orden_var)
        
        if showmedia == True:
            ax.scatter(np.array(range(len(orden_var))), [np.mean(datos1[d]) for d in orden_var],
                       c='red', s=400, alpha=1,linewidth=1, marker='_')
        if showstd == True:
            ax.scatter(np.array(range(len(orden_var))), np.array([np.mean(datos1[d]) for d in orden_var]) - np.array([np.std(datos1[d], ddof = 1) for d in orden_var]),
                       marker='_', c='blue', s=100, alpha=1,linewidth=1)
            ax.scatter(np.array(range(len(orden_var))), np.array([np.mean(datos1[d]) for d in orden_var]) + np.array([np.std(datos1[d], ddof = 1) for d in orden_var]),
                       marker='_', c='blue', s=100, alpha=1,linewidth=1)
            ax.plot([np.array(range(len(orden_var))), np.array(range(len(orden_var)))],
                    [np.array([np.mean(datos1[d]) for d in orden_var]) - np.array([np.std(datos1[d], ddof = 1) for d in orden_var]),
                                np.array([np.mean(datos1[d]) for d in orden_var]) + np.array([np.std(datos1[d], ddof = 1) for d in orden_var])],
                    '-', linewidth=1, color='blue')
        
        if SIG == True:
            if len(significativos) > 0:
                correspondencia_order = dict(zip(orden_var, list(ax.get_xticks())))
                limy = ax.get_ylim()[1] + (ax.get_ylim()[1] * 0.05)
                al = limy * 0.02
                distancia = limy * 0.05
                for i in significativos:
                    uno = correspondencia_order[i[0][0]]
                    dos = correspondencia_order[i[0][1]]
                    #print(uno, dos, i[1])
                     # mas el 10%
                    ax.plot([uno, dos], [limy, limy], '-', linewidth=0.5, color='black')
                    ax.plot([uno, uno], [limy, limy-al], '-', linewidth=0.5, color='black')
                    ax.plot([dos, dos], [limy, limy-al], '-', linewidth=0.5, color='black')
                    ax.text(dos+0.1, limy, i[1], ha = 'left', va = 'center', fontsize = LETRA2-1)
                    limy += distancia
                
        plt.yticks(np.round(np.array(labels_Y)[0:-1].astype(float), 3), np.round(np.array(labels_Y)[0:-1].astype(float), 3), size=LETRA2)
        plt.xticks(size=LETRA)
        plt.gca().spines['right'].set_color(None)
        plt.gca().spines['top'].set_color(None)
        plt.gca().set_xlabel(var, fontsize=LETRA+1)
        plt.gca().set_ylabel(tipo_indice2, fontsize=LETRA+1)
        
        if png == True:
            plt.savefig('Plots16S/'+data_ind2+'_Diversity_Index_'+var+'_'+tipo_indice2+'_'+str(SIG_val)+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
        if svg == True:
            plt.savefig('Plots16S/'+data_ind2+'_Diversity_Index_'+var+'_'+tipo_indice2+'_'+str(SIG_val)+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
        plt.show()





data_caja = widgets.ToggleButtons(options=['ASV', 'OTU'])
data_caja.layout.width = '10%'
#data.style.font_weight = 'bold'
data_caja.style.button_width = '60px'
boton_data_caja = Box(children=[data_caja], layout= Layout(border='1px solid pink', width='69px', height='63px'))


tipo_indice2 = widgets.Dropdown(options=indices_ASVs['Unnamed: 0'].tolist(),value='Taxa_S',disabled=False,
                         layout=Layout(width='88px', height='25px'))


show_sigbar = widgets.ToggleButtons(options=[True, False])
show_sigbar.layout.width = '10%'
show_sigbar.style.button_width = '60px'
boton_sigbar = Box(children=[show_sigbar], layout= Layout(border='1px solid pink', width='69px', height='63px'))


show_media = widgets.ToggleButtons(options=[True, False])
show_media.layout.width = '10%'
show_media.style.button_width = '60px'
boton_media = Box(children=[show_media], layout= Layout(border='1px solid pink', width='69px', height='63px'))


show_std = widgets.ToggleButtons(options=[True, False])
show_std.layout.width = '10%'
show_std.style.button_width = '60px'
boton_std = Box(children=[show_std], layout= Layout(border='1px solid pink', width='69px', height='63px'))


significancia = widgets.SelectionSlider(options=[0.01, 0.02, 0.03, 0.04, 0.05],value=0.05,disabled=False,
                                    layout=Layout(width='80px', height='25px'),
                                      #description = 'Size Plot:',
                                continuous_update=False,orientation='horizontal',readout=False)
def sig(significancia):
    print('  '+str(significancia)+' ')
OUTsignifi = widgets.interactive_output(sig, {'significancia':significancia})


VariableS = widgets.Dropdown(options=variables,value='Coffee_Variety',disabled=False,
                         layout=Layout(width='100px', height='25px'))



cor_size3 = dict(zip(([1, 2, 3, 4, 5, 6]), ([4, 4.2, 4.4, 4.6, 4.8, 5])))



tamano_plot_sig = widgets.SelectionSlider(options=list(cor_size3.keys()),value=4,disabled=False,
                                      description = 'Size Plot:',
                                continuous_update=False,orientation='horizontal',readout=True)



Rampas4 = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='Paired',
                          disabled=False,
                         layout=Layout(width='98px', height='25px'))
def boxrampas4(Rampas4):
    barcolor(lista = QUALITATIVE_colors[Rampas4])
OUTboxrampas4 = widgets.interactive_output(boxrampas4, {'Rampas4':Rampas4})



but_ind_sig = widgets.Button(description=" UPDATE PLOT ", icon = 'fa-refresh')
but_ind_sig.style.button_color = 'lime' #'deepskyblue'
but_ind_sig.style.font_weight = 'bold'
out_ind_sig = widgets.Output()

def button_ind_sig(b):
    with out_ind_sig:
        clear_output(True)
        index_plots(data_ind2 = data_caja.value, tipo_indice2 = tipo_indice2.value, ColormaP2 = Rampas4.value, var = VariableS.value,
               SIG_val = significancia.value, SIG = show_sigbar.value, TAMANO = cor_size3[tamano_plot_sig.value],
                   showmedia = show_media.value, showstd = show_std.value, png = False, svg = False)

but_ind_sig.on_click(button_ind_sig)
#------------------------------------------------------
but_ind_sig_save = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
but_ind_sig_save.style.button_color = 'gold' #'deepskyblue'
out_ind_sig_save = widgets.Output()

def button_ind_sig_save(b):
    with out_ind_sig_save:
        clear_output(True)
        index_plots(data_ind2 = data_caja.value, tipo_indice2 = tipo_indice2.value, ColormaP2 = Rampas4.value, var = VariableS.value,
               SIG_val = significancia.value, SIG = show_sigbar.value, TAMANO = cor_size3[tamano_plot_sig.value],
                   showmedia = show_media.value, showstd = show_std.value, png = True, svg = False)


but_ind_sig_save.on_click(button_ind_sig_save)
#------------------------------------------------------
but_ind_sig_save2 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
but_ind_sig_save2.style.button_color = 'gold' #'deepskyblue'
out_ind_sig_save2 = widgets.Output()

def button_ind_sig_save2(b):
    with out_ind_sig_save2:
        clear_output(True)
        index_plots(data_ind2 = data_caja.value, tipo_indice2 = tipo_indice2.value, ColormaP2 = Rampas4.value, var = VariableS.value,
               SIG_val = significancia.value, SIG = show_sigbar.value, TAMANO = cor_size3[tamano_plot_sig.value],
                   showmedia = show_media.value, showstd = show_std.value, png = False, svg = True)


but_ind_sig_save2.on_click(button_ind_sig_save2)


box_layout2222 = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='1px solid lime',
                    width='950px',
                   height='500px')



SIGNIFICATIVO = VBox([Box([HBox([VBox([widgets.Label('Clustering:'), boton_data_caja]),
                                 VBox([widgets.Label('Index:'), tipo_indice2]),
                                 VBox([widgets.Label('Variable:'), VariableS]),
                                 VBox([widgets.Label('Significance:'), boton_sigbar]),
                                 VBox([widgets.Label('Pvalue:'), significancia, OUTsignifi]),
                                 VBox([widgets.Label('Media:'), boton_media]),
                                 VBox([widgets.Label('STD:'), boton_std]),
                                 VBox([widgets.Label('Box Color:'), Rampas4, OUTboxrampas4])
                               ])], layout = Layout(border='1px solid silver', width='950px')),
                      Box([HBox([but_ind_sig, tamano_plot_sig, but_ind_sig_save, but_ind_sig_save2])], 
                          layout = Layout(border='1px solid silver', width='950px')),
                      HBox([Box([out_ind_sig], layout=box_layout2222)])
                    ])



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



from numpy import zeros
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import squareform




###############  obtenido a partir de SKBIO
from functools import partial
def _preprocess_input(distance_matrix, grouping, column):
    """Compute intermediate results not affected by permutations.

    These intermediate results can be computed a single time for efficiency,
    regardless of grouping vector permutations (i.e., when calculating the
    p-value). These intermediate results are used by both ANOSIM and PERMANOVA.

    Also validates and normalizes input (e.g., converting ``DataFrame`` column
    into grouping vector).

    """
    #if not isinstance(distance_matrix, DistanceMatrix):
    #    raise TypeError("Input must be a DistanceMatrix.")

    if isinstance(grouping, pd.DataFrame):
        if column is None:
            raise ValueError(
                "Must provide a column name if supplying a DataFrame.")
        else:
            grouping = _df_to_vector(distance_matrix, grouping, column)
    elif column is not None:
        raise ValueError(
            "Must provide a DataFrame if supplying a column name.")

    sample_size = distance_matrix.shape[0]
    if len(grouping) != sample_size:
        raise ValueError(
            "Grouping vector size must match the number of IDs in the "
            "distance matrix.")

    # Find the group labels and convert grouping to an integer vector
    # (factor).
    groups, grouping = np.unique(grouping, return_inverse=True)
    num_groups = len(groups)

    if num_groups == len(grouping):
        raise ValueError(
            "All values in the grouping vector are unique. This method cannot "
            "operate on a grouping vector with only unique values (e.g., "
            "there are no 'within' distances because each group of objects "
            "contains only a single object).")
    if num_groups == 1:
        raise ValueError(
            "All values in the grouping vector are the same. This method "
            "cannot operate on a grouping vector with only a single group of "
            "objects (e.g., there are no 'between' distances because there is "
            "only a single group).")

    tri_idxs = np.triu_indices(sample_size, k=1)
    distances = squareform(distance_matrix, force='tovector', checks=False) # CAMBIOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    return sample_size, num_groups, grouping, tri_idxs, distances

    
    
def _df_to_vector(distance_matrix, df, column):
    """Return a grouping vector from a ``DataFrame`` column.

    Parameters
    ----------
    distance_marix : DistanceMatrix
        Distance matrix whose IDs will be mapped to group labels.
    df : pandas.DataFrame
        ``DataFrame`` (indexed by distance matrix ID).
    column : str
        Column name in `df` containing group labels.

    Returns
    -------
    list
        Grouping vector (vector of labels) based on the IDs in
        `distance_matrix`. Each ID's label is looked up in the ``DataFrame``
        under the column specified by `column`.

    Raises
    ------
    ValueError
        If `column` is not in the ``DataFrame``, or a distance matrix ID is
        not in the ``DataFrame``.

    """
    if column not in df:
        raise ValueError("Column '%s' not in DataFrame." % column)

    grouping = df.reindex(df.index, axis=0).loc[:, column] # CAMBIOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    if grouping.isnull().any():
        raise ValueError(
            "One or more IDs in the distance matrix are not in the data "
            "frame.")
    return grouping.tolist()

def _compute_f_stat(sample_size, num_groups, tri_idxs, distances, group_sizes,
                    s_T, grouping):
    """Compute PERMANOVA pseudo-F statistic."""
    # Create a matrix where objects in the same group are marked with the group
    # index (e.g. 0, 1, 2, etc.). objects that are not in the same group are
    # marked with -1.
    grouping_matrix = -1 * np.ones((sample_size, sample_size), dtype=int)
    for group_idx in range(num_groups):
        within_indices = _index_combinations(
            np.where(grouping == group_idx)[0])
        grouping_matrix[within_indices] = group_idx

    # Extract upper triangle (in same order as distances were extracted
    # from full distance matrix).
    grouping_tri = grouping_matrix[tri_idxs]

    # Calculate s_W for each group, accounting for different group sizes.
    s_W = 0
    for i in range(num_groups):
        s_W += (distances[grouping_tri == i] ** 2).sum() / group_sizes[i]

    s_A = s_T - s_W
    return (s_A / (num_groups - 1)) / (s_W / (sample_size - num_groups))

def _index_combinations(indices):
    # Modified from http://stackoverflow.com/a/11144716
    return np.tile(indices, len(indices)), np.repeat(indices, len(indices))
def _run_monte_carlo_stats(test_stat_function, grouping, permutations):
    """Run stat test and compute significance with Monte Carlo permutations."""
    if permutations < 0:
        raise ValueError(
            "Number of permutations must be greater than or equal to zero.")

    stat = test_stat_function(grouping)

    p_value = np.nan
    if permutations > 0:
        perm_stats = np.empty(permutations, dtype=np.float64)

        for i in range(permutations):
            perm_grouping = np.random.permutation(grouping)
            perm_stats[i] = test_stat_function(perm_grouping)

        p_value = ((perm_stats >= stat).sum() + 1) / (permutations + 1)

    return stat, p_value
def _build_results(method_name, test_stat_name, sample_size, num_groups, stat,
                   p_value, permutations):
    """Return ``pandas.Series`` containing results of statistical test."""
    return pd.Series(
        data=[method_name, test_stat_name, sample_size, num_groups, stat,
              p_value, permutations],
        index=['method name', 'test statistic name', 'sample size',
               'number of groups', 'test statistic', 'p-value',
               'number of permutations'],
        name='%s results' % method_name)
def permanova(distance_matrix, grouping, column=None, permutations=999):
    """Test for significant differences between groups using PERMANOVA.

    Permutational Multivariate Analysis of Variance (PERMANOVA) is a
    non-parametric method that tests whether two or more groups of objects
    (e.g., samples) are significantly different based on a categorical factor.
    It is conceptually similar to ANOVA except that it operates on a distance
    matrix, which allows for multivariate analysis. PERMANOVA computes a
    pseudo-F statistic.

    Statistical significance is assessed via a permutation test. The assignment
    of objects to groups (`grouping`) is randomly permuted a number of times
    (controlled via `permutations`). A pseudo-F statistic is computed for each
    permutation and the p-value is the proportion of permuted pseudo-F
    statisics that are equal to or greater than the original (unpermuted)
    pseudo-F statistic.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    grouping : 1-D array_like or pandas.DataFrame
        Vector indicating the assignment of objects to groups. For example,
        these could be strings or integers denoting which group an object
        belongs to. If `grouping` is 1-D ``array_like``, it must be the same
        length and in the same order as the objects in `distance_matrix`. If
        `grouping` is a ``DataFrame``, the column specified by `column` will be
        used as the grouping vector. The ``DataFrame`` must be indexed by the
        IDs in `distance_matrix` (i.e., the row labels must be distance matrix
        IDs), but the order of IDs between `distance_matrix` and the
        ``DataFrame`` need not be the same. All IDs in the distance matrix must
        be present in the ``DataFrame``. Extra IDs in the ``DataFrame`` are
        allowed (they are ignored in the calculations).
    column : str, optional
        Column name to use as the grouping vector if `grouping` is a
        ``DataFrame``. Must be provided if `grouping` is a ``DataFrame``.
        Cannot be provided if `grouping` is 1-D ``array_like``.
    permutations : int, optional
        Number of permutations to use when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and the p-value
        will be ``np.nan``.

    Returns
    -------
    pandas.Series
        Results of the statistical test, including ``test statistic`` and
        ``p-value``.

    See Also
    --------
    anosim

    Notes
    -----
    See [1]_ for the original method reference, as well as ``vegan::adonis``,
    available in R's vegan package [2]_.

    The p-value will be ``np.nan`` if `permutations` is zero.

    References
    ----------
    .. [1] Anderson, Marti J. "A new method for non-parametric multivariate
       analysis of variance." Austral Ecology 26.1 (2001): 32-46.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    Examples
    --------
    See :mod:`skbio.stats.distance.anosim` for usage examples (both functions
    provide similar interfaces).

    """
    sample_size, num_groups, grouping, tri_idxs, distances = _preprocess_input(
        distance_matrix, grouping, column)

    # Calculate number of objects in each group.
    group_sizes = np.bincount(grouping)
    s_T = (distances ** 2).sum() / sample_size

    test_stat_function = partial(_compute_f_stat, sample_size, num_groups,
                                 tri_idxs, distances, group_sizes, s_T)
    stat, p_value = _run_monte_carlo_stats(test_stat_function, grouping,
                                           permutations)

    return _build_results('PERMANOVA', 'pseudo-F', sample_size, num_groups,
                          stat, p_value, permutations)
###########################################




def bray_curtis_distance(table, sample1_id, sample2_id):
    numerator = 0
    denominator = 0
    sample1_counts = table[sample1_id]
    sample2_counts = table[sample2_id]
    for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
        numerator += abs(sample1_count - sample2_count)
        denominator += sample1_count + sample2_count
    return numerator / denominator
def table_to_distances(table, pairwise_distance_fn):
    sample_ids = table.columns
    num_samples = len(sample_ids)
    data = zeros((num_samples, num_samples))
    for i, sample1_id in enumerate(sample_ids):
        for j, sample2_id in enumerate(sample_ids[:i]):
            data[i,j] = data[j,i] = pairwise_distance_fn(table, sample1_id, sample2_id)
    return data, list(sample_ids)



#### 'ASV':
DATAFRAME_PHYLOGENETIC_ASVs = pd.read_csv('Anexos16S/ASVs_Phylogenetic_Tree.tsv', sep = '\t')
DATAFRAME_PHYLOGENETIC_ASVs = DATAFRAME_PHYLOGENETIC_ASVs.set_index('#OTU ID')
N_Cols_ASVs = list(DATAFRAME_PHYLOGENETIC_ASVs.columns)

if os.path.exists('Anexos16S/unifrac_matrix_16S_ASVs.txt') == True:
    frame_matrix_otus = pd.read_csv('Anexos16S/unifrac_matrix_16S_ASVs.txt', sep = '\t')
    id_samples_ASVs = frame_matrix_otus.Sample.tolist()
    unifrac_matrix_ASVs = frame_matrix_otus.iloc[:, 1:].values
else:
    A=urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/unifrac_matrix_16S_ASVs.txt',
                            'Anexos16S/unifrac_matrix_16S_ASVs.txt')

    frame_matrix_otus = pd.read_csv('Anexos16S/unifrac_matrix_16S_ASVs.txt', sep = '\t')
    id_samples_ASVs = frame_matrix_otus.Sample.tolist()
    unifrac_matrix_ASVs = frame_matrix_otus.iloc[:, 1:].values


#### 'OTU':
DATAFRAME_PHYLOGENETIC_OTUs = pd.read_csv('Anexos16S/OTUs_Phylogenetic_Tree.tsv', sep = '\t')
DATAFRAME_PHYLOGENETIC_OTUs = DATAFRAME_PHYLOGENETIC_OTUs.set_index('#OTU ID')
N_Cols_OTUs = list(DATAFRAME_PHYLOGENETIC_OTUs.columns)

if os.path.exists('Anexos16S/unifrac_matrix_16S_OTUs.txt') == True:
    frame_matrix_otus = pd.read_csv('Anexos16S/unifrac_matrix_16S_OTUs.txt', sep = '\t')
    id_samples_OTUs = frame_matrix_otus.Sample.tolist()
    unifrac_matrix_OTUs = frame_matrix_otus.iloc[:, 1:].values
else:
    A=urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/unifrac_matrix_16S_OTUs.txt',
                            'Anexos16S/unifrac_matrix_16S_OTUs.txt')

    frame_matrix_otus = pd.read_csv('Anexos16S/unifrac_matrix_16S_OTUs.txt', sep = '\t')
    id_samples_OTUs = frame_matrix_otus.Sample.tolist()
    unifrac_matrix_OTUs = frame_matrix_otus.iloc[:, 1:].values



def plot_pca3(dpca3 = 'ASV', variable = 'Coffee_Variety', CoLoR = 'tab10',
              mostrar_den_pca = True, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '3D', png = False, svg = False,
             enfoque = 'No Phylogenetic'):
    
    if enfoque == 'No Phylogenetic':

        if dpca3 == 'ASV':
            DATAFRAME_SUMARY = TAX_DFS_ASVs['Species']
        elif dpca3 == 'OTU':
            DATAFRAME_SUMARY = TAX_DFS_OTUs['Species']

        #print(dpca3, len(DATAFRAME_SUMARY))

        N_Cols = list(DATAFRAME_SUMARY.iloc[:, 1:].columns)
        ONE = DATAFRAME_SUMARY.set_index('Species')
        bc_matrix, id_samples = table_to_distances(ONE, bray_curtis_distance)

        SAMple = DataFrame(id_samples, columns = ['Sample']).merge(Sampledata, on = 'Sample', how = 'left')
        SAMple =SAMple.set_index('Sample')

        bc_matrix_condensed_ASVs = squareform(bc_matrix, force='tovector', checks=False)
        ZZZ = average(bc_matrix_condensed_ASVs)

        DfDf = DataFrame(bc_matrix, columns = id_samples)
        DfDf.insert(loc = 0, column='Sample', value= id_samples)
        DfDf = DfDf.merge(Sampledata, on = 'Sample', how = 'left')
        DfDf = DfDf.set_index('Sample')

        SS = StandardScaler()
        DfDf[id_samples] = SS.fit_transform(DfDf[id_samples].astype('float'))

        pca3 = PCA(n_components=3)

        try:
            pca_3 = pca3.fit_transform(DfDf[id_samples].astype('float'))
        except np.linalg.LinAlgError:
            pca_3 = pca3.fit_transform(DfDf[id_samples].astype('float'))

        DF_3 = pd.DataFrame({'PCA1':pca_3[:,0], 'PCA2':pca_3[:, 1], 'PCA3':pca_3[:, 2], 'clase':DfDf[variable]})


        Permanova = permanova(bc_matrix, SAMple, variable)['p-value']


        mpl.rcParams.update(mpl.rcParamsDefault)

        CorreS = {}
        for i, j in zip(DF_3.index, DF_3.clase):
            CorreS[i] = j

        XXX1=DF_3['PCA1']
        XXX2=DF_3['PCA2']
        XXX3=DF_3['PCA3']

        colorr = [matplotlib.colors.to_hex(i) for i in plt.get_cmap(CoLoR)(np.arange(len(set(DF_3.clase.tolist()))))]
        
        if variable == 'Time_Dry':
            category = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                category[j] = colorr[i]
            unicos = [str(x) for x in sorted([int(i) for i in category])]
            category = dict(zip(unicos, [category[j] for j in unicos]))
            marcadores = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                marcadores[j] = mark[i]
            unicos = [str(x) for x in sorted([int(i) for i in marcadores])]
            marcadores = dict(zip(unicos, [marcadores[j] for j in unicos]))
        elif variable == 'OTA':
            category = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                category[j] = colorr[i]
            unicos = [str(x) for x in sorted([float(i) for i in category])]
            category = dict(zip(unicos, [category[j] for j in unicos]))
            marcadores = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                marcadores[j] = mark[i]
            unicos = [str(x) for x in sorted([float(i) for i in marcadores])]
            marcadores = dict(zip(unicos, [marcadores[j] for j in unicos]))
        else:
            category = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                category[j] = colorr[i]
            marcadores = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                marcadores[j] = mark[i]

        sam_marker = {}
        for i, j in zip(DF_3.index, DF_3.clase):
            sam_marker[i] = [marcadores[j], category[j]]

        col = DF_3['clase'].map(category)

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        if seleccion == '3D':

            fig = plt.figure(figsize=(plot_alto_pca3, plot_alto_pca3-2))
            ax = fig.add_subplot(111, projection = '3d')



            for x, y, z, label, c  in zip(XXX1, XXX2, XXX3, N_Cols, col):
                ax.scatter(x, y, z, c= c, s = 80, marker = marcadores[CorreS[label]], alpha = 0.7,
                           edgecolors='white', linewidths = 1, zorder=2)

            ax.text(ax.get_xlim()[1] * 0.5, ax.get_ylim()[1], ax.get_zlim()[1], 'PERMANOVA: p='+str(Permanova)+'\n\n\n\n',
                    fontsize=9, ha='left', va = 'center')

            proxies = []
            labels = []
            for cat in category:
                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                         c=category[cat], marker=marcadores[cat], alpha = 0.7,
                                         markersize = 7.5,
                                         markeredgecolor=category[cat], linewidth = 0)
                proxies.append(proxy)
                labels.append(cat)


            ax.legend(proxies, list(category.keys()), numpoints=1, loc=2,borderpad = 0.2,
                      handletextpad=-0.2,prop={'size':9},
                              bbox_to_anchor=(1.01, 0.8))


            ax.set_xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_3pca, weight="bold")
            ax.set_ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_3pca, weight="bold")
            ax.set_zlabel('PCA 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)', size=LETRA_3pca, weight="bold")

            plt.xticks(size=9)
            plt.yticks(size=9)
            ax.zaxis.set_tick_params(labelsize=9)


            for cl in DF_3.clase.drop_duplicates():
                a = DF_3[DF_3.clase == cl].PCA1.values
                b = DF_3[DF_3.clase == cl].PCA2.values
                c = DF_3[DF_3.clase == cl].PCA3.values
                ax.scatter(np.mean(a), np.mean(b), np.mean(c), zorder=0, c= 'white', edgecolors = 'black', s = 30, marker = 'o', alpha = 1)
                for A, B, C in zip(a, b, c):
                    plt.plot([np.mean(a), A], [np.mean(b), B], [np.mean(c), C],
                             linestyle='--',
                            markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=1, alpha = 0.3)

            if mostrar_den_pca == True:
                alto = 0.7
                ax2 = plt.axes([1.05, 0.07, 0.2, alto - 0.038])
                de = dendrogram(ZZZ,orientation='left',
                                labels=N_Cols)
                Z = np.asarray(de['dcoord'], order='c')
                mh = max(Z[:, 2])

                n = 5
                for i in de['ivl']:
                    ax2.scatter(-0.12, n, s = 30, c = sam_marker[i][1], marker = sam_marker[i][0])
                    ax2.text(-0.3, n, i, color = 'black', fontsize=8, ha='left', va = 'center')
                    n += 10

                ax2.text(0.2, n + 2, variable, color = 'black', fontsize=9, ha='left', va = 'center', weight="bold")
                ax2.set_xlim((mh + (mh*0.15), -0.22))
                ax2.axis('off')

            if png == True:
                plt.savefig('Plots16S/PCA3_No_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
            if svg == True:
                plt.savefig('Plots16S/PCA3_No_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
            plt.show()

        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        if seleccion == '2D':

            show_den_2pca = True
            LETRA_2pca = 9
            plot_alto_pca2 = 5.5


            fig, ax = plt.subplots(figsize=(plot_alto_pca2, plot_alto_pca2-1))


            for xx, yy, label, c in zip(DF_3.PCA1, DF_3.PCA2, N_Cols, col):
                ax.scatter(xx, yy, c= c, s = 80, marker = marcadores[CorreS[label]], alpha = 0.7,
                          edgecolors='white', linewidths = 0.5, zorder=2)#, edgecolors='black'



            ax.text(ax.get_xlim()[1]*0.2, ax.get_ylim()[1], 'PERMANOVA: p='+str(Permanova)+'\n',
                     fontsize=9, ha='left', va = 'center')

            plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
            plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

            plt.xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_2pca, weight="bold")
            plt.ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_2pca, weight="bold")

            proxies = []
            labels = []
            for cat in category:
                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                         c=category[cat], marker=marcadores[cat], alpha = 0.7,
                                        markersize = 7.5,
                                         markeredgecolor=category[cat], linewidth = 0)
                proxies.append(proxy)
                labels.append(cat)

            ax.legend(proxies, labels, numpoints=1, loc=2,borderpad = 0.15,
                       handletextpad=-0.2,prop={'size':9},
                              bbox_to_anchor=(1.01, 1.02))


            POINTS = []
            for cl in DF_3.clase.drop_duplicates():
                points = []
                a = DF_3[DF_3.clase == cl].PCA1.values
                b = DF_3[DF_3.clase == cl].PCA2.values
                ax.scatter(np.mean(a), np.mean(b), zorder=0, c= 'white', edgecolors = 'black', s = 30, marker = 'o', alpha = 1)
                for A, B in zip(a, b):
                    points.append([A, B])
                    plt.plot([np.mean(a), A], [np.mean(b), B],
                             linestyle='--',
                            markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=1, alpha = 0.3)
                POINTS.append(points)

            if mostrar_den_pca == True:
                alto = 0.8
                ax2 = plt.axes([1.05, 0.07, 0.2, alto - 0.038])
                de = dendrogram(ZZZ,orientation='left',
                                labels=N_Cols)
                Z = np.asarray(de['dcoord'], order='c')
                mh = max(Z[:, 2])

                n = 5
                for i in de['ivl']:
                    ax2.scatter(-0.12, n, s = 30, c = sam_marker[i][1], marker = sam_marker[i][0])
                    ax2.text(-0.3, n, i, color = 'black', fontsize=8, ha='left', va = 'center')
                    n += 10

                ax2.text(0.2, n + 2, variable, color = 'black', fontsize=9, ha='left', va = 'center', weight="bold")
                ax2.set_xlim((mh + (mh*0.15), -0.22))
                ax2.axis('off')

            if png == True:
                plt.savefig('Plots16S/PCA2_No_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
            if svg == True:
                plt.savefig('Plots16S/PCA2_No_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
            plt.show()
            
        if seleccion == '3Di':
            
            if dpca3 == 'ASV':
                DATAFRAME_SUMARY = TAX_DFS_ASVs['Species']
            elif dpca3 == 'OTU':
                DATAFRAME_SUMARY = TAX_DFS_OTUs['Species']

            #print(dpca3, len(DATAFRAME_SUMARY))

            N_Cols = list(DATAFRAME_SUMARY.iloc[:, 1:].columns)
            ONE = DATAFRAME_SUMARY.set_index('Species')
            bc_matrix, id_samples = table_to_distances(ONE, bray_curtis_distance)


            DfDf = DataFrame(bc_matrix, columns = id_samples)
            DfDf.insert(loc = 0, column='Sample', value= id_samples)
            DfDf = DfDf.merge(Sampledata, on = 'Sample', how = 'left')
            DfDf = DfDf.set_index('Sample')

            SS = StandardScaler()
            DfDf[id_samples] = SS.fit_transform(DfDf[id_samples].astype('float'))

            pca3 = PCA(n_components=3)

            try:
                pca_3 = pca3.fit_transform(DfDf[id_samples].astype('float'))
            except np.linalg.LinAlgError:
                pca_3 = pca3.fit_transform(DfDf[id_samples].astype('float'))

            DF_3 = pd.DataFrame({'PCA1':pca_3[:,0], 'PCA2':pca_3[:, 1], 'PCA3':pca_3[:, 2], 'clase':DfDf[variable]})
            
            DF3D = DF_3
            DF3D['Size'] = 2
            DF3D[variable] = DF3D.clase
            DF3D['Sample'] = DF3D.index
            DF3D = DF3D.reset_index(drop = True)
            
            import plotly.express as px
            
            fig = px.scatter_3d(DF3D, x = 'PCA1',  
                                y = 'PCA2',  
                                z = 'PCA3', 
                                color = variable,  
                                size='Size',
                                text = 'Sample',
                                size_max = 10,  
                                opacity = 1,
                                color_discrete_sequence=colorr,
                               labels={'PCA1': 'PC 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)',
                        'PCA2': 'PC 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)',
                        'PCA3': 'PC 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)'})
            def add_trace_copy(trace):
                fig.add_traces(trace)
                new_trace = fig.data[-1]
                new_trace.update(textfont_color=trace.marker.color, textposition='top center', 
                                 mode="text", showlegend=False)
                trace.update(mode="markers")
            fig.for_each_trace(add_trace_copy)  
            fig.show()
            
    
    if enfoque == 'Phylogenetic':
        
        if dpca3 == 'ASV':
            UniFracDF = DataFrame(unifrac_matrix_ASVs, columns = id_samples_ASVs)
            UniFracDF.insert(loc = 0, column='Sample', value= id_samples_ASVs)
            UniFracDF = UniFracDF.merge(Sampledata, on = 'Sample', how = 'left')
            UniFracDF = UniFracDF.set_index('Sample')
            id_samples = id_samples_ASVs
            unifrac_matrix_condensed_ASVs = squareform(unifrac_matrix_ASVs, force='tovector', checks=False)
            WWW = average(unifrac_matrix_condensed_ASVs)
            DisMat2 = unifrac_matrix_ASVs
            
        elif dpca3 == 'OTU':
            UniFracDF = DataFrame(unifrac_matrix_OTUs, columns = id_samples_OTUs)
            UniFracDF.insert(loc = 0, column='Sample', value= id_samples_OTUs)
            UniFracDF = UniFracDF.merge(Sampledata, on = 'Sample', how = 'left')
            UniFracDF = UniFracDF.set_index('Sample')
            id_samples = id_samples_OTUs
            unifrac_matrix_condensed_OTUs = squareform(unifrac_matrix_OTUs, force='tovector', checks=False)
            WWW = average(unifrac_matrix_condensed_OTUs)
            DisMat2 = unifrac_matrix_OTUs
        

        SAMple = DataFrame(id_samples, columns = ['Sample']).merge(Sampledata, on = 'Sample', how = 'left')
        SAMple =SAMple.set_index('Sample')

        

        UU = StandardScaler()
        UniFracDF[id_samples] = UU.fit_transform(UniFracDF[id_samples].astype('float'))

        pca3 = PCA(n_components=3)

        try:
            pca_3 = pca3.fit_transform(UniFracDF[id_samples].astype('float'))
        except np.linalg.LinAlgError:
            pca_3 = pca3.fit_transform(UniFracDF[id_samples].astype('float'))

        DF_3 = pd.DataFrame({'PCA1':pca_3[:,0], 'PCA2':pca_3[:, 1], 'PCA3':pca_3[:, 2], 'clase':UniFracDF[variable]})


        Permanova = permanova(DisMat2, SAMple, variable)['p-value']

        CorreS = {}
        for i, j in zip(DF_3.index, DF_3.clase):
            CorreS[i] = j
        

        XXX1=DF_3['PCA1']
        XXX2=DF_3['PCA2']
        XXX3=DF_3['PCA3']

        colorr = [matplotlib.colors.to_hex(i) for i in plt.get_cmap(CoLoR)(np.arange(len(set(DF_3.clase.tolist()))))]
        
        if variable == 'Time_Dry':
            category = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                category[j] = colorr[i]
            unicos = [str(x) for x in sorted([int(i) for i in category])]
            category = dict(zip(unicos, [category[j] for j in unicos]))
            marcadores = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                marcadores[j] = mark[i]
            unicos = [str(x) for x in sorted([int(i) for i in marcadores])]
            marcadores = dict(zip(unicos, [marcadores[j] for j in unicos]))
        elif variable == 'OTA':
            category = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                category[j] = colorr[i]
            unicos = [str(x) for x in sorted([float(i) for i in category])]
            category = dict(zip(unicos, [category[j] for j in unicos]))
            marcadores = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                marcadores[j] = mark[i]
            unicos = [str(x) for x in sorted([float(i) for i in marcadores])]
            marcadores = dict(zip(unicos, [marcadores[j] for j in unicos]))
        else:
            category = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                category[j] = colorr[i]
            marcadores = {}
            for i, j in enumerate(set(DF_3.clase.tolist())):
                marcadores[j] = mark[i]

        sam_marker = {}
        for i, j in zip(DF_3.index, DF_3.clase):
            sam_marker[i] = [marcadores[j], category[j]]

        col = DF_3['clase'].map(category)
    
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        if seleccion == '3D':
            mpl.rcParams.update(mpl.rcParamsDefault)
            fig = plt.figure(figsize=(plot_alto_pca3, plot_alto_pca3-2))
            ax = fig.add_subplot(111, projection = '3d')



            for x, y, z, label, c  in zip(XXX1, XXX2, XXX3, id_samples, col):
                ax.scatter(x, y, z, c= c, s = 80, marker = marcadores[CorreS[label]], alpha = 0.7,
                           edgecolors='white', linewidths = 1, zorder=2)

            ax.text(ax.get_xlim()[1] * 0.5, ax.get_ylim()[1], ax.get_zlim()[1], 'PERMANOVA: p='+str(Permanova)+'\n\n\n\n',
                    fontsize=9, ha='left', va = 'center')

            proxies = []
            labels = []
            for cat in category:
                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                         c=category[cat], marker=marcadores[cat], alpha = 0.7,
                                         markersize = 7.5,
                                         markeredgecolor=category[cat], linewidth = 0)
                proxies.append(proxy)
                labels.append(cat)


            ax.legend(proxies, list(category.keys()), numpoints=1, loc=2,borderpad = 0.2,
                      handletextpad=-0.2,prop={'size':9},
                              bbox_to_anchor=(1.01, 0.8))


            ax.set_xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_3pca, weight="bold")
            ax.set_ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_3pca, weight="bold")
            ax.set_zlabel('PCA 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)', size=LETRA_3pca, weight="bold")

            plt.xticks(size=9)
            plt.yticks(size=9)
            ax.zaxis.set_tick_params(labelsize=9)


            for cl in DF_3.clase.drop_duplicates():
                a = DF_3[DF_3.clase == cl].PCA1.values
                b = DF_3[DF_3.clase == cl].PCA2.values
                c = DF_3[DF_3.clase == cl].PCA3.values
                ax.scatter(np.mean(a), np.mean(b), np.mean(c), zorder=0, c= 'white', edgecolors = 'black', s = 30, marker = 'o', alpha = 1)
                for A, B, C in zip(a, b, c):
                    plt.plot([np.mean(a), A], [np.mean(b), B], [np.mean(c), C],
                             linestyle='--',
                            markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=1, alpha = 0.3)

            if mostrar_den_pca == True:
                alto = 0.7
                ax2 = plt.axes([1.05, 0.07, 0.2, alto - 0.038])
                de = dendrogram(WWW,orientation='left',
                                labels=id_samples)
                Z = np.asarray(de['dcoord'], order='c')
                mh = max(Z[:, 2])

                n = 5
                for i in de['ivl']:
                    ax2.scatter(-0.07, n, s = 30, c = sam_marker[i][1], marker = sam_marker[i][0])
                    ax2.text(-0.15, n, i, color = 'black', fontsize=8, ha='left', va = 'center')
                    n += 10

                ax2.text(0.2, n + 2, variable, color = 'black', fontsize=9, ha='left', va = 'center', weight="bold")
                ax2.set_xlim((mh + (mh*0.15), -0.22))
                ax2.axis('off')

            if png == True:
                plt.savefig('Plots16S/PCA3_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
            if svg == True:
                plt.savefig('Plots16S/PCA3_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')

            plt.show()

        #%%%%%%%%%%%%%%%%%%

        if seleccion == '2D':

            LETRA_2pca = 9
            plot_alto_pca2 = 5.5


            fig, ax = plt.subplots(figsize=(plot_alto_pca2, plot_alto_pca2-1))


            for xx, yy, label, c in zip(DF_3.PCA1, DF_3.PCA2, id_samples, col):
                ax.scatter(xx, yy, c= c, s = 80, marker = marcadores[CorreS[label]], alpha = 0.7,
                          edgecolors='white', linewidths = 0.5, zorder=2)#, edgecolors='black'



            ax.text(ax.get_xlim()[1]*0.2, ax.get_ylim()[1], 'PERMANOVA: p='+str(Permanova)+'\n',
                     fontsize=9, ha='left', va = 'center')

            plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
            plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

            plt.xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_2pca, weight="bold")
            plt.ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_2pca, weight="bold")

            proxies = []
            labels = []
            for cat in category:
                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                         c=category[cat], marker=marcadores[cat], alpha = 0.7,
                                        markersize = 7.5,
                                         markeredgecolor=category[cat], linewidth = 0)
                proxies.append(proxy)
                labels.append(cat)

            ax.legend(proxies, labels, numpoints=1, loc=2,borderpad = 0.15,
                       handletextpad=-0.2,prop={'size':9},
                              bbox_to_anchor=(1.01, 1.02))



            POINTS = []
            for cl in DF_3.clase.drop_duplicates():
                points = []
                a = DF_3[DF_3.clase == cl].PCA1.values
                b = DF_3[DF_3.clase == cl].PCA2.values
                ax.scatter(np.mean(a), np.mean(b), zorder=0, c= 'white', edgecolors = 'black', s = 30, marker = 'o', alpha = 1)
                for A, B in zip(a, b):
                    points.append([A, B])
                    plt.plot([np.mean(a), A], [np.mean(b), B],
                             linestyle='--',
                            markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=1, alpha = 0.3)
                POINTS.append(points)

            if mostrar_den_pca == True:
                alto = 0.8
                ax2 = plt.axes([1.05, 0.07, 0.2, alto - 0.038])
                de = dendrogram(WWW,orientation='left',
                                labels=id_samples)
                Z = np.asarray(de['dcoord'], order='c')
                mh = max(Z[:, 2])

                n = 5
                for i in de['ivl']:
                    ax2.scatter(-0.07, n, s = 30, c = sam_marker[i][1], marker = sam_marker[i][0])
                    ax2.text(-0.15, n, i, color = 'black', fontsize=8, ha='left', va = 'center')
                    n += 10

                ax2.text(0.2, n + 2, variable, color = 'black', fontsize=9, ha='left', va = 'center', weight="bold")
                ax2.set_xlim((mh + (mh*0.15), -0.22))
                ax2.axis('off')

            if png == True:
                plt.savefig('Plots16S/PCA2_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
            if svg == True:
                plt.savefig('Plots16S/PCA2_Phylogenetic_'+dpca3+'_'+variable+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')

            plt.show()
        
        if seleccion == '3Di':
            
            if dpca3 == 'ASV':
                UniFracDF = DataFrame(unifrac_matrix_ASVs, columns = id_samples_ASVs)
                UniFracDF.insert(loc = 0, column='Sample', value= id_samples_ASVs)
                UniFracDF = UniFracDF.merge(Sampledata, on = 'Sample', how = 'left')
                UniFracDF = UniFracDF.set_index('Sample')
                id_samples = id_samples_ASVs

            elif dpca3 == 'OTU':
                UniFracDF = DataFrame(unifrac_matrix_OTUs, columns = id_samples_OTUs)
                UniFracDF.insert(loc = 0, column='Sample', value= id_samples_OTUs)
                UniFracDF = UniFracDF.merge(Sampledata, on = 'Sample', how = 'left')
                UniFracDF = UniFracDF.set_index('Sample')
                id_samples = id_samples_OTUs


            UU = StandardScaler()
            UniFracDF[id_samples] = UU.fit_transform(UniFracDF[id_samples].astype('float'))

            pca3 = PCA(n_components=3)

            try:
                pca_3 = pca3.fit_transform(UniFracDF[id_samples].astype('float'))
            except np.linalg.LinAlgError:
                pca_3 = pca3.fit_transform(UniFracDF[id_samples].astype('float'))

            DF_3 = pd.DataFrame({'PCA1':pca_3[:,0], 'PCA2':pca_3[:, 1], 'PCA3':pca_3[:, 2], 'clase':UniFracDF[variable]})
            
            DF3D = DF_3
            DF3D['Size'] = 2
            DF3D[variable] = DF3D.clase
            DF3D['Sample'] = DF3D.index
            DF3D = DF3D.reset_index(drop = True)
            
            import plotly.express as px
            
            fig = px.scatter_3d(DF3D, x = 'PCA1',  
                                y = 'PCA2',  
                                z = 'PCA3', 
                                color = variable,  
                                size='Size',
                                text = 'Sample',
                                size_max = 10,  
                                opacity = 1,
                                color_discrete_sequence=colorr,
                               labels={'PCA1': 'PC 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)',
                        'PCA2': 'PC 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)',
                        'PCA3': 'PC 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)'})
            def add_trace_copy(trace):
                fig.add_traces(trace)
                new_trace = fig.data[-1]
                new_trace.update(textfont_color=trace.marker.color, textposition='top center', 
                                 mode="text", showlegend=False)
                trace.update(mode="markers")
            fig.for_each_trace(add_trace_copy)  
            fig.show()


data_pca = widgets.ToggleButtons(options=['ASV', 'OTU'])
data_pca.layout.width = '10%'
#data.style.font_weight = 'bold'
data_pca.style.button_width = '60px'
boton_data_pca = Box(children=[data_pca], layout= Layout(border='1px solid pink', width='69px', height='63px'))


show_den_pca = widgets.ToggleButtons(options=[True, False])
show_den_pca.layout.width = '10%'
show_den_pca.style.button_width = '60px'
boton_den_pca = Box(children=[show_den_pca], layout= Layout(border='1px solid pink', width='69px', height='63px'))

#select_PCA = widgets.ToggleButtons(options=['2D', '3D', '3Di'])
#select_PCA.layout.width = '10%'
#select_PCA.style.button_width = '40px'
#boton_select_PCA = Box(children=[select_PCA], layout= Layout(border='1px solid pink', width='49px', height='63px'))


select_PCA = widgets.ToggleButtons(options=['2D', '3D', '3Di'])
select_PCA.style.button_width = '50px'
boton_select_PCA = Box(children=[select_PCA], layout= Layout(border='1px solid pink', width='169px', height='35px'))


select_beta = widgets.ToggleButtons(options=['No Phylogenetic', 'Phylogenetic'])
select_beta.layout.width = '10%'
select_beta.style.button_width = '120px'
boton_select_beta = Box(children=[select_beta], layout= Layout(border='1px solid pink', width='129px', height='63px'))



Rampas5 = widgets.Dropdown(options=list(QUALITATIVE_colors.keys()),value='tab10',
                          disabled=False,
                         layout=Layout(width='98px', height='25px'))

def boxrampas5(Rampas5):
    barcolor(lista = QUALITATIVE_colors[Rampas5])
OUTboxrampas5 = widgets.interactive_output(boxrampas5, {'Rampas5':Rampas5})



VariableS_pca = widgets.Dropdown(options=variables,value='Coffee_Variety',disabled=False,
                         layout=Layout(width='100px', height='25px'))



but_PCA3 = widgets.Button(description=" UPDATE PLOT ", icon = 'fa-refresh')
but_PCA3.style.button_color = 'lime' #'deepskyblue'
but_PCA3.style.font_weight = 'bold'
out_PCA3 = widgets.Output()

def button_PCA3(b):
    with out_PCA3:
        clear_output(True)
        if select_PCA.value == '3D':
            plot_pca3(dpca3 = data_pca.value, variable = VariableS_pca.value, CoLoR = Rampas5.value,
                  mostrar_den_pca = show_den_pca.value, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '3D',
                     png = False, svg = False, enfoque = select_beta.value)
        if select_PCA.value == '2D':
            plot_pca3(dpca3 = data_pca.value, variable = VariableS_pca.value, CoLoR = Rampas5.value,
                  mostrar_den_pca = show_den_pca.value, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '2D',
                     png = False, svg = False, enfoque = select_beta.value)
        if select_PCA.value == '3Di':
            plot_pca3(dpca3 = data_pca.value, variable = VariableS_pca.value, CoLoR = Rampas5.value,
                  mostrar_den_pca = show_den_pca.value, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '3Di',
                     png = False, svg = False, enfoque = select_beta.value)
                

but_PCA3.on_click(button_PCA3)


but_PCA3_save = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
but_PCA3_save.style.button_color = 'gold' #'deepskyblue'
out_PCA3_save = widgets.Output()

def button_PCA3_save(b):
    with out_PCA3_save:
        clear_output(True)
        if select_PCA.value == '3D':
            plot_pca3(dpca3 = data_pca.value, variable = VariableS_pca.value, CoLoR = Rampas5.value,
                  mostrar_den_pca = show_den_pca.value, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '3D',
                     png = True, svg = False, enfoque = select_beta.value)
        if select_PCA.value == '2D':
            plot_pca3(dpca3 = data_pca.value, variable = VariableS_pca.value, CoLoR = Rampas5.value,
                  mostrar_den_pca = show_den_pca.value, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '2D',
                     png = True, svg = False, enfoque = select_beta.value)
        
        


but_PCA3_save.on_click(button_PCA3_save)
#----------------------------------------------------------
but_PCA3_save2 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width='68px'))
but_PCA3_save2.style.button_color = 'gold' #'deepskyblue'
out_PCA3_save2 = widgets.Output()

def button_PCA3_save2(b):
    with out_PCA3_save2:
        clear_output(True)
        if select_PCA.value == '3D':
            plot_pca3(dpca3 = data_pca.value, variable = VariableS_pca.value, CoLoR = Rampas5.value,
                  mostrar_den_pca = show_den_pca.value, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '3D',
                     png = False, svg = True, enfoque = select_beta.value)
        if select_PCA.value == '2D':
            plot_pca3(dpca3 = data_pca.value, variable = VariableS_pca.value, CoLoR = Rampas5.value,
                  mostrar_den_pca = show_den_pca.value, LETRA_3pca = 9, plot_alto_pca3 = 7, seleccion = '2D',
                     png = False, svg = True, enfoque = select_beta.value)
        
        


but_PCA3_save2.on_click(button_PCA3_save2)




box_layout3333 = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='1px solid lime',
                    width='950px',
                   height='470px')



COMPONENTES = VBox([Box([HBox([VBox([widgets.Label('Clustering:'), boton_data_pca]),
                               VBox([widgets.Label('Beta Diversity:'), boton_select_beta]),
                               VBox([widgets.Label('Variable:'), VariableS_pca]),
                               VBox([widgets.Label('Dendrog.:'), boton_den_pca]),
                               VBox([widgets.Label('PCA:'), boton_select_PCA]),
                               VBox([widgets.Label('Box Color:'), Rampas5, OUTboxrampas5])
                               ])], layout = Layout(border='1px solid silver', width='950px')),
                    Box([HBox([but_PCA3, gris, gris, gris, gris, but_PCA3_save, but_PCA3_save2])], 
                        layout = Layout(border='1px solid silver', width='950px')),
                    HBox([Box([out_PCA3], layout=box_layout3333)])
                    ])



graficas = {'GLOBAL': GLOBAL_TAX,
            'TAXONOMY': TAXONOMI_IND,
            'INDEX':INDICES_16S,
            'RICHNESS':SIGNIFICATIVO,
            'COMPONENTS':COMPONENTES}
exe = widgets.ToggleButtons(options=list(graficas.keys()),disabled=False,button_style='')
data.layout.width = '1000%'
exe.style.button_width = '190px'


box_exe = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='5px solid deepskyblue',
                    width='984px',
                   height='45px')



def ss(exe):
    display(graficas[exe])
OU = widgets.interactive_output(ss, {'exe':exe})


RESULTS_16S = VBox([blanca, Box([exe], layout = box_exe), blanca, OU])
