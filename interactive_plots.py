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
from matplotlib import cm
import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, interact_manual, Button, HBox, VBox, IntSlider, Label, IntRangeSlider
from ipywidgets import Checkbox, RadioButtons
import random
from matplotlib import cm


import funciones
from funciones import *


GO_Slim = pd.read_csv('GO_Slim.tsv', sep = '\t')
uniprotkb = pd.read_csv('uniprotkb.tsv', sep = '\t')



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
escala_uno = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
common_colors = ['red','orange','yellow','green','blue','purple','brown','magenta','tan',
                 'cyan','olive','maroon','navy','aquamarine','turquoise','silver',
                 'lime','teal','indigo','violet','pink','black','white']
common_colors.sort()
label_options = {'GO':0, 'Term':1, 'Short_Term':2}
layouts = {'Circular':nx.circular_layout,
          'Random':nx.random_layout,
          'Shell':nx.shell_layout,
          'Spectral':nx.spectral_layout,
          'Spring':nx.spring_layout}
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
v = widgets.Dropdown(options=sizes,description='size_plot:',value='5x5',disabled=False,
                    layout=Layout(width='95%', height='25px'))
#z = widgets.Checkbox(value=False,description='Save',disabled=False)
z = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                         tooltip='Description')
w = widgets.VBox([v, n, b, c, d, e, f, g, h, z])

def box(a, b, c, d, e, f, g, h, n, z, v):
    #print((a, b, c, d, e, f, g, h))
    frame1 = GO_Slim[columnas].groupby(['ASPECT']).get_group(a).drop_duplicates()
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
    frame_uniprot = frame1[['Entry', n]].merge(uniprotkb[['Entry', 'Protein_name', 'Gene']],
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
           size_bartxt = g)
    if z == True:
        plt.savefig('barras_'+a+'.png', dpi = 600, bbox_inches='tight')
    else:
        pass
out = widgets.interactive_output(box, {'a':a, 'b': b, 'c': c, 'd':d, 'e':e, 'f':f, 'g':g,
                                       'h':h, 'n':n, 'z':z, 'v':v})
def barras_interactive_plot():
    display(HBox([a]),HBox([w, out]))
    
#..........................................................................................
#..........................................................................................
#..........................................................................................

a1 = widgets.ToggleButtons(options=words,description='',disabled=False,button_style='')
b1 = widgets.Dropdown(options=list(range(0,101)),value=5,description='rows_num:',disabled=False,
                     layout=Layout(width='95%', height='25px'))
c1 = widgets.Dropdown(options=list(range(0,21)),value=10,description='size_text:',disabled=False,
                     layout=Layout(width='95%', height='25px'))
d1 = widgets.Dropdown(options=list(range(0,361)),value=0,description='angle:',disabled=False,
                     layout=Layout(width='95%', height='25px'))
e1 = widgets.Dropdown(options=list(pie_colors.keys()),value='Spectral',description='colors:',disabled=False,
                     layout=Layout(width='95%', height='25px'))
f1 = widgets.Dropdown(options=list(range(0,21)),value=10,description='size_title:',disabled=False,
                     layout=Layout(width='95%', height='25px'))
g1 = widgets.Dropdown(options=sizes,description='size_plot:',value='5x5',disabled=False,
                     layout=Layout(width='95%', height='25px'))
n1 = widgets.Dropdown(options=['GO', 'Term', 'Short_Term'],value='GO',description='label_option:',
                      disabled=False,layout=Layout(width='95%', height='25px'))

z1 = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                         tooltip='Description')
w1 = widgets.VBox([g1, n1, b1, c1, d1, f1, e1, z1])
def box1(a1, b1, c1, d1, e1, f1, g1, n1, z1):
    #print((a, b, c, d, e, f, g, h))
    frame1 = GO_Slim[columnas].groupby(['ASPECT']).get_group(a1).drop_duplicates()
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
    frame_uniprot = frame1[['Entry', n1]].merge(uniprotkb[['Entry', 'Protein_name', 'Gene']],
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
           size_text = c1,
           size_title = f1,
           color = pie_colors[e1])
    if z1 == True:
        plt.savefig('pastel_'+a1+'.png', dpi = 600, bbox_inches='tight')
    else:
        pass
out1 = widgets.interactive_output(box1, {'a1':a1, 'b1': b1, 'c1': c1, 'd1':d1, 'e1':e1,
                                         'f1':f1, 'g1':g1, 'n1':n1, 'z1':z1})
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

v2 = widgets.Dropdown(options=sizes,description='size_plot:',value='9x9',disabled=False,
                     layout=Layout(width='95%', height='25px'))

z2 = widgets.ToggleButton(value=False,description='Save',disabled=False,button_style='',
                         tooltip='Description')
w2 = widgets.VBox([v2, o2, f2, n2, m2, b2, c2, d2, e2, g2, h2, i2, j2, k2, l2, z2])

def box(a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2, l2, m2, n2, o2, z2, v2):
    #print((a, b, c, d, e, f, g, h))
    frame1 = GO_Slim[columnas].groupby(['ASPECT']).get_group(a2).drop_duplicates()
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
    frame_uniprot = frame1[['Entry', n2]].merge(uniprotkb[['Entry', 'Protein_name', 'Gene']],
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
             node_alpha = l2)
    if z2 == True:
        plt.savefig('net_'+a2+'.png', dpi = 600, bbox_inches='tight')
    else:
        pass
out2 = widgets.interactive_output(box, {'a2':a2, 'b2': b2, 'c2': c2, 'd2':d2, 'e2':e2, 'f2':f2,
                                        'g2':g2, 'h2':h2, 'i2':i2, 'j2':j2, 'k2':k2, 'l2':l2, 'm2':m2,
                                        'n2':n2, 'o2':o2, 'z2':z2, 'v2':v2, 'z2':z2})
def network_interactive_plot():
    display(HBox([a2]),HBox([w2, out2]))
