from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))

import pandas as pd
import numpy as np
import re
from functools import reduce
version = pd.__version__
if float(re.sub('[.]$', '', version[0:4])) >= 0.25:
    from io import StringIO
elif float(re.sub('[.]$', '', version[0:4])) < 0.25:
    from pd.compat import StringIO
from pandas import DataFrame
import warnings
warnings.filterwarnings("ignore")
from tabulate import tabulate
from IPython.display import clear_output, display 


import matplotlib as mpl
import ipywidgets as widgets
from ipywidgets import Layout, HBox, VBox, Box
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, OrderedDict
#import matplotlib.path as mpath
#import matplotlib.patches as mpatches

NCBI_RDP_SILVA_summary_ASVs = pd.read_csv('AnexosITS/ASVs_Taxonomy_Counts.tab', sep = '\t')
muestras_ASVs = {}
for i in list(NCBI_RDP_SILVA_summary_ASVs.columns[15:]):
    co = i
    df = NCBI_RDP_SILVA_summary_ASVs[['#OTU ID', co, 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    df = df[df[co] > 0]
    df = df.drop_duplicates()
    df['Per'] = ((np.array(df[co].tolist()) / np.sum(df[co].tolist())) * 100)
    muestras_ASVs[co] = df
NCBI_RDP_SILVA_summary_OTUs = pd.read_csv('AnexosITS/OTUs_Taxonomy_Counts.tab', sep = '\t')
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
valsota = []
for i in Sampledata.OTA:
    if float(i) < 5:
        valsota.append('Under 5')
    else:
        valsota.append('Above 5')
Sampledata['OTA'] = valsota

Sampledata = DataFrame(list(ordenado3.values()), columns = ['Sample']).merge(Sampledata, on = 'Sample', how = 'left')

variables = list(Sampledata.iloc[:, 1:].columns)


muestras_asv = TAX_DFS_ASVs['Species'].iloc[:, 1:].columns
especies_por_muestra_asv = {}
for c in muestras_asv:
    df = TAX_DFS_ASVs['Species'][['Species', c]]
    df = df[df[c] > 0].sort_values(by =c,ascending=False).reset_index(drop=True)
    especies_por_muestra_asv[c] = df

muestras_otu = TAX_DFS_OTUs['Species'].iloc[:, 1:].columns
especies_por_muestra_otu = {}
for c in muestras_otu:
    df = TAX_DFS_OTUs['Species'][['Species', c]]
    df = df[df[c] > 0].sort_values(by =c,ascending=False).reset_index(drop=True)
    especies_por_muestra_otu[c] = df


TAX_DFS_asv = TAX_DFS_ASVs
database_level_1_asv = {}
for t in TAX_DFS_asv.keys():
    sp = TAX_DFS_asv[t]
    sp['suma'] = np.sum(sp.iloc[:, 1:].values, axis = 1)
    sp = sp.sort_values(by ='suma',ascending=False).reset_index(drop=True)
    sp = sp.drop(columns = ['suma'])
    database_level_2_asv = {}
    for i in sp[t]:
        u = sp[sp[t] == i]
        u = u.set_index(t)
        x = u.T[u.T[i] > 0]
        x.insert(loc = 0, column='Sample', value=x.index) 
        x = x.reset_index(drop = True)
        y = x.merge(Sampledata, on = 'Sample', how = 'left').sort_values(by =i,ascending=False).reset_index(drop=True)
        database_level_2_asv[i] = y
    database_level_1_asv[t] = database_level_2_asv
XX_asv = NCBI_RDP_SILVA_summary_ASVs[['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].drop_duplicates()
XX_asv = TAX_DFS_asv['Species'][['Species']].merge(XX_asv, on = 'Species', how = 'left').drop_duplicates()
nivel_specie_1_asv = {}
for f, h in zip(['Phylum', 'Class', 'Order', 'Family', 'Genus'], ['Class', 'Order', 'Family', 'Genus', 'Species']):
    ef = XX_asv[[f, h]].drop_duplicates()
    nivel_specie_2_asv = {}
    for g in ef[f].tolist():
        rf = ef[ef[f] == g]
        if len(rf[h].tolist()):
            nivel_specie_2_asv[g] = rf[h].tolist()
        else:
            for v in rf[h].tolist():
                if v in list(database_level_1_asv[h].keys()):
                    nivel_specie_2_asv[g] = v
    nivel_specie_1_asv[f] = nivel_specie_2_asv



TAX_DFS_otu = TAX_DFS_OTUs
database_level_1_otu = {}
for t in TAX_DFS_otu.keys():
    sp = TAX_DFS_otu[t]
    sp['suma'] = np.sum(sp.iloc[:, 1:].values, axis = 1)
    sp = sp.sort_values(by ='suma',ascending=False).reset_index(drop=True)
    sp = sp.drop(columns = ['suma'])
    database_level_2_otu = {}
    for i in sp[t]:
        u = sp[sp[t] == i]
        u = u.set_index(t)
        x = u.T[u.T[i] > 0]
        x.insert(loc = 0, column='Sample', value=x.index) 
        x = x.reset_index(drop = True)
        y = x.merge(Sampledata, on = 'Sample', how = 'left').sort_values(by =i,ascending=False).reset_index(drop=True)
        database_level_2_otu[i] = y
    database_level_1_otu[t] = database_level_2_otu
XX_otu = NCBI_RDP_SILVA_summary_OTUs[['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].drop_duplicates()
XX_otu = TAX_DFS_otu['Species'][['Species']].merge(XX_otu, on = 'Species', how = 'left').drop_duplicates()
nivel_specie_1_otu = {}
for f, h in zip(['Phylum', 'Class', 'Order', 'Family', 'Genus'], ['Class', 'Order', 'Family', 'Genus', 'Species']):
    ef = XX_otu[[f, h]].drop_duplicates()
    nivel_specie_2_otu = {}
    for g in ef[f].tolist():
        rf = ef[ef[f] == g]
        if len(rf[h].tolist()):
            nivel_specie_2_otu[g] = rf[h].tolist()
        else:
            for v in rf[h].tolist():
                if v in list(database_level_1_otu[h].keys()):
                    nivel_specie_2_otu[g] = v
    nivel_specie_1_otu[f] = nivel_specie_2_otu



LinajE_asv = {}
for l in XX_asv.Species:
    df = XX_asv[XX_asv.Species == l].reset_index(drop=True)
    LinajE_asv[l] = '/'.join(list(df[['Phylum', 'Class', 'Order', 'Family', 'Genus']].values[0]))
    
LinajE_otu = {}
for l in XX_otu.Species:
    df = XX_otu[XX_otu.Species == l].reset_index(drop=True)
    LinajE_otu[l] = '/'.join(list(df[['Phylum', 'Class', 'Order', 'Family', 'Genus']].values[0]))


pie_colors = {"Arabica": '#1f77b4', "Robusta": '#aec7e8',
              "Wet": '#ff7f0e', "Dry": '#ffbb78',
              "DUCM":'blue', "DPS":'red',
              "Conventional": '#2ca02c', "Organic": '#98df8a',
              "3": '#d62728', "0": '#ff9896', "15": '#9467bd', "8": '#c5b0d5', "6": '#8c564b', "7": '#c49c94', "10": '#e377c2',
              "Above 5":"#de9ed6", "Under 5":"#ad494a"}
              #"0.0": '#de9ed6', "0.86": '#ce6dbd', "1.25": '#a55194', "0.73": '#7b4173', "15.83": '#e7969c', "3.19": '#d6616b', "0.7": '#ad494a',
    #"0.74": '#843c39', "0.59": '#e7cb94', "0.64": '#e7ba52', "3.11": '#bd9e39'}


def highlight_dps(x):
    return ['background-color: red' if re.search('^DPS$', str(y)) else 'black' for y in x]
def highlight_ducm(x):
    return ['background-color: blue' if re.search('^DUCM$', str(y)) else 'black' for y in x]
def parametros(s):
    return ['background-color: '+pie_colors[v] if v in list(pie_colors.keys()) else '' for v in s]
def color_negative(val):
    color = 'white' if val in ['DPS', 'DUCM'] else ''
    return 'color: %s' % color


from plotly.offline import iplot, init_notebook_mode
import plotly.graph_objects as go
from plotly.subplots import make_subplots


tipo = widgets.ToggleButtons(options=['ASV', 'OTU'])
tipo.style.button_width = '60px'
tipo.style.font_weight = 'bold'
boton_data = Box(children=[tipo], layout= Layout(border='1px solid pink', width='69px', height='63px'))

limite = widgets.Dropdown(options= range(0, 30050, 50), value = 100,
                         description='Threshold:', layout=Layout(width='200px', height='25px'))
EspecieS_asv = widgets.Dropdown(options= list(database_level_1_asv['Species'].keys()),
                         description='', layout=Layout(width='250px', height='28px'))
EspecieS_otu = widgets.Dropdown(options= list(database_level_1_otu['Species'].keys()),
                         description='', layout=Layout(width='250px', height='28px'))

bg_color = {'White':"plotly_white", 'Black':"plotly_dark"}
tema = widgets.SelectionSlider(options=list(bg_color.keys()), value='White', description='Theme:',
                               disabled=False, continuous_update=False, orientation='horizontal', readout=True)

centro = widgets.IntSlider(value=5, min=3, max=10, description='Center:', disabled=False, continuous_update=False,
                  orientation='horizontal', readout=True, readout_format='d', step=1)
numeros = {0:1,1:0.9,2:0.8,3:0.7,4:0.6,5:0.5,6:0.4,7:0.3,8:0.2,9:0.1,10:0}

all_samples = widgets.ToggleButtons(options=list(ordenado3.keys()))
all_samples.style.button_width = '40px'
all_samples.style.font_weight = 'bold'


import webbrowser
def info_pubmed(org = ''):
    org = re.sub('_', '+', org)
    url = 'https://pubmed.ncbi.nlm.nih.gov/?term='+org+'+coffee'
    webbrowser.open_new_tab(url)
    
def info_google(org = ''):
    org = re.sub('_', '+', org)
    url = 'https://www.google.com/search?client=firefox-b-d&q='+org+'+coffee'
    webbrowser.open_new_tab(url)
def show_frame(dict_db = {}, org = ''):
    z = dict_db['Species'][org]
    dps = z[z.Kit == 'DPS'].reset_index(drop=True)
    ducm = z[z.Kit == 'DUCM'].reset_index(drop=True)
    c = pd.concat([dps, ducm])
    c = c.reset_index()
    c['Reads'] = [int(i) for i in c[org]]
    c = c[['Kit', 'Sample', 'Coffee_Variety', 'Processing', 'Cultivation', 'Time_Dry', 'OTA', 'Reads']]
    u = c.style.bar(subset=['Reads'], align='mid', color='salmon').apply(parametros).applymap(color_negative, subset=['Kit'])
    return u
def linaje(EspecieS_asv, EspecieS_otu, tipo):
    if tipo == 'ASV':
        display(XX_asv[XX_asv.Species == EspecieS_asv][['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].drop_duplicates(subset = 'Species', keep = 'first').reset_index(drop = True))
    if tipo == 'OTU':
        display(XX_otu[XX_otu.Species == EspecieS_otu][['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].drop_duplicates(subset = 'Species', keep = 'first').reset_index(drop = True))
OUT2 = widgets.interactive_output(linaje, {'EspecieS_asv':EspecieS_asv, 'EspecieS_otu':EspecieS_otu, 'tipo':tipo})


mostratdf = widgets.Checkbox(value=False,description='Show Metadata?',disabled=False)
def frame(EspecieS_asv, EspecieS_otu, tipo, mostratdf):
    if tipo == 'ASV':
        if mostratdf == True:
            display(show_frame(dict_db = database_level_1_asv, org = EspecieS_asv))
        if mostratdf == False:
            print('')
    if tipo == 'OTU':
        if mostratdf == True:
            display(show_frame(dict_db = database_level_1_otu, org = EspecieS_otu))
        if mostratdf == False:
            print('')

OUT3 = widgets.interactive_output(frame, {'EspecieS_asv':EspecieS_asv, 'EspecieS_otu':EspecieS_otu, 'tipo':tipo, 'mostratdf':mostratdf})


    
    


def espe(EspecieS_asv, EspecieS_otu, tipo):
    if tipo == 'ASV':
        button = widgets.Button(description='Pubmed: '+re.sub('_', ' ', EspecieS_asv), layout=Layout(width='250px', height='25px'))
        button.style.button_color = 'moccasin'
        button.style.font_weight = 'bold'
        output = widgets.Output()
        def on_button_clicked(b):
            with output:
                info_pubmed(org = EspecieS_asv)
        button.on_click(on_button_clicked)

        button2 = widgets.Button(description='Google: '+re.sub('_', ' ', EspecieS_asv), layout=Layout(width='250px', height='25px'))
        button2.style.button_color = 'lightblue'
        button2.style.font_weight = 'bold'
        output2 = widgets.Output()
        def on_button_clicked2(b):
            with output2:
                info_google(org = EspecieS_asv)
        button2.on_click(on_button_clicked2)
        display(VBox([button, button2]))
    if tipo == 'OTU':
        button = widgets.Button(description='Pubmed: '+re.sub('_', ' ', EspecieS_otu), layout=Layout(width='250px', height='25px'))
        button.style.button_color = 'moccasin'
        button.style.font_weight = 'bold'
        output = widgets.Output()
        def on_button_clicked(b):
            with output:
                info_pubmed(org = EspecieS_otu)
        button.on_click(on_button_clicked)

        button2 = widgets.Button(description='Google: '+re.sub('_', ' ', EspecieS_otu), layout=Layout(width='250px', height='25px'))
        button2.style.button_color = 'lightblue'
        button2.style.font_weight = 'bold'
        output2 = widgets.Output()
        def on_button_clicked2(b):
            with output2:
                info_google(org = EspecieS_otu)
        button2.on_click(on_button_clicked2)
        display(VBox([button, button2]))
    
OUT = widgets.interactive_output(espe, {'EspecieS_asv':EspecieS_asv, 'EspecieS_otu':EspecieS_otu, 'tipo':tipo})


especie = list(database_level_1_asv['Species'].keys())[0]
ancho = 0.4
umbral = 100

z = database_level_1_asv['Species'][especie]
dps = z[z.Kit == 'DPS'].reset_index(drop=True)
dps = dps.sort_values(by = especie,ascending=True).reset_index(drop=True)
ducm = z[z.Kit == 'DUCM'].reset_index(drop=True)
ducm = ducm.sort_values(by = especie,ascending=True).reset_index(drop=True)


specs=[[{"rowspan": 3},{"rowspan": 3},{"rowspan": 3}, {'type':'xy'},{'type':'domain'}, {'type':'domain'}, {'type':'domain'}, {'type':'domain'}, {'type':'domain'}],
       [None, None, None, {'type':'xy'},{'type':'domain'}, {'type':'domain'}, {'type':'domain'}, {'type':'domain'}, {'type':'domain'}],
      [None, None, None, None, {}, {}, {}, {}, {}],
      [{"rowspan": 2, "colspan": 9}, None, None, None, None, None, None, None, None],
      [None, None, None, None, None, None, None, None, None]]

fig = make_subplots(rows=5, cols=9, specs=specs, column_widths = [6,8,8,1,4,4,4,4,4],
                    row_heights=[0.02, 0.02, 0.05, 0.05, 0.05], row_width=[0.01, 0.013, 0.02, 0.01, 0.01],
                    horizontal_spacing = 0.005, vertical_spacing=0.035,
                   subplot_titles = ['','<b>DPS</b>','<b>DUCM</b>'])
## Box
fig.add_trace(go.Box(y=np.array(dps[especie]), name="<b>DPS</b>", boxpoints='all', marker_color  = 'red'), 1, 1)
fig.add_trace(go.Box(y=np.array(ducm[especie]), name="<b>DUCM</b>", boxpoints='all', marker_color  = 'blue'), 1, 1)
fig.update_xaxes(tickangle=0, row = 1, col = 1, showline=True, linewidth=0.5)
fig.update_yaxes(tickangle=0, row = 1, col = 1, showline=True, linewidth=0.5)


# Scatter dps
fig.add_trace(go.Scatter(x=dps[dps[especie] < umbral][especie].tolist(), y=dps[dps[especie] < umbral].Sample.tolist(), name="<b>DPS</b>", #mode='markers',
                         line = dict(color='silver', width=1), # fill='toself',
                         mode="lines+markers",
                         #mode="lines+markers+text", text = ['<b>'+i+'</b>' for i in dps[dps[especie] < umbral].Sample.tolist()], textposition="middle right",
                        marker=dict(color='silver', size=np.log2(dps[dps[especie] < umbral][especie].tolist())*1.7, line=dict(width=1))), 1, 2)
fig.add_trace(go.Scatter(x=dps[dps[especie] >= umbral][especie].tolist(), y=dps[dps[especie] >= umbral].Sample.tolist(), name="<b>DPS</b>", #mode='markers',
                         line = dict(color='red', width=1), # fill='toself',
                         mode="lines+markers",
                         #mode="lines+markers+text", text = ['<b>'+i+'</b>' for i in dps[dps[especie] >= umbral].Sample.tolist()], textposition="middle right",
                        marker=dict(color='red', size=np.log2(dps[dps[especie] >= umbral][especie].tolist())*1.7, line=dict(width=1))), 1, 2)
fig.update_xaxes(tickangle=0, zeroline=True, row = 1, col = 2, showline=True, linewidth=0.5,
                 #range=[-np.array(dps[especie]).max()/5, np.array(dps[especie]).max()*1.5]
                )
fig.update_yaxes(tickangle=0, zeroline=True, showline=True, linewidth=0.5, showticklabels=False, row = 1, col = 2)

# Scatter ducm
fig.add_trace(go.Scatter(x=ducm[ducm[especie] < umbral][especie].tolist(), y=ducm[ducm[especie] < umbral].Sample.tolist(), name="<b>DUCM</b>",
                         line = dict(color='silver', width=1), #fill='toself',
                         mode="lines+markers",
                        #mode="lines+markers+text", text = ['<b>'+i+'</b>' for i in ducm[ducm[especie] < umbral].Sample.tolist()], textposition="middle right",
                        marker=dict(color='silver', size=np.log2(ducm[ducm[especie] < umbral][especie].tolist())*1.7, line=dict(width=1))), 1, 3)
fig.add_trace(go.Scatter(x=ducm[ducm[especie] >= umbral][especie].tolist(), y=ducm[ducm[especie] >= umbral].Sample.tolist(), name="<b>DUCM</b>",
                         line = dict(color='blue', width=1), #fill='toself',
                         mode="lines+markers",
                        #mode="lines+markers+text", text = ['<b>'+i+'</b>' for i in ducm[ducm[especie] >= umbral].Sample.tolist()], textposition="middle right",
                        marker=dict(color='blue', size=np.log2(ducm[ducm[especie] >= umbral][especie].tolist())*1.7, line=dict(width=1))), 1, 3)
fig.update_xaxes(tickangle=0, zeroline=True, row = 1, col = 3, showline=True, linewidth=0.5,
                 #range=[-np.array(ducm[especie]).max()/5, np.array(ducm[especie]).max()*1.5]
                )
fig.update_yaxes(tickangle=0, zeroline=True, row = 1, col = 3, showline=True, linewidth=0.5, showticklabels=False)




# Etiquetas
fig.add_annotation(dict(font=dict(size=13)), textangle = -90,
                   x=0, y=0, text='<b>DPS</b>', showarrow=False, yshift=0, row = 1, col = 4)
fig.update_xaxes(showgrid=False, zeroline=False, row = 1, col = 4, showticklabels=False)
fig.update_yaxes(showgrid=False, zeroline=False, row = 1, col = 4, showticklabels=False)

fig.add_annotation(dict(font=dict(size=13)), textangle = -90,
                   x=0, y=0, text="<b>DUCM</b>", showarrow=False, yshift=0, row = 2, col = 4)
fig.update_xaxes(showgrid=False, zeroline=False, row = 2, col = 4, showticklabels=False)
fig.update_yaxes(showgrid=False, zeroline=False, row = 2, col = 4, showticklabels=False)

# Pie charts
if len(dps[dps[especie] >= umbral]) > 0:
    columnas = 5
    for v in ['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA']:
        df = pd.pivot_table(dps[dps[especie] >= umbral][[especie, v]].drop_duplicates(), values=especie, index=[v], aggfunc=sum).reset_index()
        #print(df[v].tolist(), df[especie].tolist())
        fig.add_trace(go.Pie(labels=df[v].tolist(), values=df[especie].tolist(), hole=0.4, textinfo='none', name = '<b>'+v+'</b>', rotation = 0,
                        marker_colors = [pie_colors[k] for k in df[v].tolist()]), 1, columnas)
        columnas += 1
else:
    columnas = 5
    for v in ['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA']:
        fig.add_trace(go.Pie(labels=[''], values=[100], hole=ancho, textinfo='none', name = '<b>'+v+'</b>', rotation = 0,
                        marker_colors = ['silver']), 1, columnas)
        columnas += 1


        
if len(ducm[ducm[especie] >= umbral]) > 0:
    columnas = 5
    for v in ['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA']:
        df = pd.pivot_table(ducm[ducm[especie] >= umbral][[especie, v]].drop_duplicates(), values=especie, index=[v], aggfunc=sum).reset_index()
        #print(df[v].tolist(), df[especie].tolist())
        fig.add_trace(go.Pie(labels=df[v].tolist(), values=df[especie].tolist(), hole=0.4, textinfo='none', name = '<b>'+v+'</b>', rotation = 0,
                        marker_colors = [pie_colors[k] for k in df[v].tolist()]), 2, columnas)
        columnas += 1
else:
    columnas = 5
    for v in ['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA']:
        fig.add_trace(go.Pie(labels=[''], values=[100], hole=ancho, textinfo='none', name = '<b>'+v+'</b>', rotation = 0,
                        marker_colors = ['silver']), 2, columnas)
        columnas += 1

# Variables
columnas = 5
for v in ['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA']:
    fig.add_annotation(dict(font=dict(size=12)), 
                       x=0, y=0, text=v, showarrow=False, yshift=85, row = 3, col = columnas)
    fig.update_xaxes(showgrid=False, zeroline=False, row = 3, col = columnas, showticklabels=False)
    fig.update_yaxes(showgrid=False, zeroline=False, row = 3, col = columnas, showticklabels=False)
    columnas += 1
#
muestra_sel = '23A15'
df = especies_por_muestra_asv[muestra_sel]
df2 = df[df[muestra_sel] >= umbral]
pasa = len(df2)
no_pasa = len(df) - len(df2)
x = df.Species.tolist()
y = np.log2(df[muestra_sel].tolist())
colores = ['limegreen']*pasa + ['silver']*no_pasa
if especie in df.Species.tolist():
    for e, d in enumerate(df.Species.tolist()):
        if especie == d:
            posicion_sp = e
            break
    colores[posicion_sp] = 'purple'
else:
    pass



fig.add_trace(go.Bar(x=x, y=y, hovertext=[('Reads='+str(int(e))+' || '+LinajE_asv[i]) for e, i in zip(df[muestra_sel], df.Species)], marker_color=colores, name = '<b>'+muestra_sel+'</b>'), row = 4, col = 1)
fig.update_xaxes(tickangle=0, zeroline=True, showticklabels=False, showline=True, linewidth=0.5, row = 4, col = 1)
fig.update_yaxes(tickangle=0, zeroline=True, showline=True, linewidth=0.5, row = 4, col = 1)
fig.add_annotation(x=int(len(df)/2), y=np.log2(df2[muestra_sel][0])/2, text='<b>'+muestra_sel+'</b>', showarrow=False, font=dict(size=20), yshift=10, row = 4, col = 1)
    
    
fig.update_layout(autosize=True, showlegend=False, width=1000, height=620,margin=dict(l=1, r=1, b=1, t=50, pad=0),
                 font_size=11, template="plotly_white")
F = go.FigureWidget(fig)


def res(f):
    if tipo.value == 'ASV':
        database_level_1 = database_level_1_asv
        especie = EspecieS_asv.value
    if tipo.value == 'OTU':
        database_level_1 = database_level_1_otu
        especie = EspecieS_otu.value
    
    
    z = database_level_1['Species'][especie]
    dps = z[z.Kit == 'DPS'].reset_index(drop=True)
    dps = dps.sort_values(by = especie, ascending=True).reset_index(drop=True)
    ducm = z[z.Kit == 'DUCM'].reset_index(drop=True)
    ducm = ducm.sort_values(by = especie,ascending=True).reset_index(drop=True)
    
    with F.batch_update():
        ## Box
        F.data[0].y = np.array(dps[especie])
        F.data[1].y = np.array(ducm[especie])
        # Scatter dps
        F.data[2].x = dps[dps[especie] < limite.value][especie].tolist()
        F.data[2].y = dps[dps[especie] < limite.value].Sample.tolist()
        F.data[2].text = ['<b>'+i+'</b>' for i in dps[dps[especie] < limite.value].Sample.tolist()]
        F.data[2].marker.size = np.log2(dps[dps[especie] < limite.value][especie].tolist())*1.7
        #---
        F.data[3].x = dps[dps[especie] >= limite.value][especie].tolist()
        F.data[3].y = dps[dps[especie] >= limite.value].Sample.tolist()
        F.data[3].text = ['<b>'+i+'</b>' for i in dps[dps[especie] >= limite.value].Sample.tolist()]
        F.data[3].marker.size = np.log2(dps[dps[especie] >= limite.value][especie].tolist())*1.7
        # Scatter ducm
        F.data[4].x = ducm[ducm[especie] < limite.value][especie].tolist()
        F.data[4].y = ducm[ducm[especie] < limite.value].Sample.tolist()
        F.data[4].text = ['<b>'+i+'</b>' for i in ducm[ducm[especie] < limite.value].Sample.tolist()]
        F.data[4].marker.size = np.log2(ducm[ducm[especie] < limite.value][especie].tolist())*1.7
        #---
        F.data[5].x = ducm[ducm[especie] >= limite.value][especie].tolist()
        F.data[5].y = ducm[ducm[especie] >= limite.value].Sample.tolist()
        F.data[5].text = ['<b>'+i+'</b>' for i in ducm[ducm[especie] >= limite.value].Sample.tolist()]
        F.data[5].marker.size = np.log2(ducm[ducm[especie] >= limite.value][especie].tolist())*1.7
        ## Pie charts
        if len(dps[dps[especie] >= limite.value]) > 0:
            for v, p in zip(['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA'], [6, 7, 8, 9, 10]):
                df = pd.pivot_table(dps[dps[especie] >= limite.value][[especie, v]].drop_duplicates(), values=especie, index=[v], aggfunc=sum).reset_index()
                F.data[p].labels = df[v].tolist()
                F.data[p].values = df[especie].tolist()
                F.data[p].hole = numeros[centro.value]
                F.data[p].name = '<b>'+v+'</b>'
                F.data[p].marker.colors = [pie_colors[k] for k in df[v].tolist()]
        else:
            for v, p in zip(['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA'], [6, 7, 8, 9, 10]):
                F.data[p].labels = ['']
                F.data[p].values = [1]
                F.data[p].hole = numeros[centro.value]
                F.data[p].name = '<b>'+v+'</b>'
                F.data[p].marker.colors = ['silver']
        ## Pie charts
        if len(ducm[ducm[especie] >= limite.value]) > 0:
            for v, p in zip(['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA'], [11, 12, 13, 14, 15]):
                df = pd.pivot_table(ducm[ducm[especie] >= limite.value][[especie, v]].drop_duplicates(), values=especie, index=[v], aggfunc=sum).reset_index()
                F.data[p].labels = df[v].tolist()
                F.data[p].values = df[especie].tolist()
                F.data[p].hole = numeros[centro.value]
                F.data[p].name = '<b>'+v+'</b>'
                F.data[p].marker.colors = [pie_colors[k] for k in df[v].tolist()]
        else:
            for v, p in zip(['Coffee_Variety','Processing','Cultivation','Time_Dry','OTA'], [11, 12, 13, 14, 15]):
                F.data[p].labels = ['']
                F.data[p].values = [1]
                F.data[p].hole = numeros[centro.value]
                F.data[p].name = '<b>'+v+'</b>'
                F.data[p].marker.colors = ['silver']
        
        
        
        
        if tipo.value == 'ASV':
            muestra_sel = ordenado3[all_samples.value]
            dff = especies_por_muestra_asv[muestra_sel]
            
            dff2 = dff[dff[muestra_sel] >= limite.value]
            pasa = len(dff2)
            no_pasa = len(dff) - len(dff2)
            x = dff.Species.tolist()
            y = np.log2(dff[muestra_sel].tolist())
            colores = ['limegreen']*pasa + ['silver']*no_pasa

            if especie in dff.Species.tolist():
                for e, d in enumerate(dff.Species.tolist()):
                    if especie == d:
                        posicion_sp = e
                        break
                colores[posicion_sp] = 'purple'
            else:
                pass

            F.data[16].x = x
            F.data[16].y = y
            F.data[16].hovertext = [('Reads='+str(int(e))+' || '+LinajE_asv[i]) for e, i in zip(dff[muestra_sel], dff.Species)]

            F.data[16].marker.color = colores
            F.data[16].name = '<b>'+muestra_sel+'</b>'

            F.layout.annotations[-1].text = '<b>'+muestra_sel+'</b>'
            F.layout.annotations[-1].x = int(len(dff)/2)
            F.layout.annotations[-1].y = np.log2(dff2[muestra_sel][0])/2
            
            
        if tipo.value == 'OTU':
            muestra_sel = ordenado3[all_samples.value]
            dff = especies_por_muestra_otu[muestra_sel]
            
            dff2 = dff[dff[muestra_sel] >= limite.value]
            pasa = len(dff2)
            no_pasa = len(dff) - len(dff2)
            x = dff.Species.tolist()
            y = np.log2(dff[muestra_sel].tolist())
            colores = ['limegreen']*pasa + ['silver']*no_pasa

            if especie in dff.Species.tolist():
                for e, d in enumerate(dff.Species.tolist()):
                    if especie == d:
                        posicion_sp = e
                        break
                colores[posicion_sp] = 'purple'
            else:
                pass

            F.data[16].x = x
            F.data[16].y = y
            F.data[16].hovertext = [('Reads='+str(int(e))+' || '+LinajE_otu[i]) for e, i in zip(dff[muestra_sel], dff.Species)]

            F.data[16].marker.color = colores
            F.data[16].name = '<b>'+muestra_sel+'</b>'

            F.layout.annotations[-1].text = '<b>'+muestra_sel+'</b>'
            F.layout.annotations[-1].x = int(len(dff)/2)
            F.layout.annotations[-1].y = np.log2(dff2[muestra_sel][0])/2
        
        F.layout.template = bg_color[tema.value]


EspecieS_asv.observe(res, names="value")
EspecieS_otu.observe(res, names="value")
limite.observe(res, names="value")
tema.observe(res, names="value")
centro.observe(res, names="value")
tipo.observe(res, names="value")
all_samples.observe(res, names="value")

interactiveITS = VBox([HBox([widgets.Label('Samples:'), all_samples]), HBox([VBox([HBox([boton_data, VBox([EspecieS_asv, EspecieS_otu])]), limite, tema, centro, OUT, mostratdf]), VBox([F, OUT3])])])