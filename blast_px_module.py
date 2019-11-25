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
import webbrowser
import subprocess
from datetime import datetime


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
        blast = re.sub('uuuuu',"'"+str(evalue)+"'", blast)
        blast = re.sub('#wwwww', 'print("")\nprint("Finalizado")\ntime.sleep(60)', blast)
        save= open('run.py','w')
        save.write(blast)
        save.close()
        os.system('gnome-terminal -e "python3 run.py"')
###
def run_blast_jupyter(x_p = '', file = '', db = '', evalue = 0.1, max_target_seqs = 0, max_hsps = 0):
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
        "6 qacc sacc qlen slen length qstart qend sstart send score bitscore evalue pident nident\n",
        "mismatch positive gaps gapopen stitle",'-max_target_seqs',str(max_target_seqs),'\n',
        '-max_hsps',str(max_hsps),'-out', name+'\n')
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
                                "6 qacc sacc qlen slen length qstart qend sstart send score bitscore evalue pident nident mismatch positive gaps gapopen stitle",
                                '-max_target_seqs',str(max_target_seqs),'-max_hsps',str(max_hsps),'-out', identifier+'.txt'])
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
            header = ('qacc','Entry','qlen','slen','length', 'qstart','qend','sstart','send','score','bitscore','evalue','pident','nident',
                              'mismatch','positive','gaps','gapopen','stitle')
            resultados_blast.columns = header
            resultados_blast.to_csv(name, sep = '\t',index = None)
            print('\n\nTiempo: {}'.format(tim).split('.')[0], '(h:m:s)')
            print('\nResultado: ', name, '\n')
            print('Secuencias sin resultados: '+str(len(sin_resultados)))
        #for i in sin_resultados:
        #    print(i)
