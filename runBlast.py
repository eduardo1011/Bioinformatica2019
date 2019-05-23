import pandas as pd
from pandas.compat import StringIO
import subprocess
from datetime import datetime
import sys
import re
import os

print('\nCorriendo '+yyyyy+'\n')
# 
fasta1 = open(xxxxx,'r')
filename = xxxxx.split('.')[0]
name = yyyyy+'_'+filename.split('/')[-1]+'.tab'
#
print(yyyyy,'-db', zzzzz,'-query', xxxxx,'-evalue','1E-6','-outfmt\n',
"6 qacc sacc qlen slen length score bitscore evalue pident nident\n \
mismatch positive gaps gapopen stitle",'-max_target_seqs','10\n',
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
    blast = subprocess.call([yyyyy,'-db',zzzzz,'-query', identifier,'-evalue','1E-6','-outfmt',
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
        dif = max([i for i in range(1, out_bytes+1,1)]) - out_bytes
        total_length = int(out_bytes)
        dld += out_bytes
        dl = 0
        for dat in [i for i in range(1, out_bytes+1,1)]:
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
os.remove('run.py')
#for i in sin_resultados:
#    print(i)
