from subprocess import call
import shlex, subprocess
import subprocess

print('\n*****************   jupyter    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'jupyter'])
print('\n*****************   matplotlib    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'matplotlib'])
print('\n*****************   wordcloud    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'wordcloud'])
print('\n*****************   networkx    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'networkx'])
print('\n*****************   seaborn    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'seaborn'])
print('\n*****************   matplotlib-venn    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'matplotlib-venn'])
print('\n*****************   requests    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'requests'])
print('\n*****************   pandas    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'pandas'])
print('\n*****************   scipy    *****************\n')
subprocess.call(['python', '-mpip', 'install', 'scipy'])

print('\n**********   Proceso finalizado   **********\n')
