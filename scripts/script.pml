
python

import pymol
import re
import sys
import os

os.chdir(os.getcwd()+'/'+sys.argv[2]+'/pdb')

arg1 = 'starts.dat'
arg2 = 'folds.dat'

def buildDict(data):
    dataDict = {}
    for line in data:
        if (line.startswith('>')):
            pdb = line.replace('>','').replace('\n','')
            dataDict[pdb] = ''
            continue
        dataDict[pdb] += line.replace('\n','')
    return(dataDict)

def selectFold(foldDict):
    foldDict = re.split("\t|\n",foldDict)
    foldDict.remove('')
    res = ''
    for i in range(0,len(foldDict)-1,2):
        if (foldDict[i] == 'None'):
            continue
        s1 = int(foldDict[i]) + shift
        s2 = int(foldDict[i+1]) + shift
        res += '(resi '+str(s1)+':'+str(s2)+' & '+pdb+') + '
    return (res)

def colorFold(select,ss,color):   
    cmd.select(ss,select)
    cmd.color(color,ss)

with open(arg1,'r') as f:
    pdbStart = buildDict(f.readlines())

with open(arg2,'r') as f:
    foldDict = buildDict(f.readlines())

listPDB = pdbStart.keys()[:]
listPDB.sort()

for pdb in listPDB:
    if os.path.isfile(pdb+'.pdb'):
        cmd.load(pdb + '.pdb')

cmd.hide('all')
cmd.show('cartoon')

helices, strands = '', ''
for pdb,item in foldDict.items():
    cmd.color('green', pdb)
    shift = int(pdbStart[pdb])
    res1 =re.search('Helices: ',item)
    res2 =re.search('Strands: ',item)
    helices += selectFold(item[res1.end():res2.start()])
    strands += selectFold(item[res2.end():])

colorFold(helices[:-3],'helix','red')
colorFold(strands[:-3],'strand','yellow')

python end
