#! /usr/bin/env python3

import sys
import re
import os
from subprocess import call
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa

# COLORS
W = '\33[1;39m'
R = '\33[1;31m'
G = '\33[1;32m'
Y = '\33[33m'
B = '\33[1;94m'
WH = '\33[7m'
RH = '\33[41m'
GH = '\33[42m'
YH = '\33[43m'
BH = '\33[44m'
NC = '\33[0m'

# CATH SUPERFAMILIES
#################################################################################################################

FAMILY_1  = '2.30.29.30.4_2_0'
MODEL_1 = '1dynA00'
SCAFFOLD_1 = [ {'name':'B1','type':'Strand','number':  0}, {'name':'B2','type':'Strand','number':  2},
               {'name':'B3','type':'Strand','number':  3}, {'name':'B4','type':'Strand','number':  4},
               {'name':'B5','type':'Strand','number':  5}, {'name':'B6','type':'Strand','number':  6},
               {'name':'B7','type':'Strand','number':  7}, {'name':'A','type':'Helix','number':    8} ]

FAMILY_2  = '3.20.20.190.4_1_0'
MODEL_2 = '2plcA00'
SCAFFOLD_2 = [ { 'name':'B1','type':'Strand','number': 0 }, { 'name':'A1','type':'Helix','number':  3 },
               { 'name':'B2','type':'Strand','number': 4 }, { 'name':'A2','type':'Helix','number':  7 },
               { 'name':'B3','type':'Strand','number': 8 }, { 'name':'A3','type':'Helix','number':  9 },
               { 'name':'B4','type':'Strand','number': 12 } ]


print(B+'\n Which CATH superfamily would you like to study ?')
print('1. '+FAMILY_1+' (model:',MODEL_1+')')
print('2. '+FAMILY_2+' (model:',MODEL_2+')')

while(True):
    choice =  input(R+'\nPlease, choose a menu item : '+NC)
    if (choice == '1'):
        FAMILY = FAMILY_1
        MODEL = MODEL_1
        SCAFFOLD = SCAFFOLD_1
        print()
        break
    elif (choice == '2'):
        FAMILY = FAMILY_2
        MODEL = MODEL_2
        SCAFFOLD = SCAFFOLD_2
        print()
        break

#################################################################################################################

# PATHS
PATH = 'superfamilies/'+FAMILY+'/'
CATH = PATH+'fasta/cath/'
NEW = PATH+'fasta/new/'
PDB = PATH+'pdb/'
DSSP = PATH+'dssp/'
TABLES = PATH+'tables/'

# LIST OF AA
AA = ['R','K','D','E','Q','N','H','S','T','Y','C','W','A','I','L','M','F','V','P','G']

# FUNCTIONS

def createDirectory(DIR):
    if not os.path.isdir(DIR):
        call(['mkdir '+DIR], shell=True)

def buildDict(fasta):
    '''
    Transformation of a fasta into a dictionary
    '''
    fastaDict = {}
    for line in fasta:
        if (line.startswith('>')):
            pdb = line.replace('>','').replace('\n','')
            fastaDict[pdb] = ''
            continue
        fastaDict[pdb] += line.replace('\n','')
    return(fastaDict)

def fromAlnToSeq(filename):
    '''
    Writes a fasta sequence from a fasta alignment
    Returns the dictionary corresponding to the new fasta file 
    '''
    if not os.path.isfile(filename):
        with open(filename,'w') as f:
            for pdb in ALN.keys():
                f.write('>'+pdb+'\n')
                for i in ALN[pdb]:
                    if (i != '-'):
                        f.write(i)
                f.write('\n')
            print(R+' Sequences retrieved from Alignment file :'+G+' "'+f.name+'" created.\n'+NC)
    with open(filename,'r') as f:
        SEQ = buildDict(f.readlines())
    return(SEQ)

def retrieveSS(filename,dico):
    '''
    Retrieves the secondary structures from pdb files and creates a fasta file containing these informations
    Returns the dictionary corresponding
    '''
    if not os.path.isfile(filename):
        print(R+' Secondary structures retrieved from pdb files using dssp :'+NC)
        createDirectory(DSSP)
        with open(filename,'w') as filout:
            for pdb in dico.keys():
                if not os.path.isfile(PDB+'/'+pdb+'.pdb'):
                    continue
                dsspFile = DSSP+pdb+'.dssp'
                if not os.path.isfile(dsspFile):
                    call(['softwares/mkdssp -i'+PDB+pdb+'.pdb -o'+dsspFile],shell=True)
                    print(G+'   "'+dsspFile+'" created.'+NC)
                with open(DSSP+pdb+'.dssp','r') as filin:
                    data = filin.read()
                res = re.search('#',data)
                data = data[res.start():].split('\n')[1:-1]
                filout.write('>'+pdb+'\n')
                for line in data:
                    filout.write(line[16].replace(' ','C').replace('G','H').replace('I','H').replace('B','C').replace('S','C'))
                filout.write('\n')
            print(R+'\n Summary of results: '+G+'"'+filout.name+'" created.\n'+NC)
    with open(filename,'r') as f:
        SEQ_SS = buildDict(f.readlines())  
    return(SEQ_SS)

def alignSS(filename):
    '''
    Aligment of the secondary structures according to the fasta alignment
    Returns the dictionary corresponding
    '''
    if not os.path.isfile(filename):
        with open(filename,'w') as f:
            for pdb in ALN.keys():
                if (pdb not in SEQ_SS.keys()):
                    continue
                f.write('>'+pdb+'\n')
                ss = 0
                for i in ALN[pdb]:
                    if (i == '-'):
                        f.write(i)
                    else:
                        f.write(SEQ_SS[pdb][ss])
                        ss += 1
                f.write('\n')
            print(R+' Alignment of ss :'+G+' "'+f.name+'" created.\n'+NC)
    with open(filename,'r') as f:
        ALN_SS = buildDict(f.readlines())
    return(ALN_SS)

def buildFoldDict():
    '''
    Creates a dictionary containing the fold of the model 
    '''
    foldDict = {}
    pattn = re.compile('(E(|-*)){1,}E|(H(|-*)){1,}H')
    seq = ALN_SS[MODEL]
    ssNb, k = 0, 0
    for ss in re.finditer(pattn,seq):
        if ( k <len(SCAFFOLD) and SCAFFOLD[k]['number'] == ssNb):
            newFold = {'start':ss.start(),'stop':ss.end()}
            if (ss.group(0)[0] == 'H'):
                newFold['type'] = 'Helix'
            elif (ss.group(0)[0] == 'E'):
                newFold['type'] = 'Strand'
            foldDict[SCAFFOLD[k]['name']] = newFold    
            k += 1
        ssNb += 1
    return(foldDict)

def buildSSDict():
    '''
    Creates a dictionary of all ss of all proteins
    '''
    SSDict = {}
    alnPattn = re.compile('(E(|-*)){1,}E|(H(|-*)){1,}H|(T(|-*)){1,}T|T')
    fastaPattn = re.compile('E{2,}|H{2,}|T{1,}')

    for pdb in ALN.keys():

        if pdb not in SEQ_SS:
            SSDict[pdb] = SSDict[MODEL]
            continue
        k,startLoop = 0,0
        find = False

        SSList = []
        newLoop = {'site':'loop','start_aln':0,'start':0,'SS':[]}

        alnSeq = ALN_SS[pdb]
        fastaSeq = SEQ_SS[pdb]

        ssAln = list(re.finditer(alnPattn,alnSeq))
        ssFasta = list(re.finditer(fastaPattn,fastaSeq))

        for N in range(len(ssAln)):
            if (ssAln[N].group(0)[0] == 'H'):
                SStype = 'Helix'
            elif (ssAln[N].group(0)[0] == 'E'):
                SStype = 'Strand'
            elif (ssAln[N].group(0)[0] == 'T'):
                SStype = 'Turn'

            if (k < len(SCAFFOLD) and SStype != 'Turn'):
                for s in range(ssAln[N].start(),ssAln[N].end()):
                    if (find == True):
                        break
                    for f in range(k,len(foldDict)):
                        fold = foldDict[SCAFFOLD[f]['name']]
                        foldSeq = [i for i in range(fold['start'],fold['stop'])]   
                        if (s in foldSeq and SStype == SCAFFOLD[f]['type']):
                            find = True
                            break

            if (find == True):
                newLoop['stop_aln'] = ssAln[N].start()
                newLoop['stop'] = ssFasta[N].start()
                SSList.append(newLoop)
                newLoop = {'site':'loop','start_aln':ssAln[N].end(),'start':ssFasta[N].end(),'SS':[]}

                newFold = {'site':'fold','type':SStype,'name':SCAFFOLD[f]['name'],
                'start_aln':ssAln[N].start(),'stop_aln':ssAln[N].end(),
                'start':ssFasta[N].start(),'stop':ssFasta[N].end()}
                SSList.append(newFold)

                k += 1
                find = False
            else:
                newLoop['SS'].append({'type':SStype,'start':ssFasta[N].start(),'stop':ssFasta[N].end()})

        newLoop['stop_aln'] = len(ALN_SS[pdb])
        newLoop['stop'] = len(SEQ_SS[pdb])

        SSList.append(newLoop)

        SSDict[pdb] = SSList

    return(SSDict)

def colorSS(choice):
    '''
    Prints and colors the alignment by secondary structures in the terminal
    Red: Helices, Yellow: Strands, Uncolored: Coils
    '''
    PDBlist = ALN.keys()
    if (choice != 'all'):
        PDBlist = choice

    for pdb in PDBlist:
        if (pdb not in ALN.keys()):
            print(R+'\n ERROR : '+NC+pdb+R+' not in the alignment.'+W)
            continue
        print(WH+'\n'+pdb+'\n'+NC,end='')
        st = 0
        for i in range(0,len(SSDict[pdb])):
            ss = SSDict[pdb][i]
            print(ALN[pdb][st:ss['start_aln']],end='')
            if (ss['site'] == 'fold'):
                if (ss['type'] == 'Helix'):
                    COL = RH
                elif (ss['type'] == 'Strand'):
                    COL = YH
            else:
                COL = NC
            print(COL+ALN[pdb][ss['start_aln']:ss['stop_aln']]+NC,end='')
            st = ss['stop_aln']
        print(ALN[pdb][st:len(ALN[pdb])]+'\n',end='')
    print(W+'\n  Fold elements : '+NC+RH+'HELIX'+NC+', '+YH+'STRAND'+NC+'\n')

def pdbStart():
    '''
    Creates a file containing the start from each pdb files
    Necessary for PyMOL
    '''
    filename = PDB+'starts.dat'
    with open(filename,'w') as filout:
        for pdb in ALN.keys():
            if not os.path.isfile(PDB+'/'+pdb+'.pdb'):
                continue
            filout.write('>'+pdb+'\n')
            with open(PDB+pdb+'.pdb','r') as filin:
                data = filin.readlines()
            for line in data:
                if (re.search('ATOM',line)):
                    filout.write(str(int(line[22:26])-1))
                    break
            filout.write('\n')
        print(G+' "'+filout.name+'" created.\n'+NC)

def retrieveFolds():
    '''
    Creates a file containing the start and the stop of folds
    Necessary for PyMOL
    '''
    folds = {}
    for pdb in SSDict.keys():
        if not os.path.isfile(PDB+'/'+pdb+'.pdb'):
            continue
        if (pdb not in folds.keys()):
            folds[pdb] = {}
            folds[pdb]['Helices'] = ''
            folds[pdb]['Strands'] = ''
        for i in SSDict[pdb]:
            if (i['site'] == 'fold'):
                if (i['type'] == 'Helix'):
                    folds[pdb]['Helices'] += str(i['start']+1)+'\t'+str(i['stop'])+'\t'
                elif (i['type'] == 'Strand'):
                    folds[pdb]['Strands'] += str(i['start']+1)+'\t'+str(i['stop'])+'\t'

    filename = PDB+'folds.dat'
    with open(filename,'w') as filout:
        for pdb in folds.keys():
            filout.write('>'+pdb+'\n')
            for ss,item in folds[pdb].items(): 
                filout.write(ss+': '+item+'\n')
        print(G+' "'+filout.name+'" created.'+NC)

def buildTable(loop):
    '''
    Data frame of amino acids for each pdb
    '''
    df = pd.DataFrame(index=loop.keys(), columns=AA, data=0)
    for pdb in loop.keys():
        for i in loop[pdb]:
            if (i != '-'):
                df[i][pdb] += 1
    df.index.name = 'pdb'
    return(df)

def createCSV(tableLoop,filename):
    if not os.path.isfile(filename):
        tableLoop.to_csv(filename)
        print(G+'\n "'+filename+'" created.\n'+NC)

def studyFold(foldname):
    '''
    For a chosen fold, it creates a csv file containing the number of each aa for each protein
    '''
    loopSeq = {}
    for pdb,domains in SSDict.items():
        for i in range(len(domains)):
            if (domains[i]['site'] == 'fold' and domains[i]['name'] == foldname):
                loopSeq[pdb] = ALN[pdb][domains[i]['start_aln']:domains[i]['stop_aln']].replace('-','')
                break
    filename = TABLES+'fold_'+foldname+'_aa.csv'
    createCSV(buildTable(loopSeq),filename)

def studyLoop(ele1,ele2):
    '''
    For a chosen loop (between two folds), it creates a csv file containing the number of each aa for each protein
    '''
    loop = {}
    for pdb,domains in SSDict.items():
        if (ele1 == 'N' and domains[1]['site'] == 'fold' 
            and domains[1]['name'] == ele2):
            loop[pdb] = domains[0]
        elif (ele2 == 'C' and  domains[-2]['site'] == 'fold' 
            and domains[-2]['name'] == ele1):
            loop[pdb] = domains[-1]
        elif (ele1 != 'N' and ele2 != 'C'):
            for i in range(1,len(domains)-1):
                if (domains[i]['site'] == 'fold'):
                    continue
                if (domains[i-1]['name'] == ele1 and domains[i+1]['name'] == ele2):
                    loop[pdb] = domains[i]
                    break

    loopSeq, loopSS = {}, {}
    for pdb,l in loop.items():
        loopSeq[pdb] = ALN[pdb][l['start_aln']:l['stop_aln']].replace('-','')
    filename = TABLES+'loop_'+ele1+'_'+ele2+'_aa.csv'
    createCSV(buildTable(loopSeq),filename)

def addSeqAndAlign():
    '''
    Adds sequences to the structural alignment 
    '''
    call(['softwares/hmmbuild --amino '+CATH+'hmmfile.out '+CATH+'alignment.fasta'+' > '+CATH+'hmmbuild.log'],shell=True)
    print(R+' hmm profil :'+G+' "'+CATH+'hmmfile.out" and '+CATH+'"hmmbuild.log'+'" created.\n'+NC)
    call(['softwares/hmmalign --amino --informat FASTA --outformat SELEX -o '+NEW+'alignment.slx '+CATH+'hmmfile.out '+NEW+'seq_to_align.fasta'],shell=True)
    call(['softwares/seqret -supper -sformat SELEX -osformat fasta '+NEW+'alignment.slx '+NEW+'alignment.fasta'],shell=True)
    print(R+' New alignment generated :'+G+' "'+NEW+'alignment.slx" and '+NEW+'"alignment.fasta'+'" created.\n'+NC)
    with open(NEW+'alignment.fasta','r') as f:
        return(buildDict(f.readlines()))

def delete():
    call(['rm -fr '+DSSP], shell=True)
    call(['rm -f '+PDB+'starts.dat'], shell=True)
    call(['rm -f '+PDB+'folds.dat'], shell=True)
    call(['rm -f '+CATH +'alignment_ss.fasta'], shell=True)
    call(['rm -f '+CATH +'sequences.fasta'], shell=True)
    call(['rm -f '+CATH+'hmm*'], shell=True)
    call(['rm -f '+NEW+'alignment_ss.fasta'], shell=True)
    call(['rm -f '+NEW+'sequences.fasta'], shell=True)
    call(['rm -f '+NEW+'alignment.slx'], shell=True)

    
def menu():
    print(B+' MENU :')
    print('1. Display the structural alignment colored by ss')
    print('2. Display the structural alignment in PyMOL')
    print('3. Study of a fold/loop')
    print('4. Add sequences to the cath alignment')
    print('5. Generate a plot (in /tables_old)')
    print('6. Delete unused files and Quit')

# MAIN

if not os.path.isfile(CATH+'alignment.fasta'):
    print(R+' ERROR : '+NC+CATH+'alignment.fasta'+R+' file does not exist.')
    exit()
with open(CATH+'alignment.fasta','r') as f:
    ALN = buildDict(f.readlines())

SEQ = fromAlnToSeq(CATH+'sequences.fasta')
SEQ_SS = retrieveSS(CATH+'sequences_ss.fasta',ALN)
ALN_SS = alignSS(CATH+'alignment_ss.fasta')

foldDict = buildFoldDict()
eleList = [i for i in foldDict.keys()]
eleList.insert(0,'N')
eleList.append('C')
SSDict = buildSSDict()

while (True):
    menu()
    choice1 = input(R+'\nPlease, choose a menu item : '+NC)
    print()
    if (choice1 == '1'):
        while (True):
            print(B+'1. Display all sequences')
            print('2. Display a list of sequences')
            choice2 = input('\n'+R+'Please, choose a menu item : '+NC)
            if (choice2 == '1'):
                colorSS('all')
                break 
            elif (choice2 == '2'):
                print(B+'\nPDB IDs :'+NC,end=' ')
                PDBlist = [pdb for pdb in input().split()]
                PDBlist.insert(0,MODEL)
                colorSS(PDBlist)
                break
            else:
                break
    elif (choice1 == '2'):
        pdbStart()
        retrieveFolds()
        call(['pymol '+'scripts/script.pml '+PATH+'&'], shell=True)
    elif (choice1 == '3'):
        createDirectory(TABLES)
        while (True):
            print(B+'1. Fold element study')
            print('2. Loop study')
            choice2 = input('\n'+R+'Please, choose a menu item : '+NC)
            if (choice2 == '1'):
                while (True):
                    fold = input('\n'+G+'fold element : '+NC)
                    if (fold in foldDict.keys()):
                        studyFold(fold)
                        break
                    else:
                        print(R+'\n ERROR '+NC+fold+R+' is not in the Fold.')
                        print('\n Please, choose a fold within : '+NC+str(list(foldDict.keys())))
                break
            elif (choice2 == '2'):
                while (True):
                    ele1 = input('\n'+G+'Fold element 1 : '+NC)
                    ele2 = input('\n'+G+'Fold element 2 : '+NC)
                    for i in range(len(eleList)):
                        if (eleList[i] == ele1):
                            break
                    if ((ele1 and ele2) in eleList and ele2 == eleList[i+1]):
                        studyLoop(ele1,ele2)
                        break
                    else:
                        print(R+' ERROR : Please, choose two successive fold elements within : '+NC+str(eleList))
                break
    elif (choice1 == '4'):
        with open(NEW+'seq_to_align.fasta','r') as f:
            SEQ = buildDict(f.readlines())
        SEQ_SS = retrieveSS(NEW+'sequences_ss.fasta',SEQ)
        ALN = addSeqAndAlign()
        ALN_SS = alignSS(NEW+'alignment_ss.fasta')
        foldDict = buildFoldDict()
        SSDict = buildSSDict()
    elif (choice1 == '5'):
        call(['Rscript '+'scripts/'+FAMILY+'.R'], shell=True)
    elif (choice1 == '6'):
        delete()
        break
