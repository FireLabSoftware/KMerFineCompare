#!/usr/bin/env python -i

import gzip
from time import time, strftime, localtime
import os
import operator
from itertools import chain
import sys
from copy import copy,deepcopy
import datetime
from BayesRatioAH import ConservativeFoldDifference



## 10/15/22 Version KmerFineCompare AR03


klen1 = 21
Down = 0
Up = 1
DownUp = 0
UpDown = 1
ReadsPerSample1 = 'All' ##100000 ##'All' ##1000000## 'All'  ##1000000 ##
FnList0 = ['PD3523-noaux-18x-rna_S1_L001_R1_001.fastq.gz','PD3523-g18ste-18x-rna_S3_L001_R1_001.fastq.gz'] ##,'PD35041h-3-g43ste-18x-rna_S2_L001_R1_001.fastq.gz']  ## 'TestSVsUnc13R2.fastq',
FnList1 = FnList0[:]
FnType1 = [Down,Up]
DefaultFiles1 = True
FnOrder1 = DownUp
MyCode1 = 'default'
FoldMinimum1 = 5.0
Regularization1 = 1.0
JunkFiles1 = ['OP50.fa']
TetritisFiles1 = ['illuminatetritis.fa']
JunkKlen1 =21
TetritisKlen1 = 11
JoinKLen1 = 11
ATrim1 = 0 ## Trim all A residuesof the 3' end of each read, also key in on internal strings of As (len>value of this parameter) and truncate reads at that point 
TTrim1 = 0 ## Trim all T residuesof the 3' end of each read, also key in on internal strings of Ts (len>value of this parameter) and truncate reads at that point 
FTrim1 = 0 ## if >0, Trim FTrim1 bases off of 5' end of each read
ETrim1 = 0 ## if >0, Trim ETrim1 bases off of 3' end of each read (after linker trim)
RefFileList1 = ['ce10tRNAs.fa','cerDNA.fa']
NormalizationMode = 0
## Modes
## 0 Normalize by total K-mer counts for each sample
## Additional Normalizations could be added but the author of this script is not convinced that any of these
## would provide additional utility
FDR1 = 0.05
ReportGranularity1 = 200000 ## How frequently to report progress in converting FASTQ/FASTA for comparisons
Linker1 = ''
RequireLinker1 = False



t0 = time()
now=datetime.datetime.now()
def strpad2(n):
    if type(n)==str:
        return n
    s = str(n)
    if len(s)==1:
        return '0'+s
    if len(s)>2:
        return s[-2:]
    return s
timestamp1 = ''.join(map(strpad2,(now.month,now.day,now.year,'_',now.hour,now.minute,now.second)))
def LogNote1(note,LogFileName):
    LogFile=open(LogFileName,mode='a')
    LogFile.write(note+'\t'+'; Time='+"{0:.2f}".format(time()-t0)+' sec'+'\t'+strftime("%m_%d_%y__%H_%M_%S",localtime())+' \r')
    LogFile.close()
    print(note.split('#')[0].replace('\t',' ') + '; Time='+"{0:.2f}".format(time()-t0)+' sec')

def FileInfo1(FileID):
    if type(FileID)==str:
        Fn1=FileID
    else:
        Fn1=FileID.name
    s11=os.stat(Fn1)
    return '\t'.join([Fn1,
                      '#',
                      'Path='+os.path.abspath(Fn1),
                        'Size='+str(s11[6]),
                        'Accessed='+strftime("%m_%d_%y__%H_%M_%S",localtime(s11[7])),
                        'Modified='+strftime("%m_%d_%y__%H_%M_%S",localtime(s11[8])),
                        'Created='+strftime("%m_%d_%y__%H_%M_%S",localtime(s11[9])),
                        'FullInfo='+str(s11)])
LogFile1 = 'KFC_LogFile'+timestamp1+'.txt'
def LogOpeningFile1(FileID):
    LogNote1("Opening "+FileInfo1(FileID),LogFile1)
def LogClosingFile1(FileID):
    LogNote1("Closed "+FileInfo1(FileID),LogFile1)
def LogRunningFile1(FileID):
    LogNote1("Running "+FileInfo1(FileID),LogFile1)



ai1 = 0
FullDescriptor1 = ';; KmerFineCompare with parameters as follows'+'\r;;  Script Version: '+'\r;;  '.join(sys.argv)+'\r;;  Script Run Time: '+'_'.join(map(str,(now.month,now.day,now.year,'',now.hour,now.minute,now.second)))+'\r'
OutputFileBase1 = "DiffEx"+timestamp1
LogNote1(FullDescriptor1,LogFile1)
DirectionD1 = {0:'Down', 1:'Up'}

while ai1<len(sys.argv):
    ArgMode1 = "replace"
    a1=sys.argv[ai1]
    if '+' in a1.split('=')[0]:
        ArgMode1 = 'extend'
    ai1+=1
    a1=a1.replace('"','').replace("'","")
    if a1[0]=='-':
        a11=a1.strip()[1:].lower()
        if len(sys.argv)>ai1:
            a22 = sys.argv[ai1].strip().replace('"','').replace("'","")
        else:
            a22 = ''
        ai1+=1
    else:
        a11=a1.split('=')[0].strip().lower()
        a22=a1.split('=')[-1].strip()
    if a11.startswith('down'):
        if ArgMode1 == 'replace' and DefaultFiles1:
            FnList1 = []
            FnType1 = []
            FnOrder1 = DownUp
            DefaultFiles1 = False
        FnList1.extend(a22.split(','))
        FnType1.extend([Down,]*len(a22.split(',')))
    elif a11.startswith('up'):
        if ArgMode1 == 'replace' and DefaultFiles1:
            FnList1 = []
            FnType1 = []
            FnOrder1 = UpDown        
            DefaultFiles1 = False
        FnList1.extend(a22.split(','))
        FnType1.extend([Up,]*len(a22.split(',')))
    elif a11.startswith('junkf'):
        if ArgMode1 == 'replace':
            JunkFiles1=a22.split(',')
        elif ArgMode1 == 'extend':
            JunkFiles1+=a22.split(',')
    elif a11.startswith('tetritisf'):
        if ArgMode1 == 'replace':
            TetritisFiles1=a22.split(',')
        elif ArgMode1 == 'extend':
            TetritisFiles1+=a22.split(',')
    elif a11.startswith('reff'):
        if ArgMode1 == 'replace':
            RefFileList1=a22.split(',')
        elif ArgMode1 == 'extend':
            RefFileList1+=a22.split(',')
    elif a11.startswith('junkk'):
        JunkKlen1=int(a22)
    elif a11.startswith('report'):
        ReportGranularity1=int(a22)
    elif a11.startswith('tetritisk'):
        TetritisKlen1=int(a22)
    elif a11.startswith('klen'):
        klen1 = int(a22)
    elif a11.startswith('joink'):
        JoinKLen1 = int(a22)
    elif a11.startswith('reads'):
        ReadsPerSample1=int(a22)
    elif a11.startswith('ftrim'):
        FTrim1=int(a22)
    elif a11.startswith('etrim'):
        ETrim1=int(a22)
    elif a11.startswith('output'):
        OutputFileBase1=a22
    elif a11.startswith('fdr'):
        FDR1=float(a22)
    elif a11.startswith('fold'):
        FoldMinimum1=float(a22)
    elif a11.startswith('reg'):
        Regularization1=float(a22)
    elif a11.startswith('link'):
        Linker1=a22.upper().strip()
    elif a11.startswith('atrim'):
        if a22.lower().startswith('f'):
            ATrim1 = 0
        else:
            try:
                ATrim1 = int(a22)
            except:
                ATrim1 = -1
    elif a11.startswith('ttrim'):
        if a22.lower().startswith('f'):
            TTrim1 = 0
        else:
            try:
                TTrim1 = int(a22)
            except:
                TTrim1 = -1
    elif a11.startswith('requirel'):
        if a22.lower().startswith('f'):
            RequireLinker1 = False
        else:
            RequireLinker1 = True
    elif a1.lower().endswith('.py'):
        MyCode1=a1.strip()
    continue
nF0 = len(FnList1)

if FnList1 == FnList0:
    LogNote1("Running KFC in Debug Mode with default files.  If these files aren't in the current folder, program will stop: "+str(FnList0),LogFile1)

LogRunningFile1(MyCode1)
LogNote1("Python Flavor/Version: "+sys.version,LogFile1)

AllBase1=['G','A','T','C']
Numbase1={AllBase1[0]:0,AllBase1[1]:1,AllBase1[2]:2,AllBase1[3]:3,'N':0}
## e4 is a set of 4**x exponents
e4=[]
ne4=[]
for i in range(32):
    e4.append(4**i)
    ne4.append(4**i-1)
ex2=[2**i for i in range(64)]
AntisenseB1={'A':'T','T':'A','G':'C','C':'G'}
Tr2=''  ## Tr2 allows a quick translation of sequence to a numerical array.
for i in range(256):
    if chr(i) in 'AGCTagct':
        Tr2+=chr(Numbase1[chr(i).upper()])
    else:  ##in case there are unusual characters in the sequence, they will be converted to Gs
        Tr2+=chr(Numbase1['N'])
TrN1=''  ## Tr2 allows a quick translation of sequence to a numerical array.
for i in range(256):
    if chr(i) in 'AGCTagct':
        TrN1+=chr(1)
    else:  ##in case there are unusual characters in the sequence, they will be converted to Gs
        TrN1+=chr(0)

## antisense-- returns the reverse complement of argument string 
filterminus=''
ASB11="AaCcNn*nNgGtT"
for i in range(256):
    if chr(i) in ASB11:
        filterminus+=ASB11[12-ASB11.find(chr(i))]
    else:
        filterminus+=' '
def antisense(s):
    '''return an antisense and filtered version of any sequence)'''
    return s.translate(filterminus)[::-1].replace(' ','')
GATC = 'GATC'
Mask1 = (4**klen1)-1
def kSeqToNum(s,klen):
    n = 0
    sA1 = []
    for i,c in enumerate(s):
        n = ((n<<2)&Mask1) + Numbase1[c.upper()]
        if i>=klen-1:
            sA1.append(n)
    return sA1

def decodebin(n,k):    
    s = ''
    for h in range(k):
        s+=GATC[n&3]
        n >>= 2
    return s[::-1]

if ATrim1>0:
    An1 = 'A' * ATrim1
if TTrim1>0:
    Tn1 = 'T' * TTrim1
    
def ATrim(s):
    s1 = s.rstrip('A')
    if ATrim1>0:        
        p1 = s1.find(An1)
        if p1>0:
            s1=s1[:p1]
        elif p1==0:
            return ''
    return s1

def TTrim(s):
    s1 = s.rstrip('T')
    if TTrim1>0:        
        p1 = s1.find(Tn1)
        if p1>0:
            s1=s1[:p1]
        elif p1==0:
            return ''
    return s1

def FTrim(s,l):
    return s[l:]

def ETrim(s,l):
    return s[:-l]

def FilesToKDict1(Fn,reads):
    D = {}
    if type(Fn)==str:
        Fn = [Fn,]
    nF = len(Fn)
    Totals = [0,] * nF
    for n1,Fn0 in enumerate(Fn):
        if Fn0.endswith('.gz'):
            F0 = gzip.open(Fn0,mode='rt')
        else:
            F0 = open(Fn0,mode='rt')
        LogOpeningFile1(F0)
        if 'fastq' in Fn0.lower():
            FileType = 0 ##'FastQ'
            F00 = F0
        elif '.fa' in Fn0.lower():
            FileType = 1 ##'FastA'
            F00 = chain(F0,['>',])
        else:
            FileType = 2 ##'ReadList' ## line format read <tab> counts for each read
            F00 = F0                
        L1 =''
        Multiplicity = 1
        Tags = 0
        for i,L in enumerate(F00):
            Ls = L.strip()
            if len(Ls)==0: continue
            if FileType == 0:
                if i%4 != 1:
                    continue
                else:
                    L1 = Ls
            elif FileType == 1 and L[0]!='>':
                L1 += Ls
                continue
            if FileType == 2:
                L1 = Ls.split()[0]
                try:
                    Multiplicity = int(L.strip().split()[-1])
                except:
                    Multiplicity = 1
            if Linker1:
                Ps = L1.rfind(Linker1)
                if Ps<0 and RequireLinker1:
                    continue
                L1 = L1[:Ps]
            if L1 and L1[0] in 'AGCT':                                       
                if ATrim1 != 0:
                    L1 = ATrim(L1)
                if TTrim1 != 0:
                    L1 = TTrim(L1)
                if FTrim1>0:
                    L1 = FTrim(L1,FTrim1)
                if ETrim1>0:
                    L1 = ETrim(L1,ETrim1)
                if len(L1)<klen1:
                    continue
                s = kSeqToNum(L1,klen1)
                for s1 in s:
                    if not(s1 in D):
                        D[s1]=[0,] * nF
                    D[s1][n1]+=Multiplicity
                    Totals[n1]+=Multiplicity
                Tags += 1
                if Tags%ReportGranularity1==0 and Tags>0:
                    LogNote1(Fn0.split('.')[0]+"; Lines:"+str(Tags),LogFile1)
            L1 = ''
            if Tags>reads>0:
                break
        F0.close()
    return D

def KDictMetrics(D,nF):
    '''input is a dictionary with keys being anything and values that are lists of integer values, output is a listlet of metrics per sample,
        Number of species (NS) 0,
        Total Instances (TI) 1,
        Nonzero items (NI) 2,
        Total squared instances (TS) 3
        Total cubed instances (TC) 4
        AverageCountNumber (AC) 5
        CoincidenceProbability (CP) 6
        TRipleCoincidenceProbability (TP) 7
        RAtioTotaltoUnique (TU) 8
        CoincidenceMatrixRaw (CM) 9
        InterProbabilities (IP) 10'''
    NS = len(D)
    TI = [0,]*nF
    NI = [0,]*nF
    TS = [0,]*nF
    TC = [0,]*nF
    AC = [0.0,]*nF
    CP = [0.0,]*nF
    TP = [0.0,]*nF
    TU = [0.0,]*nF
    CM = [[0,]*nF for ww in range(nF)]
    IP = [[0.0,]*nF for ww in range(nF)]
    for k in D:
        v = D[k]
        for i in range(nF):
            vv = D[k][i]
            TI[i] += vv
            if vv>0:
                NI[i] += 1
            TS[i] += vv**2
            TC[i] += vv**3
            for j in range(i):
                CM[i][j] +=vv*D[k][j]
                    
    for i in range(nF):
        if NS>0:
            AC[i] = TI[i]*1.0/NS
            TU[i] = TI[i]*1.0/NI[i]
            for j in range(i):
                IP[i][j] = CM[i][j]*1.0/(TI[i]*TI[j])
                IP[j][i] = IP[i][j]
        if TI[i]>1:    
            CP[i] = (TS[i]-TI[i])*1.0/(TI[i]**2-TI[i])
            IP[i][i] = CP[i]
        if TI[i]>2:    
            TP[i] = (TC[i]-3*TS[i]+2*TI[i])*1.0/(TI[i]*(TI[i]-1)*(TI[i]-2))
        
    return NS,TI,NI,TS,TC,AC,CP,TP,TU,CM,IP

def kmerindex(Fn,klen):
    '''input is a fastA type file or list (often junk genomes), klen, output is a dictionary of k-mers and their incidence (both strands)'''  
    if type(Fn)==str:
        if ',' in Fn and not(os.path.isfile(Fn)):
            Fn = Fn.split(',')
        else:
            Fn = [Fn,]
    for Fn2 in Fn:
        if type(Fn2)==str:
            if Fn2.endswith('.gz'):
                F1 = gzip.open(Fn,mode='rt')
            else:
                F1 = open(Fn2,mode='rt')
        else:
            F1 = Fn2
        LogOpeningFile1(F1)
        S1 = []
        D = {}
        for L1 in F1:
            if L1[0]=='>':
                if len(S1)>0:
                    S1=''.join(S1)
                    A1 = antisense(S1)
                    for i1 in range(len(S1)-klen+1):
                        k1 = S1[i1:i1+klen]
                        a1 = A1[i1:i1+klen]
                        if not(k1 in D):
                            D[k1] = 0
                        if not(a1 in D):
                            D[a1] = 0
                        D[k1]+=1
                        D[a1]+=1
                    S1 = []
            else:
                S1.append(L1.strip())
        S1=''.join(S1)
        A1 = antisense(S1)
        for i1 in range(len(S1)-klen+1):
            k1 = S1[i1:i1+klen]
            a1 = A1[i1:i1+klen]
            D[k1] = 0
            D[a1] = 0
            D[k1]+=1
            D[a1]+=1
        F1.close()
    return D
def ConsensusSeq(ConsensusD,classS):
    'returns a consensus sequence; with "N" at a position if only a minority match each other'
    pmin = min(ConsensusD.keys())
    pmax = max(ConsensusD.keys())
    Consensus = []
    for m in range(pmin,pmax+1):
        TotalB = sum(ConsensusD[m].values())
        nextb = 'n'
        for b in ['G','A','T','C']:
            if (b in ConsensusD[m]) and (ConsensusD[m][b]*2>TotalB):
                nextb = b.lower()
        if m>=0 and m<=len(classS):
            nextb = nextb.upper()
        Consensus.append(nextb)
    return ''.join(Consensus),-pmin ### consensus sequence and offset for display
        
DRef1 = {}
NRef1 = ['']
def RefIndex(Fn,klen):
    ## reads file Fn as a FastA file and takes k-mers into a dictionary D (input k-mer output seqID, position), and a list N (input seqID, output seqID text line [>xxx line in FASTA file)
    if Fn.endswith('.gz'):
        F1 = gzip.open(Fn,mode='rt')
    else:
        F1 = open(Fn,mode='rt')
    LogOpeningFile1(F1)
    S1 = []
    n1 = len(NRef1)-1
    for L1 in F1:
        if L1[0]=='>':
            if len(S1)>0:
                S1=''.join(S1)
                A1 = antisense(S1)
                for i1 in range(len(S1)-klen+1):
                    k1 = S1[i1:i1+klen]
                    a1 = A1[i1:i1+klen]
                    if not(k1) in DRef1:
                        DRef1[k1] = (n1,i1)
                        DRef1[a1] = (n1,-i1)
                S1 = []
            n1+=1
            NRef1.append(L1[1:].strip().split()[0])
        else:
            S1.append(L1.strip().upper())
    S1=''.join(S1)
    A1 = antisense(S1)
    if len(S1)>0:
        S1=''.join(S1)
        A1 = antisense(S1)
        for i1 in range(len(S1)-klen+1):
            k1 = S1[i1:i1+klen]
            a1 = A1[i1:i1+klen]
            if not(k1 in DRef1):
                DRef1[k1] = (n1,i1+1)  ## 1 based first position of match on sense strand
                DRef1[a1] = (n1,i1-len(S1)) ## -DRef1 is 1 based last position of match on antisense strand
    F1.close()
    return 0

for Rn1 in RefFileList1:
    RefIndex(Rn1,klen1)

def RefFind(S1,klen):
    RRef = {}
    for i1 in range(len(S1)-klen+1):
        k1 = S1[i1:i1+klen]
        if k1 in DRef1:
            n1,p1 = DRef1[k1]
            if not(n1 in RRef):
                RRef[n1] = {}
            if not(p1 in RRef[n1]):
                RRef[n1][p1]=0
    desc1 = ''
    for n1 in sorted(RRef.keys()):
        desc1+=NRef1[n1]
        ulist1 = sorted(RRef[n1].keys())
        start1 = 0
        for g,i in enumerate(ulist1):
            if i-1 in RRef[n1]:
                span1 += 1
            if not(i-1 in RRef[n1]) or g==len(ulist1)-1:
                if start1<0:
                    desc1 += '_'+str(-start1) + '-'
                    desc1 += str(-start1-klen+1-span1)
                    desc1 += '(-)_'
                if start1>0:
                    desc1 += '_'+str(start1) + '-'
                    desc1 += str(start1+klen-1+span1)
                    desc1 += '(+)_'
                start1 = i
                span1 = klen

    return desc1
            
if type(ReadsPerSample1)==str and ReadsPerSample1.lower().startswith('all'):
    ReadsPerSample1 = -1
D1 = FilesToKDict1(FnList1,ReadsPerSample1)
NS1,TI1,NI1,TS1,TC1,AC1,CP1,TP1,TU1,CM1,IP1 = KDictMetrics(D1,nF0)



##    '''input is a dictionary with keys being anything and values that are lists of integer values, output is a listlet of metrics per sample,
##        Number of species (NS) 0,
##        Total Instances (TI) 1,
##        Nonzero items (NI) 2,
##        Total squared instances (TS) 3
##        Total cubed instances (TC) 4
##        AverageCountNumber (AC) 5
##        CoincidenceProbability (CP) 6
##        TRipleCoincidenceProbability (TP) 7
##        RAtioTotaltoUnique (TU) 8
##        CoincidenceMatrixRaw (CM) 9
##        InterProbabilities (IP) 10'''

## Figure out which pairs to compare
rR1 = range(nF0)

cL1 = []
for j in rR1:
    if j==0 or (FnOrder1==DownUp and FnType1[j]==Down and FnType1[j-1]==Up) or (FnOrder1==UpDown and FnType1[j]==Up and FnType1[j-1]==Down):
        Ld = []
        Lu = []
    if FnOrder1==UpDown and FnType1[j]==Down:            
        for k in Lu:
            cL1.append((j,k))
    elif FnOrder1==DownUp and FnType1[j]==Up:
        for k in Ld:
            cL1.append((k,j))
    elif  FnOrder1==UpDown and FnType1[j]==Up:
        Lu.append(j)
    elif  FnOrder1==DownUp and FnType1[j]==Down:
        Ld.append(j)
cfdD1 = {}
d0 = {}
d1 = {}
rR1 = range(nF0)
lD1 = len(D1)
FDRa1 = min(FDR1,FDR1*klen1/len(D1))  ## a bonferroni-like correction, but realizig that not every k-mer is independent
for n1 in D1:
    x1 = D1[n1]
    Reject1 = False
    for jd,ju in cL1:
        xd = (x1[jd]+Regularization1)/(TI1[jd]+Regularization1)
        xu = (x1[ju]+Regularization1)/(TI1[ju]+Regularization1)
        if xd*FoldMinimum1 > xu:
            Reject1 = True
            break
    if Reject1:
        continue
    d0[n1] = D1[n1]
for n1 in d0:
    x1 = D1[n1]
    Reject1 = False
    for jd,ju in cL1:
        if (x1[jd],x1[ju],jd,ju) in cfdD1:
            cfd = cfdD1[(x1[ju],x1[jd],ju,jd)]
        else:
            cfd = ConservativeFoldDifference(x1[ju],x1[jd],TI1[ju],TI1[jd],FDRa1)
            cfdD1[(x1[ju],x1[jd],ju,jd)] = cfd
        if cfd < FoldMinimum1:
            Reject1 = True
            break
    if Reject1:
        continue
    s1 = decodebin(n1,klen1)
    d1[s1] = d0[n1]

ShortestOverlap1=JoinKLen1  ## The shortest overlap that will be used to join segments

if len(d1) == 0:
    print("No differences found meeting criteria-- Exiting Program")
    aardvark
m0=max([len(v1) for v1 in d1])


StillMerging3=True
while StillMerging3:
    StillMerging3=False
    for m1 in range(m0-1,ShortestOverlap1-1,-1):
        d5={}  ## keys are subsequences from the 5' end, values are full sequence
        StillMerging2=True
        while StillMerging2:
            StillMerging2=False
            for s1 in sorted(d1.keys(),key=lambda q:sum(d1[q])):
                if len(s1)>m1:
                    d5[s1[:m1]]=s1
            StillMerging1=True
            while StillMerging1:
                StillMerging1=False
                for s2 in sorted(d1.keys(),key=lambda q:-sum(d1[q])):
                    if not s2 in d1:
                        continue
                    if len(s2)>m1:
                        s3=s2[-m1:]
                        if s3 in d5:
                            s4=d5[s3]
                            if s2==s4:
                                continue
                            StillMerging1=True
                            StillMerging2=True
                            StillMerging3=True
                            s5=s2[:-m1]+s4
                            count1=[d1[s4][p]+d1[s2][p] for p in rR1]
                            if s2 in d1:
                                del d1[s2]
                            if s4 in d1:
                                del d1[s4]
                            if s3 in d5:
                                del d5[s3]
                            if s2[:m1] in d5:
                                del d5[s2[:m1]]
                            d1[s5]=count1
                            d5[s5[:m1]]=s5
B=sorted(d1.keys(),key=lambda q:-len(q))
F=open('KmerAssemblies_'+OutputFileBase1+'.fa',mode='w')
LogOpeningFile1(F)

Junk1 = kmerindex(JunkFiles1,JunkKlen1)
Tetritis1 = kmerindex(TetritisFiles1,TetritisKlen1)

KmerToClassD1 = {}
ClassToReadsD1 ={}
ClassToReadsD2 ={}
ClassToSequence1 = {}
ClassConnections1 = {}
ClassConsensusD1 = {}
AntisenseConnections1 = {}

i=0
for b in B:
    FilterMeOut1 = False
    for j in range(len(b)-JunkKlen1+1):
        if b[j:j+JunkKlen1] in Junk1:
            FilterMeOut1 = True
    for j in range(len(b)-TetritisKlen1+1):
        if b[j:j+TetritisKlen1] in Tetritis1:
            FilterMeOut1 = True
    if not(FilterMeOut1):
        i+=1
        F.write('>s_'+str(i)+'_DiffSeq_'+OutputFileBase1+'\r')
        F.write(b+'\r')
        ClassToSequence1[i] = b
        for j in range(len(b)-klen1+1):
            KmerToClassD1[b[j:j+klen1]] = i
        ClassToReadsD1[i] = {}  ## counts
        ClassToReadsD2[i] = {}  ## Offsets (bases upstream of the consensus in read)
        ClassConsensusD1[i] = {}
        ClassConnections1[i] = {}
F.close()
for iF1,Fn0 in enumerate(FnList1):
    if Fn0.endswith('.gz'):
        F0 = gzip.open(Fn0,mode='rt')
    else:
        F0 = open(Fn0,mode='rt')
    LogOpeningFile1(F0)
    if 'fastq' in Fn0.lower():
        FileType = 0 ##'FastQ'
        F00 = F0
    elif '.fa' in Fn0.lower():
        FileType = 1 ##'FastA'
        F00 = chain(F0,['>',])
    else:
        FileType = 2 ##'ReadList' ## line format read <tab> counts for each read
        F00 = F0                
    L1 =''
    Multiplicity = 1
    Tags = 0
    for i,L in enumerate(F00):
        Ls = L.strip()
        if len(Ls)==0: continue
        if FileType == 0:
            if i%4 != 1:
                continue
            else:
                L1 = Ls
        elif FileType == 1 and L[0]!='>':
            L1 += Ls
            continue
        if FileType == 2:
            L1 = Ls.split()[0]
            try:
                Multiplicity = int(L.strip().split()[-1])
            except:
                Multiplicity = 1
        if Linker1:
            Ps = L1.rfind(Linker1)
            if Ps<0 and RequireLinker1:
                continue
            L1 = L1[:Ps]
        if ATrim1!=0:
            L1 = ATrim(L1)
        if TTrim1!=0:
            L1 = TTrim(L1)
        if FTrim1:
            L1 = FTrim(L1,FTrim1)
        if ETrim1:
            L1 = ETrim(L1,ETrim1)
        c1 = 0
        for j in range(len(L1)-klen1+1):
            su1 = L1[j:j+klen1]
            if su1 in KmerToClassD1:
                if c1 == 0:
                    cD1 ={} ## dictionary of all hit Classes                    
                c1 = KmerToClassD1[su1]
                if not(L1 in ClassToReadsD1[c1]):
                    ClassToReadsD1[c1][L1] = [0,] * nF0
                    Offset1 = j-ClassToSequence1[c1].find(su1)  ##number to add to position in class to get position in read
                    ClassToReadsD2[c1][L1] = Offset1
                if not c1 in cD1:
                    ClassToReadsD1[c1][L1][iF1]+=Multiplicity
                    for c2 in cD1:
                        if not(c2 in ClassConnections1[c1]):
                            ClassConnections1[c1][c2] = 0
                        if not(c1 in ClassConnections1[c2]):
                            ClassConnections1[c2][c1] = 0
                        ClassConnections1[c1][c2] += 1
                        ClassConnections1[c2][c1] += 1
                    cD1[c1] = 0
                    for m in range(len(L1)):
                        md1 = m-ClassToReadsD2[c1][L1]
                        if not md1 in ClassConsensusD1[c1]:
                            ClassConsensusD1[c1][md1] = {}
                        if not L1[m] in ClassConsensusD1[c1][md1]:
                            ClassConsensusD1[c1][md1][L1[m]] = 0
                        ClassConsensusD1[c1][md1][L1[m]] += 1
                cD1[c1]+=1                                                    
        L1 = ''
        Tags += 1
        if Tags>ReadsPerSample1>0:
            break
        
FReadsOut1 = open('CaughtReads_'+OutputFileBase1+'.txt',mode='w')
FReadsOut1.write(FullDescriptor1+'\r')
LogOpeningFile1(FReadsOut1)
FReadsOut1.write('FileName\t'+'\t'.join(FnList1)+'\r')
FReadsOut1.write('Direction\t'+'\t'.join([DirectionD1[x] for x in FnType1])+'\r')
FConsensusOut1 = open('Consensi_'+OutputFileBase1+'.txt',mode='w')
FConsensusOut1.write(FullDescriptor1+'\r')
LogOpeningFile1(FConsensusOut1)

for c1 in ClassToReadsD1:
    Consensus1,Offset0 = ConsensusSeq(ClassConsensusD1[c1],ClassToSequence1[c1])
    Refs1 = RefFind(Consensus1.upper(),klen1)
    FConsensusOut1.write('>sE_'+str(c1)+'_'+ClassToSequence1[c1]+'_'+Refs1+'\r')##+'_DiffSeq_'+OutputFileBase1+'\r')
    FConsensusOut1.write(Consensus1+'\r')
    FReadsOut1.write('>sE_'+str(c1)+'_'+ClassToSequence1[c1]+'_'+Refs1+'\r')
    FReadsOut1.write(Consensus1+'\r')
    MyKeys=sorted(ClassToReadsD1[c1], key= lambda x: (ClassToReadsD2[c1][x], -sum(ClassToReadsD1[c1][x])))
    for cKey in MyKeys:
        FReadsOut1.write(' '*(Offset0-ClassToReadsD2[c1][cKey])+cKey+'\t'+'\t'.join(map(str,ClassToReadsD1[c1][cKey]))+'\r')
    FReadsOut1.write('\r')
    
FReadsOut1.close()
FConsensusOut1.close()
FSummaryOut1 = open('SummaryData_'+OutputFileBase1+'.txt',mode='w')
LogOpeningFile1(FSummaryOut1)
FSummaryOut1.write(FullDescriptor1+'\r')
FSummaryOut1.write('\r')
FSummaryOut1.write('Number of k-mer species:'+str(NS1)+'\r')
FSummaryOut1.write('\r')
FSummaryOut1.write('FileName\t'+'\t'.join(FnList1)+'\r')
FSummaryOut1.write('Direction\t'+'\t'.join([DirectionD1[x] for x in FnType1])+'\r')
FSummaryOut1.write('Total Instances\t'+'\t'.join(map(str,TI1))+'\r')
FSummaryOut1.write('Nonzero items\t'+'\t'.join(map(str,NI1))+'\r')
FSummaryOut1.write('Total squared instances\t'+'\t'.join(map(str,TS1))+'\r')
FSummaryOut1.write('Total cubed instances\t'+'\t'.join(map(str,TC1))+'\r')
FSummaryOut1.write('AverageCountNumber\t'+'\t'.join(map(str,AC1))+'\r')
FSummaryOut1.write('CoincidenceProbability\t'+'\t'.join(map(str,CP1))+'\r')
FSummaryOut1.write('TripleCoincidenceProbability\t'+'\t'.join(map(str,TP1))+'\r')
FSummaryOut1.write('RatioTotaltoUnique\t'+'\t'.join(map(str,TU1))+'\r')
FSummaryOut1.write('\r')
FSummaryOut1.write('CoincidenceMatrixRaw\r')
for i in range(nF0):
    FSummaryOut1.write(FnList1[i]+'\t'+'\t'.join(map(str,CM1[i]))+'\r')
FSummaryOut1.write('\r')
FSummaryOut1.write('InterProbabilities\r')
for i in range(nF0):
    IP1[i][i] = CP1[i]  ## corrects the value of self comparisons for simple me-to-me coincidences
    FSummaryOut1.write(FnList1[i]+'\t'+'\t'.join(map(str,IP1[i]))+'\r')
FSummaryOut1.close()

LogNote1('Finishing KmerFineCompare',LogFile1)
try:
    LogNote1(open(MyCode1,mode='r').read(),LogFile1)
except:
    pass


    '''input is a dictionary with keys being anything and values that are lists of integer values, output is a listlet of metrics per sample,
        Number of species (NS) 0,
        Total Instances (TI) 1,
        Nonzero items (NI) 2,
        Total squared instances (TS) 3
        Total cubed instances (TC) 4
        AverageCountNumber (AC) 5
        CoincidenceProbability (CP) 6
        TRipleCoincidenceProbability (TP) 7
        RAtioTotaltoUnique (TU) 8
        CoincidenceMatrixRaw (CM) 9
        InterProbabilities (IP) 10'''






        




