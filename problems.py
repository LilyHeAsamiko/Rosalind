# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 17:08:30 2024

@author: QinLilyHe
"""
from itertools import permutations
from collections import Counter

string = 'ATATGATACCTGTTTGACTATCTATGCTGACGTAGCTATTACAGCGCCACAACCGAACTGTGTGAGCCCTTGAAGGCATTCCGCAGAACCACCTCATTAAGATTGGGAATTGCAGTCACAGTTACGTGGCGGAACCACATGTGTACAAAACACATAGGGCCGTGGGGCCTGGTGGTCATTCTGTTGACGAGCTTTCGTTACAGCGAGCGGCCGGATAGGAGAGTCCTTTGGGCGACTATGCGCCGCAATACCTCCATTATCCCACCCGTAAGAGTTTCAAGTCGGACAACGAATGTATCGTTACTCTACGGGCGCGCCCCTATTCGTTAGATGTGGTAAACCAGGGTCCTCAACTACCGCGATTCCGCAACCTTCCCTGTATGATTGTGGAGGGAGACTTAGGGTAGGTTCCCCTCTGTAAGCCCGTGCCCCAATTAGTTATGTCACCCTCCTAATCAGACGGCTACTTATGTCGCTAAGTGCCCACTGTAAGAACGATGAAGAGCCACGGCAGTTGTTCTAAAGCGGGTTTTTTTATTTCACCCCCAGTAACTTCATGACCAACTCTTGTACGCCCCGCATGGCACACGTAGAAGCTCGTACCCCTCAAATTCAGCGTTAGATGGCTGATCAGCCCCACCCAGCTGAGCAGCCGCACCTCCATAAGGTCGGTTCCGCGTTCGTTCACTGCAACCTCTCGTAAAGGGAGGTGACCGGTTTCTATGCGGCAAATACTCTTGTAGATACGTAGAACGAAGGCTGCCGGTGCTACCTACCAAAGGCGATAATGTAAGCCACTACCGTACATAGACTGACAGGGGTTCTTTCTGTCCGGGGGAAGTTCTACATCAAAAAGCGCTCCAGGAGCAGAGTTGGTCTGATTAGCTACACGCTGCCCATACGGAATCATTGATCGATGTCGGCGGTCCGTCGTTCGATAGGAGCATCCCCATAGTCACTCATCTCCTTCGCGTTGAGCACACCACATCCCGAAT'
print(str(dict(Counter(string))['A'])+' '+str(dict(Counter(string))['C'])+' '+str(dict(Counter(string))['G'])+' '+str(dict(Counter(string))['T']))

n = 3
print(len(list(permutations(range(1,n+1)))),'\n')
[print(f'{p} \n') for p in list(permutations(range(1,n+1)))]
#Transcribing DNA into RNA
string = 'GTGCATCGTCCTGGCAGAAGTTAATATGAGGTTAGGCGTGTCCAGTAGGCACTCGATCGACACTGGGTAAAATACGTCTCCGGGCCGTTCAGGAATAAATTCGGTGGCGGGAGTCACGCAAAGCGATATGGCTTGACGCACGCTAAACTCACCAGCTGAACCTCCGTTATGACTTCCCCCTCCACACACTTGCCGATTAAGACGCGCTCTAAAAACCGCTCGGCGACCGGTGCAGAAAGGAGGATCGAACATGTACTCCCGAATGGTGATTAATCCGGCTTAAGGGTGCAGCTGATGCAGAGTGAGAAGAAGCGCGTGGCCAGCCTTATAGAATAAGCAGTCAACTTTGTTTCCTAGTTTATACAAAAGCCACATCAGTCTGCCCGGTACAGTGTGCGTTCTAGGTTAGCCAGTGCTCGACCTAGTCATGAACCTTTATCCGAAGTCTGGCCGTGTGGTCACTCGACATGAAGAGGATATAAAGAGTATGATAACGCGCCGAAGAAAAAGGCCCAACCCTTGTGATACCAATAGAGTGTAGTGCGGTCCAAGCTCCTATGCGTGTGCGCGAATTGGCGACTCCACGTTTAGGTCTCCAAAAAAAAGCCGGCCGTCTTGGATGCATTTCAGTAACGACTAGATGACTCCTGTGATCCAGAAGGCGTAAGCCAGAAGTGAAGGAATACGACGCTCAAGCACAAACTTGTACTTCGAAAGACGTCCTAAGCGAGGGCTTGAATCCAAATTGGCGTGTGTCACGGTGCTGCCGATAGCCTGAACTTGTGATCACCAGTGGGCGTGGGGGTTTCCTCCCATCACTAGGGAACCCTTGTCAGTTCTGGGCAGATGAATGGAACCTCGTACCATATCAACTCGGGGGGCCTTTGGGCTCGGTAAACGGTGGGATGCTCCAAATCCCAGGGTCCCATGTGCGACGAAATCCTAGATTG'
print(string.replace('T', 'U'))
#Complementing a Strand of DNA(REVERSE COPLEMENT)
string = 'GAGTCTCGGGAAAACCAATGCACTAAGAGTGCTAGACAACGGTTCTTGCGTTGAGACGCTTCACCGGTGCCCTTTATCTTGATAGGACCTGGCTGTCCATTGAGATCGTCAGATTAAGCTCAGATAAAGCCACCCTACGCCCTTACGACGCTTTTTCCCCCGGTGCTGACGGTTGGTTCCTTTCTCGATATGCCGCACTGGCACATGCCCGGAAGCGCTTTTTCTTCATAGCAACGACGTTGGCTCCCAATCTAATAAGGGACGTATCTTGTGCGAGCCATGGTGCGCGTCTCGGTCCCCAGCTCCATTGTACGTTCTAAGTAGTGCCGCGAACCCTCAGACTAACATCCCGCGCCATTAGACCCCATAATTCACAATCCACTGCAGACTGTCTCTTGGCTGAATGTCGCAGAACTCTTTATGAGTGTGTTAATTTACGGCGACCAAGGGGTAGCACGTTCCTTTTTCGTCGCTAGACGTCACAGTATCCTGGGGAGGTAAAGTCCGGCGGCGTCACAATAGCTCCAAAGGTCCAGGAGTCCAATGATATCCTCAGAGGAATGCTGATCGCTACATAGAGTACATCAGCCCCGTTTCAAACAGGACAGCAACTCAAACTCTGCCGTCAGTCCTAGATTCGGGATTGCAAAGCCGCACGGGGTCGAACTTGACGTTACCTACATAGGCTGGATATTACTTTGAAAGGCTCCCAGCGCAACCTTTCGACCGCTTCCAGTGGGTCGGTTGCCTTCGTGTACTAGAGTTGATCACGGGCGGGCGCGCACAACGCCGTTGACCTCTGCGGCCTGACTGGATGGTCTCTCGATCCCTCGTTCTCTAAGATACGCAACCCGGCTCGCCGATTCAGCACGCGTTCGCATAAGAGGGCTTCTCGCGCTAAT'
dictR = {'A':'T','T':'A','C':'G','G':'C'} 
print(''.join([dictR[s] for s in string])[::-1]) 
#Rabbits and Recurrence Relations n<=40, k<=5
n = 5
k = 3
n-1 + (n-1)*k
#1                             1
#1                             1
#1 k                           1+k
#k 1 k                         2*k+1
#k(1+k) k 1 k                  k+(1+k)**2
#k(k 1 k) k(1+k) k 1 k         (k+1)(k 1 k k)
c = 1
p = 0
k = 2
m = 36
for m in range(1,m):
    c,p = k*p,c+p 
    print(c,p)
print(c+p)
#Computing GC Content
#test: 
#file = r"E:\Rosalind\Rosalind_t.txt"
file = r"C:\Users\Admin\Downloads\rosalind_gc.txt"
with open(file,'r') as f:
    rf = f.readlines()
#[print(rf[i],round(100*(''.join(rf[i+1:i+3]).count('C')+''.join(rf[i+1:i+3]).count('G'))/len(''.join(rf[i+1:i+3])),6),''.join(rf[i+1:i+3]).count('C'),''.join(rf[i+1:i+3]).count('G'),len(''.join(rf[i+1:i+3])),''.join(rf[i+1:i+3])) for i in range(0,len(rf),3) if '>' in rf[i]]
CDict = {}
GDict = {}
CGDict = {}
#tmp = {}
tmpS = []
max = 0
counts = 0
C= 0
for i in range(1,len(rf)+1):
    if i == len(rf):
        counts = len(tmpS)
        C += 1
        CDict['C'+str(C)]=''.join(tmpS).count('C')
        print(CDict['C'+str(C)])
        GDict['G'+str(C)]=''.join(tmpS).count('G')    
        print(GDict['G'+str(C)])
        CGDict['CG'+str(C)]= 100*(CDict['C'+str(C)]+GDict['G'+str(C)])/len(''.join(tmpS))
        if CGDict['CG'+str(C)] > max:
            max = CGDict['CG'+str(C)]
            print('！！！',rf[i-1-counts],tmpS,round(CGDict['CG'+str(C)],6))
    else:
        if '>' not in rf[i]:
            tmpS.append(''.join(rf[i]).replace('\n',''))
        else:
            counts = len(tmpS)
            C += 1
            CDict['C'+str(C)]=''.join(tmpS).count('C')
            print(CDict['C'+str(C)])
            GDict['G'+str(C)]=''.join(tmpS).count('G')    
            print(GDict['G'+str(C)])
            CGDict['CG'+str(C)]= 100*(CDict['C'+str(C)]+GDict['G'+str(C)])/len(''.join(tmpS))
            if CGDict['CG'+str(C)] > max:
                max = CGDict['CG'+str(C)]
                print('！！！',rf[i-1-counts],tmpS,round(CGDict['CG'+str(C)],6))
            tmpS = []            
#Counting Point Mutations
s = 'CTTTCCGGGGTTCATAAAAGCCTTACCAATGTTATTGAGTGTACAAGAGCCGCATTCACCATCAGTACCACTATACTGCGCCCGGAGGGCCTATGAATATCTGACGTTTGCACCTCTAGCCCGTAACGCCTGATAACTTACGTGTCGAGCTGCTGAATAGCTTTAAGAGTTCCTCTGGGGACGATTAGCCTCTTGAGAAGCTTTGCGGGTGCGACTGGCATTTCACTGGACGAGACTTGGTGAGCTAGATCTATTCACAGCCGGAATAAAGCATATGTAGAGCAGCCATTACTAGTGGACGACCGATGCAGACATTGCGCTGCGAGATAAGTGCTTCACTATGGCATGACCGTACCAGAGTTTCATGAGCCAAGACAGCGGTTCACCGCCATGAATTGACTTATAAAGGGCTTAAGTGTACATTTGACCACCATAACTATGACTGACATAGTGGGTGCCAATTCAAAATTGCGAAAGGCTACTCCCGCGAGGGAAAACTGACCTGACGGTCACGATATCCCCCGGCACGTCCCTTACATGATTGTCCGGATACTTTAGCCAGAAGAAGGGGAGCCCTCGACGTCTGAACCTGAGCGTTACTGTGGGATGGACATGAACAGGCAATCACCAGATCTCTACTACGTCTGAATTGCCTGAGTGAATGTTGCCCGCCCCGGTAAGTCGAAACTCATATAAAATTCTTCTTATATAAAGGGGCTATCGCCACACCACCTTATGACCGGCCCGGACTTGGCATGAGTTACCCGGCCATGCTTCAACCCTAGCATTAGTAGGGACACAATGCCACCACGTTTATTGTCAGGTGGCTGGCCGGTCCGAGTAGTCAGTCCCCGAGCCGGTCCCGGTGGTCGACCAGGTGACTCTACGGCAGCGTGTTCA'
t = 'TATAGCGGATTGCTGGACTAACACACGTAATTAGACAAGTGTGCAACGGTCCTTTACCTCATCAGTAACAGCAGACGGCGCAGCGGCGTCTACCGAATGTGTTAAATCTCCAACTTTCACCCGTAAAAGCTTCCAGCTTACGGCCCGCCCAGTGGGGGTGGTTGAGCCGCCCATATGTAGCCTTGGCCAGATAGCTCATTATTTCTGGCCTACACTGGTCTATGGGTGGACATAAGTGGGCGAGATAACGCTCGTGAAAGGCTGCGTAAACCATATCAGGTTAAGTCCTAAGAACGGGAAAAATGACGCTTACATATGGATGAGACTTTAGACATTCACTACCTTTCGCACGAACCTCAGAATCAGTGTTCAAAACAGCGGCACGCCCAGTTGTTTCGAATTATTAAGCCACTCACTTTACAACTGACGCCGCCAAATATCACTGCAATATATGGATTCAATGCAATCCGGCGAGTAGTTTCTCCTCCGAAGAAACAGGTGCCCTACAATGTGGTGGTTTTAATGTACGCCCCTTTCTACATTGAACTTTGACCGGTTGCTGAATCACAACAGCCCTCTACGGTTGAGAGAGTGCGAGAATATGCTAGTCAGGTGAGCGGGCCGCAGCGATAACATTGATCCATCTGATTTGCAAGTTTGAATGTTGCCATGCACGCGGAGACTCACTATAGCCTAGACCCATATTAAGACACACGGGTAAGGCACTGCCTCGCTAGGACCCGCACACACTTGATATAACTTCCTCACACTCGCATTCGGGCCTGGTTCGCTCGCACGTGACCTCGTACACAGTCATGGCAATGTCGCTGTTGTCTCCGTGGGGTGAGTCCGTGAGCGTCGCCCGCTGGGCCCCAACGCTAATCTTCTGAAGGGAGCACC'
#print(set(s).difference(set(t)))
#dif = 0
#for i in range(len(s)):
#    if s[i] != t[i]:
#        dif += 1
D= []
[D.append(s[i]) for i in range(len(s)) if s[i] != t[i]] 
print(len(D))         
#Mendel’s First Law（Punnet Squares）
#2HD 2HT 2HR
k=24
m=17
n=26
s = k+m+n
#1-(n/s*(n-1)/(s-1)+m/s*(m-1)/(s-1)/4+m/s*(m-1)/(s-1))
from scipy.special import comb                                                   
(comb(k,2)+k*m+k*n+0.5*m*n+0.75*comb(m,2))/comb(s,2)
#Translating RNA into Protein
import pandas as pd
import numpy as np
def RCSVF(File):
    '''Read a csv file and output its content as a dictionary.Odd columns as keys while even columns as values.'''
    dat= pd.read_csv(File,header = None)
    codon = []
    amino = []
    dictD = {}
    for c in range(np.shape(dat)[1]):
        if c % 2 == 0:
            codon.append(dat.iloc[:,c])
        else:
            amino.append(dat.iloc[:,c])
            for r in range(np.shape(dat)[0]):
                dictD[str(dat.iloc[r,c-1])] = str(dat.iloc[r,c])
    return dictD
inputTAB =r"E:\extdata\DNA codon Std table.csv"
dictCodon = RCSVF(inputTAB)
#'UUU':'F','CUU':'L','AUU':'I','GUU':'V','UUC':'F','CUC':'L','AUC':'I','GUC':'V',
#'UUA':'L', 'CUA':'L','AUA':'I','GUA':'V','UUG':'L','CUG':'L','AUG':'M','GUG':'V',
#'UCU':'S','CCU':'P','ACU':'T', 'GCU':'A','UCC':'S','CCC'：'P','ACC':'T','GCC':'A',
#UCA S      CCA P      ACA T      GCA A
#UCG S      CCG P      ACG T      GCG A
#UAU Y      CAU H      AAU N      GAU D
#UAC Y      CAC H      AAC N      GAC D
#UAA Stop   CAA Q      AAA K      GAA E
#UAG Stop   CAG Q      AAG K      GAG E
#UGU C      CGU R      AGU S      GGU G
#UGC C      CGC R      AGC S      GGC G
#UGA Stop   CGA R      AGA R      GGA G
#UGG W      CGG R      AGG R      GGG G 
file = r"C:\Users\Admin\Downloads\rosalind_prot.txt"
with open(file,'r') as f:
    rf = ''.join(f.readlines()).replace('\n','')
string = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
''.join([dictCodon[rf[i:i+3].replace('U','T')] for i in range(0,len(rf),3)])
#Finding a Motif in DNA
file = r"C:\Users\Admin\Downloads\rosalind_subs.txt"
with open(file,'r') as f:
    subS= f.readlines()[-1].replace('\n','')
with open(file,'r') as f:
    rf = ''.join(f.readlines()).replace('\n','')[0:-len(subS)]
string = 'GATATATGCATATACTT'
subs = 'ATAT'
result = []
for i in range(len(rf)-len(subS)+1):
    if string[i:i+len(subS)] == subS:
        result.append(i+1)
#Consensus and Profile
'''MSA: In “Counting Point Mutations”, we calculated the minimum number of symbol mismatches between two strings of equal length to model the problem of finding the minimum number of point mutations occurring on the evolutionary path between two homologous strands of DNA. If we instead have several homologous strands that we wish to analyze simultaneously, then the natural problem is to find an average-case strand to represent the most likely common ancestor of the given strands.'''
import numpy as np
file = r"C:\Users\Admin\Downloads\rosalind_cons.txt"
with open(file,'r') as f:
    rf = f.readlines()
#[print(rf[i],round(100*(''.join(rf[i+1:i+3]).count('C')+''.join(rf[i+1:i+3]).count('G'))/len(''.join(rf[i+1:i+3])),6),''.join(rf[i+1:i+3]).count('C'),''.join(rf[i+1:i+3]).count('G'),len(''.join(rf[i+1:i+3])),''.join(rf[i+1:i+3])) for i in range(0,len(rf),3) if '>' in rf[i]]
CMtrix = []
GMtrix = []
AMtrix = []
TMtrix = []
tmpS = []
DNA = []
Consensus = []
counts = 0
for i in range(1,len(rf)+1):
    if i == len(rf):
        counts = len(''.join(tmpS))
        DNA.append(list(''.join(tmpS)))
#        DNAM = np.array(DNA).reshape(len(DNA),counts)
        DNAM = np.array(DNA)
        for j in range(np.shape(DNAM)[1]):
            #CMtrix = ''.join(str(DNAM[:,j])).Counter['C']
            CMtrix.append(''.join(str(DNAM[:,j])).count('C'))
            GMtrix.append(''.join(str(DNAM[:,j])).count('G'))
            TMtrix.append(''.join(str(DNAM[:,j])).count('T'))
            AMtrix.append(''.join(str(DNAM[:,j])).count('A'))
            tmp =['C','G','T','A']
            Consensus.append(tmp[np.argmax(np.array([CMtrix[j],GMtrix[j],TMtrix[j],AMtrix[j]]))])
    else:
        if '>' not in rf[i]:
            tmpS.append(''.join(rf[i]).replace('\n',''))
        else:
            DNA.append(list(''.join(tmpS))) 
            tmpS = []            
print(''.join(Consensus))
print('A: '+' '.join([str(s) for s in AMtrix]))
print('C: '+' '.join([str(s) for s in CMtrix]))
print('G: '+' '.join([str(s) for s in GMtrix]))
print('T: '+' '.join([str(s) for s in TMtrix]))
#Mortal Fibonacci Rabbits
'''Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation Fn=Fn−1+Fn−2'''
'''and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.'''
'''Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months.'''
m = 89
dm = 17
c = 1
p = np.repeat(0, dm-1)
#k = 1 birth rate 
for n in range(m-1):
    c,p[1:],p[0] = sum(p),p[0:dm-2],c #!!!!!!!order
#    print(n+2, c, p, sum(p), c+sum(p))
    print(n+2, c, sum(p), c+sum(p))

 
 


#Overlap Graphs !!! surfix and prefix instead of prefix then surfix 
import numpy as np
#file = r"E:\Rosalind\Rosalind_t3.txt"
file = r"C:\Users\Admin\Downloads\rosalind_grph.txt"
with open(file,'r') as f:
    rf = f.readlines()
#[print(rf[i],round(100*(''.join(rf[i+1:i+3]).count('C')+''.join(rf[i+1:i+3]).count('G'))/len(''.join(rf[i+1:i+3])),6),''.join(rf[i+1:i+3]).count('C'),''.join(rf[i+1:i+3]).count('G'),len(''.join(rf[i+1:i+3])),''.join(rf[i+1:i+3])) for i in range(0,len(rf),3) if '>' in rf[i]]
startc = []
endc = []
tmpS = [] 
# method1 
 
def is_k_overlap(s1, s2, k):
    return s1[-k:] == s2[:k]
  
import itertools
  
def k_edges(k):
    edges = []
    for u, v in itertools.combinations(range(len(startc)), 2):  # data 里面任意取两个比较
        u_dna, v_dna = endc[u], startc[v]
        if is_k_overlap(u_dna, v_dna, k):
#            print(u_dna+'=='+v_dna)
            print(rf[u*(rows+1)].replace('\n','').replace('>',''),rf[v*(rows+1)].replace('\n','').replace('>',''))
            edges.append((rf[u*(rows+1)].replace('\n','').replace('>',''),rf[v*(rows+1)].replace('\n','').replace('>',''))) 
        u_dna, v_dna = startc[u], endc[v]
        if is_k_overlap(v_dna, u_dna, k):
#            print(v_dna+'=='+u_dna)
            print(rf[v*(rows+1)].replace('\n','').replace('>',''),rf[u*(rows+1)].replace('\n','').replace('>',''))
            edges.append((rf[v*(rows+1)].replace('\n','').replace('>',''),rf[u*(rows+1)].replace('\n','').replace('>','')))
    return edges
 

for i in range(1,len(rf)+1):
    if i == len(rf):
        rows = len(tmpS)
        startc.append(''.join(tmpS)[0:3])
        endc.append(''.join(tmpS)[-3:])
#        DNAM = np.array(DNA).reshape(len(DNA),counts)
#        [print(rf[startc.index(s)*(rows+1)],rf[endc.index(s)*(rows+1)]) for s in startc if (s in endc and startc.index(s)!=endc.index(s))]
#        p =  [print(rf[endc.index(startc[j])*(rows+1)].replace('\n','').replace('>',''),rf[j*(rows+1)].replace('\n','').replace('>','')) for j in range(len(startc)) if (startc[j] in endc and j!=endc.index(startc[j]))]
        p = [(rf[endc.index(startc[j])*(rows+1)].replace('\n','').replace('>',''),rf[j*(rows+1)].replace('\n','').replace('>','')) for j in range(len(startc)) if (startc[j] in endc and j!=endc.index(startc[j]))]                              
        q =  [(rf[j*(rows+1)].replace('\n','').replace('>',''),rf[startc.index(startc[j])*(rows+1)].replace('\n','').replace('>','')) for j in range(len(endc)) if (endc[j] in startc and j!=startc.index(endc[j]) and (rf[j*(rows+1)].replace('\n','').replace('>',''),rf[startc.index(startc[j])*(rows+1)].replace('\n','').replace('>','')) not in p)]
#        [print(rf[j*(rows+1)].replace('\n','').replace('>',''),rf[startc.index(startc[j])*(rows+1)].replace('\n','').replace('>','')) for j in range(len(endc)) if (endc[j] in startc and j!=startc.index(endc[j]) and (rf[j*(rows+1)].replace('\n','').replace('>',''),rf[startc.index(startc[j])*(rows+1)].replace('\n','').replace('>','')) not in p)]
        print(len(p)+len(q))
        print(list(enumerate(k_edges(3)))) #method1
    else:
        if '>' not in rf[i]:
            tmpS.append(''.join(rf[i]).replace('\n',''))
        else:
            startc.append(''.join(tmpS)[0:3])
            endc.append(''.join(tmpS)[-3:])
            tmpS = []
 
#Calculating Expected Offspring
#AA-AA
#AA-Aa
#AA-aa
#Aa-Aa
#Aa-aa
#aa-aa
def offspringA(a,b,c,d,e,f):
    return 2*(a*1+b*1+c*1+d*3/4+0.5*e+0*f)
print(offspringA(19381,16346,19892,18229,19936,19172))

#Finding a Shared Motif
import numpy as np
#file = r"C:\Users\Admin\Downloads\rosalind_cons.txt"
file = r'E:\Rosalind\Rosalind_t4.txt'
with open(file,'r') as f:
    rf = f.readlines()
#[print(rf[i],round(10;0*(''.join(rf[i+1:i+3]).count('C')+''.join(rf[i+1:i+3]).count('G'))/len(''.join(rf[i+1:i+3])),6),''.join(rf[i+1:i+3]).count('C'),''.join(rf[i+1:i+3]).count('G'),len(''.join(rf[i+1:i+3])),''.join(rf[i+1:i+3])) for i in range(0,len(rf),3) if '>' in rf[i]]
#SMotif={}
SMotif=[]
S2 = []
tmpS = []
DNA = []
counts = 0
for i in range(1,len(rf)+1):
    if i == len(rf):
        DNA.append(''.join(tmpS))
        counts = len(''.join(tmpS))
        S2.append(''.join(tmpS))
        for s in SMotif:
        #    flag = False
            for p in range(len(str(s))):
                ttt = []
                ttti = p
                for q in range(counts):
                    if ttti == len(str(s)):
                        break
                    if ''.join(tmpS)[q] == str(s)[ttti]:
                        ttt.append(''.join(tmpS)[ttti])
                        ttti += 1
                #           flag = True
                if q == counts-1 or ttti == len(str(s))-1 and len(ttt)>1:
#                    print('-----------------'+str(p)+' '+str(q)+' '+str(len(S2))+'-----------------')                                                                
                #           SMotif[''.join(ttt)]=len(''.join(ttt))
                    SMotif.append([''.join(ttt)])# = len(''.join(ttt)
                        
                        #if flag == False:
                        #    S2[s] = 0
                        #    maxN = 0
        for j in range(len(SMotif)):
                    #somehow?
                    #if j==len(SMotif):
                    #    break
            if len(str(SMotif[j])) == counts:
                SMotif.pop(j)
           # if S2 != SMotif:
           #     SMotif = S2
        tmpS = []            
#        [print('Result:'+keys) for keys,vals in SMotif.items() if vals == max(SMotif.values())]           
        m = 0
        for s in SMotif: 
            if len(s[:]) > m:
                m = len(str(s[:]))
                print('Result:'+str(s[:]))
    else:
        if '>' not in rf[i]:
            tmpS.append(''.join(rf[i]).replace('\n',''))
        else:
#            print(tmpS)
            DNA.append(''.join(tmpS))
            counts = len(''.join(tmpS))
            if SMotif == [] and tmpS:
#                SMotif[''.join(tmpS)]=0
                SMotif.append(''.join(tmpS))
#                S2=SMotif
#            for s in SMotif.keys():
            else:
                S2.append(''.join(tmpS))
                for s in SMotif:
                #    flag = False
                    for p in range(len(str(s))):
                        ttt = []
                        ttti = p
                        for q in range(counts):
                            if ttti == len(str(s)):
                                break
                            if ''.join(tmpS)[q] == str(s)[ttti]:
                                ttt.append(''.join(tmpS)[ttti])
                                ttti += 1
            #                    flag = True
                            if q == counts-1 or ttti == len(str(s))-1 and len(ttt)>1:
#                                print('-----------------'+str(p)+' '+str(q)+' '+str(len(S2))+'-----------------')                                                                
            #                        SMotif[''.join(ttt)]=len(''.join(ttt))
                                SMotif.append([''.join(ttt)])# = len(''.join(ttt)
                            
                        #if flag == False:
                        #    S2[s] = 0
                        #    maxN = 0
                for j in range(len(SMotif)):
                    #somehow?
                    #if j==len(SMotif):
                    #    break
                    if len(str(SMotif[j])) == counts:
                        SMotif.pop(j)
           # if S2 != SMotif:
           #     SMotif = S2
            tmpS = []            

#method2
def readfasta(filename, sample):
    fa = open(filename, 'r')
    fo = open(sample, 'w')
    res = {}
    Res = []
    ID = ''
    for line in fa:
        if line.startswith('>'):
            ID = line.strip('\n')
            res[ID] = ''
        else:
            res[ID] += line.strip('\n')
 
    for key in res.values():
        Res.append(key)
        fo.write(key + '\n')
    return Res
 
 
def fragement(seq_list):
    res = []
    seq = seq_list[0]
    for i in range(len(seq)):
        s_seq = seq[i:]
        #print s_seq
        for j in range(len(s_seq)):
            res.append(s_seq[:(len(s_seq) - j)])
            #print res
 
    return res
 
 
def main(infile, sample):
    seq_list = readfasta(infile, sample)   #['TAGACCA','ATACA','GATTACA']
    frags = fragement(seq_list)

main(r"C:\Users\Admin\Downloads\rosalind_lcsm.txt",r'E:\Rosalind\Rosalind_ttt.txt')
#Independent Alleles with Aa-Bb 
import itertools
#Punnett Square sees: 
#https://study.com/academy/lesson/mendels-second-law-the-law-of-independent-assortment.html#:~:text=Mendel's%20Second%20Law%20also%20called,separately%20from%20other%20genetic%20traits.
from math import *

def C(n,i):
    for j in range(1,i):
        n = n*(n-j)
    return n

def MendelSecondIA(k,N):
    p = []    
    child_num = 2**k
    for i in range(N,child_num+1):
#        p.append(len(list(itertools.combinations([x for x in range(child_num)],i)))*(0.25**i)*(0.75**(child_num-i)))    
        p.append(C(child_num,i)/math.factorial(child_num)*(0.25**i)*(0.75**(child_num-i)))
    return sum(p)
print(MendelSecondIA(7,33))

#out of memory problem
from math import *
def C(n,i):
    for j in range(1,i):
        n = n*(n-j)
    return n

def MendelSecondIA(k,N):
    p = []    
    child_num = 2**k
    for i in range(N):
#        p.append(len(list(itertools.combinations([x for x in range(child_num)],i)))*(0.25**i)*(0.75**(child_num-i)))    
        p.append(comb(child_num,i)*(0.25**i)*(0.75**(child_num-i)))
    return 1-sum(p)
print(MendelSecondIA(7,36))

#Finding a Protein Motif
#A component of a domain essential for its function is called a motif, a term that in general has the same meaning as it does in nucleic acids, although many other terms are also used (blocks, signatures, fingerprints, etc.) Usually protein motifs are evolutionarily conservative, meaning that they appear without much change in different species.
#https://rosalind.info/problems/mprt/
#To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows: [XY] means "either X or Y" and {X} means "any amino acid except X." For example, the N-glycosylation motif is written as N{P}[ST]{P}.
#!!!!!!Notice: extra spaces at end of line will not be passed

#input = r"E:\Rosalind\BSZC00.fasta.txt"
#input = r"E:\Rosalind\P07204.fasta.txt"
import requests
#input = r"E:\Rosalind\Rosalind_t5.txt"
input = r"C:\Users\Admin\Downloads\rosalind_mprt.txt"

rf = []
with open(input,'r') as f:
    rf = f.readlines()
downloadList = []
names = []
for s in rf:
    names.append(s[0:s.index('\n')])
    if '_' in s:
        downloadList.append(''.join(s[0:s.index('_')]))
    else:
        downloadList.append(''.join(s[0:s.index('\n')]))

localDir = 'E:\Rosalind\\'
header = {"User-Agent":"Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400"}
for urls in downloadList:
    req = requests.get(url ='https://rest.uniprot.org/uniprotkb/'+urls+'.fasta',headers=header)    
    req.encoding = 'gbk'
    html = req.text
    with open(localDir+urls+'.txt','w') as file:
        file.write(html) 

tmpProtein = []
for urls in downloadList:
    with open(localDir+urls+'.txt','r') as f:
        rfs=f.readlines()
    tmp=[]
    for lines in rfs:
        if '>' not in lines:
            tmp.append(''.join(lines[0:lines.index('\n')]))
    tmpProtein.append(''.join(tmp))
output = []
IDs = []
#answ = [[85,118,142,306,395],[47,115,116,382,409],[79,109,135,248,306,358,364,402,485,501,614]]
for i in range(len(tmpProtein)):
    if tmpProtein[i].find('N') >=0:
        tmp = []
        for j in range(len(tmpProtein[i])-3):    
            if tmpProtein[i][j] == 'N'  and tmpProtein[i][j+1] != 'P' and (tmpProtein[i][j+2] == 'S' or tmpProtein[i][j+2] == 'T') and tmpProtein[i][j+3] != 'P':
                tmp.append(j+1)
        if len(tmp)>0:
            IDs.append(downloadList[i])
            output.append(' '.join([str(s) for s in tmp]))
#            print('Difference to Answer and surrounding:')
            print(names[i])
#            [print(s) for s in tmp]
            print(' '.join([str(s) for s in tmp]))
#                print('Answ: ',answ[i-1])
#                print('Mine: ',tmp)
#                [print('Included: ',tmpProtein[i][k-2],tmpProtein[i][k-1],tmpProtein[i][k],tmpProtein[i][k+1]) for k in range(len(tmp)) if k in answ[i-1]]
#                [print('Excluded: ',tmpProtein[i][k-2],tmpProtein[i][k-1],tmpProtein[i][k],tmpProtein[i][k+1]) for k in range(len(tmp)) if k not in answ[i-1]]
print(output)
print(answ)        

#Inferring mRNA from Protein
import pandas as pd
import numpy as np
def RCSVF(File):
    '''Read a csv file and output its content as a dictionary.Odd columns as keys while even columns as values.'''
    dat= pd.read_csv(File,header = None)
    codon = []
    amino = []
    dictD = {}
    for c in range(np.shape(dat)[1]):
        if c % 2 == 0:
            codon.append(dat.iloc[:,c])
        else:
            amino.append(dat.iloc[:,c])
            for r in range(np.shape(dat)[0]):
                dictD[str(dat.iloc[r,c-1])] = str(dat.iloc[r,c])
    return dictD
inputTAB =r"E:\extdata\DNA codon Std table.csv"
dictCodon = RCSVF(inputTAB)
string = 'MVWVPVMKVWNEWYQAIVKEVPRATHNFYIRGIDKQGTYSWLQKAPQCLTTRQQPGFNGKNRTAMPHRDGQQNIFQWELQNLCCTEMMAIVIHGFWGWFMMMWNWCEIAGWDYRQISCPQFNLGVHKNWNDGWCNMPYTHTDQPQLVFVKLTEVEYMHQREIFSMWEIDEERWLWRVYEHMMFYKWMPQKWSLCMQDSSYKMCHNGEYAEEHHWCHNAMVCECMYNPDQLMYSNPIRECNLMGNSKYNSIHHSSSLIPSFDDTDRNWVRKSLWWDQLTWFLYDCHRQLFSCNNGHSAQFMFFPNPSQPFNCHASILYRAPMNLEKYAKKTRVMVYFCVPVPFCFRHDQHVRYRQKCKWETDAHHFIDLVQRTHWDKIRPPCELPVWNMRFEPYMPQENTGSQIPMRRWGAREFHFAHFYYKNEVTIKCTGMFVRAAMCWNQICSRSFSKHMHWKKCAVSNAWMWREYKGVHSKRDASKYWKVWMFKKHVHTQGCTTKPRQYWTMYYIRCARDHHHFGPHYDEENPNAWEGIHIKAWVREFHAWIWCMTSCGAVYQNIYMRNKEACNIQEVPEWNHHIHSQYNMSYTANPFKITWGAYQLPKVAVPSDEIFNKAHKTAQPKWCTIICSFYIHAPLSASSEWEPTKEYWPMRWGDEACILGDVPVMAHFVEWSHQLQMGANYFYIHKPIDMFFVGWEWQYSFVRNVITFTEDKKPDHLVNRVIKPDMFSPNNIEMAKSDNPPEFFPIMCHQRIDMSRGRKEKESNTKVAGQTDTSWGIPLLEHGNAEHRAAKPHEERTATTSCLDQWAISVFMNMYEYRQGKKLQAVAAVPVDSMLDPKKWSWIIGPLWKDKHPPDIEQRQCPEMMLTVRECNARHYRGCIEIMATWRNQTFLWKECEKFMFVSDKLGAMNYPRNMCWVSMLCIVCTSYSQFGYFATCIWNNPFHNVYSFCHFPRMFKERASISHNKSAASCLDTYKSHWEMIRSGLEGRFAWYRGVDHML'
outputN = 1
for s in string:
    outputN *= list(dictCodon.values()).count(s) 
outputN *=  list(dictCodon.values()).count('*')
print(np.mod(outputN,1000000))

#Open Reading Frames
#Notice the one passed is method3 which basically searches on DNA level instead of interval of three to translate to amino acids first. 
import pandas as pd
import numpy as np
string = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
def RCSVF(File):
    '''Read a csv file and output its content as a dictionary.Odd columns as keys while even columns as values.'''
    dat= pd.read_csv(File,header = None)
    codon = []
    amino = []
    dictD = {}
    for c in range(np.shape(dat)[1]):
        if c % 2 == 0:
            codon.append(dat.iloc[:,c])
        else:
            amino.append(dat.iloc[:,c])
            for r in range(np.shape(dat)[0]):
                dictD[str(dat.iloc[r,c-1])] = str(dat.iloc[r,c])
    return dictD

inputTAB =r"E:\extdata\DNA codon Std table.csv"
dictCodon = RCSVF(inputTAB)
RC = {'A':'T','T':'A','C':'G','G':'C'}
inputf = r"C:\Users\Admin\Downloads\rosalind_orf.txt"
with open(inputf,'r') as f:
    rf = f.readlines()
string = ''.join(rf[1:]).replace('\n','')
if np.mod(len(string),3) >0:
    string = string[0:len(string)-np.mod(len(string),3)]
def ORFF(string,startP,endP):
    ORFP = []
    tmpORF = []
    for i in range(startP,endP-3,3):
        if dictCodon[string[i:i+3]] == 'M' and tmpORF == []:
            tmpORF.append(dictCodon[string[i:i+3]])
        if dictCodon[string[i:i+3]] != 'M' and dictCodon[string[i:i+3]] != '*' and tmpORF != []:
            tmpORF.append(dictCodon[string[i:i+3]])
#        if dictCodon[string[i:i+3]] == 'M' and tmpORF != []:
#            tmpORF=[]
#            tmpORF.append(dictCodon[string[i:i+3]])
        if dictCodon[string[i:i+3]] == '*' and tmpORF != []:
            ORFP.append(''.join(tmpORF))
            tmpORF=[]
 #       if i == endP-1 and tmpOEF != []:
 #           tmpORF = []
    return ORFP
print(ORFF(string,0,len(string)))
print(ORFF(string,1,len(string)-2))
print(ORFF(string,2,len(string)-1))
print(ORFF(''.join([RC[s] for s in string[::-1]]),0,len(string)))
print(ORFF(''.join([RC[s] for s in string[::-1]]),1,len(string)-2))
print(ORFF(''.join([RC[s] for s in string[::-1]]),2,len(string)-1))

#method2
def ORFF2(seq):
    try:
        stti = seq.index('M')
        stpi = seq.index('*')
        if stti < stpi:
            print(''.join(seq[stti:stpi]))
#            print(ps)
            ORFF2(seq[stpi+1:])
        if stti > stpi:
            ORFF2(seq[stti:])
    except:
        return ''
AA1=''.join([dictCodon[string[i:i+3]] for i in range(0,len(string)-3,3)])
AA2=''.join([dictCodon[string[i:i+3]] for i in range(1,len(string)-3-2,3)])
AA3=''.join([dictCodon[string[i:i+3]] for i in range(2,len(string)-3-1,3)])
AA4=''.join([dictCodon[''.join([RC[s] for s in string[::-1]])[i:i+3]] for i in range(0,len(string)-3,3)])
AA5=''.join([dictCodon[''.join([RC[s] for s in string[::-1]])[i:i+3]] for i in range(1,len(string)-3-2,3)])
AA6=''.join([dictCodon[''.join([RC[s] for s in string[::-1]])[i:i+3]] for i in range(2,len(string)-3-1,3)])

ORFF2(AA1)
ORFF2(AA2)
ORFF2(AA3)
ORFF2(AA4)
ORFF2(AA5)
ORFF2(AA6)

#method3
def ORFF3(seq):
    pSeq = []
    for i in range(len(seq)-3):
        if dictCodon[seq[i:i+3]]=='M':
            tmpSeq = [dictCodon[seq[ii:ii+3]] for ii in range(i,len(seq)-3,3)]
            try:
                tmpk = ''.join(tmpSeq).index('*')
                pSeq.append(''.join(tmpSeq[0:tmpk]))
            except:
                print('End of searching')
                break
    return pSeq

O3 = ORFF3(string)
O3RC = ORFF3(''.join([RC[s] for s in string[::-1]]))
[print(s) for s in set(O3).union(set(O3RC))]
#Enumerating Gene Orders
import itertools
import math
n = 7
#[print(p) for p in list(itertools.permutations(range(1,n+1)))]
print(math.perm(n))
L = []
for i in list(itertools.permutations(range(1,n+1))):
    tmp = []
    for ii in i:
        tmp.append(str(ii))
    print(' '.join(tmp))
    L.append(' '.join(tmp))
    
#
MassTbl = {'A':71.03711,
'C':103.00919,
'D':115.02694,
'E':129.04259,
'F':147.06841,
'G':57.02146,
'H':137.05891,
'I':113.08406,
'K':128.09496,
'L':113.08406,
'M':131.04049,
'N':114.04293,
'P':97.05276,
'Q':128.05858,
'R':156.10111,                                                                                                                                          
'S':87.03203,
'T':101.04768,
'V':99.06841,
'W':186.07931,
'Y':163.06333}
string = 'KMAKGDQFNCSLDGHTKVRLFKKMHPSENKAPFEEMEQLGFILWSCKMLHSVQVILTQTMCDCWDYMLTIITRAVKWFFRREIFSFRIDTQPLATCMHTDHLVQPSQFTAEIQWQSIGTEFSMVWTSMPNYPLMMQLKHVEWTYIMFPCCCALICVCWWYSDHQIVDAPPYQNQLWSGSGHEIFDYHVWQRPAYNSKMYFSDYYACRNPYTKWALADNRRVIDRPQAWPMLTNAVYHRFTFPMQKDQSWDSHGAYELEANCLCVRNCHPNLPMSFTQPITYYFGTPKGYIEDVDTQNSTDPTADAKVECRASPCVNCTRDIKTHGQQRTGIAKPTMKCFNAWEPLIKLFVGRKFGKMDVFFNYYVTNGIVNTTECSTIVWNMTSLTYEHAWCVPTNGSAARYRQQNKCMRQDCKPVMECVHWVNGVEGREMDHCYDMVHYHCHATKSADFLFKVPWPRNDLTNVFYMKNIFGMMSDDQWEHNYRHYRNVSEHPRNAITCYYTKPGHEQIIITFERPPWEFHCQFTTNQCFWPFGIFRDKLKMIKVWGEYVWCNSCATSTNRVDKAATKLNEPKGRREYEVMFKNSQITHCVSMDCAQIRCCPGDPWYWGCRLKFLQSDTQDEVYDNTGAFDIKWVWYDLTTWFLWRISCAMPNNLQLPTYFIKNQIHHMNPPWSHIAGAVALPWFGHQPRDKTHQNYCERYWWHWIMEPYLPGSMNHNHKHSWTATQEEQWDPFLKRCREACGLMPIPAQFDLGTRTSIASWQSMVRSRVYTCMCPLTMFWETGNIDLAAVMHERMTQYGGLGFTYG'
m = sum([MassTbl[s] for s in string])
    
#Locating Restriction Sites
string = 'TCAATGCATGCGGGTCTATATGCAT'
string = 'GCGCCTCGCCCAGAGGTAAGTCAAGCTCAAGCATACGGATCACGTTAATCTACCTCCCCTGTAGATTAGGTGAAGTGAACCACCAGCTTTTCGTAGACTGTTGTCCCTCGTCATCAGTTACAAGTAGTGTTACGAGGCTTACCAAGCATCAGTACAAGGGTGTGGTCTGGGGAAGAACTGACGGACTCAAATGAATATAGTAGCGCGGAGAAGCAGAGCGGAATATTTCGGACCAGTCGGGTGACTGTAGCAACGATGAGGTGCGAGGCGGGGCCAAAGCTATCATCGAAATGCGAATCTATCTGGTTTAGAGCCATCTCGATTTTGCGGCAATACAACGTAAGGCGGAGCATTTTTTCATTGCGGCTCATTATTGTATCCTCGAGACTCCGGCTTTGGCGTCGTTTGCGAAGCACCAATTGAAATGCTTACTACCCTTGCGTCTTCAAGGTCTGCATGCAGAGCTCCTTGGGAAGCTTATCTTCAGAAGGGATGAGTAACACATACGCGGAAAGACCTAAAAACAAGTCAGTAATTGTGCCCGTTCACCTCTACCCAAGAGACCTGTAGTTCCGCCTTGGTTAACGTACTCTAACCTTTTTGTCCCACCTACGACCCTGTACGGCGCTTCGGTATCTGTTAGACTGCCAATGACTTGTGCAACAGAGAAACGTGGCCGCGGAACTGTGCTGCGCCGTTAGACTTTCCCATGGTTGACTTCCTAACTGGGATGGATCACATGTATCAAGCACAGATAGCTAGTGGCTGTAGTGCCACGAAAAGCAGCAAAAACTTGAGACGCTTGTAACCACAAGCGGGGACATTTAACCTCACGAATATATGACACTGATTTGAGAGTCCATAGATCACCGAGTTGGGTCGAGGGGGAGTAAATACTACACGACGTGCGCGCG'
TblDNA = {'A':'T','T':'A','C':'G','G':'C'}
stringCR = ''.join([TblDNA[s] for s in string])[::-1]
minL = 4
maxL = 12

P = []
L = []
asw = []
for i in range(len(string)-minL):
    for j in range(len(stringCR)-minL):
        if string[i:i+minL] == stringCR[j:j+minL]:
            P.append(i+1)
            for k in range(j,min(j+maxL,len(stringCR))):
                temp = k-j
                if (i+temp < len(string) and string[i+temp] != stringCR[k]) or (i+temp == len(string) and string[i+temp-1] == stringCR[k-1]): 
                    L.append(temp)
                    asw.append((i+1,temp))
                    print(i+1,temp)
                    print(i+1,j+1,temp)
                    print(string[i:i+temp],stringCR[j:k])
                    break
asw = sorted(list(set(asw)))
og = [(4,6),(5,4),(6,6),(7,4),(17,4),(18,4),(20,6),(21,4)] 
 
#should not check reverce palindrom in the whole string. instead, check every each pieces 
#method for passing (Note here it doesn't search the position in reserve palindrome different then in origin string)
for i in range(len(string)):
    for j in range(minL,maxL):
        if string[i:i+j]==''.join([TblDNA[s] for s in string[i:i+j]][::-1]) and i+j<=len(string):
            print(i+1, " ",j)
#        if (i,j) in og:
#            print(string[i-1:i-1+j],stringCR[i-1:i-1+j])

# Online Python - IDE, Editor, Compiler, Interpreter
string = 'TCAATGCATGCGGGTCTATATGCAT'
string = 'TATAAACCCTTAATCGGTTTACTATAGATTAATGTGCTTTTGCACTCACAGTAATTCCGTGTGGATAGTTCAGTTCTGACACGTAACGTCGGAGAGTACACCGTCGTTGAGCGATACTGCGCAGGTAAAATTCTTAAGTACCGCACTCGGAGGGACCTGAAGCATCCGGCATGTGATCTTGACCCTGTGTCTCGATTCGTTAGCAGACGACTTGTATATGTGCGCTTCTGACGGAAAGATAATAGAGTCCTAACGATATGACTAAAATACCTCCATCAGTTGCCCACGCCGCAATACCTGGGGCGCTTCCCGGTTGCAAACTCTCATACCCGGCACGAAGGCCTGTCACCAAAAGGCATGGATGCAACCCGTCCATTGGCCCCCTATGTACCGAACTCGCTAGAATTTGCAGGGGTATAGGTAAGCACAAAACCTAGGGTACATGTACCCCGTCGGGTATTACTTGGTAACAAATGCCGTCACGGTTCTAAAGTATAGGTTACTCGTAGGGGCGTTTGCAAAGACCAAATTCGGTGGCCGTGAAAAAAAAATCGCCGGCCGATGCCTATAAGGGGCTTGAAACCCAGGACTCGCGAGTGGCGTCTAATTATCAGTATTACAGCAGGGAGGCGAGCTGCTCAAGCAGCGAATGGCCCCCTCCACTCTCGTCATACGCTATCAGCCTCGATGTACAGGGAGAATGTTTTTAAATTCAAGGAGCCAAATTCCGCTGATAAGGCTTAGGTGGATGCTCACTGATAACGTCCCAGACACCCAACCAAGCGGTAGTGGTAAACCAGGACGCGGAAAGCCAGCGTCACAGGACAGCGACTTTTCGCACCCTGACAGAGGGTCAATGCAAGCTCCCTATTGGTGCGTTAAGCGC'
TblDNA = {'A':'T','T':'A','C':'G','G':'C'}
stringCR = ''.join([TblDNA[s] for s in string])[::-1]
minL = 4
maxL = 12
for i in range(len(string)):
    for j in range(minL,maxL+1):
        if string[i:i+j]==''.join([TblDNA[s] for s in string[i:i+j]][::-1]) and i+j<=len(string):
            print(i+1,j,sep = '\t')
            
#RNA Splicing
#Genes are Discontiguousclick to expand;  Problem:After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.
import pandas as pd
import numpy as np
def RCSVF(File):
    '''Read a csv file and output its content as a dictionary.Odd columns as keys while even columns as values.'''
    dat= pd.read_csv(File,header = None)
    codon = []
    amino = []
    dictD = {}
    for c in range(np.shape(dat)[1]):
        if c % 2 == 0:
            codon.append(dat.iloc[:,c])
        else:
            amino.append(dat.iloc[:,c])
            for r in range(np.shape(dat)[0]):
                dictD[str(dat.iloc[r,c-1])] = str(dat.iloc[r,c])
    return dictD
inputTAB =r"E:\extdata\DNA codon Std table.csv"
dictCodon = RCSVF(inputTAB)
inputf = r"E:\Rosalind\Rosalind_t6.txt"
inputf = r"C:\Users\Admin\Downloads\rosalind_splc.txt"
with open(inputf,'r') as f:
    rf = f.readlines()
exons = []
introns = []
i = 1
while i <len(rf):
    if '>' not in rf[i] and len(introns)==0: 
        exons.append(rf[i].replace('\n',''))
        i += 1 
    elif '>' in rf[i]:
        temp = rf[i+1].replace('\n','')
        for j in range(i+2,len(rf)):
            if '>' not in rf[j]:
                temp += rf[j].replace('\n','')
            else:
                introns.append(temp)
                i = j
                temp = ''
        if j == len(rf) - 1 or i == len(rf) - 2:
            introns.append(temp)
            break

                
Exons = ''.join(exons).replace('\n','')
#introns = sorted(introns, reverse=True)
def delete_introns(exons, introns):
    newexons = ''
    for s in introns:
        if exons.find(s)>0:
            print(len(exons),s,"found............")
            newexons = exons.replace(s,'')
            exons = newexons
            print(len(exons))
        else:
            print(len(exons),s," not found-----------")
            newexons = exons
    return newexons
old_exons = Exons
newExons = delete_introns(old_exons, introns)
len(newExons)
len(delete_introns(newExons, introns))    
    
protein = []
if np.mod(len(newExons),3) >0:
    newExons = newExons[0:len(newExons)-np.mod(len(newExons),3)]
for i in range(0,len(newExons),3):
    if dictCodon[newExons[i:i+3]] == '*':
        break
    else:
        protein.append(dictCodon[newExons[i:i+3]])
print(''.join(protein))

#        Enumerating k-mers Lexicographically
#len(letters)**level
string = 'A B C D E '
levels = 4
letters = sorted(string.replace(' ',''))
#m1
orderedS = []
def orderS(letter,lev):
    D = ''.join(letter)
    for i in range(len(letter)):        
        D.replace(letter[i],letter[i]*(len(letter)**lev))
        print(D,letter[i],'replaced with',letter[i]*(len(letter)**lev))
    print('D After replace is',D, lev,'as in lev')
    return D            

for level in range(levels):
    print('letters',letters)
    print('unilen',len(letters),'1stdigit',letters[0])
#    print(orderS(letters,level,np.array(digits)))
#    orderedS.append(''.join(orderS(letters,level,np.array(digits,dtype = str))*len(letters)))
    orderedS.append(orderS(letters,levels-1-level)*(len(letters)**level))
    
#np.transpose(np.array(orderedS,dtype = str))
for i in range(len(orderedS[0])):
    print(orderedS[0][i],orderedS[1][i])
#m2
llll= []
for level in range(levels):
    cstring = string
#    print(cstring)
    for l in letters:
#        print(l,cstring.replace(l+' ',l*(len(letters)**(levels-1-level))))
        cstring = cstring.replace(l+' ',l*(len(letters)**(levels-1-level)))
#    print(level,cstring*((len(letters)**level)))
    llll.append(cstring*(len(letters)**level))
#np.transpose(np.array(llll,dtype = str))
result = []
for p in range(len(llll[0])):
    temp = []
    for q in range(levels): 
#        print(llll[q][p],end='')
        temp.append(llll[q][p])
    result.append(''.join(temp))
    print(result[-1])

#Longest Increasing Subsequence
string = '51423'
#54123 541 123
#51243 543 124
#!!!!!!54213 5421 23
n=5
li = []
ld = []

#m1
def sort_iter(s,orderRev = False):
    lt = []
    IDS = []
    sorted_s = sorted(enumerate(s),key=lambda x:x[1],reverse =orderRev)
    idx =[i[0] for i in sorted_s]
    lt.append(string[idx[0]])
    IDS.append(idx[0])
    temp = idx[0]        
    for i in range(len(idx)-1):
        if idx[i] < temp:
            break
        lt.append(string[idx[i]])
        temp = idx[i]
        IDS.append(idx[i])
    return lt

sorted_string = sorted(enumerate(string),key=lambda x:x[1])
idx =[i[0] for i in sorted_string]        
for i in range(len(idx)-1):
    li.append(string[idx[i]])
    if idx[i+1] - idx[i]>=3:
        temps = range(idx[i+1],len(idx)-1)
        stemps = sort_iter(temps,orderRev = False)
        tempsin = range(idx[i],idx[i+1])
        stempsin = sort_iter(tempsin,orderRev = False)        
        if len(stemps) > len(stempsin):
            [li.append(string[s]) for s in stemps]
        else:
            [li.append(string[s]) for s in stempsin]
    elif idx[i+1] - idx[i]==2:
        temps = idx[i]+1
        if string[temps] > string[idx[i+1]]:
            li.append(string[temps])
    elif idx[i+1] < idx[i]:   
        break
#if li[-1]<string[idx[-1]]:
#    li.append(string[idx[-1]])
dsorted_string = sorted(enumerate(string),key=lambda x:x[1], reverse=True)
didx =[i[0] for i in dsorted_string]        
for i in range(len(didx)-1):
    ld.append(string[didx[i]])
    if didx[i+1] - didx[i]>=3:
        dtemps = range(didx[i+1],len(didx)-1)
        dstemps = sort_iter(dtemps,orderRev = True)
        dtempsin = range(didx[i],didx[i+1])
        dstempsin = sort_iter(dtempsin,orderRev = True)        
        if len(dstemps) > len(dstempsin):
            [ld.append(dstring[s]) for s in dstemps]
        else:
            [ld.append(dstring[s]) for s in dstempsin]
    elif didx[i+1] - didx[i]==2:
        dtemps = didx[i]+1
        if i < n-2:
            if string[dtemps] > string[didx[i+2]]:
                ld.append(string[dtemps])
        elif i == len(didx)-2 and string[idx[i]] < string[idx[i-1]]:
            ld.append(string[didx[i]])            
        else:
            ld.append(string[didx[i+1]])
    elif didx[i+1] < didx[i]:   
        break
    
#m2
#string = '51423'
#54123 541 123
#51243 543 124
#!!!!!!54213 5421 23


#string='8581 8400 7967 3634 8504 6466 7803 503 1023 2600 7282 6788 3008 5888 7214 7600 4847 3750 851 6402 6188 7269 8410 2921 1318 65 7842 8431 914 6647 2077 4752 1929 2841 2892 3972 3689 7736 760 8221 5327 3772 6277 4461 7010 190 2141 4728 1560 5880 2501 5640 4008 490 7455 4088 4311 4851 4423 6835 4184 4287 1322 5452 4993 6684 5944 2674 3446 4695 1895 3820 7957 2663 4522 1225 4657 4058 3707 8618 5082 6286 7057 7856 2922 2631 5519 4785 8524 2294 206 1404 2833 8130 8580 8751 8378 7642 6142 5243 316 5033 1860 4007 7780 5046 178 6809 531 996 2123 3447 8289 6543 1573 8672 2118 4619 6374 3653 6583 75 8113 8655 5570 6594 2946 356 5904 466 7920 1178 2149 6677 7923 6690 4509 1394 4827 3373 5657 4017 7593 4825 5442 1256 3022 2939 7176 3529 1602 1165 8637 3187 4331 8452 8339 5864 2612 4022 4194 2241 1551 4792 2176 1753 5296 5122 3746 5990 7636 6810 4334 920 7107 8680 721 516 4212 3094 2788 6629 6859 6046 1212 7625 262 822 2052 4955 6729 6934 5667 3007 8478 1885 4178 4760 1837 6390 1812 5064 7813 6691 5492 4567 4233 2353 1746 4798 5699 3626 7941 3070 5352 2370 421 2728 1979 5437 1223 3489 2450 6185 3235 1660 4251 861 3083 5311 2004 4568 7663 6476 8307 5351 6077 2228 4316 3776 7065 3996 1160 4021 8091 6154 8301 8559 2995 2867 6976 5245 5130 4465 1460 4063 6443 5301 6528 6290 3630 3300 6459 2436 5890 7796 3115 6504 636 2543 6632 2620 7908 4399 2840 918 746 1766 633 4156 2882 4166 5682 4581 134 2694 3593 3583 528 2479 230 7694 7066 7510 4592 5409 773 478 4153 5858 287 3298 117 4163 56 1236 758 4648 3941 3108 881 8629 1765 2207 3420 3876 5213 4140 84 1525 5160 4239 2388 1243 7178 3068 2572 1606 8532 4776 8750 754 2131 8203 4857 2243 842 5892 8030 3839 5440 4099 1828 6054 4534 181 1978 8741 167 8455 4378 1467 99 3711 2487 1923 8367 3943 7028 4220 1625 5680 1305 8015 240 4020 5597 6280 4852 8613 7160 1411 3552 7806 6720 8105 3358 5186 1390 5416 3635 3383 4143 5656 8681 1081 6529 2716 4292 7599 7114 1886 8644 7834 1964 6569 2054 1543 5431 6930 5180 5251 7657 3491 4307 2732 3586 8724 5743 5407 2132 6173 4940 7607 6146 2902 2768 2049 4535 1686 2366 2281 3225 1077 1137 8398 741 5403 7613 7091 6781 6361 3718 6806 7475 5259 7542 6223 3119 2603 5943 2925 7406 663 6408 4053 784 5368 5698 6530 6412 8359 1824 5745 685 3807 184 5718 5659 5417 1450 1502 6431 8661 6337 2984 6398 4831 7596 5143 659 5636 795 7684 6416 6379 7188 4783 7388 8257 3174 7339 6203 5179 7071 5154 3934 7962 7948 4734 2815 6218 1291 5511 392 2853 7207 3200 4104 1222 2404 7183 1313 789 4641 5856 7762 8281 561 1706 6169 8391 620 1586 5860 7562 3616 4259 7935 386 5003 6982 7173 6566 1548 557 3500 2764 8601 791 6657 7844 4272 4717 2447 2920 679 4747 6658 2325 5111 1172 6151 4476 8201 7231 5030 2026 2472 649 3632 5635 2280 7060 2817 7271 5843 3104 3436 6607 3294 2511 4971 5993 5773 2927 3241 2544 2960 7482 4131 3207 556 4968 4846 8240 6086 3782 1068 5095 1150 5099 400 1342 688 4172 8016 4631 1011 4105 4446 1351 7932 2763 5755 205 7773 2669 8366 7237 3496 5465 4511 3470 2143 566 7143 5930 3303 684 3617 3655 2040 4429 24 2657 2625 3860 2887 2293 1855 4384 843 1166 5123 5881 480 8021 1581 1535 1661 277 8323 4874 994 286 5847 8123 2399 3548 7980 838 1445 7827 1741 4637 1410 5615 495 469 2394 4880 412 5262 5542 1611 6602 4667 6576 4314 8330 2642 254 8603 5867 1181 854 3998 604 7748 3487 2238 3811 8156 7658 370 4180 4655 7744 7627 5664 6922 7868 6178 6266 8515 7548 4822 2318 4816 8256 2266 1078 2248 4398 4542 8052 6822 7573 730 3415 4206 886 8496 2328 4256 3965 2180 7767 434 2990 1009 3081 3289 7172 8703 1184 3193 5447 6332 7316 4864 3886 5520 2159 6448 6348 5072 5443 874 4506 325 443 2438 3631 3387 3829 7751 7692 5250 7702 7726 588 7059 212 1447 4609 3378 7093 3971 4917 3956 759 8474 7678 3482 6324 5148 1050 8392 1054 3757 8389 1319 1042 3824 4949 3116 1986 7303 5571 5187 6888 4967 3680 1601 6513 893 4698 6610 6021 8318 6299 1127 2838 2397 4530 5135 3641 6353 1186 3948 632 8311 7392 4963 1699 905 2339 7251 8587 7777 5011 2519 7246 1716 1375 892 590 3380 5158 5389 2806 2548 4799 8401 8454 2978 8005 5546 7735 2110 3263 4539 2618 7294 2126 2108 7689 5241 6932 2466 4964 4913 7727 6750 8472 5874 3419 3788 5037 2352 1352 977 762 7703 4672 2545 6608 3444 5795 137 7315 6208 2601 7371 2222 2807 6517 3670 1708 444 5706 6474 661 6913 1048 1448 2152 6966 5889 4586 1277 1615 3205 5728 3643 5226 3130 1405 2182 1113 3452 777 1610 1171 3168 7204 6106 6599 1295 5796 6224 3360 3079 2596 7041 3243 1003 1209 5954 3227 39 3431 8291 5061 1135 1426 4787 863 1019 4678 719 4309 2839 3967 7846 1126 6239 944 8551 2996 5444 1288 6122 7595 2931 6320 4388 5015 7550 7732 6924 5610 4558 1015 8471 7472 8068 735 1010 6422 8100 5600 52 6445 3363 605 7217 637 5921 8336 2101 1413 7296 8652 3004 3647 506 3230 7454 208 6580 4486 6325 6037 3602 3350 2851 6184 5598 5508 4425 3356 4407 4074 5969 4870 2453 477 6156 1070 7199 3779 5811 6275 5693 1871 1787 6407 8090 7150 6577 2863 2491 7487 7381 8412 7640 7743 2267 4750 8555 4380 6267 6500 7227 554 2009 2028 5507 7281 2340 8387 363 53 5420 3493 8451 7567 8157 5653 456 6462 5621 2416 1455 6247 2365 4449 6321 49 3687 1237 305 6070 6850 1312 2778 3765 382 3030 5785 8196 3897 2395 871 2137 1515 7737 6970 4269 5277 1061 5225 7511 1961 3084 7606 6831 8186 8523 2606 5543 7289 1858 8723 16 907 1461 2079 6206 6912 2608 104 5505 143 6111 2880 7982 8435 8514 3077 7975 2816 4753 1376 3260 2592 7340 6802 4930 3766 8385 3871 1520 2380 7100 586 4842 3073 6873 6177 2785 5317 4624 4797 8383 4018 4512 3468 3517 3242 2634 6747 1642 8593 8479 2756 3366 3218 4714 7367 3014 8752 3912 704 6945 2148 3939 389 950 7272 3146 7528 3463 7073 5032 7470 8348 610 8158 6651 5604 8783 390 1781 2321 2288 3966 2080 468 1400 3938 3345 248 6782 1861 5613 271 6000 8404 5392 535 8346 1826 5068 4093 6013 8725 3611 6994 505 2345 3142 45 1080 1838 2913 6161 6427 5462 7720 2987 1629 2691 5895 7274 7437 4972 8162 379 4356 7236 2142 5707 8574 63 1898 2720 7974 4814 2033 2762 4741 2367 6716 510 165 3659 1255 2523 3238 3283 4800 1499 6815 6989 3511 2499 3595 1761 763 4394 1079 8610 5566 5035 7681 1933 5614 8605 217 3453 8145 6063 877 7084 4077 5387 6795 6823 5518 5521 750 3137 1114 7177 7045 4770 546 888 4501 8408 3636 3620 4647 866 7784 4409 971 775 5239 1852 8026 6293 218 4520 2312 1419 5996 6592 2283 2356 4690 1105 2731 1468 565 3074 1074 3797 496 499 3900 8187 5658 2336 1863 86 2733 5826 8357 6350 3244 2962 2861 6972 4408 8302 5363 7109 8571 112 7983 7537 5829 1422 6191 3936 8153 1213 5287 939 6798 4389 3890 5800 852 8155 2937 1737 66 5724 1480 4666 3830 2704 3425 7046 7805 8361 4386 2359 6967 5964 1873 5371 1199 2727 4284 6933 5356 8351 8582 4126 8416 4830 8490 2250 3460 4023 1388 5183 8224 4231 3234 7870 5914 3794 4905 6819 1934 7023 7874 6391 4447 3767 3833 1672 1273 4505 5050 896 4405 5816 1360 2967 3854 7839 5850 3017 1901 4644 1179 3440 7668 1406 4343 5063 4289 2185 3006 855 3978 341 4095 8098 3533 716 6113 7697 1071 2678 981 7683 7205 8248 4266 1779 3282 5430 6296 1870 2869 8092 3740 91 5265 572 2121 8087 6149 8170 8108 3393 8499 727 3190 4713 8512 839 8218 1829 1572 6832 1380 1122 708 3747 2240 6558 5898 5622 8343 3915 1190 3589 4911 6109 8126 3347 7140 449 3178 3157 293 7555 5016 6957 3945 3578 4551 8154 1051 3841 1251 1920 6679 3713 6570 2820 5048 5298 6480 4906 5034 654 6302 8238 6210 5394 2116 564 4794 1849 6360 4733 8174 6319 8484 6751 7145 8627 8488 8660 6380 7042 4737 2767 7166 17 6757 5876 1132 5515 8470 4485 2757 3154 2739 5863 7739 1336 4727 3451 1609 5323 3598 5175 3048 5244 2645 7836 3615 5335 1205 935 97 8112 4793 7080 3210 4249 4980 7861 3183 5254 5539 7157 2287 4428 1129 7496 5059 5831 5692 7241 7452 3177 1429 3792 7730 1473 508 5019 6572 3706 4199 8591 1890 8696 2090 3591 5612 2426 3755 7259 6817 131 2744 2872 3188 8427 8653 8483 609 7162 3539 7544 8217 124 5341 2361 5758 5100 699 8564 6712 7588 3435 6848 6251 7181 5793 8663 3969 3134 123 320 8597 8134 7098 2060 2304 5907 2810 5537 6172 7821 1269 1817 6490 8042 5948 2688 6340 4910 4617 7170 7096 578 3273 602 1924 4404 2926 8589 2384 5675 5406 1622 7649 5107 6503 176 366 2308 4082 734 7417 8139 3921 1550 4490 6202 2859 6433 8576 7208 514 4494 7440 1645 4927 2559 1522 5110 8730 5257 7391 7582 7990 3935 5067 3879 60 33 6018 8640 130 1156 603 6508 215 6115 6537 646 8063 2391 4583 2752 7135 2524 319 4841 3066 1710 3025 4638 6231 2562 4033 5488 6273 5390 2334 5771 3475 8228 4265 22 4336 317 4585 2755 7687 878 6004 3132 4218 693 2970 5970 4244 8225 7515 8704 4346 2793 4313 7759 4196 576 6056 1650 7394 5694 767 6990 2610 5093 2577 3228 2301 7902 5267 2031 4651 592 2446 5014 4253 2753 7991 597 2344 433 6301 7106 432 3214 2761 3481 4032 8447 255 8553 6950 8242 668 946 6437 7783 2605 8495 5538 4061 1330 1637 1880 8277 933 534 7614 8368 4173 6902 3292 1821 2053 1131 2055 7136 4211 1887 4056 8373 3232 7998 7997 243 8425 5005 5883 7127 7691 8019 5151 7781 414 4357 974 897 4642 8402 2246 4369 5263 7950 3874 3118 1139 2455 7988 8720 1325 4261 7256 4650 4450 864 6801 5322 8671 5939 1123 3164 1101 906 5343 4116 7352 1267 2933 7704 991 7130 1939 5810 2676 4123 6436 2956 1972 1600 5463 8308 544 6278 7163 2553 6698 1813 5581 122 7155 7559 5347 2599 3834 196 4254 1153 6968 8599 7500 8585 8444 2002 5671 5219 7427 948 4696 5794 1083 48 1864 4855 6904 832 965 7978 15 1463 6917 1835 5112 6573 5172 5382 6428 562 8150 7701 978 6118 4756 4854 8237 2000 6499 3769 6993 6383 2747 5483 3962 807 61 8674 3801 8287 2044 5375 7750 5940 8097 3538 214 138 1872 1067 2905 8667 615 3975 3323 1662 8606 1528 3625 7232 7698 420 8670 1516 5947 2087 2566 7873 6367 4957 228 5482 4813 4894 2643 7075 5761 211 3011 7623 1102 7146 6147 8011 1612 4114 7597 4636 1420 2451 2609 739 2204 3652 613 8448 42 3763 125 1742 1673 3182 4438 1693 5525 7297 8377 2377 2623 6220 3385 6667 3712 6057 703 3818 4109 7709 3297 5328 4222 6033 889 7087 6916 1935 8189 2578 4012 6294 3863 5029 6979 656 8714 7981 7794 1549 4459 1604 5922 207 7311 7292 7038 4242 3248 5206 6175 7220 3554 2818 7039 5138 538 6813 7579 5873 4177 6025 6242 7414 919 465 7557 1386 4576 1982 7985 4188 364 8411 2917 3236 1090 8552 4130 6743 1801 3861 2801 7424 3449 6674 1682 4759 3924 1762 4975 959 6366 3036 2270 1347 2661 781 5306 5934 3669 4523 1950 7363 8536 7219 411 4167 8069 8659 5608 5825 4939 413 3117 8586 4103 7431 4271 8294 4757 560 8246 2082 2229 7504 8055 1585 2038 4102 1187 4010 2748 7672 6762 3565 6878 6630 7147 8476 7053 8437 8295 6318 8468 2976 1664 1486 7078 8436 6038 6840 8626 4135 1331 4297 3495 3490 8635 5408 1900 6252 7462 5611 6089 5381 2684 7638 2860 7402 1680 5489 2113 8409 7332 8219 3010 6737 8230 6160 1082 1553 1476 7036 1760 7845 6958 6971 6023 3981 8554 1208 1144 1170 2724 1668 1911 18 6209 8114 3067 8356 3105 6344 6315 4614 6469 1493 1932 472 1658 6114 81 3219 921 890 6333 1158 5139 6777 1348 2457 4871 3888 768 3745 8340 5528 6792 7643 7186 7947 7092 5802 5678 3528 8650 7746 3781 8159 4 5435 5592 1369 2208 820 3855 7904 4849 3046 2973 7330 4182 185 1994 5834 2015 1671 4046 4436 3009 6098 227 1002 7448 5983 7262 1644 6557 8206 8502 3785 5504 4610 5174 5333 7949 7506 7052 1430 6742 58 5647 676 7050 6741 5772 3171 3247 1485 860 4451 7342 1038 7878 2547 7919 2161 8009 381 6126 3980 2500 3031 7906 1031 1215 3730 5575 8609 6588 8306 4322 6615 4868 1648 8006 6110 1176 82 7574 2633 7666 5302 8649 1106 6642 2682 8232 4367 7215 5924 6328 4912 1534 7492 291 712 2751 2324 3173 8244 6808 8517 8579 2153 898 1355 1275 100 5261 2842 7655 8658 2907 3870 8768 2567 3255 7139 3139 8756 7804 3374 2792 3683 3307 3034 7238 3053 4418 6534 3704 7667 8556 7724 6664 6821 3176 2538 7964 7989 5886 7788 4154 2078 459 8081 1511 2286 450 3800 8675 159 7749 4096 3161 6331 8633 2951 5410 520 2018 5624 6921 4016 2307 4570 2076 7444 5784 8394 5076 8227 3280 7764 4633 7899 4031 4552 3012 7273 6655 1967 6282 3787 3410 5498 8687 8560 689 8252 7494 1777 335 6907 4823 6493 6758 3919 3357 6749 3951 2396 3546 7611 8760 4843 1362 4969 7458 6911 6759 938 6984 5806 4544 2448 7490 2430 3535 3364 233 6622 6641 5149 2051 4674 7375 5080 3725 1147 1919 6627 1000 2557 8176 6042 4563 8347 5683 8058 148 2702 7714 8003 8457 3982 6928 5240 3903 6481 1960 5545 2247 486 2653 7318 8285 5454 3128 8508 2348 2743 2736 1063 5607 4566 4365 3790 7524 6010 2011 7020 875 7944 7465 6796 8381 6981 7505 7254 2712 3314 6951 7433 5071 2224 3976 3405 3278 1663 7425 6134 3209 7326 7717 4639 4430 2593 6059 3095 8430 8202 311 4988 4353 484 8441 2422 5517 2198 1437 575 5920 8578 579 8397 3892 1224 281 4884 4754 1541 6467 5713 2765 6313 7710 5248 4080 6189 3329 8254 5038 2166 652 5284 7083 1808 6406 1250 7705 8519 5088 6964 488 2870 4138 3103 6241 2111 3266 4622 6525 2938 1773 5024 5630 7603 4531 7011 136 7365 2630 5865 8717 6694 2158 601 7099 1287 6123 8616 2421 4818 4362 3693 518 275 5677 5467 155 8710 5941 3650 6886 5401 4948 1332 1304 5499 5862 1442 3509 1109 7713 7393 2497 8743 3369 5307 653 2168 1006 161 5457 7206 2878 7169 4454 2136 1228 7076 3661 3361 879 5584 4473 166 7532 1452 4652 2230 1449 865 790 4379 1758 3590 8214 4456 3272 2919 8432 7531 135 794 3760 4840 7929 729 7044 5753 3125 611 3377 5992 6551 8611 6937 2949 8719 149 3286 3842 1559 2770 1151 5009 8755 2837 152 2530 8073 3251 8084 160 3752 7261 7000 7426 2329 481 3 4176 1857 4745 2784 4493 8679 3375 5450 6409 5820 1554 3045 4521 3610 3544 7022 8520 5210 6816 5936 1034 8449 6228 2786 3927 3336 3551 7074 2714 4306 6636 147 3930 7415 216 8706 5181 614 5344 3734 1814 1029 4232 6410 1715 5268 5027 813 3107 3954 1562 3437 817 7353 4504 5238 6346 7849 1044 4437 2540 1418 3476 3840 5877 3940 1943 8235 2439 829 4139 7195 2489 5733 7892 806 6693 7622 7180 574 3208 6547 4603 2952 553 7469 2986 8590 7968 2383 8365 5185 5984 2302 1538 5372 5650 6956 6128 3317 815 2802 3141 1357 2814 7264 1667 72 4092 1321 5270 2195 5017 3563 3644 1359 5247 3147 3885 5997 5557 156 8136 1227 8684 3501 5472 6051 4646 6637 1402 1944 4158 5726 5276 835 7583 6375 7590 8329 4236 2213 5116 2120 4555 439 14 1544 3762 2941 3889 606 3719 6345 2565 2729 3716 3099 6212 8518 2349 5972 70 5998 5096 2825 3812 2656 6139 885 3970 3096 1283 7894 6593 4134 4772 3614 4097 1738 4328 8636 3135 3574 5377 6724 153 4396 1415 4127 5448 5436 580 3562 8612 4685 6893 6874 812 7194 3258 269 6652 6096 6959 8195 4066 7097 5258 5790 8775 922 8642 3629 7865 7385 6137 3682 6727 2836 256 6991 2890 7072 5491 5809 7766 6314 1464 7249 1183 1593 1007 8364 8048 7575 8380 3955 5526 776 127 6703 6145 6058 173 4824 8271 4298 4054 8093 5221 2005 2197 7740 47 6510 7168 6953 7653 8712 8530 841 2886 2210 5681 2888 7212 5218 1417 1633 5985 4112 655 8185 5460 6889 4255 4011 4224 8247 5012 1059 949 6342 7706 2327 1384 3510 990 779 6533 492 2898 5925 667 6006 7651 4786 3213 3277 696 2401 1965 1877 7132 6281 8043 199 1475 4804 3051 3853 1866 5052 7526 1626 8125 8777 3072 3080 647 6362 1869 4553 3997 1909 1587 8429 1906 3789 5177 869 2284 7556 6415 5036 742 3076 6766 1374 6812 7123 5524 6883 7461 333 2813 7491 4455 4162 5479 4401 8099 3319 6162 3250 6841 1764 3541 2081 7572 749 3582 8465 3895 4260 6775 2431 4703 4898 2030 4300 7069 4152 6811 195 7329 7871 7715 8166 7723 2371 5326 7646 6472 5286 5398 416 3751 7733 3163 6559 7928 6715 4922 2776 1239 2517 1182 8199 7467 1782 57 2665 2073 3622 5140 5509 6235 6683 2408 6804 2249 3206 354 3584 8654 8708 3720 4257 2193 1142 4089 2862 7233 4106 5580 764 1557 4735 1204 384 440 5159 3852 4879 3109 6200 2097 682 8250 1062 7539 7191 2822 2649 7921 2169 2598 4941 2476 1569 930 3302 2092 1953 1207 3339 3851 345 183 5128 2417 3883 2835 5766 3127 4811 1210 5935 3822 8624 589 5632 6929 3905 6935 7167 6531 4669 747 1750 8562 501 1853 7344 7615 4995 6205 1914 7322 5201 313 630 6107 2654 5176 5931 2427 285 3609 5380 290 6895 7477 7117 8116 5846 811 8664 3950 7228 5541 8466 6072 301 761 2999 6221 5967 5805 7257 6856 778 2828 4208 6776 7569 7358 6470 5902 972 2881 3445 7934 625 3052 744 4492 8018 1513 4881 6818 1488 4593 5449 1736 2737 4726 5484 6505 6732 1580 6104 7451 3396 6015 8727 4419 2503 1260 6183 3160 536 4662 7520 1451 645 1913 7409 4240 5202 8749 475 3426 1268 4107 4928 6121 6919 1624 4124 111 1203 106 6591 4395 805 6682 1333 1145 5157 4543 8561 629 3929 6174 3069 5163 6259 2734 5549 4128 2415 3507 4774 1494 7346 5480 5256 5384 6983 1307 4197 7103 5134 5125 5884 6347 3390 5477 6008 4044 3087 6358 7979 3902 600 569 1723 1521 4965 6564 71 7359 5170 5354 6725 4653 6071 8713 3585 7395 2047 3111 8128 3432 2541 4435 6866 1104 1315 4782 8461 5763 4246 3804 692 1799 6606 792 2502 894 6279 7589 5646 1311 395 8296 1907 6170 5649 857 3229 6803 4470 5903 11 5212 2803 1729 4755 1191 146 970 6858 2539 5049 5108 2798 6192 8758 515 5441 8450 1289 6837 7718 1970 825 5400 7418 7786 4013 4676 5283 7808 7043 3877 4997 618 7533 6065 6306 1174 37 2355 4716 3743 2429 1905 5103 419 4464 1185 7082 8641 7268 6166 5561 796 4835 509 3488 3691 7791 4697 330 142 8608 7891 1640 4340 6387 94 4115 6552 2725 8563 5617 3441 5899 7645 4499 5468 4608 3112 5841 8529 1576 1220 6052 674 4649 4777 7253 3907 5321 6708 7247 7302 2379 28 1407 447 8316 8282 3352 2443 5252 5579 3409 2094 8622 1820 3592 8771 3368 4934 5968 7062 6186 8071 6327 7931 1769 5901 5946 3411 1807 8594 4620 8261 294 1408 5184 4482 1804 59 6828 3698 7034 3532 1373 787 4683 4406 5167 2709 2199 3850 8355 1774 2982 6357 1955 4890 2521 559 5503 6401 2824 5734 7832 3795 5273 1425 7543 8666 3987 5386 7240 6613 4141 3223 485 3191 1867 3013 1649 2212 258 6101 1641 7319 4720 4117 7644 3473 2536 5383 6713 8565 6876 8424 6894 5869 6440 4858 1789 8074 5230 6596 7270 5278 3262 6554 5079 2563 4466 4802 3803 7820 4832 5560 6833 7610 4643 6369 7628 1385 5194 5601 2514 3627 209 563 4702 2292 1810 4064 740 3914 4875 4943 5648 3553 8267 2037 1377 6544 3953 4040 5926 709 3456 7965 1409 3813 1482 1456 5842 6270 7245 7129 5136 2036 5845 4252 3828 303 105 1922 511 2306 1192 3274 5145 4280 3114 5818 3798 7474 8191 3471 3330 1730 352 2462 5788 2930 7911 2184 8673 5303 2555 8211 8757 3974 2070 7761 8509 7509 7977 2102 5817 4108 6432 8386 8415 5778 6165 2492 7434 4807 4215 5062 4645 1219 7677 2375 2961 7343 1498 3684 1087 1537 5236 8628 3334 8148 3986 4157 7442 7999 2209 7823 2954 5320 551 1428 5803 3408 3256 373 4628 1483 7390 4998 2319 8699 8773 2812 3654 4836 8507 5069 4015 4839 5654 6310 234 8700 4281 4275 5715 7937 2188 5603 570 6014 337 997 3961 5374 8434 6656 7671 537 1454 3040 8516 7225 7008 6864 5777 102 3846 7519 954 7897 3167 7032 1546 6473 2958 3881 2177 8322 3252 4288 1674 7013 7335 2908 6027 6779 5106 7327 8008 399 5124 4670 3666 7151 6985 7193 1300 4361 4457 2469 4679 7210 1338 5551 819 1570 7416 2483 4179 2414 5609 76 3253 8567 4795 2172 8575 8255 6386 1937 8728 1510 1341 3249 7525 993 1317 8584 438 2823 8110 4366 8686 3703 7633 529 3697 6016 2885 1334 8776 6885 4479 6914 3738 8172 4026 2174 4440 8572 1436 8731 2478 5075 5643 4931 3537 1639 3287 5434 7594 7153 2042 5471 8017 1725 8721 8779 4286 3320 5227 3029 4892 2381 4549 7209 2256 9 6442 5937 3397 2695 2909 617 1775 4838 3457 6790 1306 1140 2227 1865 8738 1173 8292 2086 7973 4914 4904 3564 4764 5663 7068 5022 7012 7336 6941 780 8140 7516 6002 4262 2460 5875 5945 5523 6119 7357 8406 6155 4452 6304 7460 2985 6260 1579 2706 2718 504 1678 343 7277 6195 3958 8423 296 6939 3649 326 2437 793 1285 7660 4000 2611 4086 6542 8259 6041 594 928 7026 5919 404 3389 21 8497 8620 8305 8662 4122 4350 4758 3148 2023 5373 6489 2506 3960 7126 422 2277 7859 3891 3097 2445 1 5429 4513 7187 2940 2465 7497 6905 4241 6579 8012 7341 8056 2335 830 6194 4664 967 1583 1271 4472 7604 5723 4149 2966 986 5872 915 7223 4002 8310 4865 1149 103 2105 7061 1654 8733 2534 2338 567 3753 1242 5744 3821 2551 6336 4478 8161 4039 2424 1607 4047 7115 6405 1398 7970 4332 1839 6159 6190 7328 5001 5458 4815 3773 1743 4048 3427 4374 5910 278 5739 612 4708 1121 7441 873 5792 5209 5703 3664 3172 1878 5191 1980 8326 7355 3808 2681 3044 6322 2671 4740 7213 1757 7522 6181 1588 4590 6666 4264 4775 327 6516 7755 7089 1975 3403 8062 4956 5558 5908 6653 3648 5531 5563 245 4829 911 2896 1298 4739 1973 4463 2449 1772 1244 7411 435 3305 7875 956 7924 3686 4951 4556 8625 7782 1495 4918 7312 6256 8469 4781 6600 5021 2454 2274 1200 4308 6253 1524 848 5353 2660 3156 8370 7768 6003 5897 8475 2508 7174 7200 441 7665 3434 5731 2003 6881 6587 3038 5281 2214 8639 7561 1700 5233 4315 8107 5552 4632 3835 4982 8481 4722 6555 6649 6225 3581 2041 5205 2221 2906 5854 512 521 8141 4502 7662 1655 4616 2583 2275 3524 4498 4228 4986 8178 1574 4085 4076 6700 3222 5536 6050 8067 5097 7635 5735 782 7421 8521 1756 3657 1721 3454 8243 2061 5415 8543 2147 6639 5953 2088 6396 436 3267 5290 4217 4337 2708 1556 4169 5274 6575 3658 5349 7447 3534 1308 2772 4723 7438 6207 5787 2400 8312 1146 1987 7895 4872 6619 4320 5197 951 7369 41 2865 1233 437 6152 1936 8190 452 1603 391 1290 2855 6778 1097 3543 340 7121 7581 3991 8290 2673 2251 1501 5362 3558 4414 3878 1739 1833 4629 1566 1214 593 4923 5736 6055 2459 736 4925 3295 1770 3246 3016 7669 2095 953 1889 5705 3732 3082 4193 4347 2787 3005 2997 3392 4686 3651 6060 192 5652 2192 4071 962 966 683 5668 3318 711 6080 2481 4546 69 2233 1912 2244 2942 901 6925 908 3143 1630 5234 5857 2866 6863 7578 5690 4488 5023 1438 3783 927 5242 3413 519 4519 912 3696 6197 7141 1657 5586 7374 324 3102 1676 6787 816 1951 489 1533 5962 5497 7676 1457 383 8774 5909 2742 1925 7035 5131 2883 3217 6565 6901 473 7972 1189 7443 4004 1297 5933 695 2943 193 591 5506 3338 5006 4577 872 7275 6920 7601 1466 1976 6062 1354 3433 1124 62 4527 1748 426 3571 7661 917 8705 7619 1882 8013 3572 752 517 8249 800 672 225 5870 6359 7027 6182 4230 6744 7618 3310 2945 4710 3749 5091 6728 8047 4769 3197 4903 5527 2834 3275 2585 6854 6701 461 2685 109 6464 2944 4243 4611 2809 7370 1014 4019 2808 5229 4936 4640 3478 328 6140 8362 8615 1840 8179 3333 8198 3165 8375 2412 6351 4604 2317 7652 7840 351 4889 3063 3959 1248 8337 2658 2721 5807 7552 999 6706 2124 6838 840 2119 1471 7632 2680 4210 2775 1591 6633 3525 6158 3918 7216 5737 4926 5688 5132 5196 4295 5894 6604 3313 4003 5742 1701 585 1358 314 717 8309 769 6226 4809 3092 6196 7848 3486 631 5725 6439 1272 7570 297 7621 5237 3999 6786 7756 220 2007 1202 1432 5361 4659 4027 6908 3523 479 1651 1282 2444 1884 2223 6581 1474 402 6739 194 8600 5418 453 5478 4812 7828 4014 8023 7534 5031 2556 2564 5651 6309 2117 1930 2261 5606 7841 3226 5975 4888 5117 6977 6179 4121 691 1697 6947 7324 1252 2225 6736 6900 1506 7309 3002 3422 5264 1112 2062 4119 4763 5200 3864 2854 3774 6450 2797 4371 7858 6678 4507 7656 3709 141 587 430 1632 8676 5721 8279 4305 2129 6549 425 8698 7943 7814 1792 6068 4381 5000 1041 4767 5752 3923 2420 7384 4680 8119 7700 743 5981 2690 2343 5544 2893 2100 5203 2983 6890 2268 8766 2549 2527 2923 2008 3494 7984 304 6668 6418 6509 8759 4051 3359 6520 3875 8215 984 5764 154 2581 2231 2899 6389 1618 3676 6138 3662 3758 7056 2918 4500 3990 644 4327 5595 6585 3836 5388 1981 7670 1107 6257 7102 6906 4489 3513 7149 2759 4730 2989 1695 4790 2012 6955 3047 5469 7910 4084 670 2085 7527 4148 6323 876 7250 4533 8020 4958 8280 1984 2146 5288 5587 7654 8352 7966 1595 5232 4762 5404 5929 7412 8121 4321 6465 2130 3605 6129 5906 4411 493 7125 7529 1992 961 982 2722 4748 8688 7540 6685 5661 3414 380 1500 6834 4987 1631 3645 46 5085 1012 8678 7077 2112 2857 4110 2558 7624 7818 7760 267 4445 2032 6996 1477 3331 5487 3791 4538 3335 1489 8165 3515 8643 8115 3768 8265 1188 396 2382 2646 6371 4036 8482 1453 368 6373 8718 677 2190 7422 5228 6839 6568 2777 5338 8583 3928 1164 7298 8338 4146 6623 3296 6399 8746 1647 3071 3866 6404 3867 3603 5092 3049 7101 2289 4947 6556 7267 3739 5813 2470 4151 6903 4329 1138 7301 2201 1065 3873 6532 525 3638 3381 2259 1045 7197 362 584 1805 7112 513 2994 8677 95 5325 6805 7331 140 1382 2296 1320 1459 7763 1881 385 705 8463 802 3985 1005 8031 2640 7014 5644 4282 6814 2435 4202 627 231 8748 5445 1093 5056 6085 1052 1230 6849 6403 6946 8050 2873 4189 4682 1327 6308 5455 7498 5142 2912 1412 8288 8342 3264 7741 7133 6082 4946 5989 706 3308 1263 7853 2440 8258 251 5913 1815 6271 884 8458 8669 2532 3569 1024 8573 3671 6131 4344 2387 1656 5982 707 1111 6261 7847 3677 6643 2346 6487 5848 6141 442 3472 2425 4293 4712 6291 3681 4185 3727 5514 7830 8117 3020 4810 2626 7850 7320 6730 3868 7707 5628 8001 338 5495 4594 7175 3573 3261 6494 8036 7138 7585 6987 7922 5915 4989 5355 1033 7862 8480 1555 8143 5304 2677 8440 8632 1991 1691 77 8 3786 8354 4325 5513 8142 7158 2098 2717 5776 1266 666 6163 8004 44 6456 3799 7546 1387 8537 7776 6164 7285 6395 1596 7901 2879 4205 6918 6352 5626 2285 5973 7105 6659 6538 1985 6312 8503 1168 809 1558 1759 2783 5275 2474 7104 4876 3061 2254 7351 4136 5337 3391 5564 1117 5791 4623 2884 3922 3894 6330 913 6861 1096 4481 5741 6740 7501 4587 673 1512 2950 7754 8550 7771 7203 7015 2013 2904 2621 2856 7439 6011 35 2936 2819 2010 5446 3695 6755 6784 4175 6844 2741 3911 4929 731 8604 4715 7551 5634 8439 1025 8220 3916 5717 1698 7024 7142 6868 7837 3621 3315 2456 3618 1809 7258 4778 170 7058 423 2666 2683 239 1465 4656 5282 6820 1253 8557 5808 3019 8169 8369 2791 5774 6234 2796 2392 5573 6029 1825 2072 6705 4828 8183 3701 5193 7807 4426 4273 2350 6829 2071 6180 3979 4701 6453 5740 1845 5309 4203 828 5058 1399 6274 4263 2594 5399 6754 5164 3324 6012 2637 8349 3151 6199 410 3504 7939 4475 2675 5004 4310 8135 7779 963 7857 6237 8421 3065 8320 690 3367 5313 808 582 4442 2826 648 312 1962 8769 3899 568 2205 964 2458 5166 5094 4974 4893 1518 8129 6047 5078 1592 4859 5662 8245 4706 5554 3884 7952 375 1389 1545 577 6227 67 2710 8734 3458 2364 2877 6609 3133 454 8184 4575 1383 3994 6485 4618 409 2262 2696 5070 1749 6434 7605 4621 1118 308 2525 7753 6400 1504 973 4060 1589 289 532 900 3155 7019 318 5961 3042 8607 358 5840 7789 2407 6562 3354 3973 6392 6099 8054 378 2526 3276 7685 6616 2979 8494 651 6845 5995 756 7288 715 5026 6429 6696 177 1724 5113 5220 1926 5168 1908 8414 4579 2386 2034 5885 5081 8236 6826 6998 7907 272 7626 3179 774 5439 2582 4826 5674 7946 1841 3778 3090 8527 5535 5852 8193 8486 88 6244 5308 7002 1790 1568 8127 952 5748 698 6825 7580 5977 5222 5428 583 3372 542 1361 1157 8621 2504 3869 6144 5669 1324 3536 2965 7508 6770 8399 6414 7403 4354 2981 3268 1198 7338 6061 3288 7833 4417 150 4113 6093 4873 6960 6617 4375 5582 1370 6940 8209 3399 2127 7266 4415 5047 3466 4301 3337 5474 1998 3056 8711 6097 3018 8693 3136 5294 1959 8538 5476 3465 1620 315 4895 6447 7881 6962 8558 8545 3623 7484 145 8446 8094 5974 2635 5066 6413 3674 5419 2103 573 6032 3545 3386 834 2934 924 1443 1226 6120 6718 8044 7523 1677 2075 3600 2955 6635 6076 2418 54 50 6491 4160 1234 7286 6171 1803 6603 2668 6646 5779 8631 2655 8171 1163 7047 4422 8443 1945 4392 3362 6512 4887 4994 2477 2574 168 797 7883 783 4159 4518 3384 133 634 7189 1066 5655 3844 5866 1363 2935 3003 725 6634 7377 7680 2385 8000 7196 4517 2580 5316 6791 5369 2452 8539 4729 4529 6760 8083 3175 6382 883 5329 4453 8223 3462 6584 5147 1278 6567 6108 8078 3608 1057 5837 1718 6597 4979 6974 766 6454 8390 3203 6553 5350 5053 1904 5720 3526 2507 158 4382 6495 862 8066 8722 3343 5722 427 5978 6548 3530 3110 5215 90 7778 6045 926 1526 3043 3254 8194 5702 5395 882 347 7222 4732 7674 5666 1921 1284 7349 5751 7086 3865 2059 8577 3715 8709 6871 4377 8782 786 2771 8210 8051 7772 283 5781 2001 1216 1594 1646 5727 5797 8029 3480 5769 3910 8747 8459 7473 4630 2140 350 4780 3000 7476 4731 4142 188 1018 3668 7307 1076 6847 5965 5297 3796 1421 8300 7064 4326 4319 2186 2378 6482 7201 6201 1836 1717 6612 7879 2495 6711 3216 1310 4886 7252 5260 1381 2957 417 4742 157 4779 2821 8102 4953 6654 8770 5757 8420 8510 1434 2305 4689 2845 3729 7110 2518 5804 2498 6496 4276 4352 6595 260 2617 1634 1619 738 451 5007 2069 424 5783 4992 8147 3058 7882 5891 942 8438 5028 5971 5018 115 7798 4959 3428 5510 5568 4820 7005 608 6135 4491 7688 2535 8568 6807 4391 2597 1073 4805 4606 6132 2639 2975 1575 7536 1617 6168 3181 5955 4976 1229 785 5676 2795 4907 1143 642 6193 445 1931 2726 4673 3301 2781 5590 1806 3037 274 2358 1875 4915 3091 7560 8785 671 2876 4038 2864 2093 2290 4342 5625 2226 6571 3908 2219 4984 1831 6511 470 6746 3032 5516 6699 2170 120 3309 3024 7325 530 8188 5593 7398 8037 6078 7499 8491 4434 2313 4526 7745 4285 6704 8592 6541 3346 6869 8133 681 6973 2510 669 7641 3656 4341 2064 3412 6611 4168 1462 8395 7436 702 6381 6923 5760 1503 660 5359 3604 4245 5578 5714 4001 3257 252 6216 7900 5357 5988 8173 6424 7884 3814 8315 6483 5054 7067 2473 2354 5918 6269 8033 5553 1696 1681 5980 710 6081 2652 3984 2604 7956 3131 6460 1055 6066 7255 4834 6384 7507 4547 7362 3149 1578 116 3810 1134 2980 3332 932 4441 5192 7293 6709 8122 3898 8273 4724 3540 121 235 5562 2505 4942 1896 4684 1795 6219 7903 4216 3054 540 3550 6586 6246 7428 4692 6157 8740 7211 3416 7481 8028 765 8360 3376 3777 6446 1613 8544 5887 2156 3492 7122 1675 4035 5077 5 6233 3742 7566 2590 3816 732 7503 4087 3394 6650 3015 5716 4226 8780 2300 723 846 6618 7468 8299 6295 7202 2303 4471 4768 6515 1397 823 5065 1859 321 5833 2091 6036 1627 1565 4668 4460 8269 3527 238 6673 6560 1949 7693 8501 3348 5178 1026 8422 4145 1364 7116 3521 2114 5836 6640 7488 5339 6710 3057 224 7889 5129 128 8726 2099 5223 1232 3506 4402 3519 2490 2644 5559 4891 3202 5490 2390 1401 5332 3559 12 988 257 4477 4094 8333 5830 5911 2276 201 5770 4524 8278 4125 3106 3849 1509 8638 2692 5638 6944 1940 4991 253 7495 3607 3124 6799 8522 549 2513 4981 5956 1177 975 4860 6789 4536 3293 6092 6738 232 87 8270 1734 6069 6073 3667 1818 992 5040 3023 3748 3858 4550 641 4515 3240 4420 448 3351 2264 2332 1893 2362 1350 6697 6745 3995 5928 2852 5314 541 1435 6769 7094 2628 3726 3993 1989 8403 8413 7148 4238 4200 36 3404 1206 5045 1314 694 6486 26 5486 6034 1731 7866 4819 73 3327 845 6910 2697 8513 1745 3304 7152 2194 2616 6007 6167 5567 7886 657 7108 7480 1020 3395 2650 1440 3520 4625 2738 1755 1017 1547 2316 5422 4072 7379 6954 7564 7291 7134 7308 1161 6702 6870 8080 6519 969 4439 5089 4147 3166 8767 4403 4324 5426 4853 6215 7423 3721 5815 3568 8707 6827 6451 1842 162 428 3059 5759 4333 5695 2482 2560 2804 139 6721 3421 2699 2929 3736 1903 2315 4258 1340 2740 1791 5432 4738 6767 745 4869 5550 2779 8534 2239 5423 2014 7731 2974 5951 4671 1983 3913 1086 2964 5295 4279 3498 2006 1496 3078 4213 7483 7918 3233 3271 6867 292 6550 3442 2773 3588 3159 307 6334 7333 344 8645 1891 2279 8168 2488 1095 2531 7007 3917 8325 4079 5696 367 1043 628 7699 126 7774 8691 6298 6661 5216 6638 7446 5060 6335 1148 4545 5345 8374 7396 5041 6526 4977 5839 4323 1030 6660 7018 7926 1049 1874 8371 500 1599 4432 3561 8077 6507 6461 6268 6645 34 7838 5700 1563 8685 5300 4803 8010 2178 2486 2322 1598 7299 6692 3448 2471 6355 83 5672 394 8716 4541 868 1241 1690 1064 2991 6793 7938 200 4277 3371 4416 5051 7664 6243 2903 1514 6444 3761 4861 5812 2389 937 6001 1292 833 476 222 1927 539 4681 1827 2368 3685 1091 8753 1197 1623 491 1963 1958 4181 4516 1116 5548 4765 2701 8701 5709 264 1261 8200 1665 1527 4005 51 6824 93 7471 1788 25 5340 1956 916 6341 1735 720 7800 7119 6394 3265 7885 3153 7159 910 1299 1115 361 4601 2187 7711 2687 6975 4312 6103 3722 909 4278 4850 7914 8266 467 1993 1705 3640 7224 801 2829 1099 2972 5585 4195 5747 5324 5533 7184 7629 4045 4709 3438 2515 387 2969 1344 6681 5044 4484 4330 722 3126 3784 5010 4091 837 2849 598 8167 8313 7305 1286 4744 407 6857 1196 1848 5979 5687 5905 3847 7864 3100 8327 6124 8511 4250 2428 3793 7290 191 5459 8656 6245 2183 4916 7812 2257 1614 6772 1270 3328 8384 7835 79 6020 487 1060 5376 4223 2218 7260 7916 3660 8335 5679 3245 4170 2622 4274 2022 4743 7244 5711 8525 7795 2723 550 2679 6892 306 6363 6598 6843 5334 1403 931 458 132 5086 5950 2050 6232 8075 6419 1847 4877 7337 457 2832 7453 1036 5986 3400 6877 7712 7854 5775 7634 4368 8222 6999 4186 4663 5421 40 3221 2372 1711 1776 408 4990 7190 5414 1218 3138 1004 2278 929 3122 2461 2273 3672 5249 2782 8462 5665 349 202 3968 5730 3483 189 6022 7682 5121 5673 3088 8668 3039 1027 2155 1316 5424 2799 1470 6731 960 2235 4588 7725 4844 686 1221 7287 2584 2196 6978 5020 5701 4458 3992 4833 7792 5602 1727 7620 2528 1490 2686 1478 8204 6283 7925 3098 3497 7430 7945 1897 1201 985 4589 2831 4129 355 6680 7994 1892 753 8060 2780 2963 626 6875 5952 6087 4578 85 7955 7192 8079 7070 8737 3597 6229 4600 5101 6665 4784 859 1876 4132 7234 6539 8358 8096 3469 6488 3906 4582 4383 3909 2509 6039 2413 6372 3322 3904 5645 639 8007 4885 1794 4514 3837 7364 2730 3180 4921 6927 1309 4537 377 7867 3406 2107 5126 3817 8072 8353 7728 1119 7695 279 2084 6441 1423 8595 1846 5616 4508 2423 650 5591 3270 3062 5042 7493 2220 870 1638 7936 6377 7400 2341 2789 4687 8379 19 2576 1822 8263 2068 1740 4100 4355 8089 2138 1259 107 665 247 8541 6 5279 7368 3224 7399 6931 6733 4283 250 7242 8473 5293 1302 3741 376 6117 2607 6625 5494 1108 3832 7815 310 7265 8376 5596 2971 6287 2871 8032 4235 2202 1561 7686 3587 7959 6292 3831 8160 4937 64 4607 4335 6626 4572 7961 374 8138 8132 6376 4098 4030 1635 4996 2337 5576 353 8647 2237 1337 369 263 6965 1800 8745 4751 3439 3838 7869 6133 5470 8382 7350 5708 7565 7602 1684 3220 1393 4626 4219 3989 7765 5599 6774 418 6780 7445 471 2850 6942 6986 8547 5485 5750 8619 1296 6430 5855 1653 4561 29 3859 1180 3570 8207 6303 1154 2067 8038 804 7051 6349 498 101 5461 242 2163 5851 8683 5266 8045 1747 249 4897 2154 3692 2552 3848 1751 8213 96 3705 8735 3819 7040 2529 1323 8057 4718 4961 6484 2320 2537 6590 3474 8331 298 3212 4705 5500 1519 4299 1056 1666 3596 5589 2179 5102 5087 365 2916 4267 2160 7953 2157 2139 8182 6899 2494 2236 7822 113 3699 847 5133 6397 2846 4302 8570 8467 803 3370 2928 622 2074 7617 2398 4410 5358 858 7037 3809 3988 2619 700 2915 4487 7137 3503 2800 4304 7124 824 3050 236 5627 903 2331 7986 7757 4540 6605 6717 827 8505 5127 3735 3085 4034 5732 4661 1058 265 7785 2627 8233 2265 6249 1918 6385 5204 2299 6475 7719 2433 8500 2058 6563 581 8533 4294 6289 8205 8648 5165 2468 7376 6862 887 2144 8146 6024 4448 3925 798 276 5481 6285 7479 726 2311 1974 4704 5195 5057 8053 3189 2326 4229 1103 3321 6498 8506 6830 5765 526 7690 5556 6420 118 5689 6948 4935 464 6521 7775 7419 2968 1995 8276 3932 1235 2255 3920 244 6672 3770 7819 1850 941 1888 4118 2910 4554 7348 507 2323 3896 755 4749 1069 6028 7156 2260 1343 7095 8689 1392 8665 6992 8118 7898 8433 20 1643 4495 5555 8035 2713 8109 5451 8341 1022 7829 3064 1367 1687 1670 6582 5927 7382 455 1303 1714 552 5618 8732 945 1726 7896 1902 664 2651 1692 1084 4571 5379 5008 7888 6417 6915 728 3467 1013 2749 6492 1469 4848 6855 5039 2025 1072 3599 1567 3120 7321 6601 5217 5246 6527 5214 6426 2045 6589 4443 6898 1915 360 6354 2016 2667 5464 5073 4057 4675 7787 6477 8692 523 2988 8350 6753 4204 8697 4052 2561 68 923 5391 4412 5547 1240 2171 5119 3401 3512 7996 8781 6797 675 1249 1152 5532 2165 213 4468 4237 8253 3185 1329 8264 2602 6764 2133 2282 2573 6088 2571 1531 8086 2216 6735 6896 5305 8319 3075 6040 4359 995 7280 3690 4973 1542 2554 2134 5156 5754 5105 522 2419 7383 5762 2707 6545 2096 6949 5502 6624 8464 1968 8485 397 4694 6300 2020 4190 1728 545 8456 6670 1946 5397 5987 282 4062 7790 4338 98 8192 4387 6860 7304 5188 718 2579 6250 7630 8061 3606 4654 6378 5438 4837 7722 38 6882 2035 5540 8144 6284 2546 8212 3893 4719 8426 1194 8623 331 6031 4385 1098 8103 261 6926 4349 7435 5271 1293 2632 8275 957 6449 4214 1754 1702 5289 8784 8085 2754 5224 7033 8634 5365 3041 3944 246 3344 4821 6035 6628 6365 3407 114 751 2272 7863 2662 5084 2089 2700 1844 6049 210 8418 4596 182 8324 151 6574 7696 388 6561 4373 5055 596 4431 3355 3756 1245 1492 5025 8614 6305 6094 4612 401 6535 7081 1879 3949 6726 1948 724 7055 6262 3353 7825 4469 8262 1279 4303 6687 6017 7457 446 1851 3857 7987 571 2434 7586 1771 1577 2150 3724 3566 8286 1128 1719 4006 7513 5572 5569 2769 4137 8388 3430 31 3926 1195 713 6100 3957 1047 8229 2145 8082 302 6053 3613 6689 6411 2019 7917 8284 3145 5425 4938 1136 2106 6116 8542 8762 6339 4983 1326 4564 7085 3823 1032 7553 6090 5712 1552 3505 8690 169 4599 818 2948 7360 2568 7016 5882 7817 4595 5849 497 8251 1169 3862 5141 8283 4059 2206 7679 6393 1339 8106 3825 2827 1584 3388 3805 4658 7 119 5637 3710 3101 4944 8022 372 6671 5684 4433 3211 3237 5719 7541 7992 7356 5291 3577 2900 1621 1819 6836 4363 5002 483 3417 7486 1458 4900 7387 3646 6969 1694 3033 4902 7334 7969 1120 5932 5255 1035 6748 7378 7708 1928 8598 7354 7405 2569 8177 3737 4901 6130 3764 5819 2175 2309 4400 4503 6211 7587 697 2191 979 5629 7239 8321 2711 2570 6471 3717 3195 5453 7521 6714 2485 74 7459 5768 7729 3196 2647 3887 4496 4187 4009 8694 2242 1942 3461 3516 4635 6388 1484 7592 5583 78 3754 348 32 3204 237 1092 6497 7887 2588 4390 1175 5378 6501 5963 1075 8445 4065 4999 1564 2897 4101 1608 5633 1685 5150 8372 6176 2291 4068 4073 2253 8239 5312 2924 5949 2467 4557 1193 1508 3089 2889 280 5991 2664 7927 1997 6842 7831 836 2441 7221 1683 329 5835 8049 4788 5827 7408 8303 1733 2746 8175 6338 7826 1704 3663 1811 1133 357 6125 3239 1262 5366 983 5959 547 7279 8460 4584 4043 295 7171 6043 5146 891 4467 1028 7463 5299 4268 7843 4602 7372 6240 6887 5976 80 300 3477 5923 3827 2480 1155 7616 5821 8546 7612 8407 1039 6478 2314 5413 635 4677 6005 6768 2298 4042 1345 3514 1988 7851 4924 1335 1780 3194 5074 1159 4424 3459 2066 2805 7243 5473 3560 2109 1605 1280 5729 5310 2493 4049 7810 4024 4806 2595 7432 3113 4952 6707 2830 4574 223 1100 8498 8778 4954 5162 7003 5917 5619 4510 6734 624 7571 4413 3567 5828 3093 6255 6468 7029 4700 2125 2442 1778 8417 4120 6084 7306 2151 7466 3279 3901 2575 1427 43 3826 8772 7995 8104 3285 259 1094 6019 4736 7366 6370 8101 2659 4634 2351 2203 6676 7413 2245 6095 5642 8715 7048 4198 5152 6578 219 7855 658 2496 3306 5879 4397 7930 3170 2135 1497 7198 7235 3579 1616 5522 770 3424 4766 1356 2614 1785 1856 1446 2189 7738 3942 198 5685 3382 7860 7361 4882 4707 8540 6102 5594 110 6879 180 1768 1395 5211 6127 1444 174 5171 8729 3150 6686 2894 2234 1689 108 8764 5144 309 7940 1823 6258 7951 7120 3423 1523 895 8137 5822 1366 6880 1372 8566 998 6105 1396 6988 7752 5198 899 8549 5942 6723 7872 8124 8268 1709 4225 2705 7165 4345 968 5280 4919 7647 7659 5501 4627 5780 7576 6452 13 4360 2043 799 5960 6083 2636 3026 2629 662 7218 524 5385 8208 2021 733 3576 1767 2402 2200 2048 955 1868 3086 7976 5173 1597 2374 1217 2591 7456 1053 3714 4771 3144 7401 332 8274 7793 4339 229 7568 7118 336 7824 2874 7609 1669 5697 5412 5475 1636 4150 4808 7313 8596 6136 7890 2393 5370 5235 6307 1125 2703 3522 4497 6518 3316 1966 4351 3464 6364 1707 7031 7637 821 1378 7128 4041 1530 1802 6030 4227 4164 1722 5115 4559 7397 6794 6479 406 2624 8163 616 1353 4421 4317 6997 2774 6079 5916 3642 1040 2533 7558 5346 2947 8602 2162 2719 7410 8216 3733 8152 8180 5315 4945 2638 1211 3628 4560 7958 8064 2263 7915 4613 904 2039 3198 3479 1910 6852 8241 3259 3484 7088 5691 6276 8059 3290 6222 5832 2411 6938 6265 5530 4862 2405 4070 3192 2410 5859 6187 5118 8493 6316 7514 1834 8024 607 3700 5114 7295 940 3933 3542 6980 7380 3963 1679 1990 7310 8682 2891 2641 6288 958 4376 7502 7598 6872 595 5367 429 3379 6423 8765 7912 4856 3845 7880 5182 7049 3806 3744 3186 3455 3619 6620 5189 1414 204 7769 6722 6952 5896 1832 3694 8588 6091 3201 7386 7161 1088 1130 3284 3759 1732 6719 989 7591 371 5782 3269 6614 5577 5393 2745 6851 226 5253 474 7144 4532 1797 2027 7345 7284 6198 1008 548 4155 6695 2715 393 6853 6761 6936 2766 1628 7563 7530 6148 7420 1712 8065 1046 788 1894 7001 826 8151 8630 8027 5958 144 4985 4111 5090 7577 5868 6143 867 8070 6075 5043 5749 3952 2848 5155 6150 4569 8298 2953 8314 3601 5786 3450 1391 1917 4370 3856 7963 221 6865 7954 4081 8396 1365 3340 543 175 4364 4878 7802 6752 6502 1328 7650 8293 1481 5208 6044 4075 6343 849 2911 3035 3937 2046 5427 4050 1433 5534 8076 1843 3291 7429 8317 129 1167 5994 3231 1424 3633 5493 6254 1952 3675 1021 5838 5336 5153 6297 1379 8272 4201 1441 342 3162 7278 6891 4037 6074 6514 1085 6458 7113 5512 1264 6236 1957 8695 4908 7025 5199 3299 8297 6669 6067 7248 3557 4598 680 3547 7314 7079 2115 8405 1798 4950 5396 6248 7276 6217 1571 5456 3880 1971 4444 5823 4318 2310 902 4427 3771 2269 8657 4083 7913 7876 6230 5853 6311 2512 1472 6326 3502 4591 1786 527 2 2360 2164 4688 2330 987 4171 5844 1783 8489 2017 7030 4933 3485 5013 1529 7816 179 7639 5207 7230 7164 5272 3152 7347 4801 6524 2181 5348 4932 7111 6009 6663 4548 7450 5801 8535 5893 3169 8736 463 1969 1416 3688 8477 714 3199 2689 7716 4562 6356 4660 6765 2173 4597 4899 2698 6329 4711 3060 5190 398 2586 2993 4761 5402 6800 4970 6455 3815 6457 6963 5799 4483 1231 2464 4372 299 5756 7512 814 856 687 6523 2363 1582 1037 6961 2215 2057 7404 7675 8744 5360 7547 623 5605 6213 10 6785 30 4615 4691 5798 4183 7852 6264 3780 403 8761 2875 5231 2056 1793 3121 1247 3977 6317 8739 6522 5285 1763 7811 6263 5104 7942 5109 1349 1505 405 2977 4025 8334 5620 4665 2357 1258 3639 3708 5318 7263 3624 4721 5789 2403 4191 3418 4247 2122 4209 5824 980 1540 1854 2484 6214 2735 1830 164 2847 3028 163 1744 241 5912 1274 7734 7009 6064 4960 1141 3882 3055 5738 4393 4067 5641 1089 3702 3594 6425 6272 7407 6688 3123 186 4525 4221 7006 4462 2914 460 2760 810 2342 3946 8040 8702 3429 5957 3723 1001 1713 1883 621 5411 4845 4966 8231 431 8234 7185 8492 3215 4029 7747 4144 771 6438 6112 947 6153 8651 4090 5292 3184 1899 4296 268 8120 8088 2794 6048 1954 2858 8332 5496 6463 415 5767 2613 7464 1479 1916 5083 482 6435 7449 6546 6540 8344 5814 2522 6884 4920 2295 2750 5098 3001 5639 1652 1784 6756 4580 3312 7283 346 266 7960 3531 4165 2648 1238 2333 3612 7518 3326 1265 5466 2844 7801 8419 533 3508 4909 8111 5120 3728 8548 3402 1862 3342 3931 5871 8528 7993 2672 4133 2029 7545 5623 640 6026 8304 6536 2542 2932 7721 187 2432 2589 8328 5565 1517 6995 8742 4248 6644 4234 1999 4866 5574 2369 5670 2406 772 5433 339 273 4796 678 4605 6763 8428 7373 6506 203 5660 1938 3398 2901 6675 1439 2693 3678 4791 4290 5588 599 3140 8002 3027 1246 3580 2167 3964 8393 2024 1977 2373 1346 4161 2587 4192 7770 2463 1941 8039 7742 89 7933 270 737 3341 3947 7017 2992 7758 3802 3843 7648 7090 2104 1507 3158 5631 2811 1487 2083 3518 619 2998 7154 5169 4174 322 7478 8014 8363 3775 2550 7317 8526 8617 8041 2895 7226 3325 5529 558 7799 701 831 1947 4207 6204 1796 6783 2520 2347 5319 4474 7549 757 7323 172 4725 5938 7179 7538 7971 5161 6648 1703 6943 1294 925 1281 5746 2271 3443 5966 6846 943 6631 844 4055 1371 6909 7673 853 8095 8442 4883 3679 4358 5900 2516 4773 7809 1491 1720 6771 3673 3555 7797 7131 1536 1431 643 1816 8149 2758 5710 4069 7608 3365 7489 850 5330 8025 7905 3556 638 4699 5861 8131 1539 5686 2868 1162 4480 2258 6421 334 3665 2232 3983 323 7535 880 7584 7389 502 23 288 7893 3311 1368 4078 748 3281 2670 3549 2615 5137 1532 1254 7229 8531 494 6238 2252 4817 5405 936 2376 5269 8569 3731 1659 2959 1257 171 3349 3021 2128 359 4746 2475 8646 2409 8453 7300 284 4573 8034 2297 4528 3637 6897 2065 1752 5999 8754 7554 6621 3872 4867 8197 4270 4291 1301 8181 4028 1996 1016 2063 6662 8046 2211 462 1688 2843 5364 197 8345 8226 1110 8487 5342 7004 6773 3499 4565 4978 3575 2217 7182 4693 8763 5878 55 7485 555 7909 7631 2790 1276 92 8260 4896 4962 7021 4789 5704 3129 27 4348 4863 6368 7877 7054 934 8164 7063 5331 7517 1590 976'
with open(r"C:\Users\Admin\Downloads\rosalind_lgis.txt",'r') as f:
    FL = f.readlines()
String = FL[1].replace('\n','').split(' ') 
string = [int(s) for s in String]   
n=len(string)
li = []
ld = []
lil = []
ldl = []
import numpy as np
import threading        #导入多线程库
#m1
'''def sort_iter(s,orderRev = False):
    lt = []
    sorted_s = sorted(enumerate(s),key=lambda x:x[1],reverse =orderRev)
    idx =[i[0] for i in sorted_s]        
    for i in range(len(idx)-1):
        lt.append(string[idx[i]])
        if orderRev == False:
            if idx[i+1] < idx[i]:
                if string[idx[i]+1] < string[idx[i]]:   
                    break
                elif idx[i] == len(idx)-2:
                    if string[idx[i]+1] > string[idx[i]]:
                        lt.append(string[idx[i]+1])
                    break
            else:
                if idx[i] == len(idx)-2:
                    if string[idx[i]+1] > string[idx[i]]:
                        lt.append(string[idx[i]+1])
                    break
        else:
            if idx[i+1] < idx[i]:
                if string[idx[i]+1] > string[idx[i]]:   
                    break
                elif idx[i] == len(idx)-2:
                    if string[idx[i]+1] < string[idx[i]]:
                        lt.append(string[idx[i]+1])
                    break
            else:
                if idx[i] == len(idx)-2:
                    if string[idx[i]+1] < string[idx[i]]:
                        lt.append(string[idx[i]+1])
                    break
    return lt'''



'''def compute_ldli(string):
    for si in range(1,len(string)):
        li.append(sort_iter(string[0:si+1],orderRev = False))
        lil.append(len(sort_iter(string[0:si+1],orderRev = False)))
        ld.append(sort_iter(string[0:si+1],orderRev = True))
        ldl.append(len(sort_iter(string[0:si+1],orderRev = True)))
    return [li,lil,ld,ldl]

#threading.Thread(target=compute_ldli(string)).start()
'''
#!!!!!!643125 6431 125     6432 125

for si in range(1,len(string)):        
    li.append(sort_iter(string[0:si+1],orderRev = False)[0])
    lil.append(len(sort_iter(string[0:si+1],orderRev = False)[0]))
    ld.append(sort_iter(string[0:si+1],orderRev = True)[0])
    ldl.append(len(sort_iter(string[0:si+1],orderRev = True)[0]))

li = []
ld = []
lil = []
ldl = []
tempi = 0
tempd = 1000000
for i in range(len(string)):
    templ = []
    templd = []
    if string[i] > tempi:
        templ.append(string[i])
        tempi = string[i]
    elif string[i] < tempd:
        templd.append(string[i])
        tempd = string[i]
        for j in range(i,len(string)):
            if string[j] > tempi:
                templ.append(string[j])
                tempi = string[j]
            elif string[j] <tempd:
                templd.append(string[j])
                tempd = string[j]
    li.append(templ)
    lil.append(len(templ))
    ld.append(templd)
    ldl.append(len(templd))
    
        

import numpy as np    
#print(' '.join(li[np.argmax(lil)]))
#print(' '.join(ld[np.argmax(ldl)]))
print(li[np.argmax(lil)])
print(ld[np.argmax(ldl)])

#recheck
for i in range(len(li)):
    if min(li[i]) != li[i][0]:
        li.remove(li[i])
        lil.remove(lil[i])
    if max(ld[i]) != ld[i][0]:
        ld.remove(ld[i])
        ldl.remove(ldl[i])
        
        
from itertools import combinations

list_i = []
list_d = []
bli = False
bld = False
for i in reversed(range(2, n+1)):
    for j in combinations(string, r=i):
        if list(j) == sorted(j) and not bli:
            list_i.append(j)
            bli = True
        if list(j) == sorted(j, reverse=True) and not bld:
            list_d.append(j)
            bld = True
        prin
    if bli and bld:
        break

print(' '.join(str(i)for i in list_a[0]))
print(' '.join(str(i)for i in list_d[0]))

#correction
with open(r"C:\Users\Admin\Downloads\rosalind_lgis.txt",'r') as f:
    FL = f.readlines()
String = FL[1].replace('\n','').split(' ') 
string = [int(s) for s in String]   
n=len(string)

def l_i_s(sequence):
    if not sequence:
        return []
 
    lengths = [1] * len(sequence)
    for i in range(1, len(sequence)):
        for j in range(i):
            if sequence[i] > sequence[j]:
                lengths[i] = max(lengths[i], lengths[j] + 1)
 
    max_length = max(lengths)
    lis = []
    for i in range(len(sequence) - 1, -1, -1):
        if lengths[i] == max_length:
            lis.append(sequence[i])
            max_length -= 1
 
    return lis[::-1]
 
def l_d_s(sequence):
    return l_i_s(sequence[::-1])
 
#sample_sequence = list(map(int, combined_sequence.split()))
 
lis = l_i_s(string)
lds = l_d_s(string)
 
print(' '.join(map(str, lis)))
print(' '.join(map(str, lds[::-1])))

#last try
li = []
ld = []
lil = []
ldl = []
mintemp = 10000
maxtemp = 0
li.append(string[0])
lil.append(1)
ld.append(string[0])
ldl.append(1)
for i in range(1,len(string)):
    litemp = []
    litemp.append(string[i])
    mintemp = string[i]
    for j in range(i-1,0,-1):
        if string[j] < mintemp:
            mintemp = string[j]
            litemp.append(mintemp)
        print(i,j,li,lil)       
    if len(litemp) >0:
        li.append(len(litemp))
        litemp = []
    
li = []
ld = []
lil = []
ldl = []

sorted_s = sorted(enumerate(string),key=lambda x:x[1],reverse =False)
idx =[i[0] for i in sorted_s]

templi = []
for i in range(len(string)):
    temp = string[i]
    templi.append(temp)
    for j in range(i):
        templi = sort_iter(string[j:i],orderRev = False)
        li.append(templi)
        templi = []
        
import numpy as np 
def lgis(path): 
    with open(path) as fp: 
        n=int(fp.readline()) 
        r=[int(i) for i in fp.readline().split()] 
        t=list(range(1,n+1)) 
    return lis(n,r,t),lis(n,r,t[::-1]) 
def lis(n,s,t): 
    a=[[0]*(n+1)for _ in range(n+1)] 
    b=[[2]*(n+1)for _ in range(n+1)] #1为斜上(来自a[i-1][j-1]) 2为上(来自a[i-1][j]) 3为左(来自a[i][j-1]) 
    for i in range(1,n+1): 
        for j in range(1,n+1): 
            score=[a[i-1][j-1]+1 if s[i-1]==t[j-1] else -1,a[i-1][j],a[i][j-1]] 
            index=np.argmax(score) 
            a[i][j]=score[index] 
            b[i][j]=index+1 
    return backtrack(b,s,n) 
def backtrack(b,s,n):
    i=n
    j=n
    res=[]
    while i>0 and j>0:
        if b[i][j]==1: 
            res.append(s[i-1])
            i-=1
            j-=1
        elif b[i][j]==2:
            i-=1 
        else:
            j-=1 
            print(i,j)
    return res[::-1] 

print(lgis(r"C:\Users\Admin\Downloads\rosalind_lgis.txt"))

#Genome Assembly as Shortest Superstring
test = r"E:\Rosalind\Rosalind_tttt.txt"
with open(r"C:\Users\Admin\Downloads\rosalind_long.txt",'r') as f:
    rfs=f.readlines()
with open(test,'r') as f:
    rfs=f.readlines()
sequence = []
tmp=[]
for lines in rfs:
    if '>' not in lines:
        #tmp.append(''.join(lines[0:lines.index('\n')]))
        print(lines.replace('\n',''))
        if '\n' in lines:
            tmp.append(lines.replace('\n',''))
        else:
            tmp.append(lines)
    else:
        if tmp:
            print(''.join(tmp))
            sequence.append(''.join(tmp))
            tmp = []
    if lines == rfs[-1]:
        sequence.append(''.join(tmp))

                       
superstring = sequence[0]
for s in sequence[1:]:
    i = superstring.find(s[0])
    if i <0:
        superstring += s
    else:
        while i < len(superstring):
            if len(superstring[i:])>= len(s):
                if superstring[i:i+len(s)] == s:
                    break
                else:
                    temp = superstring[i+1:].find(s[0])
                    if temp <0:
                        superstring += s
                        break
                    else:
                        i += temp+1
            else:
                if superstring[i:] == s[0:len(superstring[i:])]:
                    superstring += s[len(superstring[i:]):]
                    break
                else:
                    temp = superstring[i+1:].find(s[0])
                    if temp <0: 
                        superstring += s
                        break
                    else:
                        i += temp +1
print(superstring)

superstring = sequence[0]
for s in sequence[1:]:
    overlap = 0
    i = superstring.find(s[0])
    for sj in range(i,len(superstring)):
        for si in range(1,len(s)+1):    
            if superstring[sj:sj+si] == s[0:0+si] and si > overlap:
                overlap = si
            else:
                break
    if overlap ==0:
        superstring += s
    else:
        if overlap < len(s):
            superstring += s[overlap:]        
print(superstring)


import re 
def gen_cf(m, n):
    """该函数是用来计算 m序列尾部和 n序列头部重叠碱基个数
    """
    max_gen = 0
    for j in range(len(m), 3, -1):
        if m[-j:] == n[0:j]:
            max_gen = j
    return max_gen
 
 
def order(m, n):
    """该函数通过调用gen_cf函数，来判断 m序列尾部和n序列头部连接 还是
    n序列尾部和 m序列头部连接,最后返回一个列表。列表包含序列之间连接的信息
    还有重叠碱基数
    """
    order_1 = gen_cf(m, n)
    order_2 = gen_cf(n, m)
 
    if order_1 > order_2:
        return [1, order_1]
    else:
        return [2, order_2]
 
 
def add(m, n):
    """
    该函数通过调用order函数，将m n 连接起来
    """
    j = order(m, n)
    if j[0] == 1:
        return m + n[j[1]:]
    else:
        return n + m[j[1]:]
 
 
gen_begin = sequence.pop(0)
    # 
while sequence:
    l = []
    for i in sequence:
        l.append(order(gen_begin, i)[1])
    print(l)  #
    max_index = l.index(max(l))
    n = sequence.pop(max_index)    #
    gen_begin = add(gen_begin, n)        # 
 
print(gen_begin)

#Perfect Matchings and RNA Secondary Structures
with open(r"C:\Users\Admin\Downloads\rosalind_pmch.txt",'r') as f:
    FL = f.readlines()
String = ''.join(FL[1:]).replace('\n','').split(' ')[0]
string = [int(s) for s in String]   
n=len(String)
#for test
String = 'AGCUAGUCAU'
nA = String.count('A')
nU = String.count('U')
nC = String.count('C')
nG = String.count('G')

def factorial(N):
    if N == 0 or N == 1:
        return 1
    elif N >1:
        return N*factorial(N-1)

if nA == nU and nC == nG:
#perfect matching
    print(factorial(nA)*factorial(nG))
else:
#maximal matching
    if nA*nU*nC*nG == 0:
        print(0)
    else:
        result = 1
        minAU = min(nA,nU)
        maxAU = max(nA,nU)
        minCG = min(nC,nG)
        maxCG = max(nC,nG)
        while minAU > 0:
            result *= maxAU
            maxAU -= 1
            minAU -= 1
        while minCG > 0:
            result *= maxCG
            maxCG -= 1
            minCG -= 1
        print('1:',result)
        print('2:',factorial(maxAU)/factorial(maxAU-minAU))
        
#Partial Gene Orderings
'''
A partial permutation is an ordering of only k
 objects taken from a collection containing n
 objects (i.e., k≤n
). For example, one partial permutation of three of the first eight positive integers is given by (5,7,2)
.

The statistic P(n,k)
 counts the total number of partial permutations of k
 objects that can be formed from a collection of n
 objects. Note that P(n,n)
 is just the number of permutations of n
 objects, which we found to be equal to n!=n(n−1)(n−2)⋯(3)(2)
 in “Enumerating Gene Orders”.
 
 Given: Positive integers n and k
 such that 100≥n>0 and 10≥k>0.

Return: The total number of partial permutations P(n,k)
, modulo 1,000,000.
'''
N = 92 
k = 8
def factorial(N):
    if N == 0 or N == 1:
        return 1
    elif N >1:
        return N*factorial(N-1)

P = factorial(N)/factorial(N-k)
print(np.mod(P,1000000))

#Introduction to Random Strings 
'''
Modeling Random Genomesclick to collapse
We already know that the genome is not just a random strand of nucleotides; recall from “Finding a Motif in DNA” that motifs recur commonly across individuals and species. If a DNA motif occurs in many different organisms, then chances are good that it serves an important function.

At the same time, if you form a long enough DNA string, then you should theoretically be able to locate every possible short substring in the string. And genomes are very long; the human genome contains about 3.2 billion base pairs. As a result, when analyzing an unknown piece of DNA, we should try to ensure that a motif does not occur out of random chance.

To conclude whether motifs are random or not, we need to quantify the likelihood of finding a given motif randomly. If a motif occurs randomly with high probability, then how can we really compare two organisms to begin with? In other words, all very short DNA strings will appear randomly in a genome, and very few long strings will appear; what is the critical motif length at which we can throw out random chance and conclude that a motif appears in a genome for a reason?

In this problem, our first step toward understanding random occurrences of strings is to form a simple model for constructing genomes randomly. We will then apply this model to a somewhat simplified exercise: calculating the probability of a given motif occurring randomly at a fixed location in the genome.
'''
import numpy as np
s ='ACGATACAA'
p = [0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]
asw = [-5.737, -5.217, -5.263, -5.360, -5.958, -6.628, -7.009]

#exam
s = 'AGCTCGCAGTAGCTCGTGCCGGAATCTCTCACAAAGTTGGCCTATGCTACGGGATGCTCATATTTTACTACCCAGTATTGTGTCTGGT'
p = [0.096,0.148,0.200,0.223,0.288,0.373,0.405,0.455,0.538,0.569,0.641,0.711,0.743,0.806,0.883,0.915]
#method 1
result = []
for i in enumerate(p):
    tp = 0
    for si in s:
        if si in ['A','T']:
            tp +=  np.log10((1-p[i[0]])/2)
        elif si in ['C','G']:
            tp +=  np.log10(p[i[0]]/2)
    result.append(round(tp,3))

#Enumerating Oriented Gene Orderings
'''
Problem
A signed permutation of length n
 is some ordering of the positive integers {1,2,…,n}
 in which each integer is then provided with either a positive or negative sign (for the sake of simplicity, we omit the positive sign). For example, π=(5,−3,−2,1,4)
 is a signed permutation of length 5
.

Given: A positive integer n≤6
.

Return: The total number of signed permutations of length n
, followed by a list of all such permutations (you may list the signed permutations in any order).
'''
from  itertools import permutations
import numpy as np
def factorial(N):
    if N == 0 or N == 1:
        return 1
    elif N >1:
        return N*factorial(N-1)

n = 4
P = list(permutations(set([int(l) for l in np.linspace(-n,n,2*n+1)])-set([0]),n))
result = []
for item in P:
#    print(item,len(set([abs(a) for a in item])),len(item))
    if len(set([abs(a) for a in item])) == len(item):
#        P.remove(item)
        result.append(item)
print(factorial(n)*2**n)
for r in result:
    print(r[0],r[1],r[2],r[3])

#Finding a Spliced Motif
s = 'TGCCAATTGCTAGGCTGGCTAAGATGAGTGATCTGTAAATGCCATTCCTCCTAACTCAGTCTGGACTTCGTTTCAGGTGTAGGATGCAGCCGCGTTGCCATTCGCGGTACCAGACATAGGTGTCGCGAATGGCAATTTTACTGCCTAGGCTAAGACCCACAAGGTTCGTTCTGTAGAATCCCGAGCGGCGCGGACATAGTTTCATGGTCTGGTCTATGTGTTATATGCTTTGGTAAAGCAAGCTGCGAGCGGACAAGAGCTTGAACGTCCCATAAGCCGGGTATGAGTGACCCCCCGCAGGGACGCTACAAGACACACCCTAGCCCTGCCAAGGCGCCGATTTGCATATGAGACTGGAGCCGGACGCAGTAACAGCGGGTCTGGTTGTCGGTGGCAAATCTCTCTATAGTCCTGACGACCGTAGGCATTTATGTTTTAGTTCGCAGAAACAGAGCGAAGGAAAGGACCCTGGCAATGGGTCGTGATATTGCTACTGTGCACTTCCTGCGAGGGCATAACTCGTTCATCAAACATAGGGTGTCGACACGGCGGCATCCGTTTAGAGACTCATCTTATTACGGGGCCGTCTGGCCATCAGGTAATCTTTGAGATGCGGAGAGTCTTTTCCCGGAGAGAGTTGCTGCATGGGTAAATCACCACGTCCCCCCACCTTTTCGCCCTCTCAGATGCCACTGTTAGATATTAGTAGCGGGTTACCTGTTGCACACGTGAACTACAATAACTATAGCAACAGTAGCAGAGAGGAAACTAGGGGCAAAAAATAGTTAGGCGGGCTCGCCTTATAAAGGCACAAACGTAACCTGGGTCAGTTCTCACCGTCCTACAAATTAC'
p = 'GGTAGTCGAGCTCTTCTACGAAGTTC'
result = []

ci = 0
for c in p:
     fi = s.find(c, ci)
     result.append(fi + 1)
     ci = fi + 1  
print(' '.join(map(str,result)))


#Transitions and Transversions 
transitions = {'A':'G','G':'A','C':'T','T':'C'}
transversions = {'A':['T','C'],'T':['A','G'],'C':['A','G'],'G':['C','T']}
s1='ATTTGCAGGGTATGTCGCGTTATTGTGCTTACCGTTTTTCGCTGCCGGAGTGTAGGAATATGCCGTTATAGAGGGCGTAGAAGCAACGTAAATTATGTAGCTGGAGGCTTCATCCCTTGCCTTGAGAAAGACTGCTTATGGTTCTTCGCTCTTCGTCGTAAACGTATAAACACCCTGAGACCATCACACCGTGTCAACTTCTTCACCGTGTTACTGCAGCTGTTTAACAAGAGATTGAAACACGAATATAGCCTTGGGGACTTCGCTTTCCTGGCATTCCTCCTCCTTCGTCCAGGAAAGATTCCTCGTTGGAATCAGGACCAGCCTCAGGGCAGATCCAACACCTTACTTTTTTAAGGACGGCCTTCCACGTCAGTGGACCGCATGACCTAAGAACGAGCCCTGATTGAACGACAAAACGGCCTGTATGTCGAGCTTGAATAGTGTTGTACGCATATAGGCGCTTCGCATAGTGGGGCTTTTCTCTATCGATCCCTGAGTCTCGCCGTTCGAACTAGATTATCATGTCGAGTTGTTATTCCGCTGCGTGAGTGCTCAGTCCACACCCATTCGCTAGTCACGCAATGATCTGCTGCGGGGCCCAGTCAAGAATTTTTAAACCGAGCGCAGGTACACGGCACCACCGGCTCTATGACGTGGAAGGATGGCGTGAATGTCTCGTGGCACGGCAAGGCGCGTTTACATGCGTCTTCCGACTCAACAACGCCCTAGGCATAGGGATTTACACATTATCTAGCACTTCTCAACCAACCTCTGACTTCACAAATAGCGTAATAGGGGCAGCAGCCAGTGCCAAAGCTTATGCCCTATTTCTTTCAGTGCGAGGAACCAGACTGACGTATGCGGCGGTTGTTGCGAATAGTCTCACAACAGGCAGGGA'
s2='ACTCGCGGGGTACGTCACATTCTCGTGCTCACCGTTTGTCACTGCTAAAATGCAGGGACACGCTGGTATAGAATGCGTAGAAGCGACGGGAATCATGTCACGAGTAGCTTCACTTTCTACCTGGAGAACTACCGCTTACGGTTCTTGACGCTCTGTTGAAAACATAGAAGCGCCTTGAGATCATTAGGCCGCACCGACTCATTCGCCGCGGTACTGTAGCTGTTTGACAGGCGATCCCGGTACGAACAGAGCTTTTCGGACTTCCCTTCCCGGGCATCCCCCTTCGTTAGTCTAAGGAAGATGTCCCGTTGCAGCCAAGACTAGCTCCAAGACAGAGACAATTCCTAGCTTCCTCATGTCCAGACGCAGACGTCGGTCGACTGCATGCGCGAAAGACAAGTTCCGGACGAACAATAAAACGTCTTGTATGTTGGGTTTAAGTGCTCTTGTGTGCGAGAGGCCGCTTTGTATCGTAGAGATATCCCCTATCGATCCCTGAATGTGACCGTTCCAATCTGGTTATATCGCCGAATTGCCGTTCTTATGCTTGCGTGCCCAATCCATATCCACTCGCTAGTCATGTCATGATCTGTCGCAGGACTCAGCCAATAACCTTTAGACCGTGTGCTGACTCACGGTATCGCCGACTCTACGTTATGGAAGGATGGCGCGAATGTCCCGTAGCGCGACAAGGTGCAGCAGCATGCGGCCTCCGACTCAGCAACGCCCCAAGCGAAGGGATTTATACATCATATATCACTCCTCAACCAACATTCTGCTTCACAAATAGCCCCATAGGGGCAATGCTCAATGTCAAAGCTTATGTGCCACTTCTTGCGGCGCGAGAAACCCGACTGAAATATGCGACAGCCCTCGCGAATGACCCCACAACAAACACAGA'
trst = 0
trsv = 0
if len(s1) == len(s2):
    for i in enumerate(s1):
        if transitions[i[1]] == s2[i[0]]:
            trst += 1
        elif s2[i[0]] in transversions[i[1]]:
            trsv += 1
        else:
            if i[1] != s2[i[0]]:
                print('error!',i,transitions[i[1]],transversions[i[1]],s2[i[0]])
    print(round(trst/trsv,11))
    
#Completing a Tree
import numpy as np
N = 10
AdjL = [[1,2],[2,8],[4,10],[5,9],[6,10],[7,9]]
# 1 2        8
#      4      10 6
#       5 9 7 
#add  [2,4] [4,5]    
#     [2,3]  

with open(r"C:\Users\Admin\Downloads\rosalind_tree.txt",'r') as f:
    FL = f.readlines()
N = int(FL[0].replace('\n',''))
#data = [s.replace('\n','').split(' ') for s in FL[1:]]
data = [[int(s.replace('\n','').split(' ')[0]),int(s.replace('\n','').split(' ')[1])] for s in FL[1:]]
originData = np.array(data)
AdjL = data.copy()

def creatTr(L):
    treeL = []
    for t in L.copy():
        print(set(t))
        print('and',set(treeL))
        if len(treeL)==0 and t not in treeL:
            [treeL.append(titem) for titem in set(t)]
            L.remove(t)
        elif len(treeL)>0:
            if len(set(t).intersection(set(treeL))) >0:
                [treeL.append(s) for s in (set(t) - set(treeL))]
                L.remove(t)
    return L,treeL
LIST = []
if len(LIST) == 0:
    data,treeL = creatTr(AdjL)
    LIST.append(treeL)
while len(AdjL) >0:
    AdjL,treeL = creatTr(AdjL)
    LIST.append(treeL)
nodes = []
for s in LIST:
    if len(set(nodes).intersection(set(s)))== 0:
        nodes += s
n = len(set(nodes))

N-n + len(LIST)-1

N - len(originData)-1

#Catalan Numbers and RNA Secondary Structures
RNA ='UGGCUAAUCAAUAUAUUAAUUAAUUAUGCGCGCGUACGUACCCGAGCCGUCUAGCGGAGAUCUCCGGCGGACGUCGUGCACGAUCCAUGAGUUCCGAUCCGCGGGCGUGUGCAAUCUCGAUAUAUAUAUCUAUACGGGCACGAAAGAUCUAUAAUAUCGCAGCUGCGCGUAUGUACAUGCUGCAGCUAUAUUCGAACGAGCUCGUUGCAUUAUAACAUGACAUGAUACUAGUGCAACGCGCGGAUCCGAAUAGAUCUUAUACAUGUAUCGUCUUAAGC'
RNA= 'CGGCUGCUACGCGUAAGCCGGCUGCUACGCGUAAGC'
#Pair = {'A':'U','U':'A','C':'G','G':'C'}
#cPair = RNA[0]
#nPair = 0
def PRNA(string):
    Pair = {'':1,'AU':1,'UA':1,'CG':1,'GC':1}
    #if Pair[si] != string:
    if string not in Pair:
        nPair =0
        for i in range(1,len(string),2):        
            #if Pair[si] == string[i]:
                #nPair += PRNA(string[1],string[2:i]) * PRNA(string[i+1],string[i+2:]) 
#            print(string[0]+string[i])
            if (string[0]+string[i]) in Pair.keys():
                nPair += PRNA(string[1:i]) * PRNA(string[i+1:])
#                print(string[0],string[1:i],nPair)
        Pair[string] = nPair % 1000000
    return Pair[string]

print(PRNA(RNA))

with open(r"C:\Users\Admin\Downloads\rosalind_cat.txt",'r') as f:
    rf=f.readlines()
RNA = ''.join(rf[1:]).replace('\n','')
def PRNA(string,Pair):
    #if Pair[si] != string:
    nPair =0
    if string not in Result.keys():
        for i in range(1,len(string),2):        
            #if string[0]+string[i] in Pair.keys():
            print(string[0]+string[i])
            if string[0]+string[i] in Pair.keys():
                nPair += PRNA(string[1:i],Pair) * PRNA(string[i+1:],Pair) 
                print (nPair)
        Pair[string] = nPair
        Result[string] = nPair
    return Pair[string]

Pair = {'':1,'AU':1,'UA':1,'CG':1,'GC':1}
Result = {}
print(PRNA(RNA,Pair) % 1000000)


memorize = {} # 建立一个字典，存储已经出现过的字符串及不交叉完美匹配的数量
memorize[''] = 1 # 如果序列为空，说明只有这一种情况

def ismatch(c1, c2):
    """判断是否碱基配对的函数"""

    if (c1 == 'A' and c2 == 'U') or (c1 == 'U' and c2 == 'A') or \
            (c1 == 'G' and c2 == 'C') or (c1 == 'C' and c2 == 'G'):
        return 1
    else:
        return 0


def noncross(seq):
    """判断是否有不交叉的完美匹配"""

    if seq in memorize.keys(): # 如果这段序列之前已经被分析过了，直接取出结果即可
        return memorize[seq]
    if len(seq) % 2 == 1: # 如果这个序列长度是奇数，不可能存在完美匹配
        memorize[seq] = 0
        return 0

    if seq.count('A') != seq.count('U') or seq.count('G') != seq.count('G'): # 如果这个序列中配对的碱基数量不相同，不可能存在完美匹配
        memorize[seq] = 0
        return 0

    i = 1
    num = 0
    while i < len(seq): # 在序列中搜索所有可以与第一个碱基配对的碱基
        if ismatch(seq[0], seq[i]) == 1: # 如果第i个碱基配对
            num += (noncross(seq[1:i]) * noncross(seq[i+1:])) # 去检验被第k个碱基分出的两个序列
        i += 2 # 只需从第偶数个碱基中搜索
    memorize[seq] = num # 记录这个序列的不交叉完美匹配结果

    return num

res = noncross(RNA)
print(res % 1000000)

#Error Correction in Reads

with open(r"E:\Rosalind\Rosalind_52_rwsr.txt",'r+') as f:
    rfs = f.readlines()

with open(r"C:\Users\Admin\Downloads\rosalind_corr.txt",'r+') as f:
    rfs = f.readlines()
    
sequence = []
tmp=[]
for lines in rfs:
    if '>' not in lines:
        #tmp.append(''.join(lines[0:lines.index('\n')]))
        print(lines.replace('\n',''))
        if '\n' in lines:
            tmp.append(lines.replace('\n',''))
        else:
            tmp.append(lines)
    else:
        if tmp:
            print(''.join(tmp))
            sequence.append(''.join(tmp))
            tmp = []
    if lines == rfs[-1]:
        sequence.append(''.join(tmp))

def HammingD(s1,s2):
    temp = [1 for i in range(len(s1)) if s1[i] != s2[i]]
    return sum(temp)

dictC = {'A':'T','T':'A','C':'G','G':'C'}
oldr = []
newr = []
pair = []
for s in sequence:
#    [(oldr.append(s),newr.append(ss)) for ss in set(sequence)-set(s) if HammingD(s,ss)<=1]
    CR = ''.join([dictC[si] for si in s][::-1])  
    print(s,CR)                                                 
#    [(oldr.append(ss),newr.append(CR)) for ss in set(sequence)-set(s) if HammingD(ss,CR)<=1]
    pair.append([s,CR])
oldr=[k for k in sequence if pair.count(k)<2] 
newr=[k for k in sequence if pair.count(k)>=2] 
for so in oldr:
    for sn in newr:
        if HammingD(so,sn)==1:
            print(so+' -> '+sn)
            break
        if HammingD(so,''.join([dictC[si] for si in sn][::-1]))==1:
            print(so+' -> '+''.join([dictC[si] for si in sn][::-1]))
            break

ref=[k for s in sequence for k in [s,''.join([dictC[si] for si in s][::-1])]] 
old=[k for k in sequence if ref.count(k)<2] 
new=[k for k in sequence if ref.count(k)>=2] 
res=[] 
for s in old: 
    for t in new: 
        if HammingD(s,t)==1: 
            res.append(s+'->'+t) 
            print(res[-1])
            break 
        if HammingD(s,''.join([dictC[si] for si in t][::-1]))==1: 
            res.append(s+ '->'+''.join([dictC[si] for si in t][::-1])) 
            print(res[-1])
            break
        
#
n = 1306
if n>4:
    result = (2+n/2)*n/4/2
else:
    result = n/2
print(result)

n = 1306
if n>3:
    result = n-2
else:
    result = 0
print(result)

# k-Mer Composition
from itertools import *
string = 'CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGGCCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGTTTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCAAATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCGGGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGACTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTACCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG'
with open(r"C:\Users\Admin\Downloads\rosalind_kmer.txt",'r') as f:
    frl = f.readlines()
string = ''.join(frl[1:]).replace('\n','')    
k = 4
skeys = [''.join(v) for v in product('ACGT',repeat=k)] 
kmer = []
kmers = []
for s in skeys:
    kmer.append(string.count(s))
    n = 0
    for si in range(len(string)-k+1):
        if string[si:si+k].find(s)>=0:
            n +=1
    kmers.append(str(n))
print(' '.join(kmers))


#Speeding Up Motif Finding
string = 'CAGCATGGTATCACAGCAGAG'
with open(r"C:\Users\Admin\Downloads\rosalind_kmp.txt",'r') as f:
    rf=f.readlines()
string = ''.join(rf[1:]).replace('\n','')
#[print(i) for i in enumerate(string)]
#(0, 'C')
#(1, 'A')
#(2, 'G')
#(3, 'C')
#(4, 'A')
#(5, 'T')
#(6, 'G')
#(7, 'G')
#(8, 'T')
#(9, 'A')
#(10, 'T')
#(11, 'C')
#(12, 'A')
#(13, 'C')
#(14, 'A')
#(15, 'G')
#(16, 'C')
#(17, 'A')
#(18, 'G')
#(19, 'A')
#(20, 'G')
result = []
result.append(0)
i = 0
while i < len(string)-1:
    fs = string[0:i+1]
    bs = string[i+1:]
    n =  min(len(fs),len(bs))    
    commonn = 0
    for j in range(1,n+1):        
        if fs[0:j]== bs[0:j]:
            commonn = j
            result.append(commonn)
        elif commonn > 0:
            i = len(result)
            if string[0] != string[i] and string[:i].find(string[i]) >=0:
                idx=[si for si in range(len(string[:i])) if string[si]== string[i]]
                print('idx',idx)
                commonrn = 0
                for ii in idx:                    
                    for k in range(1,ii+1):
#                        print(k,(ii,string[ii:ii-k-1:-1]),(i,string[i:i-k-1:-1]))                        
                        #if string[ii:ii-k-1:-1] == string[i:i-k-1:-1]:
                        if ii-k-1 >=0:
#                            print('>',k,(ii,string[ii:ii-k-1:-1]),(i,string[i:i-k-1:-1]))                        
                            if string[ii:ii-k-1:-1] == string[i:i-k-1:-1]:
                                if k> commonrn:
                                    commonrn = k 
                            else:
                                break
#                            print(commonrn)
                        else:
#                            print('<',k,(ii,string[ii::-1]),(i,string[i:i-k-1:-1]))
                            if string[ii::-1] == string[i:i-k-1:-1]:
                                if k> commonrn:
                                    commonrn = k 
                            else:
                                break
#                            print(commonrn)
#                print(commonrn)
                if commonrn >1:
                    result.append(commonrn+1)
            else:
                i = len(result)-1
                commonn = 0
                break
        else:
            result.append(commonn)
            i = len(result)-1
            break
print(result,i,fs,bs,commonn)
   
'''
m2
'''     
string = 'CAGCATGGTATCACAGCAGAG'
with open(r"C:\Users\Admin\Downloads\rosalind_kmp.txt",'r') as f:
    rf=f.readlines()
string = ''.join(rf[1:]).replace('\n','')
        
def PRES(f,b):
    n = min(len(f),len(b))
    common = 0
    for i in range(1,n+1):
        if f[:i] == b[:i]:
            common = i
        elif f[:i] != b[:i] and common >0:
            return common
        else:
            return 0 
    return common
            
def SURFS(f,s):
    idx = [si for si in range(len(f)) if s[si] == s[len(f)]]
    common = 0
    for ii in idx:
        n = min(ii,len(f))
        for i in range(1,n+1):
            if f[ii:ii-i:-1] == s[len(f):len(f)-i:-1]:
                if i > common:
                    common = i
#                print('test1',common)
            elif f[ii:ii-i:-1] != s[len(f):len(f)-i:-1] and common>0:
                break
            else:
                break
#    print('test',common)
    return common

import time
a=time.time()  
result = [0]*len(string)
idx = [si for si in range(1,len(string)) if string[si] == string[0]]
for i in idx:
    fs = string[0:i]
    bs = string[i:]
    c = PRES(fs,bs)
    if result[i] == 0 and c > 0:
        result[i:i+c] = range(1,c+1)
        cb = SURFS(string[0:i+c],string)
        if result[i+c] == 0 and cb >1:
            result[i+c] = cb
t= time.time()-a    
#Out[99]: 11.05  


a=time.time()  
result = [0]*len(string)
idx = [si for si in range(0,len(string)) if string[si] == string[0]]
for i in range(1,len(idx)):
    fs = string[0:idx[i]]
    bs = string[idx[i]:]
    pc = c
    c = PRES(fs,bs)
#    print(len(result),idx[i],result[idx[i]],c)
    if result[idx[i]] == 0 and c > 0:        
#        print(i,c,idx[i],idx[i]+c,result[idx[i]:idx[i]+c])
        result[idx[i]:idx[i]+c] = range(1,c+1)
#        print(result)
    elif result[idx[i]] > 0 and c > 0:
#        print(i,c,idx[i],idx[i-1])
#        print(i,pc,idx[i-1]+pc)
#        print(i,idx[i]+c)
#        print(i,c,idx[i-1]+pc,idx[i]+c,result[idx[i-1]+pc:idx[i]+c],list(range(c-(idx[i-1]+pc-1-idx[i]-1),c+1)))
        result[idx[i-1]+pc:idx[i]+c] = range(c-(idx[i-1]+pc-1-idx[i]-1),c+1)
        #print(result)        
print(result)
R = ''
#[R+str(r) for r in result]
for r in result:
    R += str(r)+' '
print(R)
t= time.time()-a    
[print(r) for r in result]

a = time.time()
result = [0]*len(string)
longest_motif_length = 0 
for i in range(1, len(string)):
    for j in range(1, len(string)-i+1):
        if string[:i] == string[j:j+i]:
            result[j+i-1] = len(string[:i])
            longest_motif_length = len(string[:i])
            # print(i, j)

    if longest_motif_length < len(string[:i]):
        break
    
print(' '.join(map(str, result)))
t = time.time()-a

import threading
a = time.time()  
result = [0]*len(string)
#thread = threading.Thread(target=PRES,args=[fs,bs,])
idx = [si for si in range(1,len(string)) if string[si] == string[0]]
for i in idx:
    fs = string[0:i]
    bs = string[i:]
    c = threading.Thread(target=PRES,args=[fs,bs,]).start()
    if result[i] == 0 and c > 0:
        result[i:i+c] = range(1,c+1)
        if c >0:
            cb = threading.Thread(target=SURFS,args=[string[0:i+c],string,]).start()
        if result[i+c] == 0 and cb >1:
            result[i+c] = cb
t= time.time()-a    

 

import asyncio
async def my_coroutine():
    result = [0]*len(string)
    idx = [si for si in range(1,len(string)) if string[si] == string[0]]
    for i in idx:
        fs = string[0:i]
        bs = string[i:]
        c = PRES(fs,bs)
        if result[i] == 0 and c > 0:
            result[i:i+c] = range(1,c+1)
            if c >0:
                cb = SURFS(string[0:i+c],string)
            if result[i+c] == 0 and cb >1:
                result[i+c] = cb
a=time.time()  
asyncio.create_task(my_coroutine())
t= time.time()-a    


async def pres(f,b):
    n = min(len(f),len(b))
    common = 0
    for i in range(1,n+1):
        if f[:i] == b[:i]:
            common = i
        elif f[:i] != b[:i] and common >0:
            return common
        else:
            return 0 
    return common
            
async def surfs(f,s):
    idx = [si for si in range(len(f)) if string[si] == s[len(f)]]
    common = 0
    for ii in idx:
        n = min(ii,len(f))
        for i in range(1,n+1):
            if f[ii:ii-i:-1] == s[len(f):len(f)-i:-1]:
                if i > common:
                    common = i
#                print('test1',common)
            elif f[ii:ii-i:-1] != s[len(f):len(f)-i:-1] and common>0:
                break
            else:
                break
#    print('test',common)
    return common 

s1='AACCTTGG'
s2='ACACTGTGA' 
s1 = 'CAACATCGCAATTGCGGGCGACTTGCACAGGGATTCCGCTTGAGGTGTCATCATCCTGATAGGAGTGCCAGTGCTTGCGAAGTTAGTACTACATGGAATTTCTAATGTCGTTGTCCCAGGTGCGGATGTCCAAGCATCTTAAGTAAGAGTCGATATACGCTACAGGGTCCCTGTCGCAGTTATCTCCAGACAGTCATCGAGGCATCGGTAAAACTCTTAGTTGGAAGGCATATAGGAGTGTAAACCTCCTCCGGATGTCTCCACTGATTGGGCGATGAGCCGGCCGAAGATCTCGCACCATCCCGTGGTTGTGTGCAAAGTTGGTACGGGATCTATTGGGGCAGACGCCGGCGATGCGCGTTGTAGGGTCACGCTGACGACGAAAGCCGCTGCTCAGCCGAGTTAATGCGTTGCTACATCCTCTCATATCAGGACCCATGTGTTCTCGCTTACCAGCCGATATCTAAAATAAGTTATGGAGTTTAAGGATGGCGCAAGGATAGAATTGATCTATCTTTTCTCCTTTGGACATGGCGGGAGATGTGATTCCAGGTTTAGAAGCGGACCTGCAGAATGCCGCGTTTCTTGTAGTTCATAGTTTAATTTATACTCTCTAGGATTGACCGGGCACATCACCGAAGAGGTCACGCCAACTGACCCGATGGCTCTCCGGGTAGGCTGGTGAGGATGCGCCCCGCATTGAGGTCAACGTGTGGTGGTTAGTTGATACACAAATCCCACCTCCAGATGTAGGCTTCCACAAGAAGTGCATTTCACTCTAAGCACTAGGGAGCGTGTTATGTTGAGACTGCGGGTATCGTATATGGGTACTAGGGTTTTACTTCCGTCGACCCGATGTGAGCCTAGAGACGCTCCTTCGGATTTTAAAGATTGAGCTCAAGAAGGGCGCTACAGGCACTACAGAGGGGCGTGTCCGACCCGTTGAGCGAAGTTTGGAAGACAATTATCAGACTAGGAGGTAGCATACCTGGCAGT'
s2 = 'TTAACTACAGTTTCATCGTTCCTCGGTCTTTGTAACGGATCATGACAACTTAGGATTGTGTGAACCAGGTGATCAGCTGGGTGTGCACAACGCGCGGTAGAACTAACCTGCTTCGCAGACGCGCCTCATTAGGGTTAAGCCTAAACTTCCAATCTCGTAGAGACAGTTTTCTAGCGCCACTCCTCACCTGTGTTGGGATAGGTTACCATCCATGGGTTATTTCTCCACTGAAGATTCGTAAATGAATAGGGTAGCCCTGTCCGGTTCGCAAATCATAGCTTGAACTTGGACAAGAGTCGCTTATATTGTACTACGAGCAAAGCGGTGCACTAGGGGGCGTTGAGAAGACATTGTGCCTATACAGACTTTCTCAGGAAGACAGCCTGTAACCGGCCTTAGGGCCCGTTGCCACCTAAAATGTCGGATTATTGCAACGTTTCACGTTTGGCCATAGGCGTGACGTTTCTTGCCGGTGCCTTTGATCTTCAGGTGTCGGGATAATCGTTACCCTCCGATAGCGCCGGGAGGTCACGCTCTTCCTAAATAAGGACTGAATGTGAGCAATTCTTAAGGACCAAAAACGTTCGCCTCCCTCTCAGATGCTGAGTCAACTACTAGATGTCGGGAACATACAACGCAGTGAGTTATCGCGAAATGTAGACCCCCCCGTTAGGATGCAGTAGATTACCTAACTCCCCTTTCAGGTTTACAACCTCCTCTAATTTCCCGGCTTGTTTGCATTAAGAAGGAGTCCTTCCGGGTTTCGAGCAATTTAAGGTGATAGGCCCACTAACCTCAAATAGCGATGGCTCCCATCATCTGATTAAACGGACTCCCTGAGTACCGTGTGATGATTTTGACACAAAATTCGAGCAAGGTGTGCCGCTGGTGTAAGGGGGACGTCACGCGCCATGTCTCACGCAAATATCCGGCACTCAACCCGCGACAACATGAGA'

ss = [s for s in s2]
result = []
si = 0 
sj = 0
while si < len(s1)-1:
    #if s in ss:   
    if ss[sj:].find(s1[si])<0:
        si += 1
    while sj < len(ss):
        if s1[si] == ss[sj]:
            result.append(s1[si])
            sj += 1
            break
        sj += 1
    si += 1
print(''.join(result))

while si < len(ss)-1:
    #if s in ss:   
    if s1[sj:].find(ss[si])<0:
        si += 1
    while sj < len(s1):
        if ss[si] == s1[sj]:
            result.append(ss[si])
            sj += 1
            break
        sj += 1
    si += 1
print(''.join(result))
#jubuzuiyou
#dynamic matrix！！！！！
s = s1
t = s2

lengths = [[0 for j in range(len(t) + 1)] for i in range(len(s) + 1)]
for i, x in enumerate(s):
    for j, y in enumerate(t):
        if x == y:
            lengths[i + 1][j + 1] = lengths[i][j] + 1
        else:
            lengths[i + 1][j + 1] = max(lengths[i + 1][j], lengths[i][j + 1])

spliced_motif = ''
x, y = len(s), len(t)
while x * y != 0:
    if lengths[x][y] == lengths[x - 1][y]:
        x -= 1
    elif lengths[x][y] == lengths[x][y - 1]:
        y -= 1
    else:
        spliced_motif = s[x - 1] + spliced_motif
        x -= 1
        y -= 1
print(spliced_motif)

async def main():
    result = [0]*len(string)
    idx = [si for si in range(1,len(string)) if string[si] == string[0]]
    for i in idx:
        fs = string[0:i]
        bs = string[i:]
        #c = await PRES(fs,bs)
        ctsk = asyncio.create_task(pres(fs,bs))
        c = await ctsk
        if result[i] == 0 and c > 0:
            result[i:i+c] = range(1,c+1)
            if c >0:
                #cb = await SURFS(string[0:i+c],string)
                cbtsk = asyncio.create_task(surfs(string[0:i+c],string))
                cb = await cbtsk
            if result[i+c] == 0 and cb >1:
                result[i+c] = cb
a=time.time()  
asyncio.run(main())
t= time.time()-a  
  
#Ordering Strings of Varying Length Lexicographically
import numpy as np
seq = ['D','N','A']
n = len(seq)
result = []
def AS0(String):
    result = []
    for n in range(len(String)):
        result.append(list(String[n]*(len(String)+1)))
    return result
def AS1(String):
    result = []
    result.append('')
    for n in range(len(String)):
        result.append(String[n])
    return list(result)
def AS(String):
    result = []
    for n in range(len(String)):
        result.append(list(String[n]*(len(String)+1)))
        temp = []
        temp.append('')
        for n in range(len(String)):
            temp.append(String[n])
        result.append(temp)
    return result
def Iters(string):
    result = []
    n = len(string)
    temp = list(itertools.permutations(string, n-1))
    [temp.append((t,t)) for t in string]
    Temp = [list(t) for t in temp]
    temps = string.copy()
    temps.append('')
    for m in range(len(string)):
        t0 = list(string[m]+''*(n-1))
        result.append(t0)
    for t in Temp:
        for ti in temps:
            print(t,ti)
            #tj = t.append(ti)
            tj = []
            for titem in t:
                tj.append(titem)
            print(tj)
            tj.append(ti)
            result.append(tj)
    return result
SEQ = seq.copy()
SEQ.insert(0,'')
def compare(x,y):
    n = min(len(x),len(y))
    print('n',n)
    for i in range(n):
        if ''.join(SEQ).find(str(x[i]))<''.join(SEQ).find(str(y[i])):
            print('r',x,y,i,-1)
            return -1
        elif ''.join(SEQ).find(str(x[i]))>''.join(SEQ).find(str(y[i])):
            print('r',x,y,i,1)
            return 1
        else:
            if n==1:
                if len(x) < len(y):
                    print('r',x,y,i,-1)
                    return -1
                elif len(x) > len(y):
                    print('r',x,y,i,1)
                    return 1
            elif n>1:
                compare(x[i+1:],y[i+1:])
ii=Iters(string)   
from functools import cmp_to_key 
sorted(ii,key=cmp_to_key(compare))

def Iter(string):
    n = len(string)
    if n > 2:
        temp = list(itertools.permutations(string, n-1))
        [temp.insert(n*(n-i-1),(string[n-i-1],string[n-i-1])) for i in range(n)]
        Temp = [list(t) for t in temp]
        temps = string.copy()
        temps.insert(0,'')
        result = []    
        for t in enumerate(Temp):
            for ti in temps:
                print(t,ti)
                #tj = t.append(ti)
                tj = []
                for titem in t[1]:
                    tj.append(titem)
                print(tj)
                tj.append(ti)
                result.append(tj)
        for m in range(len(string)):
            tp = []
            tp.append(string[n-m-1])
            for mm in range(n-1):
                tp.append('')
            result.insert((n-m-1)*(n+1)*n,tp) 
    else:
        temp = list(itertools.permutations(string, n-1))
        Temp = [list(t) for t in temp]
        temps = string.copy()
        temps.insert(0,'')
        result = []
        for t in enumerate(temp):
            for ti in temps:
                print(t,ti)
                #tj = t.append(ti)
                tj = []
                tj.append(t[1][0])
                print(tj)
                tj.append(ti)
                result.append(tj)              
    return result

import itertools
itertools.combinations_with_replacement(seq, 2)
for i in range(n):
    result.append(seq[i])
    for j in range(n+1):
        result.append(seq[s]*(n+1)*n)
    temps = []
    for ss in range(n):
        temps.append(seq[ss]*(n+1))
    result.append(temps[:])
    tempss = []
    for sss in range(n):
        tempss.append(seq[sss])
        tempss.append('')
    result.append(tempss*n)
            
#
'''
sss='D N A R C'
S= sss.split(' ')
#S = ['C', 'Z', 'R', 'G', 'Q', 'T', 'N', 'H', 'B', 'V', 'I', 'O']
ids = list(enumerate(S))
k=4
Result =[]
idx =[]
stps=[]
seqs = itertools.combinations(S,k)
for s in seqs:
    result = list(itertools.product(s,repeat = 1))
    m0 = len(result) #4
    init = m0
    for j in range(2,k+1):
        temp = list(itertools.product(s,repeat = j))
        mj = len(temp) #16 64 256 
        stp = int(mj/m0) #4 fixed
        if j ==3:
           for m in range(k):
#            for m in range(j-1):
                    [idx.append(int(i+ii-m*(k+1))) for ii in range(stp)]#20212223 19202122 18192021 17181920  15161718...  10...  5678 4567 3456 2345
                    [stps.append(int(mj-(m*k+(1-k**(j-1))/(1-k)*k+1-i)*stp+ii)) for ii in range(stp)]#60616263 56575858 52... 48...   44454647 40... 36... 32...  28..
                    [result.insert(int(i+ii-m*(k+1)),temp[int(mj-(m*k+(1-k**(j-1))/(1-k)*k+1-i)*stp+ii)]) for ii in range(stp)]             
        elif j>3: 
           for n in range(j,j+k*(j-3)):
                for m in range(k):
    #            for m in range(j-1):
                    for i in range(init,init-k,-1):  
    #                [idx.append(i+ii-m*(k+1)) for ii in range(stp)]#20212223 19202122 18192021 17181920  15161718...  10...  5678 4567 3456 2345
    #                [stps.append((i-m*k)*stp-stp+ii) for ii in range(stp)]#60616263 56575858 52... 48...   44454647 40... 36... 32...  28..
    #                [result.insert(i+ii-m*(k+1),temp[(i-m*k)*stp-stp+ii]) for ii in range(stp)]            
                        [idx.append(int(i+ii-m*(k+1)-(n-j)*init/k)) for ii in range(stp)]#20212223 19202122 18192021 17181920  15161718...  10...  5678 4567 3456 2345
                        [stps.append(int(mj-((n-j)*k*k+m*k+(1-k**(j-1))/(1-k)*k+1-i)*stp+ii)) for ii in range(stp)]#60616263 56575858 52... 48...   44454647 40... 36... 32...  28..
                        [result.insert(int(i+ii-m*(k+1)-(n-j)*init/k),temp[int(mj-((n-j)*k*k+m*k+(1-k**(j-1))/(1-k)*k+1-i)*stp+ii)]) for ii in range(stp)]            
        else:
            for i in range(m0,0,-1):  
                [idx.append(i+ii) for ii in range(stp)]#4567 3456 2345 1234
                [stps.append(i*stp-stp+ii) for ii in range(stp)]#12131415 891011 4567 0123
                [result.insert(i+ii,temp[i*stp-stp+ii]) for ii in range(stp)]            
        m0 = len(temp) #9
        init = len(result) #(（1-q**n）/(1-q)) = init/k    q=k
    Result.append(result)
            
R = pd.DataFrame(np.concatenate(Result)).drop_duplicates()
for row in Result:
    print(''.join(row))
'''
#Maximum Matchings and RNA Secondary Structures

string ='AUGCUUC'
string ='AUAAAUGCCCUUUGGCGAAUCCUACAAGACUUGUGUCAGGAUCCCUCAAGAUAAAUGCCAAGGAAAUGGGCUGCUCCCCUCUAUAGACGUCAG'
def factorial(N):
    if N == 0 or N == 1:
        return 1
    elif N >1:
        return N*factorial(N-1)
nA=string.count('A')
nU=string.count('U')
nG=string.count('G')
nC=string.count('C')
result = 1
if Na*Nu*Ng*Nc >0:
    minAU = min(nA,nU)
    maxAU = max(nA,nU)
    minCG = min(nC,nG)
    maxCG = max(nC,nG)
    while minAU > 0:
        result *= maxAU
        maxAU -= 1
        minAU -= 1
    while minCG > 0:
        result *= maxCG
        maxCG -= 1
        minCG -= 1
    print('1:',result)

def MINMAX(N,n):
    result=1    
    if N*n >0:
        minN = min(N,n)
        maxN = max(N,n)
        if minN == maxN:
            result = factorial(maxN)
        else:
            result = factorial(maxN)/factorial(maxN-minN)
    print('2:',result)
    return result

MINMAX(nA,nU)*MINMAX(nC,nG)

#Creating a Distance Matrix
rfs = ['>Rosalind_9499\n','TTTCCATTTA\n','>Rosalind_0942\n','GATTCATTTC\n','>Rosalind_6568\n','TTTCCATTTT\n','>Rosalind_1833\n','GTTCCATTTA\n']


with open(r"C:\Users\Admin\Downloads\rosalind_pdst.txt",'r') as f:
    rfs=f.readlines()
sequence = []
tmp=[]
for lines in rfs[1:]:
    if '>Rosalind_' not in lines:
        #tmp.append(''.join(lines[0:lines.index('\n')]))
        print(lines.replace('\n',''))
        tmp.append(lines.replace('\n',''))
        if lines == rfs[-1]:
            sequence.append(''.join(tmp))
    else:
        if tmp:
            print(''.join(tmp))
            sequence.append(''.join(tmp))
            tmp = []

import numpy as np
n=len(sequence)
D=np.zeros((n,n))
m=len(sequence[0])
for i in range(n-1):
    for j in range(i+1,n):
        D[i][j] = round(len([(sequence[i][t],sequence[j][t]) for t in range(m) if sequence[i][t]!=sequence[j][t]])/m,5)
        D[j][i] = D[i][j]
[[format(ll,'.5f') for ll in row] for row in D]
DM = np.array([[format(ll,'.5f') for ll in row] for row in D],dtype=float).reshape(n,n)
np.set_printoptions(precision=5)
print(np.around(DM,5))



#Reversal Distance
import numpy as np
import itertools
import copy
N=[]
with open(r"C:\Users\Admin\Downloads\rosalind_test.txt",'r') as f:
    rfs=f.readlines()
A0 = []
B0 = []
for i in range(len(rfs)):
    if np.mod(i,3)==0:
        A0.append(rfs[i].replace('\n','').split(' '))
    elif np.mod(i,3)==1:
        B0.append(rfs[i].replace('\n','').split(' '))        
for i in range(len(A)):
    n=0
    tocheck=list(range(len(A[i])))
#    print('A:',A0[i],'B:',B0[i],'bc:',tocheck)
    for j in tocheck.copy():
        if A0[i][j] == B0[i][j]:
            tocheck.remove(j)
        elif A0[i][int(B0[i][j])-1] == B0[i][j] and B0[i][int(B0[i][j])-1] == A0[i][j]:
            n += 1
            tocheck.remove(j)
            tocheck.remove(int(B0[i][j])-1)            
#    print(i,tocheck,n)
    if not tocheck:
        minn = 0
    else:    
        tocheck_all = list(itertools.permutations(tocheck))
        minn = n +100
        for checks in tocheck_all:
            A = copy.deepcopy(A0)#A0.copy()
            B = copy.deepcopy(B0)#B0.copy()
            nn = n
            cl= list(checks)
#            print(cl,A,B)
#            print('bcl:',B[i],A[i])
            for k in cl.copy():                
                tmp = [ii for ii in range(len(B[i])) if B[i][ii]==A[i][k]]
                temp = B[i][k]
                B[i][k] = A[i][k]
                B[i][tmp[0]] = temp
                nn+=1
#                print('nBA',nn,B[i],A[i])
                if nn >= minn or B[i] == A[i]:
                    break
                cl.remove(k)
                #print(cl)
#            print('nn:',nn)    
            if nn < minn:
                minn = nn
#                print(cl,A[i],B[i],nn)
#    print(i,minn)
    N.append(minn)           
    
    

import collections
 
def get_all_permutations(s): # 起生成器作用的函数，作用是依次生成排列
     for i in range(len(s)):
         for j in range(i + 2, len(s) + 1):
             yield s[:i] + s[i:j][::-1] + s[j:]
 
 
def get_reversal_distance(p1, p2): # 计算翻转距离的函数
     if p1 == p2: # 如果两个排列相等，可以直接返回翻转距离为0
         return 0
 
     target = tuple(p2) # 排列p2作为翻转的目标排列
     fromfirst = {tuple(p1): 0} # 创建字典，key为排列，value为初始排列翻转几次能得到该排列，先把p1放进去
     q = collections.deque((p1,)) # deque是一种数据对象，类似列表，但插入和删除的效率更高。创建这么一个对象q，                                                          #并把p1放进去。不难发现，q中数据是按翻转次数排列的。

 
     while len(q): # 如果q中还有排列，就继续循环
         s = q.popleft() # 把q中第一个排列移除，放到s中
         c = fromfirst[s] # 从字典中查出来该排列对应的翻转次数
 
         for j in get_all_permutations(s): # j接受生成器产生的排列
             if j == target: # 如果新产生的排列和目标排列相同，意味着最小翻转次数找到了，就可以把次数直接返回输出了
                 return c + 1
 
             if not j in fromfirst: # 如果排列还不在字典里，把它添加进去
                 fromfirst[j] = c + 1
 
                 if c != 4: # 如果翻转已经超过5次了，就先停止，不再记录次数更多的排列
                     q.append(j)
 
     # 如果之前的5次翻转都没翻转出目标排列，就从另一头开始
     fromsecond = {tuple(p2): 0} # 创建新的字典，key为排列，value为初始排列翻转几次能得到该排列，把p2放进去
     target = tuple(p1) # 这次以排列p1作为翻转的目标排列
     q = collections.deque((p2,)) # 把p2排列放进q
     answer = 100000 
 
     while len(q): # 和上个翻转循环过程含义相同，这次翻转不超过4次
         s = q.popleft()
         c = fromsecond[s]
 
         if c == 4:
             break
 
         for j in get_all_permutations(s):
             if j == target:# 如果新产生的排列和目标排列相同，意味着最小翻转次数找到了，就可以把次数直接返回输出了
                 return c + 1
 
             if not j in fromsecond:
                 fromsecond[j] = c + 1
 
                 if c != 3:
                     q.append(j)
 
             if j in fromfirst: # p2经过翻转和第一次翻转循环产生的某个排列相同
                 answer = np.min([answer, fromfirst[j] + fromsecond[j]]) # 把两次翻转的次数相加即为最终结果
     return answer
 
 
if __name__ == __main__: # 决定以下代码是否需要运行
     distances = []
 
     with open(r"C:\Users\Admin\Downloads\rosalind_rear.txt",'r') as f:
         dataset = list(map(str.strip, f.readlines())) # 用map函数读入题目
 
     for i in range(0, len(dataset), 3): # 整理输入数据
         s = tuple(map(int, dataset[i].split()))
         t = tuple(map(int, dataset[i + 1].split()))
 
         distances.append(get_reversal_distance(t, s))
 
     print(''.join(map(str, distances))) # 结果
    
#More Random Motifs
N = 90000
p = 0.6
s = 'ATAGCCGA'
import numpy as np
import itertools
import math
n = int(np.floor(len(s)*p))
temp = list(itertools.permutations(range(len(s)),n))
tempCG = list(itertools.permutations(range(n),int(n/2)))
tempAT = list(itertools.permutations(range(len(s)-n),int((len(s)-n)/2)))
'''(1+2+...+N)/(len(temp)*len(tempCG)*len(tempAT)*N)'''
(1+N)*N/2/(len(temp)/len(list(itertools.permutations(range(len(s)))))*len(tempCG)/len(list(itertools.permutations(range(n))))*len(tempAT)/len(list(itertools.permutations(range(n))))*N)
m = (math.comb(len(s),n)*math.comb(int(n),int(n/2))*(2**(len(s)-n)))
mmm = math.comb(len(s),n)*math.comb(int(n),int(n/2))*2*(2**(len(s)-n))
1-math.comb(mmm-1+N-1,mmm-1-1)/math.comb(mmm+N-1,mmm-1)  
'''def randomStrings(s, A):
        n = s.count("G") + s.count("C")
        return([(len(s)-n)*log10(1-gc) + n*log10(gc) - len(s)*log10(2) for gc in A])'''
N = 89885 
p = 0.501820
s= 'TCTAACGT'

AT = s.count('A')+s.count('T')
CG = s.count('C')+s.count('G')
prob = ((1-p)/2)**AT*((p/2)**CG)
P = 1-(1-prob)**N
print(round(P,3))

#Counting Subsets
import numpy as np    
m = 907  
print(np.mod(2**m,1000000))
#Introduction to Alternative Splicing
import numpy as np
def fact(k):
    if k == 1:
        return 1
    while k>1:
        f,k = k*(k-1),k-1
    return f
def sf(n,m):
    s = 0
    while m<n:
        print(s)
        s,m = s+fact(n)/fact(n-m)/fact(m),m+1
    return s+1
a=3
n=6
print(np.mod(sf(n,a),1000000))
n=1725
a=1272
import math
ss = 0
for b in range(a,n+1):
    temp = math.comb(n,b)
    ss += temp
print(np.mod(ss,1000000))
