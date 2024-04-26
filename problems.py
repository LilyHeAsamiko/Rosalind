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


