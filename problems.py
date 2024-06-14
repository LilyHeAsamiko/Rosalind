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
