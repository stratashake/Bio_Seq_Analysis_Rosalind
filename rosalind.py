##Problem 1
##Count base pairs in a string

# with open("./Rosalind_files/rosalind_dna.txt") as f:
#     dna = f.read().strip()
# print(dna)

# A = dna.count('A')
# C = dna.count('C')
# T = dna.count('T')
# G = dna.count('G')
# print(A, C, G, T)

###############################################################################################
#Problem 2
##convert DNA to RNA

# def dna2rna(string):
#     print(string.replace("T", "U"))

# with open("./Rosalind_files/rosalind_rna.txt") as f:
#     dna = f.read().strip()
           
# dna2rna(dna)

##############################################################################################
##Problem 3
##Reverse complement DNA

# def dnarevcomp(s):
#     complement = {"A":"T", "T":"A", "C":"G", "G":"C"}
#     return "".join(complement[base] for base in reversed(s))

# with open("./Rosalind_files/rosalind_revc.txt") as f:
#     dna = f.read().strip()
# print(dnarevcomp(dna))

# #Problem 3 easier solution
# print(dna[::-1].translate(str.maketrans('ACGT', 'TGCA')))

##################################################################################################
##Problem 4
##Recursion/dynamic programming

# with open("./Rosalind_files/rosalind_fib.txt") as f:
#     line = f.readline()
#     num1, num2 = map(int, line.split())

# def rabbits(n, k):
#     if n < 3:
#         return 1
#     else:
#         return rabbits(n-1, k) + k*rabbits(n-2, k)
# print(rabbits(num1, num2))

# #better memory solution
# memo = {}
# def fib(n,k):
#     args = (n, k)
#     if args in memo:
#         return memo[args]  # Aha! We have computed this before!

#     # We haven't computed this before, so we do it now
#     if n == 1:
#         ans = 1
#     elif n == 2:
#         ans = 1
#     else:
#         ans = fib(n-1, k) + k * fib(n-2, k)
#     memo[args] = ans  # don't forget to remember the result!
#     return ans

# print(fib(num1, num2))

######################################################################################################
# #Problem 5
# #Find the highest percentage GC fasta file and report as name and percentage

# from Bio.SeqUtils import gc_fraction
# from Bio import SeqIO
# from operator import itemgetter
# records = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_gc.txt", "fasta"):
#     records.append((rec.id, gc_fraction(rec.seq)))

# max_record = (max(records, key=itemgetter(1)))
# gc_percentage = round(max_record[1] *100, 6)
# print(max_record[0])
# print(f"{gc_percentage}%")

######################################################################################################
# #Problem 6
# #Counting Point Mutations; Calculate the Hamming distance

# with open("./Rosalind_files/rosalind_hamm.txt") as f:
#     lines = f.readlines()

# s1 = lines[0]
# s2 = lines[1]
# h = 0

# length = len(s1)
# for x in range(length):
#     if s1[x] != s2[x]:
#         h += 1
# print(h)

######################################################################################################
# #Problem 7
# #Mendel's First Law; Mendelian Inheritance

# with open("./Rosalind_files/rosalind_iprb.txt") as f:
#     line = f.readline()
#     num1, num2, num3 = map(int, line.split())

# #Someone else's code
# def firstLaw(k,m,n):
#     N = float(k+m+n)
#     return 1 - ( m*n + .25*m*(m-1) + n*(n-1) ) / ( N*(N-1) )

# print(firstLaw(num1, num2, num3))

# #finally figured it out
# def mendel(D,d,r):
#     T = float(D+d+r)
#     top = ((D*(D-1))+(D*d)+(D*r)+(r*D)+(d*D)+(0.75*d*(d-1))+(0.5*d*r)+(0.5*r*d))
#     bottom = (T*(T-1))
#     ans = top/bottom
#     print(ans)

# mendel(num1, num2, num3)

######################################################################################################
# #Problem 8
# #Translating RNA into Protein

# codons = {    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
#               'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
#               'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
#               'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
#               'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
#               'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
#               'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
#               'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#               'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
#               'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
#               'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
#               'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
#               'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
#               'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
#               'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
#               'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

# with open("./Rosalind_files/rosalind_prot.txt") as f: 
#     rna = f.read().strip()

# def translate(rna):
#     protein = ""
#     start = rna.find("AUG")
#     for x in range(start, len(rna), 3):
#         codon = (rna[x:x+3])
#         if codon in codons:
#             if codons[codon] == "Stop":
#                 break
#             else:
#                 protein += codons[codon]
#     return protein
# print(translate(rna))

######################################################################################################
# #Problem 9
# #Finding a Motif in DNA; substrings
# import re

# with open("./Rosalind_files/rosalind_subs.txt") as f:
# 	lines = f.readlines()
	
# s = lines[0].strip()
# print(s)
# t = lines[1].strip()
# print(t)

# def CntSubstr(pattern, string):
# 	a = [(m.start() + 1) for m in re.finditer(
# 		'(?={0})'.format((pattern)), string)]
# 	return a

# # Calling the function
# result = (CntSubstr(t, s))
# result_str = ' '.join(map(str, result)) # convert each integer to a string and concatenate them without commas
# print(result_str) #ignore dna sequence

######################################################################################################
# #Problem 10
# #Consensus and Profile
# from Bio import SeqIO

# sequences = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_cons.txt", "fasta"):
#     str_seq = str(rec.seq)
#     sequences.append(str_seq)    

# #These sequences are just for troubleshooting, call it legacy code.
# # sequences = ["ATCCAGCT",
# #              "GGGCAACT",
# #              "ATGGATCT",
# #              "AAGCAACC",
# #              "TTGGAACT",
# #              "ATGCCATT",
# #              "ATGGCACT"]

# n = len(sequences[0])

# profile_matrix = {
#     'A': [0]*n,
#     'C': [0]*n,
#     'G': [0]*n,
#     'T': [0]*n
#     }

# for dna in sequences:
#     for position, nucleotide in enumerate(dna):
#         profile_matrix[nucleotide][position] += 1
# # print(profile_matrix[nucleotide])

# result = []

# for position in range(n):
#     max_count = 0
#     max_nucleotide = None
#     for nucleotide in ["A", "C", "G", "T"]:
#         count = profile_matrix[nucleotide][position]
#         if count > max_count:
#             max_count = count
#             max_nucleotide = nucleotide
#     result.append(max_nucleotide)

# consensus = ''.join(result)
# print(consensus)
# for key in profile_matrix:
#     print(key + ':', end=' ')
#     print(' '.join(str(i) for i in profile_matrix[key]))
# print(profile_matrix)

######################################################################################################
# #Problem 11
# #Mortal Fibonacci Rabbits

# with open("./Rosalind_files/rosalind_fibd.txt") as f:
#     line = f.readline()
#     num1, num2 = map(int, line.split())

# n = num1                                                                       
# m = num2                                                                       
# bunnies = [1, 1]                                                               
# months = 2                                                                    
# while months < n:                                                              
#     if months < m:                                                             
#         bunnies.append(bunnies[-2] + bunnies[-1])                              
#     elif months == m:                                      
#         bunnies.append(bunnies[-2] + bunnies[-1] - 1)                          
#     else:                                                                      
#         bunnies.append(bunnies[-2] + bunnies[-1] - bunnies[-(m + 1)])                                                           
#     months += 1                                                               
# print(bunnies[-1])

######################################################################################################
# #Problem 12
# #Overlap graphs
# from Bio import SeqIO
# fastas = {}
# for rec in SeqIO.parse("./Rosalind_files/rosalind_grph.txt", "fasta"):
#     fastas[(rec.seq)]=[(rec.id)]

# seq = []
# for x in fastas.keys():
#     seq.append(x)

# matches = []
# for i in range(len(seq)):
#     for j in range(len(seq)):
#         if i != j and (seq[j][:3]) == (seq[i][-3:]):
#             print(fastas[seq[i]][0], fastas[seq[j]][0])

######################################################################################################
# #Problem 13
# #Calculating Expected Offspring

# pop = []
# with open("./Rosalind_files/rosalind_iev.txt") as f:
#     for line in f:
#         numbers = line.strip().split()
#         for number in numbers:
#             pop.append(int(number))
            
# a = 1
# b = 1
# c = 1
# d = 0.75
# e = 0.5

# print(2*a*pop[0] + 2*b*pop[1] + 2*c*pop[2] + 2*d*pop[3] + 2*e*pop[4])

######################################################################################################
##Problem 14
##Finding a Shared Motif; Longest common substring

# from Bio import SeqIO
# import time

# start_time = time.time()
# fastas = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_lcsm.txt", "fasta"):
#     str_seq = str(rec.seq)
#     fastas.append(str_seq)  


# n = len(fastas)
# s = fastas[0]
# l = len(s)

# stems = []
# for i in range(l):
#     for j in range(i + 1, l + 1):
#         stem = s[i:j]
#         if len(stem) != 1:
#             all_present = True
#             for x in fastas[1:]:
#                 if stem not in x:
#                     all_present = False
#                     break
#             if all_present:
#                 stems.append(stem)

# print(max(stems, key = len))

# end_time = time.time()
# elapsed_time = end_time - start_time
# print("Elapsed time:", elapsed_time, "seconds")

######################################################################################################
# #Problem 15
# #Independent Alleles

# from math import factorial

# pop = []
# with open("./Rosalind_files/rosalind_lia.txt") as f:
#     for line in f:
#         numbers = line.strip().split()
#         for number in numbers:
#             pop.append(int(number))

# k = pop[0]
# n = pop[1]
# P = 2**k
# probability = 0

# for i in range(n, P+1):
#     prob = (factorial(P)/(factorial(i) * factorial(P-i))) * (0.25**i) * (0.75**(P - i))
#     probability += prob
# print(round(probability, 3))

######################################################################################################
# #Problem 16
# #Find protein motif

# import requests
# import re
# from Bio import SeqIO
# from io import StringIO

# proteins_full = []
# with open ("./Rosalind_files/rosalind_mprt.txt") as f:
#     for line in f:
#         numbers = line.strip()
#         proteins_full.append(numbers)

# proteins_short = []
# for x in proteins_full:
#     short = x.split('_')[0]
#     proteins_short.append(short)

# def protein_fasta(organism):
#     url = f'https://www.uniprot.org/uniprot/{organism}.fasta'
#     fastas = requests.get(url).text
#     return fastas

# sequences = []
# for x in proteins_short:
#     sequences.append(protein_fasta(x))

# l = len(proteins_short)
# test = sequences[0:l]
# strings = '\n'.join(test)
# wrapper = StringIO(strings)

# proteins_short_run = []
# for rec in SeqIO.parse(wrapper, "fasta"):
#     proteins_short_run.append(str(rec.seq))

# motif = "(?=N[^P][ST][^P])"
# m_index = []
# ans = {}

# def find_motif(fasta):
#     return [(match.start()+1) for match in re.finditer(motif, fasta)]

# for x in proteins_short_run:
#     m_index.append(find_motif(x))

# for x in proteins_full:
#     ans[x] = []

# for i, key in enumerate(ans.keys()):
#     for x in m_index[i]:
#         ans[key].append(x)

# new_dict = {k: v for k, v in ans.items() if len(v) != 0}
# for keys, values in new_dict.items():
#     print(keys) 
#     print(*values)

######################################################################################################
# #Problem 17
# #Find protein motif

# from operator import countOf

# codons = {    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
#               'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
#               'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
#               'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
#               'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
#               'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
#               'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
#               'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#               'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
#               'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
#               'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
#               'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
#               'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
#               'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
#               'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
#               'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

# with open("./Rosalind_files/rosalind_mrna.txt") as f: 
#     protein = f.read().strip()

# count = 1
# for aa in protein:
#     count *= (countOf(codons.values(), aa))
# count = count*3 #accounting for the stop codon
# print(count % 1000000)

# #Protein generator, unrelated to this problem but useful for trouble shooting
# # k = 100
# # use codons.values() to get a list of amino acids, excluding 'Stop'
# # protein = ''.join(random.choices([aa for aa in codons.values() if aa != 'Stop'], k=k))
# # print(len(protein))

######################################################################################################
# #Problem 18
# #Open Reading Frames

# from Bio import SeqIO
# import re
# DNA = ""
# for rec in SeqIO.parse("./Rosalind_files/rosalind_orf.txt", "fasta"):
#     DNA = str(rec.seq)
# DNA = DNA.replace('\n', '')

# codons = {    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
#               'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
#               'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
#               'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
#               'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
#               'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
#               'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
#               'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#               'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
#               'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
#               'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
#               'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
#               'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
#               'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
#               'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
#               'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


# RNA = ""
# RNA = DNA.replace("T", "U")

# def dnarevcomp(x):
#     complement = {"A":"U", "T":"A", "C":"G", "G":"C"}
#     return "".join(complement[base] for base in reversed(x))

# RNA_rev = dnarevcomp(DNA)

# pattern = re.compile(r'(?=(AUG(?:(?!UAA|UAG|UGA)...)*(?:UAA|UAG|UGA)))') #'forward facing' regex that gives us all the orf's for RNA

# def forward(RNA):
#     """Using our regex to find all ORFs with start and stop codon
#     and at least one amino acid in between"""
#     RNA_orfs = (pattern.findall(RNA)) #forward strand search
#     return RNA_orfs
# def reverse(rev_RNA):
#     """Same as above, but the for the reverse complement strand"""
#     Rev_RNA_orfs = (pattern.findall(rev_RNA)) #backward strand  search
#     return Rev_RNA_orfs

# RNA_orfs = (forward(RNA))
# Rev_RNA_orfs = (reverse(RNA_rev))

# def translate(rna):
#     protein = ""
#     start = rna.find("AUG")
#     for x in range(start, len(rna), 3):
#         codon = (rna[x:x+3])
#         if codon in codons:
#             if codons[codon] == "Stop":
#                 break
#             else:
#                 protein += codons[codon]
#     return protein

# final_list = []

# for x in RNA_orfs:
#     final_list.append(translate(x))
# for x in Rev_RNA_orfs:
#     final_list.append(translate(x))

# final_set = set(final_list)
# for x in final_set:
#     print(x)

######################################################################################################
# #Problem 19
# #Enumerating Gene Orders; Permutations

# from itertools import permutations

# with open("./Rosalind_files/rosalind_perm.txt") as f:
#     number = int(f.read().strip())

# perms = []

# def calc_permutations(x): 
#     global perms
#     for perm in permutations(range(1, x + 1)):
#         perms.append(' '.join(map(str, perm)))
#     return perms

# calc_permutations(number)
# print(len(perms))
# print('\n'.join(perms))

######################################################################################################
# # Problem 20
# # Calculating Protein Mass

# #from Bio.SeqUtils.ProtParam import ProteinAnalysis #doesn't work, gives slightly different values
# #protein = ProteinAnalysis(x)

# with open ("./Rosalind_files/rosalind_prtm.txt") as f:
#     x = f.read().strip()

# monoisotopic_mass_table = {
#     'A': 71.03711,
#     'C': 103.00919,
#     'D': 115.02694,
#     'E': 129.04259,
#     'F': 147.06841,
#     'G': 57.02146,
#     'H': 137.05891,
#     'I': 113.08406,
#     'K': 128.09496,
#     'L': 113.08406,
#     'M': 131.04049,
#     'N': 114.04293,
#     'P': 97.05276,
#     'Q': 128.05858,
#     'R': 156.10111,
#     'S': 87.03203,
#     'T': 101.04768,
#     'V': 99.06841,
#     'W': 186.07931,
#     'Y': 163.06333
# }

# weight = sum(monoisotopic_mass_table[aa] for aa in x)
# print("%0.3f" % weight)

######################################################################################################
# #Problem 21
# #Locating Restriction Sites

# from Bio import SeqIO
# seq = ""
# for rec in SeqIO.parse("./Rosalind_files/rosalind_revp.txt", "fasta"):
#     seq = str(rec.seq)

# complement = {"A":"T", "T":"A", "C":"G", "G":"C"}

# # Function to print all sub strings
# def subString(s, n):
#     # Pick starting point in outer loop
#     # and lengths of different strings for a given starting point
#     for i in range(n):
#         for length in range(i+1,n+1):
#             substring = s[i: length]
#             revcomp = "".join(complement[base] for base in reversed(substring))
#             if substring == revcomp:
#                 if 3 < len(substring) < 13:
#                     print(i+1, len(substring))

# subString(seq,len(seq))

######################################################################################################
# #Problem 22
# #RNA Splicing

# from Bio import SeqIO

# introns = []
# seq = ""
# sequences = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_splc.txt", "fasta"):
#     sequences.append(str(rec.seq))

# seq = sequences[0]
# introns = sequences[1:]   

# for x in introns:
#     result = seq.find(x)
#     intron = (seq[(result):(len(x) + result)])
#     seq = seq.replace(intron, "")

# spliced_rna = (seq.replace("T", "U"))

# codons = {    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
#               'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
#               'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
#               'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
#               'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
#               'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
#               'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
#               'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#               'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
#               'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
#               'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
#               'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
#               'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
#               'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
#               'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
#               'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

# def translate(rna):
#     protein = ""
#     start = rna.find("AUG")
#     for x in range(start, len(rna), 3):
#         codon = (rna[x:x+3])
#         if codon in codons:
#             if codons[codon] == "Stop":
#                 break
#             else:
#                 protein += codons[codon]
#     return protein

# print(translate(spliced_rna))

######################################################################################################
# # Problem 23
# # Enumerating k-mers Lexicographically

# from itertools import product

# with open ("./Rosalind_files/rosalind_lexf.txt") as f:
#     x = f.readlines()

# alphabet = x[0].strip()
# alphabet = alphabet.replace(" ", "")
# n = int(x[1])

# len_alpha = (len(alphabet))

# alphadict = {}

# for x,y in zip(range(len_alpha), alphabet):
#     alphadict[x] = y

# perms = product(alphadict.values(), repeat=n)

# for perm in perms:
#     print(''.join(perm))

######################################################################################################
##Problem 24
##Longest Increasing Subsequence

## I received quite a bit of help with this problem - it's still difficult for me to understand
# def longest_increasing_subsequence(sequence): 
#     length = len(sequence)
#     print(length)
#     # Initialized with 0
#     tail_indexes =[0 for _ in range(length + 1)]
#     # Initialized with -1
#     prev_indice =[-1 for _ in range(length + 1)]
#     # it will always point to empty location
#     length = 1
#     for i in range(1, len(sequence)):
#         if (sequence[i] < sequence[tail_indexes[0]]):
#             # new smallest value
#             tail_indexes[0] = i
#         elif (sequence[i] > sequence[tail_indexes[length-1]]):
#             # sequence[i] extends largest subsequence
#             prev_indice[i] = tail_indexes[length-1]
#             tail_indexes[length] = i
#             length += 1
#         else:
#             # sequence[i] is a candidate for future subsequence, replace cell value in tail_indexes
#             l = -1
#             r = length-1
#             while (r - l > 1):
#                 m = l + (r - l)//2
#                 if (sequence[tail_indexes[m]] >= sequence[i]):
#                     r = m
#                 else:
#                     l = m
#             pos = r
#             prev_indice[i] = tail_indexes[pos-1]
#             tail_indexes[pos] = i
#     i = tail_indexes[length-1]
#     result = []
#     while(i >= 0):
#         result.append(sequence[i])
#         i = prev_indice[i]
#     return result[::-1]  # reversing the result to get increasing order

# with open ("./Rosalind_files/rosalind_lgis.txt") as f:
#     x = f.readlines()

# sequence = [int(item) for line in x[1:] for item in line.split()]

# lis = longest_increasing_subsequence(sequence)
# print(" ".join(str(num) for num in lis))

# sequence = sequence[::-1]
# lds = longest_increasing_subsequence(sequence)
# print(" ".join(str(num) for num in lds[::-1]))
# #ignore the largest number that's a line above each subsequence, that's just the integer

######################################################################################################
# #Problem 25
# #Genome Assembly as Shortest Superstring

# # I received quite a bit of help with this problem, too. It's still difficult for me to understand

# import sys
# from Bio import SeqIO

# # Utility function to calculate minimum of two numbers
# def minimum(a, b):
#     return a if a < b else b

# # Function to calculate maximum overlap in two given strings
# def findOverlappingPair(str1, str2):
#     # Max will store maximum overlap i.e maximum length of the matching prefix and suffix
#     max_len = -sys.maxsize
#     len1 = len(str1)
#     len2 = len(str2)
#     str_ = ""

#     # Check suffix of str1 matches with prefix of str2
#     for i in range(1, minimum(len1, len2)+1):
#         # Compare last i characters in str1 with first i characters in str2
#         if str1[len1-i:] == str2[:i]:
#             if max_len < i:
#                 # Update max and str_
#                 max_len = i
#                 str_ = str1 + str2[i:]
#     # print(str_)

#     # Check prefix of str1 matches with suffix of str2
#     for i in range(1, minimum(len1, len2)+1):
#         # compare first i characters in str1 with last i characters in str2
#         if str1[:i] == str2[len2-i:]:
#             if max_len < i:
#                 # Update max and str_
#                 max_len = i
#                 str_ = str2 + str1[i:]
#     # print(str_)

#     return max_len, str_

# # Function to calculate smallest string that contains each string in the given set as substring.
# def findShortestSuperstring(arr, n):
#     # Run n-1 times to consider every pair
#     while n != 1:
#         # To store maximum overlap
#         max_len = -sys.maxsize
#         # To store array index of strings
#         l, r = 0, 0
#         # Involved in maximum overlap
#         res_str = ""

#         # Maximum overlap
#         for i in range(n):
#             for j in range(i+1, n):
#                 str_ = ""
#                 # res will store maximum length of the matching prefix and suffix str is passed by reference and will store the resultant
#                 # string after maximumoverlap of arr[i] and arr[j], if any.
#                 res, str_ = findOverlappingPair(arr[i], arr[j])

#                 # check for maximum overlap
#                 if max_len < res:
#                     max_len = res
#                     res_str = str_
#                     l, r = i, j

#         # Ignore last element in next cycle
#         n -= 1

#         # If no overlap, append arr[n-1] to arr[0]
#         if max_len == -sys.maxsize:
#             arr[0] += arr[n]
#         else:
#             # Copy resultant string to index l
#             arr[l] = res_str
#             # Copy string at last index to index r
#             arr[r] = arr[n]

#     return arr[0]

# reads = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_long.txt", "fasta"):
#     reads.append(str(rec.seq))

# n = len(reads)

# Function Call
# print(findShortestSuperstring(reads, n))