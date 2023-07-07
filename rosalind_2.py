######################################################################################################
#Problem 26
#Perfect Matchings and RNA Secondary Structures

# from math import factorial as fac

# def perfect_matchings(rna):
#    return fac(rna.count('A')) * fac(rna.count('C'))

# from Bio import SeqIO
# reads = []
# for rec in SeqIO.parse("c:/Users/Matt/downloads/rosalind_pmch.txt", "fasta"):
#     reads.append(str(rec.seq))

# for x in reads:
#     print(x)
#     print(perfect_matchings(x))

######################################################################################################
#Problem 27
#Partial Permutations

# from math import factorial as fac

# with open("c:/Users/Matt/downloads/rosalind_pper.txt") as f: 
#     protein = f.read()

# print(protein)
# n, k = protein.split()
# n, k = int(n), int(k)

# print(n, k)
# # n = 99
# # k = 9

# print(fac(n) / fac(n-k) % 1000000)

######################################################################################################
#Problem 28
#Introduction to Random Strings

# import math

# # s = "ACGATACAA"
# # # arr = [0.129]
# # arr = [0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]

# with open("c:/Users/Matt/downloads/rosalind_prob.txt") as f: 
#     string = f.readlines()

# s = string[0]
# arr = (string[1].split())

# def rand_string_calc(s, arr):

#     GC = (arr/2)
#     AT = ((1 - arr) / 2)

#     C = (s.count('C'))
#     G = (s.count('G'))
#     A = (s.count('A'))
#     T = (s.count('T'))

#     result = (math.log10(GC**(C+G)) + (math.log10(AT ** (A+T))))
#     print("%0.3f" % result, end=' ')

# for x in arr:
#     rand_string_calc(s, float(x))

######################################################################################################
#Problem 29
#Enumerating Oriented Gene Orderings 

# from itertools import permutations, product

# def signedperms(items):
#     for p in permutations(items):
#         # print(f"p is:", p)
#         for signs in product([-1,1], repeat=len(items)):
#             # print(f"signs are:", signs)        
#             yield [a*sign for a,sign in zip(p,signs)]

# n = 4

# result = list(signedperms(range(1, n+1)))

# print(len(result))
# for item in result:
#     print(' '.join(map(str, item)))

######################################################################################################
#Problem 30
#Finding a Spliced Motif

# string = "ACGTCACGTGACG"
# sub = "GTA"

# from Bio import SeqIO

# info = []
# for rec in SeqIO.parse("c:/Users/Matt/downloads/rosalind_sseq.txt", "fasta"):
#     info.append(str(rec.seq))

# string = info[0]
# sub = info[1]

# seq = []
# skip = [0]
# for x in sub:
#     index = string.find(x) + 1
#     seq.append(index + sum(skip))
#     skip.append(index)
#     string = string[index:]
   
# for x in seq:
#     print(x, end= " ")

######################################################################################################
#Problem 31
#Finding a Spliced Motif

# from Bio import SeqIO

# info = []
# for rec in SeqIO.parse("c:/Users/Matt/downloads/rosalind_tran.txt", "fasta"):
#     info.append(str(rec.seq))

# str1 = info[0]
# str2 = info[1]

# str1 = "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT"
# str2 = "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"

# transitions, transversions = 0, 0

# for x in range(len(str1)):
#     if str1[x] == str2[x]:
#         pass
#     elif str1[x] == "A" and str2[x] == "G":
#         transitions += 1
#     elif str1[x] == "C" and str2[x] == "T":
#         transitions += 1
#     elif str1[x] == "T" and str2[x] == "C":
#         transitions += 1
#     elif str1[x] == "G" and str2[x] == "A":
#         transitions += 1
#     else:
#         transversions += 1
    
# print(transitions, transversions)
# print(transitions/transversions)

######################
#better way

# str1 = "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT"
# str2 = "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"

# # Map of transitions
# transition_map = {"AG", "GA", "CT", "TC"}

# transitions, transversions = 0, 0

# for s1, s2 in zip(str1, str2):
#     if s1 != s2:
#         if f"{s1}{s2}" in transition_map:
#             transitions += 1
#         else:
#             transversions += 1

# print(transitions, transversions)
# print(transitions/transversions)

######################################################################################################
#Problem 32
#Completing a Tree

#this is pseudo code that doesn't use an algorithm. 

# with open('c:/Users/Matt/downloads/rosalind_tree(1).txt','r') as f:
#     lines = f.readlines()

# print(int(lines[0]) - len(lines[1:]) - 1)

######################################################################################################
#Problem 33
#Catalan Numbers and RNA Secondary Structures

#catalan number formula, psuedo code

# from Bio import SeqIO

# rna = ""
# for rec in SeqIO.parse("c:/Users/Matt/downloads/rosalind_cat(3).txt", "fasta"):
#     rna = (str(rec.seq))

# def solve(rna):
#     """
#     An input RNA consisting of {A, U, C, G}
#     The number of non-overlapping perfect 
#     matchings.
#     """
#     return helper(rna, 0, len(rna) - 1, {})


# def helper(rna, lo, hi, dp):
    
#     mapping = {
#     "A" : "U",
#     "U" : "A",
#     "G" : "C",
#     "C" : "G"
#     }
#     characters = hi - lo + 1
    
#     # if there are an odd number of nucleotides, 
#     # this is an invalid matching.
#     if characters % 2 == 1:
#         return 0

#     # handles tricky edge cases.
#     if lo >= hi or lo >= len(rna) or hi < 0:
#         return 1

#     # return answer if it is memoized.    
#     if (lo, hi) in dp:
#         return dp[(lo, hi)]
#     else:
#         curr = rna[lo]
#         target = mapping[curr]
#         acc = 0
#         for i in range(lo + 1, hi + 1, 2):
#             if rna[i] == target:
#                 left = helper(rna, lo + 1, i - 1, dp)
#                 right = helper(rna, i + 1, hi, dp)
#                 acc += (left * right) % 1000000
#         dp[(lo, hi)] = acc
#         return acc

# print(solve(rna.strip()) % 1000000)

######################################################################################################
#Problem 34
#Catalan Numbers and RNA Secondary Structures

from Bio import SeqIO
import itertools

def revcomp(s):
    return s.translate(str.maketrans('ACTG','TGAC'))[::-1]

def compare_matches(s1, s2):
    return sum(s1[i] != s2[i] for i in range(len(s1))) == 1

rna = []
revrna = []
for rec in SeqIO.parse("c:/Users/Matt/downloads/rosalind_cat.txt", "fasta"):
    rna.append(str(rec.seq))
    revrna.append(revcomp(str(rec.seq)))

matches=[]
revmatches=[]
hammings = []

for a, b in itertools.combinations(rna, 2):
    if a == revcomp(b):
        revmatch = (a,b)
        revmatches.append(revmatch)
    if a == b:
        match = (a)
        matches.append(match)

print(f"revmatches", revmatches) 

for x in revmatches:
    for y in x:
        for a in rna:
            if compare_matches(y, a) == True:
                hammings_match = (y,a)
                hammings.append(hammings_match)

# print(f"hammings_match", hammings)
# print(f"revmatches", revmatches)        
# print(f"matches", matches)

correct = []
final = []

for x in revmatches:
    if y in x:
        if y in hammings_match:
            correct.append(y)

for x in correct:
    if x in hammings_match:
        # Check if x is the first element
        if hammings_match[0] == x:
            # Reverse the tuple
            reversed_hammings_match = hammings_match[::-1]
        else:
            # Keep the tuple as it is
            reversed_hammings_match = hammings_match

        final.append("->".join(reversed_hammings_match))

for x in matches:
    x = revcomp(x)
    for a in rna:
        if compare_matches(x,a) == True:
            temp_match = (a,x)
            final.append("->".join(temp_match))

for x in revmatches:
    for y in x:
        for a in rna:
            if compare_matches(y,a) == True:
                temp_match = (a,y)
                final.append("->".join(temp_match))
                
for x in set(final):
    print(x)

#############################################################
# #corrected via chatgpt:
# from Bio import SeqIO
# import itertools

# rna = []
# for rec in SeqIO.parse("c:/Users/Matt/downloads/rosalind_cat.txt", "fasta"):
#     rna.append(str(rec.seq))

# matches = []
# revmatches = []
# hammings = {}

# def revcomp(s):
#     return s.translate(str.maketrans('ACTG','TGAC'))[::-1]

# def compare_matches(s1, s2):
#     return sum(s1[i] != s2[i] for i in range(len(s1))) == 1

# for a, b in itertools.combinations(rna, 2):
#     if a == revcomp(b):
#         revmatch = (a, b)
#         revmatches.append(revmatch)
#     elif a == b:
#         match = (a,)
#         matches.append(match)

# for x in revmatches:
#     for y in x:
#         for a in rna:
#             if compare_matches(y, a):
#                 if y in hammings:
#                     hammings[y].append(a)
#                 else:
#                     hammings[y] = [a]

# corrections = []

# for x in revmatches:
#     if x[0] in hammings and x[1] in hammings[x[0]]:
#         correction = f"{x[0]}->{x[1]}"
#         corrections.append(correction)

# for x in matches:
#     x = revcomp(x[0])
#     if x in hammings:
#         for a in hammings[x]:
#             correction = f"{x}->{a}"
#             corrections.append(correction)

# for correction in set(corrections):
#     print(correction)




