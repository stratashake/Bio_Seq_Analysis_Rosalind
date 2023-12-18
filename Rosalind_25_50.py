######################################################################################################
# #Problem 26
# #Perfect Matchings and RNA Secondary Structures

# from math import factorial as fac
# from Bio import SeqIO
# reads = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_pmch.txt", "fasta"):
#     reads.append(str(rec.seq))

# def perfect_matchings(rna):
#    return fac(rna.count('A')) * fac(rna.count('C'))

# for x in reads:
#     print(perfect_matchings(x))

######################################################################################################
# #Problem 27
# #Partial Permutations

# # Another way I've imported variables from files
# # with open("./Rosalind_files/rosalind_pper.txt") as f: 
# #     protein = f.read()
# # n, k = protein.split()
# # n, k = int(n), int(k)

# from math import factorial as fac

# with open("./Rosalind_files/rosalind_pper.txt") as f:
#     line = f.readline()
#     n, k = map(int, line.split())

# # Calculate partial permutation using modulo at the end
# ans = int(fac(n) / fac(n-k)) % 1000000
# print(ans)

######################################################################################################
# #Problem 28
# #Introduction to Random Strings

# import math

# with open("./Rosalind_files/rosalind_prob.txt") as f: 
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
## Problem 29
## Enumerating Oriented Gene Orderings 

# from itertools import permutations, product

# def signedperms(items):
#     for p in permutations(items):
#         for signs in product([-1,1], repeat=len(items)):       
#             yield [a*sign for a,sign in zip(p,signs)]

# with open("./Rosalind_files/rosalind_sign.txt") as f:
#     n = int(f.read().strip())

# result = list(signedperms(range(1, n+1)))

# print(len(result))
# for item in result:
#     print(' '.join(map(str, item)))

######################################################################################################
# #Problem 30
# #Finding a Spliced Motif

# from Bio import SeqIO

# info = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_sseq.txt", "fasta"):
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
# #Problem 31
# #Finding a Spliced Motif

# from Bio import SeqIO

# info = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_tran.txt", "fasta"):
#     info.append(str(rec.seq))

# str1 = info[0]
# str2 = info[1]

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
    
# # print(transitions, transversions)
# print(transitions/transversions)

# #####################
# #better way

# # Map of transitions
# transition_map = {"AG", "GA", "CT", "TC"}

# transitions, transversions = 0, 0

# for s1, s2 in zip(str1, str2):
#     if s1 != s2:
#         if f"{s1}{s2}" in transition_map:
#             transitions += 1
#         else:
#             transversions += 1

# # print(transitions, transversions)
# print(transitions/transversions)

######################################################################################################
# #Problem 32
# #Completing a Tree

# #this is pseudo code that doesn't use an algorithm. 

# with open("./Rosalind_files/rosalind_tree.txt",'r') as f:
#     lines = f.readlines()

# print(int(lines[0]) - len(lines[1:]) - 1)

######################################################################################################
##Problem 33
##Catalan Numbers and RNA Secondary Structures

# from Bio import SeqIO

# rna = ""
# for rec in SeqIO.parse("./Rosalind_files/rosalind_cat.txt", "fasta"):
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
## Problem 34
## Error Correction in Reads

# """
# #This was my original best attempt, leaving here just for my benefit

# from Bio import SeqIO
# import itertools

# def revcomp(s):
#     return s.translate(str.maketrans('ACTG','TGAC'))[::-1]

# def compare_matches(s1, s2):
#     return sum(s1[i] != s2[i] for i in range(len(s1))) == 1

# rna = []
# revrna = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_corr.txt", "fasta"):
#     rna.append(str(rec.seq))
#     revrna.append(revcomp(str(rec.seq)))

# matches=[]
# revmatches=[]
# hammings = []

# for a, b in itertools.combinations(rna, 2):
#     if a == b:
#         match = (a)
#         matches.append(match)

# # Convert both lists to sets
# rna_set = set(rna)
# revrna_set = set(revrna)

# # Find the intersection of the two sets (i.e., the sequences that are in both)
# revmatches = rna_set & revrna_set

# # If you need the matches as a list, you can convert the set back to a list
# revmatches = list(revmatches)

# for x in revmatches:
#     for a in rna:
#         if compare_matches(x, a) == True:
#             hammings_match = (x,a)
#             hammings.append(hammings_match)

# correct = []
# final = []

# for x in revmatches:
#     if x in hammings_match:
#         correct.append(x)

# for x in correct:
#     if x in hammings_match:
#         # Check if x is the first element
#         if hammings_match[0] == x:
#             # Reverse the tuple
#             reversed_hammings_match = hammings_match[::-1]
#         else:
#             # Keep the tuple as it is
#             reversed_hammings_match = hammings_match

#         final.append("->".join(reversed_hammings_match))

# for x in matches:
#     x = revcomp(x)
#     for a in rna:
#         if compare_matches(x,a) == True:
#             temp_match = (a,x)
#             final.append("->".join(temp_match))

# for y in revmatches:
#     for a in rna:
#         if compare_matches(y,a) == True:
#             temp_match = (a,y)
#             final.append("->".join(temp_match))
                
# for x in set(final):
#     print(x)
# """

# # Much Corrected code

# from Bio import SeqIO
# import itertools

# def revcomp(s):
#     return s.translate(str.maketrans('ACTG','TGAC'))[::-1]

# def compare_matches(s1, s2):
#     return sum(s1[i] != s2[i] for i in range(len(s1))) == 1

# rna = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_corr.txt", "fasta"):
#     rna.append(str(rec.seq))

# matches=[]
# hammings = []

# for a, b in itertools.combinations(rna, 2):
#     if a == b:
#         matches.append(a)

# for x in matches:
#     for a in rna:
#         if compare_matches(x, a):
#             hammings.append((x, a))

# final = []

# for x in hammings:
#     # Check if x[0] is the first element
#     if x[0] == x[0]:
#         # Reverse the tuple
#         reversed_hammings_match = x[::-1]
#     else:
#         # Keep the tuple as it is
#         reversed_hammings_match = x

#     final.append("->".join(reversed_hammings_match))

# for x in matches:
#     x = revcomp(x)
#     for a in rna:
#         if compare_matches(x,a):
#             final.append("->".join((a, x)))

# for x in set(final):
#     print(x)


######################################################################################################
# #Problem 35
# #Counting Phylogenetic Ancestors

# # More psuedocode

# with open("./Rosalind_files/rosalind_inod.txt") as f:
#     n = int(f.read().strip())

# print(n-2)

######################################################################################################
# #Problem 36
# #k-Mer Composition

# from Bio import SeqIO
# from collections import defaultdict

# s = ""
# for rec in SeqIO.parse("./Rosalind_files/rosalind_kmer.txt", "fasta"):
#     s = (rec.seq)

# sdict = defaultdict(int) #makes all values 0

# import itertools
# x = ["A", "C", "G", "T"]
# for y in itertools.product(x, repeat=4): #cartisean product is what we want
#     sdict[("".join(y))] #allows us to create mutable and connected keys like AAAC instead of "A", "A", "A", "C"
    
# for x in range(len(s)):
#     kmer = (s[x:x+4]) #go back 4 letter increments
#     if len(kmer) == 4:
#         sdict[kmer] += 1

# print(' '.join(map(str, sdict.values()))) #print with no commas or brackets, just spaces 

######################################################################################################
# #Problem 37
# #Speeding Up Motif Finding

from Bio import SeqIO

<<<<<<< HEAD
# for rec in SeqIO.parse("./Rosalind_files/rosalind_kmp.txt", "fasta"):
#     s = (rec.seq)

# y = len(s)
# array = [0 for _ in range(y)]

# j = 0
# for i in range(1, y):
#     while j > 0 and s[i] != s[j]:
#         j = array[j-1]
#     if s[i] == s[j]:
#         j += 1
#         array[i] = j
# print(' '.join(map(str, array))) 
=======
for rec in SeqIO.parse("./Rosalind_files/rosalind_kmp.txt", "fasta"):
    s = (rec.seq)

y = len(s)
array = [0 for _ in range(y)]

j = 0
for i in range(1, y):
    while j > 0 and s[i] != s[j]:
        j = array[j-1]
    if s[i] == s[j]:
        j += 1
        array[i] = j

print(' '.join(map(str, array)))
>>>>>>> 8bb2c93 (Updated up to probelm 37 with all .txt files)
