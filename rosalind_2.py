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

from itertools import permutations, product

def signedperms(items):
    for p in permutations(items):
        for signs in product([-1,1], repeat=len(items)):
            yield [a*sign for a,sign in zip(p,signs)]

n = 4

result = list(signedperms(range(1, n+1)))

print(len(result))
for item in result:
    print(' '.join(map(str, item)))

