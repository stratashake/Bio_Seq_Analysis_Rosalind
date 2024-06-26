######################################################################################################
#Problem 26
#Perfect Matchings and RNA Secondary Structures

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
# #Transitions and Transversions 

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

# from Bio import SeqIO

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

######################################################################################################
# #Problem 38
# #Finding a Shared Spliced Motif

# from Bio import SeqIO

# fastas = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_lcsq.txt", "fasta"):
#     str_seq = str(rec.seq)
#     fastas.append(str_seq)  

# def lcs(X, Y):
#     """Function is an iterative based longest common subsequence algorithm using dynamic programming"""
#     m = len(X)
#     n = len(Y)
#     L = [[0] * (n + 1) for _ in range(m + 1)] #creates the empty 2D array, add 1 because arrays start at 0
 
#     for i in range(1, m + 1): #iterates over the "rows"
#         for j in range(1, n + 1): #iterates over the "columns"
#             if X[i - 1] == Y[j - 1]:
#                 L[i][j] = L[i - 1][j - 1] + 1 #adds 1 to the cell in the array where the subsequence matches 
#                                               #also, this line is utilized if the last character of X and Y match
#             else:
#                 L[i][j] = max(L[i - 1][j], L[i][j - 1]) #If the last character of the two strings don't match
#                                                         # we find which direction to go in to the find the highest value 
#                                                         # relative to it being the x or y string 
 
#     lcs = ""    # Initialize an empty string to hold the longest common subsequence
#     i = m       # Set i to the length of X, which is the number of rows in L
#     j = n       # Set j to the length of Y, which is the number of columns in L
#     while i > 0 and j > 0:  # While we have not reached the beginning of either string
#         if X[i - 1] == Y[j - 1]: #if we don't subtract one, we get out of range error because python starts indexing at 0
#             lcs = X[i-1] + lcs  # Add the current character to the beginning of the lcs string each time, which "reverses" the 2D array
#                                 #with regard to how the dynamic programing is searching it, which is bottom right to top left
#             i -= 1
#             j -= 1
#         elif L[i][j - 1] > L[i - 1][j]: #these lines move us from the bottom right to the 
#             j -= 1                      #upper left where we get to the 0 corner of the array
#         else:                           #the direction we move depends on which string is longer at the moment, so either up or left
#             i -= 1
#     return lcs


# output = (lcs(*fastas))
# print(output)

######################################################################################################
# # #Problem 39
# # #Ordering Strings of Varying Length Lexicographically 

# from itertools import product

# with open("./Rosalind_files/rosalind_lexv.txt") as f:
#     lines = f.readlines()


# alpha = (lines[0].strip())
# alpha = "".join(alpha.split())
# num = int(lines[1].strip())

# lst = []

# for x in range((num)):
#     lst.append((list(product(alpha, repeat = x+1))))

# flatlst = [item for sublist in lst for item in sublist]

# alphadict = {}
# for i, alph in enumerate(alpha):
#     alphadict[alph] = int(i)

# # Function to get the sorting key for each product
# def sort_key(product):
#     # Convert the product tuple into a list of ranks based on the custom alphabet
#     return [alphadict.get(letter, float('inf')) for letter in product]

# # Sort the list of products
# flatlst.sort(key=sort_key)

# for x in flatlst:
#     x = ''.join(map(str, x))
#     print(x)

# # I suggest using this command to view the results instead of printing to the terminal: "python .\Rosalind_25_50.py > output.txt"

######################################################################################################
# #Problem 40
# #Maximum Matchings and RNA Secondary Structures

# from Bio import SeqIO
# from math import factorial

# for rec in SeqIO.parse("./Rosalind_files/rosalind_mmch.txt", "fasta"):
#     rna = str(rec.seq)

# A = rna.count('A')
# C = rna.count('C')
# U = rna.count('U')
# G = rna.count('G')

# #Have to use '//' instead of '/' to remove memory inprecision that can arise with large number
# #calculations in Python
# AU = factorial(max(A, U)) // factorial(max(A, U) - min(A, U))
# GC = factorial(max(G, C)) // factorial(max(G, C) - min(G, C))
# print(int(AU * GC))

######################################################################################################
# #Problem 41
# #Creating a Distance Matrix

# from Bio import SeqIO

# fastas = []
# for rec in SeqIO.parse("./Rosalind_files/rosalind_pdst.txt", "fasta"):
#     str_seq = str(rec.seq)
#     fastas.append(str_seq) 

# # Initialize a matrix to store the distances
# distance_matrix = [[0 for x in range(len(fastas))] for x in range(len(fastas))]

# for i, fasta1 in enumerate(fastas):
#     for j, fasta2 in enumerate(fastas):
#         if i != j:
#             mismatches = sum(1 for a, b in zip(fasta1, fasta2) if a != b)
#             distance = mismatches / len(fasta1)
#             distance_matrix[i][j] = distance

# for row in distance_matrix:
#     print(" ".join(f"{dist:.5f}" for dist in row))

######################################################################################################
# #Problem 42
# #Reversal Distance

# #This was as far as I could get on my own, it could find a reversal distance, but not the minimum. 
# with open("./Rosalind_files/rosalind_rear.txt") as f:
# 	lines = [line.strip().split() for line in f.readlines() if line.strip()]

# def brute_reversals(s1, s2, reversals):
#     # Check for consecutive matching elements from the beginning
#     match_from_start = 0
#     while match_from_start < len(s1) and match_from_start < len(s2) and s1[match_from_start] == s2[match_from_start]:
#         match_from_start += 1

#     # Check for consecutive matching elements from the end
#     match_from_end = 0
#     while match_from_end < len(s1) and match_from_end < len(s2) and s1[-(match_from_end + 1)] == s2[-(match_from_end + 1)]:
#         match_from_end += 1

#     # Remove the matching elements from both ends
#     if match_from_start > 0 or match_from_end > 0:
#         s1 = s1[match_from_start:len(s1) - match_from_end]
#         s2 = s2[match_from_start:len(s2) - match_from_end]

#     # Special case: s1 and s2 are in reverse order
#     if s1 == s2[::-1]:
#         s2.reverse()
#         reversals += 1
#         return s1, s2, reversals

#     # Handling other cases
#     if s1[-1] == s2[-1] and s1[0] == s2[0]:
#         s1 = s1[1:-1]
#         s2 = s2[1:-1]

#     elif s1[-1] != s2[-1] and s1[0] == s2[0]:
#         end_index = s2.index(s1[-1])
#         s2[end_index:] = s2[end_index:][::-1]
#         reversals += 1

#     elif s1[-1] == s2[-1] and s1[0] != s2[0]:
#         begin_index = s2.index(s1[0])
#         s2[:begin_index + 1] = s2[:begin_index + 1][::-1]
#         reversals += 1

#     elif s1[-1] != s2[-1] and s1[0] != s2[0]:
#         end_index = s2.index(s1[-1])
#         begin_index = s2.index(s1[0])
#         s2[end_index:] = s2[end_index:][::-1]
#         s2[:begin_index + 1] = s2[:begin_index + 1][::-1]
#         reversals += 2

#     return s1, s2, reversals

# for i in range(0, len(lines), 2):
#     s1 = list(map(int, lines[i]))
#     s2 = list(map(int, lines[i + 1]))
#     total_reversals = 0

#     while s1 and s2 and s1 != s2:
#         s1_before, s2_before = s1[:], s2[:]  # Make a copy of s1 and s2 before the operation
#         s1, s2, total_reversals = brute_reversals(s1, s2, total_reversals)

#         # Check if no progress is made
#         if s1 == s1_before and s2 == s2_before:
#             break

#     print(total_reversals)

# """This was a script I was able to create with quite a bit of help with chatgpt. I learned more about 
# data structure and algorithms by learning how this code works and feel more confident in implementing solutions like 
# this in the future. This is a Bidirectional Breadth-First Search (BFS) algorithm.
# The algorithm efficiently explores the space of permutations by generating all 
# possible reversals of each permutation from both ends of the search. It progressively works towards finding a common 
# permutation that can be reached from both the original and the target permutations, effectively meeting in the middle.
# """ 

# def generate_reversed_permutations(perm):
#     """Generate all permutations of 'perm' by reversing its contiguous subsequences."""
#     for i in range(len(perm)):
#         for j in range(i+2, len(perm)+1):
#             yield perm[:i] + perm[i:j][::-1] + perm[j:]

# def find_min_reversal_distance(s1, s2):
#     """Find the minimum reversal distance between two permutations."""
#     if s1 == s2:
#         return 0

#     visited = {tuple(s1): 0, tuple(s2): 0}
#     queue = [(tuple(s1), 1), (tuple(s2), -1)]  # Positive for s1, negative for s2

#     while queue:
#         current, direction = queue.pop(0) # Explore the queue in the order it's made
#         current_distance = visited[current]
#         next_distance = current_distance + (1 if direction > 0 else -1)

#         for next_perm in generate_reversed_permutations(current):
#             next_perm = tuple(next_perm)
#             if next_perm in visited:
#                 if visited[next_perm] * direction < 0:  # Opposite direction found, only want opposite direction
#                     # Add 1 to account for the final connecting reversal
#                     return abs(visited[next_perm]) + abs(current_distance) + 1
#                 continue

#             visited[next_perm] = next_distance # Appends each permuation with corresponding number of reversions to that point, if not already apended
#             queue.append((next_perm, direction))

#     return ("No solution found")

# # Read data and calculate reversal distances
# with open("./Rosalind_files/rosalind_rear.txt", "r") as f:
#     lines = [list(map(int, line.strip().split())) for line in f.readlines() if line.strip()]
#     for i in range(0, len(lines), 2):
#         distance = find_min_reversal_distance(lines[i], lines[i + 1])
#         print(distance, end=' ') # Gives correct formated output

######################################################################################################
# #Problem 43
# #Matching Random Motifs
   
# with open("./Rosalind_files/rosalind_rstr.txt") as f:
#     lines = [line for line in f]
#     n, x = lines[0].split()
#     n = int(n)
#     x = float(x)
#     s = lines[1].strip() 

# def rand_string_calc(x, s, n):

#     GC = (x/2)
#     AT = ((1 - x) / 2)

#     C = (s.count('C'))
#     G = (s.count('G'))
#     A = (s.count('A'))
#     T = (s.count('T'))

#     odds = (GC ** G) * (GC ** C) * (AT ** A) * (AT ** T)
#     odds_appear = 1 - ((1 - odds) ** n)
#     return(odds_appear)

# print("%0.3f" % rand_string_calc(x, s, n))

######################################################################################################
# #Problem 44
# #Counting Subsets

# with open("./Rosalind_files/rosalind_sset.txt") as f:
#     n = int(f.read().strip())
    
# print((2**n % 1000000))

## This is an unwieldly and bad attempt that uses too much memory, but it works
# # from itertools import combinations

# # n = 10
# # numbers = list(range(1,n+1))
# # list = []

# # for x in range(1,n+1):
# #     for p in combinations(numbers, x):
# #         list.append(p)

# # answer = (len(list)) % 1000000 + 1
# # print(answer)

######################################################################################################
# #Problem 45
# #Introduction to Alternative Splicing

# with open("./Rosalind_files/rosalind_aspc.txt") as f:
#     line = f.readline()
#     n, m = map(int, line.split())

# from math import factorial as fac
# sum = 0

# for k in range(m, n + 1):
#     sum += (fac(n))//(fac(k)*(fac(n-k)))

# print(sum % 1000000)

######################################################################################################
# #Problem 46
# #Edit Distance

# from Bio import SeqIO

# sequences = [str(rec.seq) for rec in SeqIO.parse("./Rosalind_files/rosalind_edit.txt", "fasta")]

# s, t = sequences[0], sequences[1]

# def edit_distance(s, t): 
#     dp = [[0] * (len(t)+1) for _ in range(len(s)+1)]

#     for i in range(len(s)+1):
#         for j in range(len(t)+1):
#             if i == 0:
#                 dp[i][j] = j
#             elif j == 0:
#                 dp[i][j] = i
#             elif s[i-1] == t[j-1]:
#                 dp[i][j] = dp[i-1][j-1]
#             else:
#                 dp[i][j] = 1 + min(dp[i - 1][j],    # Insertion
#                                    dp[i][j - 1],    # Deletion
#                                    dp[i - 1][j - 1]) # Substitution

#     return dp[len(s)][len(t)]

# print(edit_distance(s,t))

######################################################################################################
# #Problem 47
# #Expected Number of Restriction Sites

# with open("./Rosalind_files/rosalind_eval.txt") as f: 
#     string = f.readlines()

# n = int(string[0])
# s = string[1]
# length_s = len(s)
# array = [float(x) for x in string[2].split()]

# answer = []

# for x in array:
#     G = x/2
#     C = x/2
#     A = (1-x)/2
#     T = (1-x)/2

#     prob = 1
#     for char in s:
#         if char == "G":
#             prob *= G
#         elif char == 'C':
#             prob *= C
#         elif char == 'A':
#             prob *= A
#         elif char == 'T':
#             prob *= T

#     combos = n - length_s + 1 #number of possible combinations according to the length of "s"
#     chance = (prob * combos) #additive property of expected number
#     answer.append(round(chance, 3)) #round to 3 decimals

# print(' '.join(map(str, answer)))

######################################################################################################
# #Problem 48
# #Motzkin Numbers and RNA Secondary Structures

# """This is entirely over my head and I do not understand most of this. I copied this from someone
# else's codebase. I'll come back to it later but dynamic programming and this kind of combinatorics 
# is over my head"""

# from Bio import SeqIO

# # Reading the sequence from a FASTA file
# for rec in SeqIO.parse("./Rosalind_files/rosalind_motz.txt", "fasta"):
#     s = str(rec.seq)

# # Initializing the dictionary for memoization with base cases
# c = {
#     '': 1, 'A': 1, 'C': 1, 'G': 1, 'U': 1,
#     'AA': 1, 'AC': 1, 'AG': 1, 'AU': 2,
#     'CA': 1, 'CC': 1, 'CG': 2, 'CU': 1,
#     'GA': 1, 'GC': 2, 'GG': 1, 'GU': 1,
#     'UA': 2, 'UC': 1, 'UG': 1, 'UU': 1
# }

# def motzkin(s):
#     if s not in c:
#         # Initialize the sum for the recursive formula
#         total = motzkin(s[1:])
#         # Replace the list comprehension with a loop
#         for k in range(1, len(s)):
#             factor = c[s[0]+s[k]] - 1
#             result = motzkin(s[1:k]) * factor * motzkin(s[k+1:])
#             total += result
#         c[s] = total
#     return c[s]

# # Print the result modulo 10**6
# print(motzkin(s) % 10**6)

######################################################################################################
# #Problem 49
# #Distances in Trees

# from Bio import Phylo
# from io import StringIO

# with open("./Rosalind_files/rosalind_nwck.txt") as f: 
#     trees = f.read().splitlines() 

# phylo = [] #the newick tree
# animals = [] #the two nodes of interest
# ans_list = []

# for x in range(0,len(trees),3):
#     phylo.append(trees[x])
#     animals.append(trees[x+1].split(" "))


# def newick_distance_finder(phy, ani):
#     tree = Phylo.read(StringIO(phy), "newick")
#     clades = tree.find_clades()
#     for clade in clades: #biopython was reading each branch/clade with a length of 0, originally. 
#         clade.branch_length = 1
#     ans_list.append(tree.distance(ani[0], ani[1]))

# for x,y in zip(phylo, animals):
#     newick_distance_finder(x, y)

# print(' '.join(map(str, ans_list)))

######################################################################################################
# #Problem 50
# #Interleaving Two Motifs

# with open("./Rosalind_files/rosalind_scsp.txt") as f: 
#     string = f.readlines()

# def scs_length(X, Y):
#     m, n = len(X), len(Y)
#     L = [[0] * (n + 1) for _ in range(m + 1)]

#     # Initialize base cases
#     for i in range(m+1):
#         L[i][0] = i
#     for j in range(n+1):
#         L[0][j] = j

#     # Fill the DP table
#     for i in range(1, m+1):
#         for j in range(1, n+1):
#             if X[i-1] == Y[j-1]:
#                 L[i][j] = L[i-1][j-1] + 1
#             else:
#                 L[i][j] = min(L[i-1][j], L[i][j-1]) + 1
#     return L

# def reconstruct_scs(L, X, Y):
#     scs = ""
#     m, n = len(X), len(Y)

#     while m > 0 and n > 0:
#         if X[m-1] == Y[n-1]:
#             scs = X[m-1] + scs
#             m -= 1
#             n -= 1
#         elif L[m-1][n] <= L[m][n-1]:
#             scs = X[m-1] + scs
#             m -= 1
#         else:
#             scs = Y[n-1] + scs
#             n -= 1

#     # Add any remaining characters from X or Y
#     while m > 0:
#         scs = X[m-1] + scs
#         m -= 1
#     while n > 0:
#         scs = Y[n-1] + scs
#         n -= 1

#     return scs

# s = string[0].strip()
# t = string[1].strip()

# table = scs_length(s, t)
# scs = reconstruct_scs(table, s, t)
# print(scs)
