######################################################################################################
# #Problem 51
# #Introduction to Set Operations

# import ast  # Importing ast which will help safely evaluate string literals as Python expressions
# import os

# # Open the file
# with open("./Rosalind_files/rosalind_seto.txt", 'r') as file:
#     # Read the first line for the universal set size
#     n = int(file.readline().strip())
#     u_set = set(range(1, n + 1))  # Create the universal set from 1 to n

#     # Read the second line for the first set
#     a = file.readline().strip()
#     a = ast.literal_eval(a)  # Convert string representation of the set to a Python set
    
#     # Read the third line for the second set
#     b = file.readline().strip()
#     b = ast.literal_eval(b)  # Convert string representation of the set to a Python set

# # u_set = set(range(1, 11))
# # a = {1, 2, 3, 4, 5}
# # b = {2, 8, 5, 10}

# set_list =union = a | b
# intersection = a & b
# diff_a = a.difference(b)
# diff_b = b.difference(a)
# comp_a = u_set.difference(a)
# comp_b = u_set.difference(b)

# folder_path = './Rosalind_files/Seto_output'
# # Create the folder if it doesn't already exist
# os.makedirs(folder_path, exist_ok=True)
# # Define the full path for the new file within the new subfolder
# file_path = os.path.join(folder_path, 'output_seto.txt')

# with open(file_path, "w") as f:
#     print("\n".join(map(str, (union, intersection, diff_a, diff_b, u_set.difference(a), u_set.difference(b)))), file=f)

######################################################################################################
# #Problem 52
# #Sorting by Reversals

"""This is not my code. I modified it from another repository"""

# import random

# # get the all reversal array of a list "s".
# def _get_reverse_array(s):
#     reverse_arrays = []
#     for i in range(len(s)-1):
#         for j in range(i+1, len(s)):
#             r_list = s[i:j+1]
#             r_list.reverse()
#             reverse_arrays.append(s[:i] + r_list + s[j+1:])
#     return reverse_arrays

# # get the reversal_distance from a list "s1" to another list "s2".
# def _get_reversal_distance(s1, s2, distance, s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals):
#     # reverse s1 to s2, and reverse s2 to s1 at same time.
#     if s1 & s2:
#         return s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals, distance
#     # get the reveral array of s1.
#     new_s1 = set()
#     s1_s2 = {}
#     for s in s1:
#         reverse_arrays = _get_reverse_array(list(s))
#         s1_s2[s] = reverse_arrays
#         for r in reverse_arrays:
#             new_s1.add(tuple(r))
#     s1_s2_reversals_path.append(s1_s2)
#     # get the reveral array of s2.
#     new_s2 = set()
#     s2_s1 = {}
#     for s in s2:
#         reverse_arrays = _get_reverse_array(list(s))
#         s2_s1[s] = reverse_arrays
#         for r in reverse_arrays:
#             new_s2.add(tuple(r))
#     s2_s1_reversals_path.append(s2_s1)
#     # thus we reverse s1 and s2 at same time, so distance plus 2. 
#     distance += 2
#     # if s1 and the reversal array of s2 has the same array, distance substract 1.
#     if s1 & new_s2:
#         meet_reversals = list(s1 & new_s2)
#         return s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals, distance-1
#     # if s2 and the reversal array of s1 has the same array, distance substract 1.
#     if s2 & new_s1:
#         meet_reversals = list(s2 & new_s1)
#         return s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals, distance-1
#     # if reversal array of s1 and the reversal array of s2 has the same array, return distance.
#     if new_s1 & new_s2:
#         meet_reversals = list(new_s1 & new_s2)
#         return s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals, distance
    
#     s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals, distance = _get_reversal_distance(new_s1, new_s2, distance, s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals)
#     return s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals, distance

# # get the endpoints of the interval of the two indices.
# def _get_invert_endpoints(a, b):
#     a_reverse = []
#     for i in range(len(a)-1):
#         for j in range(i+1, len(a)):
#             a_reverse = a[:i] + a[i:j+1][::-1] + a[j+1:]
#             if a_reverse == b:
#                 return i+1, j+1

# # get the collections of reversals sorting π into γ
# def _get_collections_of_reversals_sorting(s1_s2_reversals_path, meet_reversals, s2_s1_reversals_path):
#     collections_of_reversals_sorting = []
#     for l in meet_reversals:
#         # print(l)
#         collection_of_reversals_sorting=[]
#         current_array = list(l)
#         for reversals_path in s1_s2_reversals_path[::-1]:
#             for k, v in reversals_path.items():
#                 if current_array in v:
#                     i, j =  _get_invert_endpoints(list(k), current_array)
#                     collection_of_reversals_sorting.append([i, j])
#                     current_array = list(k)
#                     break
#         collection_of_reversals_sorting = collection_of_reversals_sorting[::-1]
#         current_array = list(l)
#         for reversals_path in s2_s1_reversals_path[::-1]:
#             for k, v in reversals_path.items():
#                 if current_array in v:
#                     i, j =  _get_invert_endpoints(list(k), current_array)
#                     collection_of_reversals_sorting.append([i, j])
#                     current_array = list(k)
#                     break
#         collections_of_reversals_sorting.append(collection_of_reversals_sorting)
#     return collections_of_reversals_sorting
    		
# with open("./Rosalind_files/rosalind_sort.txt", "r") as f:
#     a = [int(i) for i in f.readline().strip().split(" ")]
#     b = [int(i) for i in f.readline().strip().split(" ")]

# # main solution
# distance, s1, s2 = 0, set(), set()
# s1.add(tuple(a)), s2.add(tuple(b))
# s1_s2_reversals_path = [] # the list contains reversals path (dict) from s1 to s2.
# s2_s1_reversals_path = [] # the list contains reversals path (dict) from s2 to s1.
# meet_reversals = [] # the reversals met by s1-s2 direction and s2-s1 direction.
# s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals, distance = _get_reversal_distance(s1, s2, distance, s1_s2_reversals_path, s2_s1_reversals_path, meet_reversals)
# collections_of_reversals_sorting = _get_collections_of_reversals_sorting(s1_s2_reversals_path, meet_reversals, s2_s1_reversals_path)
# a_collection = random.sample(collections_of_reversals_sorting, k=1)[0]
# print("[INFO] the reversal distance drev(π,γ): {}".format(len(a_collection)))
# for x in a_collection:
#     print(f"{x[0]} {x[1]}")

######################################################################################################
# #Problem 53
# #Inferring Protein from Spectrum

# import numpy as np

# mslist = []
# with open("./Rosalind_files/rosalind_spec.txt") as f: 
#     string = f.readlines()
#     for x in string:
#         mslist.append(float(x))
    
# amino_acid_masses = {
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

# protein_weights = []
# answers = []

# for x in range(len(mslist)-1):
#     difference = (abs(mslist[x] - mslist[x+1]))
#     r_difference = round(difference, 5)
#     protein_weights.append(r_difference)

# k = list(amino_acid_masses.keys())
# v = np.array(list(amino_acid_masses.values()))

# for x in protein_weights:
#     dist = abs(v - x)
#     arg = np.argmin(dist)
#     answer = k[arg]
#     answers.append(answer)

# print("".join(map(str, answers)))

######################################################################################################
# #Problem 54
# #Introduction to Pattern Matching

# import os

# def insert(root, text, node_count):
#     node = root
#     for char in text:
#         if char not in node:
#             node_count[0] += 1
#             node[char] = {'node_id': node_count[0]}
#         node = node[char]

# def print_trie(node, parent_id, result):
#     for char, child in node.items():
#         if char != 'node_id':  # Check to skip the node_id entry
#             result.append((parent_id, child['node_id'], char))
#             print_trie(child, child['node_id'], result)

# def build_and_print_trie(strings):
#     root = {'node_id': 1}  # Root node with node_id initialized to 1
#     node_count = [1]  # Node counter starting from 1
#     result = []

#     for s in strings:
#         insert(root, s, node_count)

#     print_trie(root, root['node_id'], result)
#     return result

# with open("./Rosalind_files/rosalind_trie.txt") as f: 
#     dna = [line.strip() for line in f.readlines()]
# trie = build_and_print_trie(dna)

# folder_path = './Rosalind_files/trie_output'
# # Create the folder if it doesn't already exist
# os.makedirs(folder_path, exist_ok=True)
# # Define the full path for the new file within the new subfolder
# file_path = os.path.join(folder_path, 'output_trie.txt')

# with open(file_path, "w") as f:
#     for edge in trie:
#         print(f"{edge[0]} {edge[1]} {edge[2]}", file=f)

######################################################################################################
# #Problem 55
# #Comparing Spectra with the Spectral Convolution

# def spectral_convolution(S1, S2):
#     from collections import defaultdict
    
#     # Dictionary to count multiplicities of each difference
#     differences = defaultdict(int)
    
#     # Compute all differences and their multiplicities
#     for s1 in S1:
#         for s2 in S2:
#             diff = s1 - s2
#             diff = round(diff, 5)
#             differences[diff] += 1

#     # Determine the maximal multiplicity and the corresponding difference
#     max_multiplicity = max(differences.values())
#     shift_value = max(differences, key=lambda r: differences[r])  # Get the difference with the max multiplicity
    
#     return max_multiplicity, abs(shift_value)

# with open("./Rosalind_files/rosalind_conv.txt", 'r') as file:
#     lines = file.readlines()
#     S1 = [float(num) for num in lines[0].strip().split()]
#     S2 = [float(num) for num in lines[1].strip().split()] 

# max_multiplicity, shift_value = spectral_convolution(S1, S2)
# print(max_multiplicity)
# print(shift_value)


# ######################################################################################################
# """The same, but with counter"""

# from collections import Counter

# # Calculate all possible differences
# differences = [round(s1 - s2, 5) for s1 in S1 for s2 in S2]

# # Use Counter to count the frequency of each difference
# difference_counts = Counter(differences)

# # Find the most common difference and its multiplicity
# most_common_difference, max_multiplicity = difference_counts.most_common(1)[0]

# print(max_multiplicity)
# print(most_common_difference)

######################################################################################################
# #Problem 56
# #Creating a Character Table

# from Bio import Phylo
# from io import StringIO

# with open("./Rosalind_files/rosalind_ctbl.txt") as f: 
#     newick_tree = f.read()

# # Read the Newick formatted tree using Bio.Phylo
# tree = Phylo.read(StringIO(newick_tree), "newick")

# def find_splits(tree):
#     splits = set()  # Using a set to avoid duplicates
#     taxa = sorted(clade.name for clade in tree.find_clades() if clade.name)

#     def collect_taxa(clade):
#         return set(node.name for node in clade.find_clades() if node.name)

#     for clade in tree.get_nonterminals():  # Only non-terminal nodes define splits
#         taxa_one_side = collect_taxa(clade)
#         taxa_other_side = set(taxa) - taxa_one_side

#         # Make sure both sides have at least one taxon
#         if taxa_one_side and taxa_other_side:
#             # Choose the smaller side for '1's to match the example output
#             if len(taxa_one_side) < len(taxa_other_side):
#                 smaller_side = taxa_one_side
#             else:
#                 smaller_side = taxa_other_side

#             # Create binary string for this split
#             binary_split = ''.join(['1' if taxon in smaller_side else '0' for taxon in taxa])
#             splits.add(binary_split)  # Add to set to keep unique

#     return splits

# # Compute and print all unique splits
# splits = find_splits(tree)
# for split in sorted(splits):
#     print(split)

######################################################################################################
# #Problem 57
# #Constructing a De Bruijn Graph

"""Skip for now"""

######################################################################################################
# #Problem 58
# #Reconstructing Edit Distance

# from Bio import SeqIO
# s, t = [str(rec.seq).strip() for rec in SeqIO.parse("./Rosalind_files/rosalind_edta.txt", "fasta")]

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
#                 L[i][j] = L[i-1][j-1]
#             else:
#                 L[i][j] = 1 + min(L[i - 1][j],    # Insertion
#                                    L[i][j - 1],    # Deletion
#                                    L[i - 1][j - 1]) # Substitution
                
#     return L, L[m][n]

# def reconstruct_alignment(L, X, Y):
#     s_edit, t_edit = "", ""
#     m, n = len(X), len(Y)

#     while m > 0 and n > 0: #character match
#         if X[m - 1] == Y[n - 1]:
#             s_edit = X[m - 1] + s_edit
#             t_edit = Y[n - 1] + t_edit
#             m -= 1
#             n -= 1

#         elif L[m][n] == L[m - 1][n] + 1: #Insertion in Y (Deletion in X)
#             s_edit = X[m - 1] + s_edit
#             t_edit = "-" + t_edit
#             m -= 1

#         elif L[m][n] == L[m][n - 1] + 1: #Insertion in X (Deletion in Y)
#             s_edit = "-" + s_edit
#             t_edit = Y[n - 1] + t_edit
#             n -= 1

#         else: #Substitution
#             s_edit = X[m - 1] + s_edit
#             t_edit = Y[n - 1] + t_edit
#             m -= 1
#             n -= 1

#     return s_edit, t_edit

# table, edit_distance_value = scs_length(s, t)
# edited_s, edited_t = reconstruct_alignment(table, s, t)

# print(edit_distance_value)
# print(edited_s)
# print(edited_t)

######################################################################################################
# #Problem 59
# #Inferring Peptide from Full Spectrum

"""This is not my code, but I understand how it works.
It is very clever in a few different ways."""

# MASS_TABLE = {
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

# def find_weight_match(current_w, w_list):
#     for weight in w_list:
#         for aa, wt in MASS_TABLE.items():
#             if abs(wt - (weight - current_w)) < 0.01: #somewhat inefficient
#                 return aa
#     return -1


# with open("./Rosalind_files/rosalind_full.txt", 'r') as file:
#     input_lines = file.read().splitlines()

# weights = [float(x) for x in input_lines]
# n = (len(weights) - 3) // 2  # Ensure integer division

# protein = ""
# current = weights[1]
# remaining_weights = weights[2:]

# while len(protein) < n:
#     aa = find_weight_match(current, remaining_weights)
#     if aa == -1:
#         break
#     else:
#         protein += aa
#         current += MASS_TABLE[aa]
#         remaining_weights = list(filter(lambda w: w - current > 0, remaining_weights)) # filter instead of sort
#         print(f"Added {aa}, new current mass: {current}, remaining weights count: {len(remaining_weights)}")

# print(protein)

######################################################################################################
# #Problem 60
# #Independent Segregation of Chromosomes

# import math
# from scipy.stats import binom

# with open("./Rosalind_files/rosalind_indc.txt") as f:
#     n = int(f.read())

# def binomial_commmon_log(n):
#     total_chromosomes = 2*n
#     probabilities = []

#     for k in range(1, total_chromosomes+1):
#         cumulative_prob = binom.sf(k-1, total_chromosomes, 0.5)
#         com_log = math.log10(cumulative_prob)
#         probabilities.append(round(com_log, 3))
#     return probabilities

# result = binomial_commmon_log(n)
# print(' '.join(map(str, result)))

######################################################################################################
# #Problem 61
# #Finding Disjoint Motifs in a Gene

# with open("./Rosalind_files/rosalind_itwv.txt", 'r') as file:
#     lines = file.readlines()
#     s = lines[0]
#     strings = [line.rstrip() for line in lines[1:]]

# def interleave(s, t, res, i, j, lis, seen):
#     # Base case: when both indices reach the end of their respective strings
#     if i == len(s) and j == len(t):
#         if res not in seen:  # Check if the result is unique
#             lis.append(res)
#             seen.add(res)  # Mark this interleaving as seen
#         return
#     # Recursive call if more characters can be taken from s
#     if i < len(s):
#         interleave(s, t, res + s[i], i + 1, j, lis, seen)
#     # Recursive call if more characters can be taken from t
#     if j < len(t):
#         interleave(s, t, res + t[j], i, j + 1, lis, seen)

# def generate_all_interleavings(strings):
#     interleaves = []  # List to hold all the lists of interleavings for each pair
#     # Generate interleavings for each possible combination
#     for s1 in strings:
#         for s2 in strings:
#             current_interleaves = []  # List for interleavings of s1 with s2
#             seen = set()  # Set to track unique interleavings
#             interleave(s1, s2, "", 0, 0, current_interleaves, seen)
#             interleaves.append(current_interleaves)
#     return interleaves

# interleaves = generate_all_interleavings(strings)

# def find_if_in_string(section, s):
#     if any(subsec in s for subsec in section):
#         return 1
#     else:
#         return 0

# n = len(strings)
# results = [[0] * n for _ in range(n)]

# count = 0
# for j in range(n):
#     for k in range(n):
#         results[j][k] = find_if_in_string(interleaves[count], s)
#         count += 1
        
# for row in results:
#     print(' '.join(map(str, row)))

######################################################################################################
# #Problem 62
# #Finding the Longest Multiple Repeat

