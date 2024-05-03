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