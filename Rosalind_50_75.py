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

import itertools

p = [1,6,7,5,3,10,4,9,8,2]
y = [4,1,8,6,3,2,9,5,10,7]

p_perm1 = []
y_perm1 = []

for i in range(len(p)-1):
    for j in range(i+1,len(p)+1):
        if i != j and len(p[i:j]) > 1:
            new = ((p[:i] + p[i:j][::-1] + p[j:]), (i+1,j))
            p_perm1.append(new)

for x in p_perm1:
    if x == y:
        print(x)

for i in range(len(y)-1):
    for j in range(i+1,len(y)+1):
        if i != j and len(y[i:j]) > 1:
            new = ((y[:i] + y[i:j][::-1] + y[j:]), (i+1,j))
            y_perm1.append(new)

for x in p_perm1:
    for y in y_perm1:
        if x[0] == y[0]:
            print(x[1], y[1])

#now we need to write code that finds which of x in p_perm1 and y in y_perm1 are closest to each other and then find the reversal distances of those, if it exists.
#else, we do as many reversals as required
