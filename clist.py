import copy
from timeit import default_timer as timer
from random import randrange

import numpy as np
import sympy as sym
from particled import particled

# number of particles
n = 12
# collision limit
c_lim = 9


# essentially does part_num C n and returns a list of these combinations
def create_base_arr(part_num):
    base = []
    for p in range(1, part_num + 1):
        i = p + 1
        while i <= part_num:
            base.append([p, i]) # if computing all lists -> [[p, i]] if random -> [p, i]
            i += 1
    # base = np.array(base)
    return base


start = timer()
base_arr = create_base_arr(n)
global_col_arr = create_base_arr(n)
end = timer()
print("time taken for base:", end - start)


# def create_col_lists(iteration, col_lim, prev_arr, global_col_arr):
#     new_arr = [[]] # <---- problem here
#     new_arr = np.array([[0, 0]])
#     print(prev_arr)
#     if iteration >= col_lim:
#         return global_col_arr
#     else:
#         for clist in prev_arr:
#             for base in base_arr:
#                 print(clist)
#                 l = np.concatenate((clist, base), axis=0)
#                 if np.array_equal(clist, base) == False:
#                     # print(np.concatenate((global_col_arr, l)))
#                     # global_col_arr = np.append(global_col_arr, l, axis=0, dtype=object)
#                     new_arr = np.insert(new_arr, l, axis=1)
#                     print(new_arr)
#         create_col_lists(iteration + 1, col_lim, new_arr, global_col_arr)
#
#
# global_col_arr = create_col_lists(1, c_lim, base_arr, global_col_arr)
# end = timer()
# print("global: ")
# for i in range(len(global_col_arr)):
#     print(global_col_arr[i])
#     print()
# print("number of collision lists: ", len(global_col_arr))
#
# print("time taken for whole program:", end - start)

#f = open("c_arrs", "w")


# delete count var later
# def create_col_lists(iteration, col_lim, prev_arr, count):
#     new_arr = []
#     if iteration >= col_lim:
#         print(count)
#         return count
#     else:
#         for clist in prev_arr:
#             for base in base_arr:
#                 l = copy.deepcopy(clist)
#                 if clist[-1] != base[0]:
#                     l.append(base[0])
#                     # global_col_arr.append(l)
#                     #global_col_arr.append([clist, base[0]])
#                     new_arr.append(l)
#                     # new_arr.append([clist, base[0]])
#         print(len(new_arr))
#         count += len(new_arr)
#         f.write(str(new_arr))
#         f.write("\n")
#         create_col_lists(iteration + 1, col_lim, new_arr, count)

clist = []

def rand_list(clim):
    while (len(clist) < clim):
        i = randrange(0, len(base_arr))
        if len(clist) == 0:
            clist.append(base_arr[i])
        elif clist[-1] != base_arr[i]:
            clist.append(base_arr[i])

def create_col_lists(iteration, col_lim, prev_arr):
    new_arr = []
    if iteration >= col_lim:
        return
    else:
        for clist in prev_arr:
            for base in base_arr:
                l = copy.deepcopy(clist)
                if clist[-1] != base[0]:
                    l.append(base[0])
                    global_col_arr.append(l)
                    new_arr.append(l)
        create_col_lists(iteration + 1, col_lim, new_arr)


#create_col_lists(1, c_lim, create_base_arr(n))
rand_list(c_lim)
end = timer()
print()
print(clist)
# for i in range(len(global_col_arr)):
#     print(global_col_arr[i])
#     print()
# print("number of collision lists: ", len(global_col_arr))

print("time taken:", end - start)
#f.close()

'''
vars we need - we always know the initial position at all collision times (these can be updated after each collision)
             - v initial
             - v final
             - x at collision
             - x at the end - we know the x after a particle's final collision (this could be checked for)
             - time of collision 
'''