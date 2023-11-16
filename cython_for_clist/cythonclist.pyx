import copy
from timeit import default_timer as timer

# python3 setup.py build_ext --inplace

# number of particles
cdef int n = 7
# collision limit
cdef int c_lim = 4


# essentially does part_num C n and returns a list of these combinations
def create_base_arr(int part_num):
    base = []
    cdef int i
    for p in range(1, part_num + 1):
        i = p + 1
        while i <= part_num:
            base.append([[p, i]])
            i += 1
    # base = np.array(base)
    return base


start = timer()
base_arr = create_base_arr(n)
global_col_arr = create_base_arr(n)
end = timer()
print("time taken for base:", end - start)

def create_col_lists(int iteration, int col_lim, prev_arr):
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


create_col_lists(1, c_lim, create_base_arr(n))
end = timer()
print()
# for i in range(len(global_col_arr)):
#    print(global_col_arr[i])
#    print()
print("number of collision lists: ", len(global_col_arr))

print("time taken:", end - start)
# f.close()

