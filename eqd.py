import numpy as np
from scipy.optimize.nonlin import NoConvergence
from sympy import sympify, lambdify, nsolve
from particled import particled
from timeit import default_timer as timer
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton_krylov, broyden1, broyden2, anderson, minimize
import warnings
from sympy.utilities.lambdify import NUMPY_TRANSLATIONS

NUMPY_TRANSLATIONS["zoo"] = "nan"
from random import randrange

# number of particles
n = 5

# number of spacial dimensions
d = 2

# collision limit
c_lim = 5

# total time of the process
T = 100

# Time array
tarr = np.empty(c_lim + 2, dtype=object)
for i in range(c_lim + 1):
    if i == 0:
        tarr[i] = 0
    else:
        tarr[i] = "t" + str(i)
tarr[-1] = T


# collision list -> create rand clist using base_arr
def create_base_arr(part_num):
    base = []
    for p in range(1, part_num + 1):
        i = p + 1
        while i <= part_num:
            base.append([p, i])
            i += 1
    # base = np.array(base)
    return base


base_arr = create_base_arr(n)

clist = []


def rand_list(clim):
    rand = np.random.default_rng()
    while len(clist) < clim:
        i = rand.integers(low=0, high=len(base_arr))
        if len(clist) == 0:
            clist.append(base_arr[i])
        elif clist[-1] != base_arr[i]:
            clist.append(base_arr[i])


rand_list(c_lim)
clist = np.array(clist)
print(clist)
# clist = np.array([[1, 2], [2, 3], [3, 4], [4, 5]])
# clist = np.array([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [7, 8],
#                   [8, 9], [9, 10], [10, 11], [11, 12], [1, 2], [2, 3], [3, 4]])
# clist = np.array([[1, 2], [1, 2], [3, 4], [12, 2], [2, 3], [4, 1], [5, 2], [2, 1], [3, 1]])

'''
part_info_list = np.array([[-12, 10, 10, 10, 17, 10, 0],
                           [20, 20, 12, 20, 20, 20, 0],
                           [30, 30, 30, 30, 30, 30, 0],
                           [40, 42, 80, 50, 60, 70, 0],
                           [50, 100, 110, 67, 80, 150, 0],
                           [110, 107, 180, 120, 140, 165, 0]])

part_info_list = np.array([[1, 1, 1, 1, 1, 1, 0],
                           [2, 2, 2, 2, 2, 2, 0],
                           [3, 3, 3, 3, 3, 3, 0],
                           [4, 4, 4, 4, 4, 4, 0],
                           [5, 5, 5, 5, 5, 5, 0],
                           ])


part_info_list = np.array([[-1.2, 1.0, 1.0, 1.0, -1.7, 10, 0],
                           [2.0, 2.0, 1.2, 2.0, 2.0, 2.0, 0],
                           [3.0, 3.0, 3.5, 3.3, 3.1, 3.0, 0],
                           [4.0, 4.2, 8.0, -5.0, 6.0, 7.0, 0],
                           [5.0, 10.0, 11.0, 6.7, 8.0, 15.0, 0],
                           # [11.0, 10.7, 18.0, 12.0, 14.0, 16.5, 0],
                           # [12, 12, 13, 16, 18, 15.9, 0],
                           # [21, 20.1, 15, 17, 16, 2, 0],
                           # [6, 15, 14, 20, 20, 20, 0],
                           # [18, 18, 18, 21, 22, 23, 0],
                           # [20.5, -1, 20.5, 19, 22, 15, 0],
                           # [20, 20.111, 20.432, 10, 20.8, 19, 0]
                           ])
'''
part_info_list = np.array([[-1.2, 1.0, -1.7, 10, 0],
                           [2.0, 2.0, 2.0, 2.0, 0],
                           [3.0, 3.0, 3.1, 3.0, 0],
                           [4.0, 4.2,6.0, 7.0, 0],
                           [5.0, 10.0, 8.0, 15.0, 0]])

# part_info_list = 100000* part_info_list

# MAYBE INC TIME SO V GOES DOWN -> THIS WAY THE RANGE OF RAND GUESSES ARE SMALL
# OR CREATE RAND IN A RANGE AND THEN TAKE THAT RANGE/T (TOTAL TIME)

guess = []


def count_col_num(part_num):
    count = 0
    for x in np.nditer(clist):
        if part_num == x:
            count += 1
    return count


p_num = 1
for r in part_info_list:
    r[-1] = count_col_num(p_num)
    p_num += 1

# particle list
part_list = np.empty(n, dtype=particled)


def times_for_part(part_num, col_num):
    part_tarr = np.empty(col_num + 2, dtype=object)
    i = 1
    j = 1
    part_tarr[0] = 0
    part_tarr[-1] = T
    for r in clist:
        for c in r:
            if c == part_num:
                part_tarr[j] = tarr[i]
                j += 1
        i += 1
    return part_tarr


i = 0
for r in part_info_list:
    xi = np.empty(d)
    xf = np.empty(d)
    for x_comp in range(d):
        xi[x_comp] = r[x_comp]
        xf[x_comp] = r[x_comp + d]
    p = particled(xi, xf, int(r[-1]), i + 1, d, times_for_part(i + 1, int(r[-1])))
    part_list[i] = p
    i += 1


def reset_index():
    for p in part_list:
        particled.reset_vindex(p)
        particled.reset_pindex(p)
        particled.reset_tindex(p)


# print statement
# for p in part_list:
#     print("num: " + str(p.pnum))
#     print(p.parr)


def reset_all_pos():
    for p in part_list:
        particled.reset_pos(p)


# xi + v1(delta t) + v2(delta t2) + ... = xf
def first_eqs():
    i = 0
    for r in clist:  # there will still be eqs1 for particles with no collisions
        for c in r:
            p = part_list[c - 1]
            for p_comp in range(d):
                p.pos[p_comp] += " + " + str(p.varr[p.vindex, p_comp]) + " * (" + str(p.tarray[p.tindex + 1]) \
                                 + " - " + str(p.tarray[p.tindex]) + ")"
                p.col_arr[p.pindex, p_comp] = p.pos[p_comp]
            particled.inc_vindex(p)
            particled.inc_pindex(p)
            particled.inc_tindex(p)
            particled.update_cindex(p, i + 1)  # might cause problems later on
        i += 1
    eqs = np.empty(d * n, dtype=object)
    i = 0
    for p in part_list:
        for p_comp in range(d):
            eqs[i + p_comp] = p.pos[p_comp] + " + " + p.varr[-1, p_comp] + " * (" + str(tarr[-1]) + \
                              " - " + str(tarr[p.cindex]) + ") - " + str(p.parr[-1, p_comp])
        i += d
    return eqs


# momentum consv
def m_cons_eqs():
    i = 0
    eqs = np.empty(d * c_lim, dtype=object)
    for r in clist:
        p1 = part_list[r[0] - 1]
        p2 = part_list[r[1] - 1]
        for v_comp in range(d):
            eq = p1.varr[p1.vindex, v_comp] + " + " + p2.varr[p2.vindex, v_comp] + \
                 " - " + p1.varr[p1.vindex + 1, v_comp] + " - " + p2.varr[p2.vindex + 1, v_comp]
            eqs[i + v_comp] = eq
        particled.inc_vindex(p1)
        particled.inc_vindex(p2)
        i += d
    return eqs


def m_cons_pos():
    i = 0
    eqs = np.empty(d * c_lim, dtype=object)
    for r in clist:
        p1 = part_list[r[0] - 1]
        p2 = part_list[r[1] - 1]
        for p_comp in range(d):
            eq = "(" + str(p1.parr[p1.pindex, p_comp]) + " / (" + str(p1.tarray[p1.tindex + 1]) \
                 + " - " + str(p1.tarray[p1.tindex]) + ")) + (" \
                 + str(p2.parr[p2.pindex, p_comp]) + " / (" + str(p2.tarray[p2.tindex + 1]) \
                 + " - " + str(p2.tarray[p2.tindex]) + "))"
            eq += " - (" + str(p1.parr[p1.pindex + 1, p_comp]) + " / (" + str(p1.tarray[p1.tindex + 2]) \
                  + " - " + str(p1.tarray[p1.tindex + 1]) + "))" + " - (" \
                  + str(p2.parr[p2.pindex + 1, p_comp]) + " / (" + str(p2.tarray[p2.tindex + 2]) \
                  + " - " + str(p2.tarray[p2.tindex + 1]) + "))"
            eqs[i + p_comp] = eq
        particled.inc_pindex(p1)
        particled.inc_pindex(p2)
        particled.inc_tindex(p1)
        particled.inc_tindex(p2)
        i += d
    return eqs


# E consv
def e_cons_eqs():
    i = 0
    eqs = np.empty(c_lim, dtype=object)
    # eqs3 = np.empty(c_lim, dtype=object)
    for r in clist:
        p1 = part_list[r[0] - 1]
        p2 = part_list[r[1] - 1]
        eq = ""
        for v_comp in range(d):
            if v_comp != 0:
                eq += " + "
            eq += p1.varr[p1.vindex, v_comp] + "**2 + " + p2.varr[p2.vindex, v_comp] + "**2"
            eq += " - " + p1.varr[p1.vindex + 1, v_comp] + "**2 - " + p2.varr[p2.vindex + 1, v_comp] + "**2"
        particled.inc_vindex(p1)
        particled.inc_vindex(p2)
        eqs[i] = eq
        i += 1
    return eqs


def e_cons_pos():
    i = 0
    eqs = np.empty(d * c_lim, dtype=object)
    for r in clist:
        p1 = part_list[r[0] - 1]
        p2 = part_list[r[1] - 1]
        for p_comp in range(d):
            eq = "(" + str(p1.parr[p1.pindex, p_comp]) + " / (" + str(p1.tarray[p1.tindex + 1]) \
                 + " - " + str(p1.tarray[p1.tindex]) + "))**2 + (" \
                 + str(p2.parr[p2.pindex, p_comp]) + " / (" + str(p2.tarray[p2.tindex + 1]) \
                 + " - " + str(p2.tarray[p2.tindex]) + "))**2"
            eq += " - (" + str(p1.parr[p1.pindex + 1, p_comp]) + " / (" + str(p1.tarray[p1.tindex + 2]) \
                  + " - " + str(p1.tarray[p1.tindex + 1]) + "))**2" + " - (" \
                  + str(p2.parr[p2.pindex + 1, p_comp]) + " / (" + str(p2.tarray[p2.tindex + 2]) \
                  + " - " + str(p2.tarray[p2.tindex + 1]) + "))**2"
            eqs[i + p_comp] = eq
        particled.inc_pindex(p1)
        particled.inc_pindex(p2)
        particled.inc_tindex(p1)
        particled.inc_tindex(p2)
        i += d
    return eqs


# col positions are equal
def fourth_eqs():
    i = 0
    eqs = np.empty(d * c_lim, dtype=object)
    for r in clist:
        p1 = part_list[r[0] - 1]
        p2 = part_list[r[1] - 1]
        for p_comp in range(d):
            eq = p1.col_arr[p1.pindex, p_comp] + " - (" + p2.col_arr[p2.pindex, p_comp] + ")"
            eqs[i + p_comp] = eq
        particled.inc_pindex(p1)
        particled.inc_pindex(p2)
        i += d
    return eqs


def fourth_eq_pos():
    i = 0
    eqs = np.empty(d * c_lim, dtype=object)
    for r in clist:
        p1 = part_list[r[0] - 1]
        p2 = part_list[r[1] - 1]
        for p_comp in range(d):
            eq = str(p1.parr[p1.pindex + 1, p_comp]) + " - " + str(p2.parr[p2.pindex + 1, p_comp])
            eqs[i + p_comp] = eq
        particled.inc_pindex(p1)
        particled.inc_pindex(p2)
        i += d
    return eqs


# for p in part_list:
#     print(p.parr)
# eq_mp = m_cons_pos()
# print(eq_mp)
# reset_index()
# eq_ep = e_cons_pos()
# print(eq_ep)
# reset_index()
# eq_four = fourth_eq_pos()
# print(eq_four)
# reset_index()
# print(len(eq_ep) +len(eq_mp) + len(eq_four))

# changed
def sym_length():
    count = 0
    count += len(tarr) - 2
    for p in part_list:
        count += len(p.varr.reshape(d * (p.ncol + 1)))
        # count += len(p.parr.reshape(d * (p.ncol + 2))) - d * 2
    return count


# changed
def sym_list():
    i = 0
    arr = np.empty(sym_length(), dtype=object)
    for j in range(1, len(tarr) - 1):
        arr[i] = (sympify(tarr[j]))
        i += 1
    for p in part_list:
        for j in range(len(p.varr)):
            for v_comp in range(d):
                arr[i] = (sympify(p.varr[j, v_comp]))
                i += 1
    # for p in part_list:
    #     for k in range(len(p.parr) - 2):
    #         for p_comp in range(d):
    #             arr[i] = (sympify(p.parr[k + 1, p_comp]))
    #             i += 1
    return arr


start = timer()


# changed
def run_eqs():
    eqs1 = first_eqs()
    reset_index()
    eqs2 = m_cons_eqs()
    reset_index()
    eqs3 = e_cons_eqs()
    reset_index()
    eqs4 = fourth_eqs()
    reset_index()
    eq_mp = m_cons_pos()
    reset_index()
    eq_ep = e_cons_pos()
    reset_index()
    eq_four_p = fourth_eq_pos()
    eqs = np.concatenate((eqs1, eqs2, eqs3, eqs4), axis=None)
    reset_index()
    reset_all_pos()
    return eqs


def convert_sym_eqs(eqs):
    eqs_vars = np.empty(len(eqs), dtype=object)
    eindex = 0
    for eq in eqs:
        eqs_vars[eindex] = sympify(eqs[eindex])
        eindex += 1
    return eqs_vars


def first_appearance(p_num):
    for r in clist:
        if p_num == r[0]:
            return True
        if p_num == r[1]:
            return False


# assuming odd = pos and odd = negative
# -   +    -   +   -   +   -   +   -  +

#  r[0] => positive
#  r[1] => negative
# go through the list and check the pos of that particle -> if statements
# first part (v intial) -> even = negative and odd = positive
# might be fine to just alternate between pos -> negative since at least 1/d will be right

# particle information list---> structure:
# xi, xf, 0 -> if d > 1 xi and xf elements are just written out next to each other

# changed
def create_guess():
    rand = np.random
    tarr[-1] = T
    for p in part_list:
        particled.update_final_time(p, T)
    guess_t = np.empty(tarr.size - 2)
    for i in range(guess_t.size):
        guess_t[i] = rand.random() * T
    guess_t.sort()
    guess_v = []
    # try no rules
    for p in part_list:
        if first_appearance(p.pnum):
            for i in range(int(p.varr.size / d)):
                if i % 2 == 0:
                    for v_comp in range(d):
                        guess_v.append(rand.random())

                elif i % 2 == 1:
                    for v_comp in range(d):
                        guess_v.append(-1 * rand.random())
        else:
            for i in range(int(p.varr.size / d)):
                if i % 2 == 0:
                    for v_comp in range(d):
                        guess_v.append(-1 * rand.random())
                elif i % 2 == 1:
                    for v_comp in range(d):
                        guess_v.append(rand.random())
    # guess_p = []
    # for p in part_list:
    #     for i in range(len(p.parr) - 2):
    #         for p_comp in range(d):
    #             low = min(p.parr[0, p_comp], p.parr[-1, p_comp])
    #             high = max(p.parr[0, p_comp], p.parr[-1, p_comp])
    #             guess_p.append(rand.randint(low, high + 1))
    guess = np.concatenate((guess_t, guess_v))
    # print("guess: " + str(len(guess)))
    return guess


vars_list = sym_list()
sol = []

warnings.filterwarnings("error")
print(run_eqs())
eqs_vars = convert_sym_eqs(run_eqs())
print(eqs_vars)
f = lambdify([vars_list], eqs_vars)
while True:
    # print("T: " + str(T))
    wrongSol = False
    try:
        # eqs_vars = convert_sym_eqs(run_eqs())
        # print("vars: " + str(len(vars_list)))
        # print("eqs: " + str(len(eqs_vars)))
        # sol = nsolve(eqs_vars, vars_list, create_guess())
        f = lambdify([vars_list], eqs_vars)
        sol = fsolve(f, create_guess())
        for i in range(len(tarr) - 1):
            if sol[i] > T or sol[i] < 0:
                wrongSol = True
                break
        if not wrongSol:
            print("rerun")
            break
    except RuntimeWarning as e:
       # e.msg
    #     # RuntimeWarning or ValueError or NoConvergence:
    #     T *= 10  # is 10 the best number?
    #     tarr[-1] = T
    #     for p in part_list:
    #         particled.update_final_time(p, T)
        print("except")



i = 0
for s in sol:
    print(str(vars_list[i]) + ": ", s)
    i += 1


def count_total_col_pos():
    count = 0
    for p in part_list:
        count += p.col_arr.size
    return count


def pos_eqs():
    eqs = np.empty(count_total_col_pos(), dtype=object)
    i = 0
    for p in part_list:
        for j in range(len(p.col_arr)):
            for p_comp in range(d):
                eq = p.col_arr[j, p_comp]  # + " - " + p.parr[j + 1, p_comp]
                eqs[i + p_comp] = eq
            i += d
    return eqs


def variable_to_value(var, solnt):
    arr = np.empty((len(var), 2), dtype=object)
    for i in range(len(var)):
        arr[i, 0] = var[i]
        arr[i, 1] = solnt[i]
    return arr


def K_check_eqs():
    eqs_K = np.empty((len(tarr) - 1), dtype=object)
    kindex = 0
    for i in range(len(tarr) - 1):
        eq = ''
        if i != 0:
            p1 = part_list[clist[i - 1][0] - 1]
            p2 = part_list[clist[i - 1][1] - 1]
            p1.inc_vindex()
            p2.inc_vindex()
        for p_comp in range(d):
            for p in part_list:
                eq += p.varr[p.vindex, p_comp] + "**2"
                if p.pnum != n or p_comp != d - 1:
                    eq += " + "
        eqs_K[kindex] = eq
        kindex += 1
    j = 0
    print(eqs_K)
    eqs_vars = convert_sym_eqs(eqs_K)
    eqs_num = np.empty(len(eqs_vars))
    for eq in eqs_vars:
        eqs_num[j] = eq.subs(variable_to_value(vars_list, sol))
        j += 1
    print("KE arr:")
    print(eqs_num)


# [[1 4]
#  [1 5]
#  [4 5]
#  [2 4]
#  [2 3]]

K_check_eqs()


def find_col_pos():
    eqs_vars = convert_sym_eqs(pos_eqs())
    eqs_pos = np.empty(len(eqs_vars))
    i = 0
    for eq in eqs_vars:
        eqs_pos[i] = eq.subs(variable_to_value(vars_list, sol))
        i += 1
    eqs_pos = eqs_pos.reshape((int(eqs_vars.size / d), d))
    # might be better to create a new array to store this data in
    j = 0
    for p in part_list:
        for x_fill in range(len(p.col_arr)):
            p.parr[x_fill + 1] = eqs_pos[j]
            j += 1
    return eqs_pos


find_col_pos()
#
for p in part_list:
     # print("pnum: " + str(p.pnum))
     # print(p.parr)
    p.print()

end = timer()
print("\ntime taken: " + str(end - start))

# for i in range(len(tarr) - 2):
#     print(str(clist[i]) + " : " + str(sol[i]))

fig = plt.figure()  # figsize=(5, 5))
axes = fig.add_subplot(projection='3d')
xlist = []
ylist = []
zlist = []
for p in part_list:
    x_pos = p.parr[:, 0]
    for x in x_pos:
        xlist.append(x)
    y_pos = p.parr[:, 1]
    for y in y_pos:
        ylist.append(y)
    if d == 3:
        z_pos = p.parr[:, 2]
        for z in z_pos:
            zlist.append(z)
    else:
        z_pos = 0
    axes.scatter(xs=x_pos, ys=y_pos, zs=z_pos, label='particle' + str(p.pnum))
    axes.plot(x_pos, y_pos, z_pos)
xlist.sort()
ylist.sort()
zlist.sort()
axes.set_xlim((xlist[0] - 1), (xlist[-1] + 1))
axes.set_ylim((ylist[0] - 0.001), (ylist[-1] + 0.001))
if d != 2:
    axes.set_zlim((zlist[0] - 0.001), (zlist[-1] + 0.001))
axes.legend()
axes.grid(False)
fig.set_size_inches(6, 6)
plt.figure()
data = np.array(sol[0: len(tarr) - 2])
if c_lim != 1:
    data -= data.min()
    data /= data.max()
plt.hist(data)
plt.show()