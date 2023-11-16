from scipy import *
from pylab import *
from scipy.optimize import minimize
import sys


def min_fun1(var, r1i, r2i, r1f, r2f):
    theta = var[0]
    tc = var[1]

    return (theta - 0.5) ** 2 + (tc - 0.4) ** 2 + 0.00001


def min_fun(var, r1i, r2i, r1f, r2f):
    global nhat
    theta = var[0]
    tc = var[1]
    nhat[0] = cos(theta)
    nhat[1] = sin(theta)
    return coll(nhat, tc, r1i, r2i, r1f, r2f)


def com_v(r1i, r2i, r1f, r2f, rc1, rc2):
    vcmi = zeros((2))
    vcmf = zeros((2))
    vcmi = 0.5 * ((rc1 - r1i) / (tc - ti1) + (rc2 - r2i) / (tc - ti2))
    vcmf = 0.5 * ((rc1 - r1f) / (tc - tf1) + (rc2 - r2f) / (tc - tf2))
    return vcmi, vcmf


def com_r(r1i, r2i, r1f, r2f, rc1, rc2):
    vcmi, vcmf = com_v(r1i, r2i, r1f, r2f, rc1, rc2)
    r1i_p = r1i - vcmi * ti1
    r2i_p = r2i - vcmi * ti2
    r1f_p = r1f - vcmi * tf1
    r2f_p = r2f - vcmi * tf2
    rc1_p = rc1 - vcmi[:] * tc
    rc2_p = rc2 - vcmi[:] * tc
    return r1i_p, r2i_p, r1f_p, r2f_p, rc1_p, rc2_p


def does_intersect(rc1, rc2, r1i, r2i, r1f, r2f):
    r1i_p, r2i_p, r1f_p, r2f_p, rc1_p, rc2_p = com_r(r1i, r2i, r1f, r2f, rc1, rc2)
    n1 = r1i_p - rc1_p
    n1 /= sqrt(sum(n1 ** 2))
    n2 = r2i - rc2
    n2 /= sqrt(sum(n2 ** 2))

    nc = rc2_p - rc1_p
    nc /= sqrt(sum(nc ** 2))
    cos1 = sum(n1 * nc)
    cos2 = -sum(n2 * nc)
    print("cos1 = ", cos1, "cos2 = ", cos2)
    acos1 = arccos(cos1)
    acos2 = arccos(cos2)
    # print("cos1 = ", cos1, " cos2 = ", cos2, "acos1 = ", acos1, "acos2 = ", acos2, "sin(1+2) =", sin(acos1+acos2))
    # print("r1i = ", r1i, "r2i = ", r2i, "rc1 = ", rc1,  "rc2 = ", rc2)
    # print("n1 = ", n1,"n2 = ", n2,"nc = ", nc)

    sign = sqrt(1 - cos1 ** 2) * cos2 + sqrt(1 - cos2 ** 2) * cos1
    return sign > 0


def coll(nhat, tc, r1i, r2i, r1f, r2f):
    global r1c_perp, r2c_perp, r1c_par, r2c_par, r1c_perp_new, ti1, ti2, tf1, tf2, vi1_perp, vi2_perp, vf1_perp, vf2_perp, vi1_par, vi2_par, vf1_par, vf2_par, r1i_perp

    nhat_perp = array([-nhat[1], nhat[0]])

    r1i_perp = dot(r1i, nhat)
    r2i_perp = dot(r2i, nhat)
    r1f_perp = dot(r1f, nhat)
    r2f_perp = dot(r2f, nhat)

    r1i_par = dot(r1i, nhat_perp)
    r2i_par = dot(r2i, nhat_perp)
    r1f_par = dot(r1f, nhat_perp)
    r2f_par = dot(r2f, nhat_perp)

    ti1 = 0.5001
    ti2 = 0.0
    tf1 = 1.0
    tf2 = 0.8

    ti1 = 0.2001
    ti2 = 0.0
    tf1 = 1.0
    tf2 = 0.8

    ti1 = 0.0
    ti2 = 0.0001
    tf1 = 1.0
    tf2 = 0.999

    ti1 = 0.13269306245562518
    ti2 = 0.0
    tf1 = 0.54545454545454541
    tf2 = 0.53587266720110616

    ti1 = 0.363636
    ti2 = 0.545455
    tf1 = 0.818182
    tf2 = 0.909091

    ti1 = 0.000000
    ti2 = 0.197845
    tf1 = 0.674864
    tf2 = 0.363636

    #   bad case:
    ti1 = 0.000000
    ti2 = 0.329561
    tf1 = 0.441895
    tf2 = 0.363636

    delta_ti1c = ti1 - tc
    delta_ti2c = ti2 - tc

    delta_tf1c = tf1 - tc
    delta_tf2c = tf2 - tc

    r1c_par = (r1f_par / delta_tf1c - r1i_par / delta_ti1c) / (1.0 / delta_tf1c - 1.0 / delta_ti1c)

    r2c_par = (r2f_par / delta_tf2c - r2i_par / delta_ti2c) / (1.0 / delta_tf2c - 1.0 / delta_ti2c)
    # add = 0.0;
    # if (r2c_par < 0 or r1c_par < 0):
    #    add =  r2c_par**2 + r1c_par**2

    r2c_perp = (r1i_perp - r1f_perp + (delta_tf1c / delta_ti2c) * r2i_perp - (delta_ti1c / delta_tf2c) * r2f_perp) / (
                (delta_tf1c / delta_ti2c) - (delta_ti1c / delta_tf2c))

    r1c_perp = r1i_perp - (delta_ti1c / delta_tf2c) * (r2f_perp - r2c_perp)

    r1c_perp_new = 0
    #    r1c_perp_new = (r2i_perp-r2f_perp+(delta_tf2c/delta_ti1c)*r1i_perp-(delta_ti2c/delta_tf1c)*r1f_perp)/((delta_tf2c/delta_ti1c)-(delta_ti2c/delta_tf1c))

    # if (r2c_perp < 0 or r1c_perp < 0):
    #    return 1e10

    delta_rc_perp = r2c_perp - r1c_perp

    delta_rc_par = r2c_par - r1c_par

    dsq_c = delta_rc_perp ** 2 + delta_rc_par ** 2

    # print("distance at collisions = ", sqrt(dsq_c))

    # print("perp diff, par diff, at collision = ", delta_rc_perp, delta_rc_par)

    vi1_perp = (r1i_perp - r1c_perp) / delta_ti1c
    vi2_perp = (r2i_perp - r2c_perp) / delta_ti2c
    vf1_perp = (r1f_perp - r1c_perp) / delta_tf1c
    vf2_perp = (r2f_perp - r2c_perp) / delta_tf2c

    vi1_par = (r1i_par - r1c_par) / delta_ti1c
    vi2_par = (r2i_par - r2c_par) / delta_ti2c
    vf1_par = (r1f_par - r1c_par) / delta_tf1c
    vf2_par = (r2f_par - r2c_par) / delta_tf2c

    return (delta_rc_perp - 0.1) ** 2 + delta_rc_par ** 2  # + add*1e-6


#seed(1007)
seed(1013)

# r1i = array([0.0,0.0])
# r2i = array([1.0,0.5])
# r1f = array([0.5,0.0])
# r2f = array([1.0,0.2])

# r1f = array([0.0,0.0])
# r2f = array([1.0,0.5])
# r1i = array([0.5,0.0])
# r2i = array([1.0,0.2])

# r1f = array([0.0,0.0])
# r2f = array([1.0,0.0])
# r1i = array([0.0,1.0])
# r2i = array([1.0,1.0])

#
r2f = array([0.0, 0.0])
r1f = array([0.0, 1.0])
r2i = array([1.0, 0.0])
r1i = array([1.0, 1.0])

r1f = array([0.0, 0.0])
r2f = array([0.0, 1.0])
r1i = array([1.0, 0.0])
r2i = array([1.0, 1.0])

r1i = array([0.100308, 0.441266])
r1f = array([0.650302, 0.611996])
r2i = array([0.732206, 0.036926])
r2f = array([0.068949, 0.960871])

x1 = array([0.975389, 0.555387, 0.461634])
y1 = array([0.885400, 0.535914, 0.321390])
x2 = array([0.472861, 0.472861, 0.134347])
y2 = array([0.479501, 0.479501, 0.166971])

r1i = array([0.975389, 0.885400])
r1f = array([0.461634, 0.321390])
r2i = array([0.472861, 479501])
r2f = array([0.134347, 0.166971])

# bad case:
x1 = array([0.213155, 0.515712, 0.508876])
y1 = array([0.492147, 0.659858, 0.612460])
x2 = array([0.588111, 0.588107, 0.610903])
y2 = array([0.728752, 0.728750, 0.755001])
r1i = array([0.213155, 0.492147])
r1f = array([0.508876, 0.612460])
r2i = array([0.588111, 0.728752])
r2f = array([0.610903, 0.755001])

# nhat = array([0.5,0.5])
# nhat = array([0.0,1.0])
nhat = array([1.0, 0.0])
nhat /= sqrt(nhat[0] ** 2 + nhat[1] ** 2)

rout = zeros((2, 2))

N = 100

out = zeros((N, N))

# tc = 0.706166
# theta = 7.308277
# tc = 0.197845
# theta = 3.741190
tc = 0.329571
theta = 7.043432

in_array = array([theta, tc])
print("min_fun (...) = ", min_fun(in_array, r1i, r2i, r1f, r2f))

# sys.exit()

for i in range(0, N):
    tc = 0.000001 + (1.0 / N) * (i + 1)
    for j in range(0, N):
        theta = (pi / N) * j
        nhat = array([cos(theta), sin(theta)])
        # out[i,j] = 1./coll(nhat,tc,r1i,r2i,r1f,r2f)
        out[i, j] = min_fun(array([theta, tc]), r1i, r2i, r1f, r2f)
        # print("out[",i,j,"]",out[i,j])

print("out[40,50] = ", out[40, 50], "50,40 = ", out[50, 40], "60,40 = ", out[60, 40])
# print("out", out[35:45,45:55])

# imshow(out)
max_i = argmax(1 / out)
# global nhat
i = max_i // N
j = max_i % N
tc = 0.000001 + (1.0 / N) * (i + 1)
theta = (pi / N) * j
nhat = array([cos(theta), sin(theta)])
print("max indices = ", i, j, "values = ", out[i, j], 1 / out[i, j], " theta = ", theta, "tc = ", tc)

found_min = False

maxti = max(ti1, ti2)
mintf = min(tf1, tf2)

theta = in_array[0] - 2 * pi
tc = in_array[1]
xin = array([theta, tc])
print("input value of dependent variables = ", xin, "opt func = ", min_fun(xin, r1i, r2i, r1f, r2f))
# sys.exit()
# for i in range(100):
for i in range(1):
    # theta = 2*pi*random()
    # tc = random()
    x0 = array([theta, tc])
    bnds = ((0.0, 2 * pi), (0.0, 1.0))
    # res = minimize(min_fun, x0, args=(r1i,r2i,r1f,r2f,),method='TNC',bounds=bnds)
    # res = minimize(min_fun, x0, args=(r1i,r2i,r1f,r2f,),method='CG',bounds=bnds)
    print("again: input value of dependent variables = ", x0, "opt func = ", min_fun(x0, r1i, r2i, r1f, r2f))
    res = minimize(min_fun, x0, args=(r1i, r2i, r1f, r2f,), method='Nelder-Mead', bounds=bnds)
    tc = res.x[1]
    # if res.fun < 1e-8:
    # if res.fun < 1e-8 and r1c_perp >= 0 and r2c_perp >= 0 and tc > maxti and tc < mintf:
    if res.fun < 1e-8 and tc > maxti and tc < mintf:
        print("in if res.fun = ", res.fun)
        sign = 1
        rout[0, 0] = nhat[0] * r1c_perp - sign * nhat[1] * r1c_par
        rout[0, 1] = sign * nhat[1] * r1c_perp + nhat[0] * r1c_par
        rout[1, 0] = nhat[0] * r2c_perp - sign * nhat[1] * r2c_par
        rout[1, 1] = sign * nhat[1] * r2c_perp + nhat[0] * r2c_par
        rc1 = rout[0, :]
        rc2 = rout[1, :]
        if does_intersect(rc1, rc2, r1i, r2i, r1f, r2f):
            continue
        print("res = ", res, "i = ", i)
        print("value of dependent variables = ", res.x, "opt func = ", min_fun(res.x, r1i, r2i, r1f, r2f), res.fun)
        found_min = True
        break

if not found_min:
    print("no minimum was found")
else:
    print("minimum was found")

theta = res.x[0]
tc = res.x[1]

x1 = array([r1i[0], rout[0][0], r1f[0]])
y1 = array([r1i[1], rout[0][1], r1f[1]])
x2 = array([r2i[0], rout[1][0], r2f[0]])
y2 = array([r2i[1], rout[1][1], r2f[1]])
print("dist = ", sqrt(sum((rout[1] - rout[0]) ** 2)))
print("r2c_perp  = ", r2c_perp, " r1c_perp  = ", r1c_perp)
print("r2c_par  = ", r2c_par, " r1c_par  = ", r1c_par)
print("r1c_perp_new = ", r1c_perp_new)
figure(0)
plot(x1, y1, x2, y2)
plot(array([0.0, nhat[0]]), array([0.0, nhat[1]]))

print("r1i = ", r1i, " r2i = ", r2i, "r1f = ", r1f, " r2f = ", r2f)
print("rc1 = ", rout[0, :], " rc2 = ", rout[1, :])

figure(1)

rc1 = rout[0, :]
rc2 = rout[1, :]
vcmi, vcmf = com_v(r1i, r2i, r1f, r2f, rc1, rc2)

print("vcmi = ", vcmi)
print("vcmf = ", vcmf)
print(" vi1_perp, vi2_perp, vf1_perp, vf2_perp = ", vi1_perp, vi2_perp, vf1_perp, vf2_perp)
print(" vi1_par, vi2_par, vf1_par, vf2_par = ", vi1_par, vi2_par, vf1_par, vf2_par)
print("does_intersect = ", does_intersect(rc1, rc2, r1i, r2i, r1f, r2f, ))

print("COM frame:\n\n")
# r1i -= vcmi*ti1
# r2i -= vcmi*ti2
# r1f -= vcmi*tf1
# r2f -= vcmi*tf2

# rout[0,:] -= vcmi[:]*tc
# rout[1,:] -= vcmi[:]*tc

r1i_p, r2i_p, r1f_p, r2f_p, rc1_p, rc2_p = com_r(r1i, r2i, r1f, r2f, rc1, rc2)

x1 = array([r1i_p[0], rc1_p[0], r1f_p[0]])
y1 = array([r1i_p[1], rc1_p[1], r1f_p[1]])
x2 = array([r2i_p[0], rc2_p[0], r2f_p[0]])
y2 = array([r2i_p[1], rc2_p[1], r2f_p[1]])

print("r1i_p = ", r1i_p, " r2i_p = ", r2i_p, "r1f_p = ", r1f_p, " r2f_p = ", r2f_p)
print("rc1_p = ", rc1_p, " rc2_p = ", rc2_p)

vcmi_p = 0.5 * ((rc1_p - r1i) / (tc - ti1) + (rc1_p - r2i) / (tc - ti2))
vcmf_p = 0.5 * ((rc1_p - r1f) / (tc - tf1) + (rc2_p - r2f) / (tc - tf2))

print("vcmi_p = ", vcmi_p)
print("vcmf_p = ", vcmf_p)

plot(x1, y1, x2, y2)

print("ti1,ti2,tf1,tf2,tc = ", ti1, ti2, tf1, tf2, tc)

print("r1c_perp-r1i_perp = ", r1c_perp - r1i_perp)

show()