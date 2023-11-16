import numpy as np

#vars are written as -> var time + particle num + component + pos in array
class particled:

    def __init__(self, xi, xf, ncol, pnum, d, tarray):
        self.d = d
        self.xi = xi
        self.xf = xf
        # pos is 1-D array that keeps track of each position component
        self.pos = np.empty(d, dtype=object)
        for x_comp in range(d):
            self.pos[x_comp] = str(xi[x_comp])
        #print(self.pos)
        self.ncol = ncol
        self.pnum = pnum

        self.tarray = np.copy(tarray)
        self.tindex = 0

        # varr is a 2-D array of shape r = num collisions and c = d
        self.varr = np.empty([(ncol + 1), d], dtype=object)
        for i in range(ncol + 1):
            for v_comp in range(d):
                if v_comp == 0:
                    self.varr[i, v_comp] = "i" + str(pnum) + str(i)
                elif v_comp == 1:
                    self.varr[i, v_comp] = "j" + str(pnum) + str(i)
                else:
                    self.varr[i, v_comp] = "k" + str(pnum) + str(i)
                #form of comp - pnum - col num
                #earlier: given in form of particle num - component number - collision number of the particle

        self.vindex = 0
        self.cindex = 0

        self.parr = np.empty([(ncol + 2), d], dtype=object)
        for i in range(ncol + 1):
            if i == 0:
                for x_comp in range(d):
                    self.parr[i, x_comp] = xi[x_comp]
            else:
                for x_comp in range(d):
                    self.parr[i, x_comp] = "x" + str(pnum) + str(i) + str(x_comp + 1)
        for x_comp in range(d):
            self.parr[-1, x_comp] = xf[x_comp]
        #print(self.parr)

        self.pindex = 0

        self.col_arr = np.empty([ncol, d], dtype=object)

    def print(self):
        print("Particle " + str(self.pnum) + ": \n")
        print("Pos arr: ")
        print(self.parr)
        print("Vel arr: ")
        print(self.varr)
        print("Col arr: ")
        print(self.col_arr)
        print("time arr: ")
        print(self.tarray)
        print("\n")


    def inc_vindex(self):
        self.vindex += 1

    def reset_vindex(self):
        self.vindex = 0

    def inc_pindex(self):
        self.pindex += 1

    def reset_pindex(self):
        self.pindex = 0

    def inc_tindex(self):
        self.tindex += 1

    def reset_tindex(self):
        self.tindex = 0

    def update_cindex(self, c_num):
        self.cindex = c_num

    def update_final_time(self, new_T):
        self.tarray[-1] = new_T

    def reset_pos(self):
        self.pos = np.empty(self.d, dtype=object)
        for x_comp in range(self.d):
            self.pos[x_comp] = str(self.xi[x_comp])