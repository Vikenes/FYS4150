#Please place all analysis here

import numpy as np
import matplotlib.pyplot as plt
import os
import plot as PLOT
import pyarma as pa
import sys



binfile_path = os.path.abspath(".") + "/../../output/binfiles/"

class Analysis:
    def __init__(self, binary_filename, title, label, num_of_slits=0, total_time=None, default=1):
        U_arma = pa.cx_cube()
        U_arma.load(binfile_path+binary_filename+".bin")
        self.U = np.asarray(U_arma)
        self.P = np.abs(self.U)**2

        self.Nt, self.M, _ = np.shape(self.P)
        self.title = title
        self.label = label
        self.num_of_slits = num_of_slits
        self.default_setup(default, total_time=total_time)
        self.set_grid()

        self.wall_y = self.get_walls_ycoords()

    def set_params(self, xc=(0.25,0.5), sigma=(0.05,0.05), p=(200,0), v0=1e10, T=0.002):
        self.xc = xc
        self.sigma = sigma
        self.p = p
        self.v0 = v0
        if self.num_of_slits == 0:
            self.v0 = 0
    
    def default_setup(self, which="first", num_of_slits=None, total_time=None):
        self.num_of_slits = num_of_slits or self.num_of_slits
        if str(which) in ["1", "first", "normal", "usual", "long"]:
            self.T = total_time or 0.008
            self.set_params()
        elif str(which) in ["2", "second", "other", "short"]:
            self.T = total_time or 0.002
            self.set_params(sigma=(0.05, 0.20))
        else:
            print("No default params set.")

    
    def get_walls_ycoords(self, num_of_slits=None):
        Ns = num_of_slits or self.num_of_slits
        if Ns == 0:
            yc = []
        else:
            Nw = Ns + 1 # num of walls
            h = 0.05 # height of wall
            a = 0.05 # aperture
            yc = np.zeros(Nw)
            sep = (a+h)
            yc_first = 0.5 - sep*(Ns/2)
            for j in range(Nw):
                yc[j] = yc_first + j*sep
        return yc

    def get_total_probability(self):
        # Wack method, but struggeled with dimesions
        # Feel free to fix
        Ptot = []
        for j in range(self.Nt):
            Ptot.append(np.sum(self.P[j]))
        return Ptot

        
    def set_grid(self):
        self.x = np.linspace(0, 1, self.M)
        self.y = np.linspace(0, 1, self.M)
        self.h = self.x[1]-self.x[0]
        self.t = np.linspace(0, self.T, self.Nt)
        self.dt = self.t[1]-self.t[0]

    def animate(self, mp4name=None, total_time=None):
        T = total_time or self.T
        mp4name =  mp4name or self.label + "_anim"
        PLOT.animate_probability_density(self.t[self.t<=T], self.P, title=self.title, wall_y=self.wall_y, mp4name=mp4name)

    def snapshots(self, time_points, pdfnames=[None, None, None]):
        num_maps = len(time_points)
        t_points = np.zeros(num_maps)
        P_points = np.zeros((num_maps, self.M, self.M))
        ReU_points = np.zeros((num_maps, self.M, self.M))
        ImU_points = np.zeros((num_maps, self.M, self.M))

        for j, time in enumerate(time_points):
            diff = np.abs(self.t-time)
            idx = np.argmin(diff)
            t_points[j] = self.t[idx]
            P_points[j] = self.P[idx]
            ReU_points[j] = np.real(self.U[idx])
            ImU_points[j] = np.imag(self.U[idx])

        yc = self.wall_y

        ### Prob. density:
        pdfname = pdfnames[0] or self.label + "_snapshots_P"
        PLOT.snapshot_probability_density(t_points, P_points, wall_y=yc, pdfname=pdfname)
        ### Real part of U:
        pdfname = pdfnames[1] or self.label + "_snapshots_ReU"
        PLOT.snapshot_real_wavefunction(t_points, ReU_points,  wall_y=yc, pdfname=pdfname)
        ### Imag. part of U:
        pdfname = pdfnames[2] or self.label + "_snapshots_ImU"
        PLOT.snapshot_imaginary_wavefunction(t_points, ImU_points, wall_y=yc, pdfname=pdfname)

    def deviation(self, pdfname=None, others=[]):
        experiments = [self] + others
        label = ""
        Ptot = []
        titles = []
        for exp in experiments:
            label += exp.label + "_"
            Ptot.append(exp.get_total_probability())
            titles.append(exp.title)
        
        pdfname = pdfname or label + "Ptot_deviation"
        PLOT.total_probability_deviation(self.t, Ptot, labels=titles, pdfname=pdfname)


    def probability_vertical_screen(self, time_point=0.002, horisontal_point=0.8, pdfname=None, label="auto"):
        idx_x = np.argmin(np.abs(self.x-horisontal_point))
        x = self.x[idx_x]
        idx_t = np.argmin(np.abs(self.t-time_point))
        t = self.t[idx_t]
        p = self.P[idx_t, idx_x, :]
        if label=="auto":
            label = self.title
        pdfname = pdfname or self.label + "_P_along_screen"
        PLOT.probability_density_along_screen(self.y, p/np.sum(p), label=label, pdfname=pdfname)

    def __str__(self):
        l = 50
        s = "-"*l + "\n" + self.title + "\n"
        s += f"(labeled: '{self.label}')" + "\n"
        s += f"\n> System configuration:"
        s += "\n   " + f" h = {self.h:5.1e}"
        s += "\n   " + f" M = {self.M:5.0f}"
        s += "\n   " + f"dt = {self.dt:5.1e}"
        s += "\n   " + f" T = {self.T:5.1e}"
        s += "\n> Slits:" 
        s += "\n   " + f"Nslits = {self.num_of_slits:5}"
        s += "\n   " + f"    v0 = {self.v0:5.1e}"
        s += "\n> Initial wave packet:" 
        s += "\n   " + f"xc = ({self.xc[0]:5.2f}, {self.xc[1]:5.2f})"
        s += "\n   " + f" ?? = ({self.sigma[0]:5.2f}, {self.sigma[1]:5.2f})"
        s += "\n   " + f" p = ({self.p[0]:5.1f}, {self.p[1]:5.1f})"
        s += "\n" + "-"*l + "\n"
        return s



### 1. No slits
NOSLITS = Analysis("NS_arma_cube", "No slits", label="NS")
### 2. Double-slit (1)
DSLIT1 = Analysis("DS1_arma_cube", "Double-slit (1)", label="DS1", num_of_slits=2)
DSLIT1.set_params(sigma=(0.05, 0.10))
### 3. Double-slit (2)
DSLIT2 = Analysis("DS2_arma_cube", "Double-slit (2)", label="DS2", num_of_slits=2, default=2)
### 4. Single-slit
SSLIT = Analysis("SS_arma_cube", "Single-slit", label="SS", num_of_slits=1, total_time=0.004, default=2)
### 5. Triple-slit
TSLIT = Analysis("TS_arma_cube", "Triple-slit", label="TS", num_of_slits=3, total_time=0.004, default=2)

print(NOSLITS)
print(DSLIT1)
print(DSLIT2)
print(SSLIT)
print(TSLIT)


which = {"NS":"NOSLITS", "DS1":"DSLIT1", "DS2":"DSLIT2", "SS":"SSLIT", "TS":"TSLIT", "ALL":"ALL"}
anim = ["animate", "anim", "video", "animations"]
fig = ["plot", "plots", "figures", "fig", "figure"]



try:
    arg1 = sys.argv[1].upper().strip()
    arg2 = sys.argv[2].lower().strip()

    ### Animate
    if arg2 in anim:
        if arg1 in which.keys():
            if arg1 != "ALL":
                eval(f"{which[arg1]}.animate()")
            else:
                NOSLITS.animate()
                DSLIT1.animate()
                DSLIT2.animate()
                SSLIT.animate()
                TSLIT.animate()
            PLOT.show_all()
        
    ### Plot
    elif arg2 in fig:
        if arg1 in which.keys():
            if arg1 in ["NS", "DS1"]:
                NOSLITS.deviation(others=[DSLIT1])
            elif arg1 in ["DS2", "SS", "TS"]:
                eval(f"{which[arg1]}.snapshots((0, 0.001, 0.002))")
                # eval(f"{which[arg1]}.probability_vertical_screen()")
            elif arg1 == "ALL":
                NOSLITS.deviation(others=[DSLIT1])
                DSLIT2.snapshots((0, 0.001, 0.002))
                SSLIT.probability_vertical_screen()
                TSLIT.probability_vertical_screen()
            PLOT.show_all()

    else:
        print("Please provide a valid command line argument.")
        print("It can be one of the following: {NS, DS1, DS2, SS, TS, ALL}")

except IndexError:
    pass



