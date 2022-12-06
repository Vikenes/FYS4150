#Please place all analysis here

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import plot as PLOT
import pyarma as pa

here = os.path.abspath(".")
binfile_path = here + "/../../output/binfiles/"


class Analysis:
    def __init__(self, filename="no_slit_arma_cube", title="No slit", total_time=0.008, time_step_size=2.5e-5, space_step_size=0.005, num_of_slits=0):
        U_arma = pa.cx_cube()
        U_arma.load(binfile_path+filename+".bin")
        self.U = np.asarray(U_arma)
        self.P = np.abs(self.U)**2
        self.title = title
        self.h = space_step_size
        self.T = total_time
        self.dt = time_step_size
        self.set_up()
       
    def set_params(self, **params):
        pass

    def set_up(self):
        x_points = np.arange(0, 1+self.h, self.h)
        y_points = np.arange(0, 1+self.h, self.h)
        self.x, self.y = np.meshgrid(x_points, y_points, sparse=True)
        self.t = np.arange(0, self.T+self.dt, self.dt)

    def animate(self, mp4name="animation", total_time=None):
        T = total_time or self.T
        PLOT.animate_probability_density(self.t[self.t<=T], self.P, title=self.title, mp4name=mp4name)

    def __str__(self):
        l = 50
        s = "-"*l + "\n" + self.title + "\n"
        s += "\n" + f" h = {self.h:5.2e}"
        s += "\n" + f"dt = {self.dt:5.2e}"
        s += "\n" + f" T = {self.T:5.2e}"
        s += "\n" + "-"*l
        return s

        


NOSLIT = Analysis()
DSLIT = Analysis("double_slit_arma_cube", "Double slit (1)", num_of_slits=2)
DSLIT2 = Analysis("double_slit_short_time_arma_cube", "Double slit (2)", num_of_slits=2)

NOSLIT.animate(mp4name="no_slit")
DSLIT.animate(mp4name="double_slit")
DSLIT2.animate(mp4name="double_slit2")
