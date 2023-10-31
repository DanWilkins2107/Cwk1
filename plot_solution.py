"""
This script plots a surface plot of the Finite Volume
Solution from solution.dat. solution.dat is produced in
the correct format by OutputSolution() from
finite_volume.cpp
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

data = np.loadtxt("solution.dat")

ax = plt.figure().add_subplot(projection='3d')

vmin = data[:,-1].min()
vmax = data[:,-1].max()



for row in data:
  x = row[0:-3:2]
  y = row[1:-3:2]
  solution = row[-1]*np.ones(3)

  ax.plot_trisurf(x, y, solution, linewidth=0, antialiased=True,cmap=cm.plasma,vmin=vmin,vmax=vmax)

  x = [row[4],row[6],row[0]]
  y = [row[5],row[7],row[1]]
  solution = row[-1]*np.ones(3)

  ax.plot_trisurf(x, y, solution, linewidth=0, antialiased=True,cmap=cm.plasma,vmin=vmin,vmax=vmax)

ax.set_proj_type('ortho')
ax.view_init(elev=90, azim=-90, roll=0)

plt.show()

