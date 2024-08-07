import csv
import math

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from py4incompact3d.postprocess.postprocess import Postprocess
from scipy.integrate import simps, trapz
from py4incompact3d.tools.misc import avg_over_axis
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker


mpl.use('TkAgg')

plt.ioff()



#####################################################################################
#up vs yp filename
upyp_c = np.loadtxt('upyp43125.dat', skiprows=1)
upyp_m = np.loadtxt('upyp52500.dat', skiprows=1)
upyp_f = np.loadtxt('upyp90000.dat', skiprows=1)

yp_c = upyp_c[:,0]
up_c = upyp_c[:,1]

yp_m = upyp_m[:,0]
up_m = upyp_m[:,1]

yp_f = upyp_f[:,0]
up_f = upyp_f[:,1]

    
# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(yp_c, up_c, label='Coarse Grid')
plt.plot(yp_m, up_m, label='Medium Grid')
plt.plot(yp_f, up_f, label='Fine Grid')

# Add labels and title
plt.xlabel('y+')
plt.xscale('log')
plt.ylabel('u+')
plt.title('Plot of u+ vs y+ for Different Grids')
plt.legend()
plt.grid(True)

# Show the plot
plt.savefig('upyp.png', bbox_inches='tight')
plt.clf()
######################################################################################

#####################################################################################
#upuprime filename
upyp_c = np.loadtxt('upup43125.dat', skiprows=1)
upyp_m = np.loadtxt('upup52500.dat', skiprows=1)
upyp_f = np.loadtxt('upup90000.dat', skiprows=1)

yp_c = upyp_c[:,0]
upup_c = upyp_c[:,1]

yp_m = upyp_m[:,0]
upup_m = upyp_m[:,1]

yp_f = upyp_f[:,0]
upup_f = upyp_f[:,1]

    
# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(yp_c, upup_c, label='Coarse Grid')
plt.plot(yp_m, upup_m, label='Medium Grid')
plt.plot(yp_f, upup_f, label='Fine Grid')

# Add labels and title
plt.xlabel('$y^+$')
plt.xscale('log')
plt.ylabel('$<upup>^+$')
plt.title('$<upup>^+$')
plt.legend()
plt.grid(True)

# Show the plot
plt.savefig('upup.png', bbox_inches='tight')
plt.clf()
########################################################################################

#####################################################################################
#vpvprime filename
upyp_c = np.loadtxt('vpvp43125.dat', skiprows=1)
upyp_m = np.loadtxt('vpvp52500.dat', skiprows=1)
upyp_f = np.loadtxt('vpvp90000.dat', skiprows=1)

yp_c = upyp_c[:,0]
vpvp_c = upyp_c[:,1]

yp_m = upyp_m[:,0]
vpvp_m = upyp_m[:,1]

yp_f = upyp_f[:,0]
vpvp_f = upyp_f[:,1]

    
# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(yp_c, vpvp_c, label='Coarse Grid')
plt.plot(yp_m, vpvp_m, label='Medium Grid')
plt.plot(yp_f, vpvp_f, label='Fine Grid')

# Add labels and title
plt.xlabel('$y^+$')
plt.xscale('log')
plt.ylabel('$<vpvp>^+$')
plt.title('$<vpvp>^+$')
plt.legend()
plt.grid(True)

# Show the plot
plt.savefig('vpvp.png', bbox_inches='tight')
plt.clf()
########################################################################################

#####################################################################################
#wpwprime filename
upyp_c = np.loadtxt('wpwp43125.dat', skiprows=1)
upyp_m = np.loadtxt('wpwp52500.dat', skiprows=1)
upyp_f = np.loadtxt('wpwp90000.dat', skiprows=1)

yp_c = upyp_c[:,0]
wpwp_c = upyp_c[:,1]

yp_m = upyp_m[:,0]
wpwp_m = upyp_m[:,1]

yp_f = upyp_f[:,0]
wpwp_f = upyp_f[:,1]

    
# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(yp_c, wpwp_c, label='Coarse Grid')
plt.plot(yp_m, wpwp_m, label='Medium Grid')
plt.plot(yp_f, wpwp_f, label='Fine Grid')

# Add labels and title
plt.xlabel('$y^+$')
plt.xscale('log')
plt.ylabel('$<wpwp>^+$')
plt.title('$<wpwp>^+$')
plt.legend()
plt.grid(True)

# Show the plot
plt.savefig('wpwp.png', bbox_inches='tight')
plt.clf()
########################################################################################


