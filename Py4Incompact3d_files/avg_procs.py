import csv
import math

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.use('TkAgg')

plt.ioff()

########################################################################################

#Calculating flow properties
print("FLOW PROPERTIES")
utau = 0.00812466799874915
nu = 1
retau = 180
reb = (np.exp(np.log(180/0.09) / 0.88) / 2) 
print(f"Friction Reynolds Number (Re_tau) {retau:.4f}")
print(f"Bulk Reynolds Number (Re_b) {reb:.4f}")
    
rho = 1  # density
nu = 1/reb  # kinematic viscosity
cf = 0.00812466799874915

tau_w = cf / 2
utau = np.sqrt(tau_w / rho)

print(f"Friction velocitty (u_tau) {utau:.4f}")
print(f"Kinematic viscosity {nu:.4f}")
print("Density (rho) is set to 1")

t = "80000"

#load data
umean = np.loadtxt('avg_umean_80000.dat', skiprows=1)
vmean = np.loadtxt('avg_vmean_80000.dat', skiprows=1)
wmean = np.loadtxt('avg_wmean_80000.dat', skiprows=1)

uumean = np.loadtxt('avg_uumean_80000.dat', skiprows=1)
vvmean = np.loadtxt('avg_vvmean_80000.dat', skiprows=1)
wwmean = np.loadtxt('avg_wwmean_80000.dat', skiprows=1)

uvmean = np.loadtxt('avg_uvmean_80000.dat', skiprows=1)


yp = np.loadtxt('yp.dat')
upyp_ref = np.loadtxt('upyp_reference.dat')
yp_ref = upyp_ref[1:,1]
up_ref = upyp_ref[1:,2]

utau_ref = 6.37309e-02
nu_ref = 3.50000e-04 
    
y_ref = yp_ref * nu_ref / utau_ref
u_ref = up_ref * utau_ref



########################################################################################


def main():
    
    uplus, yplus = velplus(umean, yp, utau, nu)

    diff = 1/(np.max(uplus)/np.max(up_ref))
    print("diff", diff)

    upyp_plot(uplus*diff, yplus, yp_ref, up_ref, t, "uplus")
    uy_plot(umean*diff, yp,  y_ref, u_ref, t)

    velvelplus(umean, umean, uumean, t, "$<u'u'>^+$", "upup")
    velvelplus(vmean, vmean, vvmean, t, "$<v'v'>^+$", "vpvp")
    velvelplus(wmean, wmean, wwmean, t, "$<w'w'>^+$", "wpwp")

    velvelplus(umean, vmean, uvmean, t, "$<u'v'>^+$", "upvp")
    

########################################################################################



def velplus(amean, yp, utau, nu):
    aplus = amean[1:(len(amean)//2 + 1)] / utau
    yplus = yp[1:(len(amean)//2 + 1)] * utau / nu

    return aplus, yplus

def velvelplus(a1_avg, a2_avg, aa_avg, t, varname, titlename):
     #Multiplication
    a1a2_avg = a1_avg*a2_avg
    
    print("For point 50...")
    print("a1", a1_avg[50])
    print("a2", a2_avg[50])
    
    print("a1a2", a1a2_avg[50])
    
    print("aa", aa_avg[50])
    
    #prime-prime calculation
    apap = (aa_avg - a1a2_avg)
    
    print ("<a1'a2'>", apap[50])
    
    apap_plus = apap / (utau**2)
    
    print("  ")
    print("u_tau", utau)
    print("   ")
    print("<a1'a2'>+ , max", max(apap_plus))
    
    yplus = yp[:(len(yp)//2 + 1)] * utau / nu
    
    
    #Plotting...
    # Plot <a1'_a2'>+ vs u+
    plt.figure()
    plt.plot(yplus, apap_plus[:(len(apap_plus)//2 + 1)])
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel(varname)
    plt.title('$y^+$ vs' + varname)
    plt.grid(True, which="both", ls="--")
    plt.savefig(titlename + t + ".png" , bbox_inches='tight')
    plt.clf()
    
    # Save data to file
    data = np.column_stack((yplus, apap_plus[:(len(apap_plus)//2 + 1)]))
    np.savetxt(titlename + t + '.dat', data, header='y+ primeprime+')

def upyp_plot(aplus, yplus, yp_ref, up_ref, t, titlename):
    plt.figure()
    plt.plot(yplus, aplus, label="data")
    plt.plot(yp_ref, up_ref, linestyle='dashed', color="red", label="reference")
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$u^+$')
    plt.title('$y^+$ vs $u^+$')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig('upyp_log'+t+'.png', bbox_inches='tight')
    plt.show()
    plt.clf()

    # Save data to file
    data = np.column_stack(aplus)
    np.savetxt(titlename + "_" + t + '.dat', data, header='values')

def uy_plot(aplus, yplus, yp_ref, up_ref, t):
    plt.figure()
    plt.plot(yplus, aplus, label="data")
    plt.plot(yp_ref, up_ref, linestyle='dashed', color="red", label="reference")
    plt.xlabel('$y^+$')
    plt.ylabel('$u^+$')
    plt.title('$y^+$ vs $u^+$')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig('upyp_real'+t+'.png', bbox_inches='tight')
    plt.clf()

if __name__ == "__main__":
    main()