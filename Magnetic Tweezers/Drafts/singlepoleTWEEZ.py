
"""
The single pole calculations give an indication of the dimensions of the required poles and beads.
"""


import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3)
plt.rcParams["figure.figsize"] = (30,20)

pico = (10**(-12))
nano = (10**(-9))
micro = (10**(-6))


"""
INPUTS
"""
bead_diameter = 2.8 *micro
viscosity = (8.90 * (10**(-4))) # Water Pa·s
velocity = (10)*micro
B_gradient = 30000 #T/m
M_bead = 11.5 #KA/m







M = 1430 #kA/m   magnitization of cobalt (pole material)
perm0 = (4*np.pi)*(10**(-7)) #vacuum permiabiltiy

Z = np.linspace(0,1000*micro,10000) #distance from pole

R = 1000*micro #
beta = Z/(R**2) #whole pole surface, described by the quadratic equation

'''
expression for the magnetic field outside the magnetic material and along the paraboloid axis (R=0):
'''
Hz = M/(4*beta*Z+1)


ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel("distance z away from parabolic axis(m)")
ax[0].set_ylabel("H field (T)")
ax[0].grid()
ax[0].plot(Z,Hz)


'''
expression for the gradient in the magnetic flux density outside of the magnetic pole (permeability μ=μ0)
'''
dBdz = (4*perm0*M*beta)/((4*beta*Z+1)**2)

'''
maximum attainable gradient for an optimal diameter relative to the distance from the tip
'''
beta_max  = (0.25)*Z  #optimum curvature for a given distance, follows from the condition ∂Fz/∂β=0 which yields
dBdz_max = (perm0*M)/4*Z




ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel("distance z away from parabolic axis(m)")
ax[1].set_ylabel("magnetic gradient (kT/m)")

ax[1].set_xlim(.1*micro,1000*micro)
ax[1].set_ylim(.01,10000)

ax[1].grid()
ax[1].plot(Z,dBdz)
ax[1].plot(Z,dBdz_max)




'''
maximum magnetic force on a uniform magnetized magnetic spherical particle (magnetized in the z-direction) 
'''


bead_diameter_lst = np.linspace(0,1000*nano,2000)

m_bead = (1/6)*np.pi*(bead_diameter_lst**3)*M_bead #saturation magnitzation  => magnetic moment = volume * magnetiziation

F_mag_single_pole = B_gradient *((1/6)*np.pi*((bead_diameter)**3)*M_bead)

print("F_mag_single_pole = ",F_mag_single_pole,"pN")

db = 30 #kT/m
F_mag =  (db * m_bead)




ax[2].set_yscale('log')
ax[2].set_xlabel("bead diameter (m)")
ax[2].set_ylabel("force(N)")

ax[2].set_xlim(0*pico,1000*nano)
ax[2].set_ylim(1*pico,10000*pico)

ax[2].grid()
ax[2].plot(bead_diameter_lst,F_mag)
plt.show()


'''
viscous force/drag force on microsphere from creeping flow/stokes flow
neglecting boynacy and gravity
'''

r = bead_diameter/2


F_drag = 6*np.pi*viscosity*r*velocity * (1/pico)

print("F_drag on sphere = ",F_drag, "pN")