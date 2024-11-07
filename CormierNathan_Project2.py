import numpy as np
import scipy as sp

# Parameters
m_bb = 0.145 # baseball mass [kg]
diam_bb = 0.074 # baseball diameter [m]
acc_g = -9.81 # gravitational acceleration [m/s^2]
atmDens = 1.2 # Air density [kg/m^3]
dragC = 0.35  # drag coefficient

# Computes cross sectional area of the baseball [m^2]
AreaX_bb = 0.25*np.pi*diam_bb**2

# part 1

# solve the projectile problem by writing a Python function that implements the Euler, Euler-Cromer and Midpoint methods 
# as described in your textbook, pp28-30. Required input should be the initial ball speed (in m/s), the launch angle (in degrees),
# the time-step used to numerically integrate the equations of motion (in seconds, called τ in your textbook), and which of the 
# three methods (Euler, Euler-Cromer or Midpoint) should be used for the integration. It would be useful to be able to turn off 
# air resistance; handle this with an argument to your function. Your routine should return the horizontal range of the ball (in metres).
# Be sure to look at Figure 2.4, p33, in your textbook in relation to computing the range.