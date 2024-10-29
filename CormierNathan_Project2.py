import numpy as np
import scipy as sp

# part 1

# solve the projectile problem by writing a Python function that implements the Euler, Euler-Cromer and Midpoint methods 
# as described in your textbook, pp28-30. Required input should be the initial ball speed (in m/s), the launch angle (in degrees),
# the time-step used to numerically integrate the equations of motion (in seconds, called τ in your textbook), and which of the 
# three methods (Euler, Euler-Cromer or Midpoint) should be used for the integration. It would be useful to be able to turn off 
# air resistance; handle this with an argument to your function. Your routine should return the horizontal range of the ball (in metres).
# Be sure to look at Figure 2.4, p33, in your textbook in relation to computing the range.

