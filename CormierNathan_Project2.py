import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import time

# part 1

# solve the projectile problem by writing a Python function that implements the Euler, Euler-Cromer and Midpoint methods 
# as described in your textbook, pp28-30. Required input should be the initial ball speed (in m/s), the launch angle (in degrees),
# the time-step used to numerically integrate the equations of motion (in seconds, called τ in your textbook), and which of the 
# three methods (Euler, Euler-Cromer or Midpoint) should be used for the integration. It would be useful to be able to turn off 
# air resistance; handle this with an argument to your function. Your routine should return the horizontal range of the ball (in metres).
# Be sure to look at Figure 2.4, p33, in your textbook in relation to computing the range.

# Projectile motion of baseball function

def projMotion(v_launch,ang_launch,tstep,method,AirResYN):
    '''Takes initial velocity[m/s], launch angle[deg], time step[s], solving method ('Euler' or 'Euler-Cromer' or 'Midpoint'), and air resistance toggle as parameters.
    Returns the horizonal range of the ball.'''

    ang0 = ang_launch
    tau = tstep
    stepLim = 6000 # massive step limit so I dont run out of position logging space (within reason)

    # baseball parameters
    # Parameters
    m_bb = 0.145 # baseball mass [kg]
    diam_bb = 0.074 # baseball diameter [m]
    acc_g = -9.81 # gravitational acceleration [m/s^2]
    atmDens = 1.2 # Air density [kg/m^3]
    dragC = 0.35  # drag coefficient
    hitHeight = 1 # ball is hit from starting height of 1m

    # Computes cross sectional area of the baseball [m^2]
    AreaX_bb = 0.25*np.pi*diam_bb**2

    # calculating the air constant for later use in a formula for cases with air resistance
    airConst = -0.5*dragC*atmDens*AreaX_bb/m_bb

    # Setting initial position and velocity
    r0 = [0,hitHeight] #x0, y0
    v0 = [np.cos(ang0*np.pi/180)*v_launch,np.sin(ang0*np.pi/180)*v_launch] #Vx0, Vy0
    
    # Setting initial values of time steppable r and v to initial values
    r = r0
    v = v0

    # x initialization for plotting
    xPos = np.empty(stepLim)

    # y initialization for plotting
    yPos = np.empty(stepLim)
    yPos[0] = hitHeight

    ######## Euler method section
    if method == 'Euler':
        # Perform numerical analysis with Euler's Method
        if AirResYN is True:
            # Do calculations with Air resistance
            acc = np.zeros(2)
            i = 1
            while r[1]>=0:
                acc[0] = airConst*abs(v[0])*v[0]           # air resistance (only acc on x)
                acc[1] = airConst*abs(v[1])*v[1] + acc_g   # air res and gravity (acc on y)

                r[0] = r[0] + tau*v[0]       # Euler's method step for position in x
                r[1] = r[1] + tau*v[1]       # Euler's method step for position in y
                v[0] = v[0] + tau*acc[0]     # Euler's method step for velocity in x
                v[1] = v[1] + tau*acc[1]     # Euler's method step for velocity in y
                xPos[i] = r[0]    # saving x pos at each step for plotting
                yPos[i] = r[1]    # saving y pos at each step for plotting
                i +=1

            return r[0],i-1,xPos,yPos
        
        elif AirResYN is False:
            # Do calculations w/o Air resistance
            i = 1
            while r[1]>=0:
                r[0] = r[0] + tau*v[0]       # Euler's method step for position in x
                r[1] = r[1] + tau*v[1]       # Euler's method step for position in y
                v[1] = v[1] + tau*acc_g      # Euler's method step for velocity in y
                xPos[i] = r[0]    # saving x pos at each step for plotting
                yPos[i] = r[1]    # saving y pos at each step for plotting
                i +=1

            return r[0],i-1,xPos,yPos
        
        else:
            return 'Input Variable for AirResYN was not True or False. Please try again.'
    


    ######## Euler-Cromer method section
    elif method == 'Euler-Cromer':
        # Perform numerical analysis with Euler-Cromers Method
        if AirResYN is True:
            # Do calculations with Air resistance
            acc = np.zeros(2)
            i = 1
            while r[1]>=0:
                acc[0] = airConst*abs(v[0])*v[0]           # air resistance (only acc on x)
                acc[1] = airConst*abs(v[1])*v[1] + acc_g   # air res and gravity (acc on y)

                v[0] = v[0] + tau*acc[0]     # Euler-Cromer step for velocity in x (same as Euler's)
                v[1] = v[1] + tau*acc[1]     # Euler-Cromer step for velocity in y (same as Euler's)
                r[0] = r[0] + tau*v[0]       # Euler-Cromer step for position in x (now using updated vx to find rx)
                r[1] = r[1] + tau*v[1]       # Euler-Cromer step for position in y (now using updated vy to find ry)
                xPos[i] = r[0]    # saving x pos at each step for plotting
                yPos[i] = r[1]    # saving y pos at each step for plotting
                i +=1
            
            return r[0],i-1,xPos,yPos
        
        elif AirResYN is False:
            # Do calculations w/o Air resistance
            i = 1
            while r[1]>=0:
                v[1] = v[1] + tau*acc_g      # Euler-Cromer step for velocity in y (same as Euler's)
                r[0] = r[0] + tau*v[0]       # Euler-Cromer step for position in x
                r[1] = r[1] + tau*v[1]       # Euler-Cromer step for position in y (now using updated vy to find ry)
                xPos[i] = r[0]    # saving x pos at each step for plotting
                yPos[i] = r[1]    # saving y pos at each step for plotting
                i +=1

            return r[0],i-1,xPos,yPos
        
        else:
            return 'Input Variable for AirResYN was not True or False. Please try again.'
    


    ######## Midpoint method section
    elif method == 'Midpoint':
        # Perform numerical analysis with Midpoint Method
        if AirResYN is True:
            # Do calculations with Air resistance
            acc = np.zeros(2)
            i = 1
            while r[1]>=0:
                acc[0] = airConst*abs(v[0])*v[0]                # air resistance (only acc on x)
                acc[1] = airConst*abs(v[1])*v[1] + acc_g        # air res and gravity (acc on y)

                r[0] = r[0] + tau*v[0] + 0.5*acc[0]*(tau**2)    # Midpoint step for position in x (now using vx and ax to find rx)
                r[1] = r[1] + tau*v[1] + 0.5*acc[1]*(tau**2)    # Midpoint step for position in y (now using vy and ay to find ry)
                v[0] = v[0] + tau*acc[0]                        # Midpoint step for velocity in x (same as Euler's)
                v[1] = v[1] + tau*acc[1]                        # Midpoint step for velocity in y (same as Euler's)
                xPos[i] = r[0]    # saving x pos at each step for plotting
                yPos[i] = r[1]    # saving y pos at each step for plotting
                i +=1
                
            return r[0],i-1,xPos,yPos
        
        elif AirResYN is False:
            # Do calculations w/o Air resistance
            i = 1
            while r[1]>=0:
                r[0] = r[0] + tau*v[0]                        # Midpoint step for position in x (same as Eulers)
                r[1] = r[1] + tau*v[1] +0.5*acc_g*(tau**2)    # Midpoint step for position in y (now using vy and ay to find ry)
                v[1] = v[1] + tau*acc_g                       # Midpoint step for velocity in y (same as Euler's)
                xPos[i] = r[0]    # saving x pos at each step for plotting
                yPos[i] = r[1]    # saving y pos at each step for plotting
                i +=1
                
            return r[0],i-1,xPos,yPos
        
        else:
            return 'Input Variable for AirResYN was not True or False. Please try again.'

# getting range, x positions, and y positions for recreation of Fig 2.3 for Euler Theory, Euler, Euler-Cromer, and Midpoint methods
rngT, stepsT, xplotT, yplotT = projMotion(50,45,0.1,'Euler',False)    # Theory (No air resistance)
rngE, stepsE, xplotE, yplotE = projMotion(50,45,0.1,'Euler',True)     # Euler with air resistance
rngEC, stepsEC, xplotEC, yplotEC = projMotion(50,45,0.1,'Euler-Cromer',True)     # Euler-Cromer with air resistance
rngMP, stepsMP, xplotMP, yplotMP = projMotion(50,45,0.1,'Midpoint',True)     # Midpoint with air resistance

# setting the plotting arrays to the values recorded for: 
xPT = xplotT[0:stepsT]     # Theoretical Euler xP
yPT = yplotT[0:stepsT]     # Theoretical Euler yP
xPE = xplotE[0:stepsE]     # Euler xP
yPE = yplotE[0:stepsE]     # Euler yP
xPEC = xplotEC[0:stepsEC]  # Euler-Cromer xP
yPEC = yplotEC[0:stepsEC]  # Euler-Cromer yP
xPMP = xplotMP[0:stepsMP]  # Midpoint xP
yPMP = yplotMP[0:stepsMP]  # Midpoint yP

# setting axis boundary lines as shown in figure 2.3
xG = np.array([0.,xplotT[stepsT]])
yG = np.array([0.,0.])

# marking ground for plotting
plt.plot(xPT,yPT,'-',xPE,yPE,'+',xPEC,yPEC,'<',xPMP,yPMP,'1',xG,yG,'r-')
plt.legend(['Theory (No AR)','Euler method (AR)','Euler-Cromer method (AR)','Midpoint method (AR)'])
plt.xlabel('Range (m)')
plt.ylabel('Height (m)')
plt.title('Projectile Motion')
plt.show()


# Part 2

# Determine the AB/HR-ratio for the proposed RDH. To do this, run your projectile program for a number of 
# simulated at-bats with the starting v and θ drawn from the random distributions given above. Remember 
# that a normally-distributed random variable of mean μ and standard deviation σ, N(μ, σ2), can be obtained 
# by scaling the output of the numpy.random.randn function, N(0, 1), via sigma * np.random.randn(...) + mu

# You can assume any projectile with a range of 400 feet or greater is a home run as this is a typical 
# distance between home plate and the centre field.

# defining mean exit speed, stdev of exit speed, mean launch angle, and stdev of launch angle
exitspd_mean_ms = 44.704
exitspd_stdev_ms = 6.7056
ang0_mean_deg = 45
ang0_stdev_deg = 10
hr_distance_ft = 400

# setting size of output random value arrays
randSize = 50

# generating array of random exit speeds and launch angles for testing
rand_v0 = np.random.normal(exitspd_mean_ms,exitspd_stdev_ms,randSize)
rand_ang0 = np.random.normal(ang0_mean_deg,ang0_stdev_deg,randSize)

# initializing empty arrays to be filled with range values in m and feet
range_out = np.zeros(randSize)
range_out_feet = np.zeros(randSize)
feet_per_m = 3.28084     # conversion value to turn range [m] into range [ft]
numHRs = 0
a,b,c = 1,2,3            # arbitrary variables to make extracting ranges easy

# determining ranges in feet for the number of pairs of randomly generated launch speed and angles
for j in range(randSize):
    # setting the range output for iteration j equal to the calculated range with
    # initial velocity from rand_v0[j], launch angle from rand_ang0[j], timestep 0.01 s,
    # using the Midpoint method with Air resistance considered
    range_out[j], a, b, c, = projMotion(rand_v0[j],rand_ang0[j],0.01,'Midpoint',True)

    #converting the calculated values of range to feet for Homerun evaluation
    range_out_feet[j] = range_out[j]*feet_per_m

    # if the range in feet is greater than 400, count 1 homerun
    if range_out_feet[j] >= hr_distance_ft:
        numHRs += 1

ABHR = randSize/numHRs
print('AB/HR ratio : ', np.round(ABHR,2))