# Import libraries relevant to the program

import numpy as np
import math
import time
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import animation
import os
from matplotlib.animation import FFMpegWriter
import ffmpy
from mpl_toolkits.mplot3d import Axes3D

# Function to replicate an array 'count' times along its x-axis


def repeat(arr, count):
    return np.repeat(arr[np.newaxis, ...], count, axis=0)# Import libraries relevant
    # to the program

import numpy as np
import math
import time
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import animation
import os
from matplotlib.animation import FFMpegWriter
import ffmpy
from mpl_toolkits.mplot3d import Axes3D

# Function to replicate an array 'count' times along its x-axis


def repeat(arr, count):
    return np.repeat(arr[np.newaxis, ...], count, axis=0)

# Function to rearrange a continuous list of data into the form
# [blackholeindex, number of bodies, number of dimensions (=3)]


def reshapelist(longlist, num):
    reshaped = np.reshape(longlist, (numblackholes, num, dimensions))
    return reshaped

# Function to create an array of the same shape as the dynamical data, whose entries are
# the relevant mass values.


def massarray(sequence, num):
    m = np.repeat(sequence, 3*num)
    n = np.reshape(m, (len(sequence), num, 3))
    return n

# Function which takes in initial conditions in a list of the form
# [x0, y0, z0, x1, ...zi, vx0, vy0, vz0, vx1, ..., vzi]
# and returns accelerations in a list of the form [ax0, ay0, az0, ..., azi]. These are
# calculated for a point mass. The intermediate array 'denominatorlimited' sets a
# minimum distance between two masses, to remove problems arising from singularities.


def accelerationcalculator(listofconds):

    reshapedconds = np.append(repeat(listofconds[0:3], bodiespergalaxy), repeat(listofconds[3:6], bodiespergalaxy))
    rj = np.reshape(reshapedconds, (numblackholes, bodiespergalaxy, 3))

    allposns = listofconds[0:3 * bodiespergalaxy]
    doublepositions = np.append(allposns, allposns)
    ri = np.reshape(doublepositions, (numblackholes, bodiespergalaxy, 3))

    difference = rj-ri
    numerator = difference * massarray(masses, bodiespergalaxy)
    sums = np.sum(difference**2, 2)
    denominator = sums**1.5
    denominatorlimited = np.where(denominator < 0.1, 0.1, denominator)
    denominatorarray = (reshapelist(np.repeat(denominatorlimited, 3), bodiespergalaxy))

    a_ik = np.sum((numerator/denominatorarray), 0)
    return np.ravel(a_ik)


# Function called by 'odeint' to evaluate the derivatives
# of the differential equation being investigated.


def derivatives(absolutevariables, t):
    positions = absolutevariables[:dimensions * bodiespergalaxy]
    velocities = absolutevariables[dimensions * bodiespergalaxy:]
    accelerations = accelerationcalculator(positions)
    dynamics = np.append(velocities, accelerations)
    return dynamics

# Function which takes in the initial conditions, as output by the
# 'initialconditions' function, and solves the equations of motion of the system using
# the odeint routine. The output is a NumPy array of the form
# [timestep, massnumber, dimension]


def solveeqn(initconds, time):
    dynamicalvars = scipy.integrate.odeint(derivatives, initconds, time, mxstep=1000000)
    return np.array(dynamicalvars)

# Alternative to solveeqn which uses the ODE solver, 'solve_ivp',  used to compare
# the odeint routine to various others- see report section 3.1.2


def alternativesolver(initconds, time):
    dynamicalvars = scipy.integrate.solve_ivp(derivatives, (0, tsteps),
                                              initconds, method ='LSODA',
                                              t_eval=t, min_step=0.01)
    return np.array(dynamicalvars.y).transpose()

# Function to determine initial conditions necessary for a parabolic orbit, as a
# function of initial separation between the galaxies and the distance of closest
# approach between them. Output values are initial positions, initial velocities.
# Calculations are made in the COM frame.


def parabolicorbit(initial_separation, closest_approach):
    startpoints = np.array([[initial_separation * masses[1]/np.sum(masses), 0, 0],
                            [-initial_separation * masses[0]/np.sum(masses), 0, 0]])
    launchangle = np.arcsin(np.sqrt(closest_approach*masses[1] / (np.linalg.norm(startpoints[0]) * np.sum(masses))))
    launchvelocity = -np.sqrt(2*masses[1]**2/(np.sum(masses) * initial_separation))
    v1 = np.array([launchvelocity * np.cos(launchangle), launchvelocity * np.sin(launchangle), 0])
    v2 = -masses[0]/ masses[1] * v1
    return startpoints, np.array([v1, v2])

# Function to determine initial conditions necessary for a circular orbit, as a
# function of separation between galaxies. Calculations are made in the COM frame.


def circularorbit(separation):
    startpoints = np.array([[-separation * masses[1]/np.sum(masses), 0, 0],
                            [separation * masses[0]/np.sum(masses), 0, 0]])
    launchvels = np.array([[0, -np.sqrt((separation * masses[1]/
                                         np.sum(masses))*masses[1]/(separation**2)), 0],
                           [0, np.sqrt((separation * masses[0]/
                                        np.sum(masses))*masses[0]/(separation**2)), 0]])
    return startpoints, launchvels


# Choose parameters for the simulation. The total number of stars per galaxy is given by
# (starsperrun * runs) (explained in report section 3.1.2.) 'radii' is the characteristic
# e-folding length of the radial distribution. 'offset' is the distance away from the star
# where the radial distribution begins (ie. no stars within r<offset.) 'direction' refers
# to the rotation direction of the galaxies- 1 for clockwise, 0 for anticlockwise.
# 'tsteps' is the number of timesteps taken, 't' is the time array for which the simulation
# runs. 'centerpositions' and 'velocitypositions' are the initial positions and velocities
# of the galaxies respectively.

starsperrun = 100
runs = 50
masses = [100, 100]
radii = [12, 12]
offset = 5
direction = [0, 0]

tsteps = 1000
t = np.array(np.linspace(0, 275, tsteps))
bodiespergalaxy = (starsperrun + 2)
numblackholes = 2
dimensions = 3
thickness = [0.1, 0.1]

centerpositions = parabolicorbit(120, 70)[0]
centervelocities = parabolicorbit(120, 70)[1]
print('number of stars per galaxy is ', starsperrun * runs)


# Calculate the initial positions necessary for test masses to execute
# circular orbits about their respective black holes. Output is of the form
# [blackholepositions, starpositions, blackholevelocities, starvelocities].


def initialconditions(blackholeindex):
    radii1 = np.random.exponential(radii[blackholeindex], starsperrun) + offset
    radii1structured = np.reshape(np.repeat(radii1, dimensions), (starsperrun, dimensions))
    table1 = np.reshape(np.random.uniform(0, 2 * np.pi, starsperrun * dimensions), (starsperrun, dimensions))
    trig1 = table1 * 0
    trig1[:, 0] = np.cos(table1[:, 0])
    trig1[:, 1] = np.sin(table1[:, 0])
    positions1 = trig1 * radii1structured
    positions1[:, 2] = np.random.normal(0, thickness[blackholeindex], starsperrun)
    velocities1 = table1 * 0
    velocities1[:, 0] = np.sqrt(masses[blackholeindex] / radii1structured[:, 0]) * -positions1[:, 1] / radii1structured[:, 0] \
                        * (-1) ** direction[blackholeindex]
    velocities1[:, 1] = np.sqrt(masses[blackholeindex] / radii1structured[:, 0]) * positions1[:, 0] / radii1structured[:, 0] \
                        * (-1) ** direction[blackholeindex]
    velocities1[:, 2] = 0
    abspositions1 = positions1 + centerpositions[blackholeindex]
    absvelocities1 = velocities1 + centervelocities[blackholeindex]
    return np.concatenate((centerpositions[0], centerpositions[1], np.ravel(abspositions1),
                           centervelocities[0], centervelocities[1], np.ravel(absvelocities1)))

# Function to act as a timer to test the time taken to solve the differential equations;
# refer to report section 3.1.2 for more details. The system is solved 'repeats' times;
# the output is the average solving time and the standard deviation of the results.


def timetest(repeats):
    times = np.zeros((repeats))
    totaltime = 0
    for counter in range(repeats):
        t0 = time.time()
        testarray1 = solveeqn(initialconditions(0), t)
        testarray2 = solveeqn(initialconditions(1), t)
        t1 = time.time()
        times[counter] = t1-t0
    return np.average(times), np.std(times)

# Code to solve the differential equations. The positions and velocities
# are stored in separate arrays.


t0 = time.time()
gal1 = solveeqn(initialconditions(0), t)
gal2 = solveeqn(initialconditions(1), t)

positions1 = np.split(gal1, 2, axis=1)[0]
velocities1 = np.split(gal1, 2, axis=1)[1]
positions2 = np.split(gal2, 2, axis=1)[0]
velocities2 = np.split(gal2, 2, axis=1)[1]
print('Finished loop 1 /', runs, '\n Total time/s =', time.time() - t0)

for counter in range(runs-1):
    gal1 = solveeqn(initialconditions(0), t)
    gal2 = solveeqn(initialconditions(1), t)
    positions1 = np.concatenate((positions1, np.split(gal1, 2, axis=1)[0]), axis=1)
    velocities1 = np.concatenate((velocities1, np.split(gal1, 2, axis=1)[1]), axis=1)
    positions2 = np.concatenate((positions2, np.split(gal2, 2, axis=1)[0]), axis=1)
    velocities2 = np.concatenate((velocities2, np.split(gal2, 2, axis=1)[1]), axis=1)
    print('Finished loop', counter + 2, '/', runs, '\n Total time/s =', time.time()-t0)

# Data is rearranged into the form necessary for animation and analysis.
# In 'coords1' and 'coords 2', position data is stored in the form
# [timestep, massindex, x/y/z](ie. the y position of the 70th mass at the 347th timestep
# is [347, 70, 1]. The format of 'vels1' and 'vels2' is entirely equivalent.


coords1 = np.reshape(np.ravel(positions1), (len(t), bodiespergalaxy * runs, 3))
coords2 = np.reshape(np.ravel(positions2), (len(t), bodiespergalaxy * runs, 3))
vels1 = np.reshape(np.ravel(velocities1), (len(t), bodiespergalaxy * runs, 3))
vels2 = np.reshape(np.ravel(velocities2), (len(t), bodiespergalaxy * runs, 3))

# Animate the motion of the masses in three dimensions.
# In the final lines of the function, the resulting plot can be saved as a .mp4
# file or plotted using matplotlib.


def animate3d():
    plt.style.use('dark_background')
    def update_graph(frame):
        line1xs = []
        line1ys = []
        line1zs = []
        line2xs = []
        line2ys = []
        line2zs = []

        for m in range(bodiespergalaxy*runs):
            line1xs.append(coords1[frame][m][0])
            line1ys.append(coords1[frame][m][1])
            line1zs.append(coords1[frame][m][2])
            line2xs.append(coords2[frame][m][0])
            line2ys.append(coords2[frame][m][1])
            line2zs.append(coords2[frame][m][2])

        galaxy1a.set_data(line1xs, line1ys)
        galaxy1a.set_3d_properties(line1zs)
        galaxy2a.set_data(line2xs, line2ys)
        galaxy2a.set_3d_properties(line2zs)
        center1.set_data(coords1[frame][0][0], coords1[frame][0][1])
        center1.set_3d_properties(coords1[frame][0][2])
        center2.set_data(coords2[frame][1][0], coords2[frame][1][1])
        center2.set_3d_properties(coords2[frame][1][2])

        return galaxy1a, galaxy2a, center1, center2


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d', xlim=(-190, 190), ylim=(-190, 190), zlim=(-50, 50))
    ax.view_init(elev=90, azim=60)

    ax.axis('off')
    galaxy1a, = ax.plot([], [], [], linestyle="", marker=",", color='w')
    galaxy2a, = ax.plot([], [], [], linestyle="", marker=",", color='y')
    center1, = ax.plot([], [], [], linestyle="", marker=",", color='r')
    center2, = ax.plot([], [], [], linestyle="", marker=",", color='r')

    ani2 = animation.FuncAnimation(fig, update_graph, tsteps, interval=40, blit=True)
    ani2.save('Massimbalance.mp4', dpi=300)
    # ax.view_init(elev=25, azim=60)
    # ani3 = animation.FuncAnimation(fig, update_graph, tmax, interval=40, blit=True)
    # ani3.save('zimpactparam2.mp4', dpi=300)
    # plt.show()
    return

############################### Analyse final data #################################

# Change style of subsequent plots


plt.style.use('classic')
plt.rcParams.update({'font.size': 32})


# Function to pick out position  or velocity coordinates of a given particle for
# all time steps. e.g. onemass(vels1, 0) returns the velocity of the first black hole
# for all time steps.


def onemass(data, massindex):
    variable = np.array([])
    if massindex >= 2 + starsperrun * runs:
        print('choose mass in appropriate range')
    else: variable = data[:, massindex, :]
    return variable


# Function to pick out speed of a given particle for all time steps.


def massspeed(velocitydata, massindex):
    speeds = np.linalg.norm(velocitydata, axis=2)
    return speeds[:, massindex]

# Function to find the positions of the masses relative to a central mass. blackholeindex may be 0 or 1 to
# denote which black hole the position is calculated relative to.


def relativepositions(simulationpositiondata, blackholeindex):
    relposns = simulationpositiondata - np.reshape(np.ravel(np.tile(onemass(simulationpositiondata, blackholeindex),
                                                                   bodiespergalaxy*runs)), (tsteps, bodiespergalaxy * runs, 3))
    return relposns


# Function to find 3D distance between all test masses and a given black hole
# for all time steps.


def relativeradii3d(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    dists = np.linalg.norm(coorddifference, axis=2)
    return dists


# Function to find 2D distance (in xy-plane) between all test masses
# and a given black hole for all time steps.


def relativeradii2d(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    coorddifference[:, :, 2]= 0
    dists = np.linalg.norm(coorddifference, axis=2)
    return dists

# Function to find 2D distance (in xz-plane) between all test masses
# and a given black hole for all time steps.


def relativeradiixz(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    coorddifference[:, :, 1]= 0
    dists = np.linalg.norm(coorddifference, axis=2)
    return dists


# Function to find azimuthal angle between all test masses
# and a given black hole for all time steps.


def azimuthalangle(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    angles = np.arctan(coorddifference[:, :, 1] / (coorddifference[:, :, 0]+0.1))
    return angles

# Function to find polar angle between all test masses
# and a given black hole for all time steps.


def polarangle(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    angles = np.arctan(coorddifference[:, :, 0]/(0.1+coorddifference[:, :, 2]))
    return angles


# Function to plot 2d radial distribution of given galaxy at timestep t


def radialdistributionplot(testmassdata, blackholeindex, time):
    histogram = plt.figure()
    ax = histogram.add_subplot(111)
    ax.set_xlabel('Radius/Arbitrary units')
    ax.set_ylabel('Frequency')
    plt.title('Radial distribution of stars')
    plt.ticklabel_format(axis='both', style='plain')
    plt.hist(relativeradii2d(testmassdata, blackholeindex)[time], 100)
    plt.show()
    return


# Function to plot angular distribution of stars in given galaxy at timestep t


def angulardistributionplot(testmassdata, blackholeindex, time):
    histogram = plt.figure()
    ax = histogram.add_subplot(111)
    ax.set_xlabel('Azimuthal angle/rad')
    ax.set_ylabel('Frequency')
    plt.title('Angular distribution of stars')
    plt.ticklabel_format(axis='both', style='plain')
    plt.hist(azimuthalangle(testmassdata, blackholeindex)[time], 100)
    plt.show()
    return


# Function to plot 2D scatter showing angular distribution
# in xy-plane against xy-radial distribution


def radialangularplot(testmassdata, blackholeindex, time):
    plot = plt.figure(figsize=(15, 10))
    ax = plot.add_subplot(111)
    ax.set_xlabel('Azimuthal angle/rad')
    ax.set_ylabel('Radius/arbitrary units')
    ax.set_xlim([-np.pi/2, np.pi/2])
    #plt.title('Radial distribution against angular distribution of stars')
    plt.ticklabel_format(axis='both', style='plain')
    ax.plot(azimuthalangle(testmassdata, blackholeindex)[time],
               relativeradii2d(testmassdata, blackholeindex)[time], 'k.', markersize=5)
    i = 0
    while os.path.exists('{}{:d}.png'.format('radang', i)):
        i += 1
    plt.savefig('{}{:d}.png'.format('radang', i), dpi=170)
    #plt.show()
    return

# Function to plot 2D scatter showing angular distribution
# in xz-plane against 3d radial distribution


def heightpolarplot(data, blackholeindex, time):
    plot = plt.figure(figsize=(15, 10))
    height = relativepositions(data, blackholeindex)[time, :, 2] - relativepositions(data, blackholeindex)[0, :, 2]
    ax = plot.add_subplot(111)
    ax.set_ylabel('Polar angle/rad')
    ax.set_xlabel('Height/arbitrary units')
    ax.set_ylim([2, 0])
    plt.ticklabel_format(axis='both', style='plain')
    ax.plot(height, polarangle(data, blackholeindex)[time],'k.', markersize=5)
    plt.savefig("polang.png", dpi=170)
    #plt.show()
    return


# Function which: if D=2, counts the number of stars within radii [R]
# of a given black hole as a function of time.
# If D=3, the function counts the number of stars within heights [R]
# of a given black hole as a function of time.

def countmasses(data, blackholeindex, L, D):
    D=2
    if D == 2:
        radiiwitht = relativeradii2d(data, blackholeindex)
        count = np.zeros(tsteps)
        for h in range(tsteps):
            count[h] = len([1 for i in radiiwitht[h] if i < L])
    if D == 3:
        height = relativepositions(data, blackholeindex)[:, :, 2] - relativepositions(data, blackholeindex)[0, :, 2]
        count = np.zeros(tsteps)
        for h in range(tsteps):
            count[h] = len([1 for i in height[h] if i < L])
    else: print('Input valid number of dimensions')

    return count

# Function which plots the radii projected into the xy- plane of all stars at a timestep t,
# against their initial radius at timestep t=0


def diffwithradius(data, blackholeindex, time):
    radiusdiff = relativeradii2d(data, blackholeindex)[time]-relativeradii2d(data, blackholeindex)[0]
    plot = plt.figure(figsize=(15, 10))
    ax = plot.add_subplot(111)
    ax.set_xlabel('Initial radius/arbitrary units')
    ax.set_xlim([0, 80])
    ax.set_ylabel('Radius difference/arbitrary units')
    #plt.title('Radial distribution against angular distribution of stars')
    plt.ticklabel_format(axis='both', style='plain')
    ax.plot(relativeradii2d(data, blackholeindex)[0], radiusdiff, 'k.', markersize=5)
    i = 0
    while os.path.exists('{}{:d}.png'.format('radrad', i)):
        i += 1
    plt.savefig('{}{:d}.png'.format('radrad', i), dpi=170)

    #plt.show()
    return

# Function which plots the height relative to a given black hole of all stars at a timestep t,
# against their initial height at timestep t=0.


def heightwithradius(data, blackholeindex, time):
    height = relativepositions(data, blackholeindex)[time, :, 2] - relativepositions(data, blackholeindex)[0, :, 2]
    plot = plt.figure(figsize=(15, 10))
    ax = plot.add_subplot(111)
    ax.set_xlabel('Initial radius/arbitrary units')
    ax.set_xlim([0, 80])
    ax.set_ylabel('Height from black hole/arbitrary units')
    #plt.title('Radial distribution against angular distribution of stars')
    plt.ticklabel_format(axis='both', style='plain')
    ax.plot(relativeradii2d(data, blackholeindex)[0], height, 'k.', markersize=5)
    plt.savefig("radheight.png", dpi=170)
    #plt.show()
    return


# Function which: if D=2, plots the number of stars within radii [R]
# of a given black hole as a function of time.
# If D=3, the function plots the number of stars within heights [R]
# of a given black hole as a function of time.


def starswithinL(testmassdata, blackholeindex, R, D=2):
    plot = plt.figure(figsize=(15, 10))
    ax = plot.add_subplot(111)
    ax.set_xlabel('Time/arbitrary units')
    if D ==2:
        ax.set_ylabel('Number of stars within radius R')
    if D ==3:
        ax.set_ylabel('Relative height/arbitrary units')
    #plt.title('Number of stars within radius R of Black Hole with time')
    plt.ticklabel_format(axis='both', style='plain', useOffset=False)
    for counter in R:
        plt.plot(t, countmasses(testmassdata, blackholeindex, counter, dimensions), label= 'R = ' + str(counter))
    plt.legend(prop={'size': 30
                     })
    i = 0
    while os.path.exists('{}{:d}.png'.format('Captured', i)):
        i += 1
    plt.savefig('{}{:d}.png'.format('Captured', i), dpi=170)
    #plt.show()
    return

# Function to plot a projection of both galaxies into the xz- plane at a timestep 'time'.


def plotsnapshotz(time):
    plot = plt.figure(figsize=(15, 15))
    ax = plot.add_subplot(111)
    ax.axis('off')
    ax.set_xlim([-230, 230])
    ax.set_ylim([-230, 230])
    # plt.title('Radial distribution against angular distribution of stars')
    plt.ticklabel_format(axis='both', style='plain')
    ax.plot(coords1[time, :, 0], coords1[time, :, 1], 'k.', markersize=5)
    ax.plot(coords2[time, :, 0], coords2[time, :, 1], 'b.', markersize=5)
    plt.savefig('{}{:d}.png'.format('time = ', time + 1), dpi=170)
    # plt.show()
    return


# Function to calculate the total energy and angular momentum (about the origin)
# of the system as a function of time. Also plots both as a function of time.
def conservedquantities():
    GPE = - masses[0]*masses[1]/abs(relativeradii3d(coords1, 0)[:, 1])
    KE1 = (0.5 * masses[0] * massspeed(vels1, 0)**2)
    KE2 = (0.5 * masses[1] * massspeed(vels2, 1)**2)
    energy = GPE + KE1 + KE2
    angmom1 = masses[0] * np.array([np.linalg.norm(x) for x in np.cross(onemass(coords1, 0),
                                                                       onemass(vels1, 0))])
    angmom2 = masses[1] * np.array([np.linalg.norm(x) for x in np.cross(onemass(coords1, 1),
                                                                       onemass(vels1, 1))])
    angmom = angmom1 + angmom2
    plot = plt.figure(figsize=(12, 8))
    ax = plot.add_subplot(211)
    plt.ticklabel_format(style='sci', useOffset=False)
    plt.subplots_adjust(hspace=0.5)
    ax.set_xlabel('Time/arbitrary units')
    ax.set_ylabel('Energy/arbitrary units')
    #ax.set_ylim([0, energy[0]+20])
    plt.ticklabel_format(axis='y', style='sci')
    plt.title('Energy of system as a function of time')
    plt.plot(t, energy)

    ax = plot.add_subplot(212)
    plt.ticklabel_format(style='sci', useOffset=False)
    ax.set_xlabel('Time/arbitrary units')
    ax.set_ylabel('Angular Momentum/arbitrary units')
    #ax.set_ylim([-0.0000001, 0.0000001])
    plt.ticklabel_format(axis='y', style='sci')
    plt.title('Angular momentum of system about the origin as a function of time')
    plt.plot(t, angmom)
    plt.savefig("ConservedQuantities.png", dpi=170)
    plt.show()
    return energy, angmom

