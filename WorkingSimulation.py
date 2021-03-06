import numpy as np
import math
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FFMpegWriter
import ffmpy
from mpl_toolkits.mplot3d import Axes3D


def repeat(arr, count):
    return np.repeat(arr[np.newaxis, ...], count, axis=0)

# in these functions, num is the total number of bodies. Same as i in previous versions of the program.


def reshapelist(longlist, num):
    reshaped = np.reshape(longlist, (j, num, k))
    return reshaped


def massarray(sequence, num):
    m = np.repeat(sequence, 3*num)
    n = np.reshape(m, (len(sequence), num, 3))
    return n

# Function which takes in initial conditions in a list of the form [x0, y0, z0, x1, ...zi, vx0, vy0, vz0, vx1, ..., vzi]
# and returns accelerations in a list of the form [ax0, ay0, az0, ..., azi]


def accelerationcalculator(listofinitialconds, numbodies):

    thingy = np.append(repeat(listofinitialconds[0:3], numbodies), repeat(listofinitialconds[3:6], numbodies))
    rj = np.reshape(thingy, (j, numbodies, 3))

    allposns = listofinitialconds[0:3*numbodies]
    thang = np.append(allposns, allposns)
    ri = np.reshape(thang, (j, numbodies, 3))

    difference = rj-ri
    numerator = difference * massarray(masses, numbodies)
    sums = np.sum(difference**2, 2)
    denominator = sums**1.5
    denominatorlimited = np.where(denominator < 0.1, 0.1, denominator)
    denominatorarray = (reshapelist(np.repeat(denominatorlimited, 3), numbodies))

    a_ik = np.sum((numerator/denominatorarray), 0)
    return np.ravel(a_ik)


def derivatives(absolutevariables, t, numbodies):
    positions = absolutevariables[:k*numbodies]
    velocities = absolutevariables[k*numbodies:]
    accelerations = accelerationcalculator(positions, numbodies)
    dynamics = np.append(velocities, accelerations)
    return dynamics


def solveeqn(initconds, time):
    numbodies = int(len(initconds)/6)
    dynamicalvars = scipy.integrate.odeint(derivatives, initconds, time, mxstep=1000000, args=(numbodies,))
    return np.array(dynamicalvars)

tmax = 1500
t = np.array(np.linspace(0, 100, tmax))
numstars = 1000
masses = [60, 150]
direction = [1, 0]
i = numstars + 2
j = 2
k = 3
h = 1
print('number of stars per galaxy is ', numstars)


centerinitialposn1 = np.array([-30, -12, 0])
centerinitialposn2 = np.array([30, 12, 0])
centerinitialvel1 = np.array([7, 0, 0])
centerinitialvel2 = np.array([-7, 0, 0])

# Calculate the initial positions necessary for test masses to execute circular orbits

centerposn1 = np.reshape(repeat(centerinitialposn1, numstars), (numstars, k))
centervel1 = np.reshape(repeat(centerinitialposn1, numstars), (numstars, k))
centerposn2 = np.reshape(repeat(centerinitialposn2, numstars), (numstars, k))
centervel2 = np.reshape(repeat(centerinitialposn2, numstars), (numstars, k))

radii1 = np.random.exponential(20, numstars)+1
radii1structured = np.reshape(np.repeat(radii1, k), (numstars, k))
table1 = np.reshape(np.random.uniform(0, 2*np.pi, numstars*k), (numstars, k))
trig1 = table1 * 0
trig1[:, 0] = np.cos(table1[:, 0])
trig1[:, 1] = np.sin(table1[:, 0])
positions1 = trig1 * radii1structured
positions1[:, 2] = np.random.normal(0, h)
velocities1 = table1 * 0
velocities1[:, 0] = np.sqrt(masses[0] / radii1structured[:, 0]) * -positions1[:, 1] / radii1structured[:, 0] \
                    * (-1) ** direction[0]
velocities1[:, 1] = np.sqrt(masses[0] / radii1structured[:, 0]) * positions1[:, 0] / radii1structured[:, 0] \
                    * (-1) ** direction[0]
velocities1[:, 2] = 0
abspositions1 = positions1 + centerinitialposn1
absvelocities1 = velocities1 + centerinitialvel1
initialconditions1a = np.concatenate((centerinitialposn1, centerinitialposn2, np.ravel(abspositions1),
                                      centerinitialvel1, centerinitialvel2, np.ravel(absvelocities1)))

radii2 = np.random.exponential(20, numstars)+1
radii2structured = np.reshape(np.repeat(radii2, k), (numstars, k))
table = np.reshape(np.random.uniform(0, 2*np.pi, numstars*k), (numstars, k))
trig2 = table * 0
trig2[:, 0] = np.cos(table[:, 0])
trig2[:, 1] = np.sin(table[:, 0])
positions2 = trig2 * radii2structured
positions2[:, 2] = np.random.normal(0, h)
velocities2 = table * 0
velocities2[:, 0] = np.sqrt(masses[1] / radii2structured[:, 0]) * -positions2[:, 1] / radii2structured[:, 0]\
                    * (-1) ** direction[1]
velocities2[:, 1] = np.sqrt(masses[1] / radii2structured[:, 0]) * positions2[:, 0] / radii2structured[:, 0]\
                    * (-1) ** direction[1]
velocities2[:, 2] = 0
abspositions2 = positions2 + centerinitialposn2
absvelocities2 = velocities2 + centerinitialvel2
initialconditions2a = np.concatenate((centerinitialposn1, centerinitialposn2, np.ravel(abspositions2),
                                      centerinitialvel1, centerinitialvel2, np.ravel(absvelocities2)))
print('finished calculating initial conditions')

galaxy1a = solveeqn(initialconditions1a, t)
print('finished solving galaxy1')
galaxy2a = solveeqn(initialconditions2a, t)
print('finished solving galaxy2')
galaxynormaldist = np.random.normal(0, h, (k*len(t)*i*2))

data1 = np.reshape(np.ravel(galaxy1a), (len(t), i * 2, 3))
data2 = np.reshape(np.ravel(galaxy2a), (len(t), i * 2, 3))

# Create 2D animation

plt.style.use('dark_background')
fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-50, 50), ylim=(-50, 50), aspect='equal')
plt.axis('off')

line1, = ax.plot([], [], '.', markersize=1.2, color='w')
line2, = ax.plot([], [], '.', markersize=1.2, color='w')


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2,


def animate(frame):
    line1xs = []
    line1ys = []
    line1zs = []
    line2xs = []
    line2ys = []
    line2zs = []

    for n in range(i):
        line1xs.append(data1[frame][n][0])
        line1ys.append(data1[frame][n][1])
        line2xs.append(data2[frame][n][0])
        line2ys.append(data2[frame][n][1])

    line1.set_data(line1xs,  line1ys)
    line2.set_data(line2xs,  line2ys)

    return line1, line2,


ani = animation.FuncAnimation(fig, animate, len(t),
                              interval=3, blit=True, init_func=init)

# Create 3D animation

def update_graph(frame):
    line1xs = []
    line1ys = []
    line1zs = []
    line2xs = []
    line2ys = []
    line2zs = []

    for m in range(i):
        line1xs.append(data1[frame][m][0])
        line1ys.append(data1[frame][m][1])
        line1zs.append(data1[frame][m][2])
        line2xs.append(data2[frame][m][0])
        line2ys.append(data2[frame][m][1])
        line2zs.append(data2[frame][m][2])

    galaxy1a.set_data(line1xs, line1ys)
    galaxy1a.set_3d_properties(line1zs)
    galaxy2a.set_data(line2xs, line2ys)
    galaxy2a.set_3d_properties(line2zs)
    center1.set_data(data1[frame][0][0], data1[frame][0][1])
    center1.set_3d_properties(data1[frame][0][2])
    center2.set_data(data2[frame][1][0], data2[frame][1][1])
    center2.set_3d_properties(data2[frame][1][2])

    return galaxy1a, galaxy2a, center1, center2


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', xlim=(-120, 120), ylim=(-120, 120), zlim=(-30, 30))
ax.view_init(elev=80, azim=60)

# ax.xaxis.pane.set_edgecolor('black')
# ax.yaxis.pane.set_edgecolor('black')
# ax.xaxis.pane.fill = False
# ax.yaxis.pane.fill = False
# ax.zaxis.pane.fill = False
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.set_zticklabels([])
# ax.xaxis._axinfo["grid"]['linewidth'] = 0.5
# ax.yaxis._axinfo["grid"]['linewidth'] = 0.5
# ax.zaxis._axinfo["grid"]['linewidth'] = 0.5
ax.axis('off')
galaxy1a, = ax.plot([], [], [], linestyle="", marker=",", color='w')
galaxy2a, = ax.plot([], [], [], linestyle="", marker=",", color='y')
center1, = ax.plot([], [], [], linestyle="", marker=",", color='r')
center2, = ax.plot([], [], [], linestyle="", marker=",", color='r')

ani2 = animation.FuncAnimation(fig, update_graph, len(t), interval=3, blit=True)

########################## Analyse final data #################################

# Rearrange data into more manageable format- store position and velocity data separately for each galaxy.

positionswithtime1 = data1[:, :2+numstars, :]
positionswithtime2 = data2[:, :2+numstars, :]
velocitieswithtime1 = data1[:, :2+numstars, :]
velocitieswithtime2 = data2[:, :2+numstars, :]

# Function to pick out position coordinates of a given particle for all time steps.

def massposition(positiondata, massindex):
    posn = np.array([])
    if massindex >= 2 + numstars:
        print('choose mass in appropriate range')
    else: posn = positiondata[:, massindex, :]
    return posn

# Function to pick out velocity of a given particle for all time steps.

def massvelocity(velocitydata, massindex):
    velty = np.array([])
    if massindex >= 2 + numstars:
        print('choose mass in appropriate range')
    else: velty = velocitydata[:, massindex, :]
    return velty

# Function to find the positions of the masses relative to a central mass. blackholeindex may be 0 or 1.

def relativepositions(simulationpositiondata, blackholeindex):
   relposns = simulationpositiondata - np.reshape(np.ravel(np.tile(massposition(simulationpositiondata, blackholeindex),
                                                                   numstars + 2)), (tmax, numstars + 2, 3))
   return relposns


# Function to find 3D distance between all test masses and a given black hole for all time steps.

def relativeradii3d(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    distances = np.zeros((tmax, numstars + 2))
    for time in range(tmax):
        for star in range(numstars + 2):
            distances[time, star] = (coorddifference[time, star, 0]**2 + coorddifference[time, star, 1]**2
                                     + coorddifference[time, star, 2]**2)**0.5
    return distances

# Function to find 2D distance (in xy-plane) between all test masses and a given black hole for all time steps.

def relativeradii2d(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    distances = np.zeros((tmax, numstars + 2))
    for time in range(tmax):
        for star in range(numstars + 2):
            distances[time, star] = (coorddifference[time, star, 0]**2 + coorddifference[time, star, 1]**2)**0.5
    return distances

# Function to find azimuthal angle between all test masses and a given black hole for all time steps.

def azimuthalangle(simulationpositiondata, blackholeindex):
    coorddifference = relativepositions(simulationpositiondata, blackholeindex)
    angles = np.zeros((tmax, numstars + 2))
    for time in range(tmax):
        for star in range(numstars + 2):
            angles[time, star] = np.arctan(coorddifference[time, star, 1] / (coorddifference[time, star, 0]+0.1))
    return angles

# Function to plot 2d radial distribution of given galaxy at timestep t

def radialdistributionplot(testmassdata, blackholeindex, time):
    histogram = plt.figure()
    ax = histogram.add_subplot(111)
    ax.set_xlabel('Radius/Arbitrary units')
    ax.set_ylabel('Frequency')
    plt.title('Radial distribution of stars at time t = ', time)
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
    plt.title('Angular distribution of stars at time t = ', time)
    plt.ticklabel_format(axis='both', style='plain')
    plt.hist(azimuthalangle(testmassdata, blackholeindex)[time], 100)
    plt.show()
    return

# Function to plot 2D histogram showing angular distribution against radial distribution

def radialangularplot(testmassdata, blackholeindex, time):
    plot = plt.figure()
    ax = plot.add_subplot(111)
    ax.set_xlabel('Azimuthal angle/rad')
    ax.set_ylabel('Radius/arbitrary units')
    plt.title('Radial distribution against angular distribution of stars at time t = ', time)
    plt.ticklabel_format(axis='both', style='plain')
    plt.hist2d(azimuthalangle(testmassdata, blackholeindex)[time],
               relativeradii2d(testmassdata, blackholeindex)[time], 200)
    plt.show()
    return

# Function to count how many masses are within a given radius R at a time t

def countmasses(testmassdata, blackholeindex, R):
    radiiwitht = relativeradii2d(testmassdata, blackholeindex)
    count = np.zeros(tmax)
    for h in range(tmax):
        for w in range(numstars + 2):
            if radiiwitht[h][w] > R:
                count[h] +=1
    return count

# Function to plot how many stars are within a radius R as a function of time

def starswithinR(testmassdata, blackholeindex, R):
    plot = plt.figure()
    ax = plot.add_subplot(111)
    ax.set_xlabel('Time/arbitrary units')
    ax.set_ylabel('Number of stars within radius R')
    plt.title('Number of stars within radius R of Black Hole with time')
    plt.ticklabel_format(axis='both', style='plain')
    for counter in R:
        plt.plot(t, countmasses(testmassdata, blackholeindex, counter))
    plt.legend()
    plt.show()
    return

radialangularplot(positionswithtime1, 0, 1000)
starswithinR(positionswithtime1, 0, [20, 25, 30, 50])
