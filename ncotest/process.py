import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime
import pysolar.solar as ps

begin = datetime.datetime(2011,1,1)

def loadFile():
	f = open("/mnt/c/dev/KNMITEST/ncotest/out.txt", 'r')
	fstring = f.read()

	columns = [7, 5, 6, 4]
	timecolumn = 7
	interpolate = [5, 6, 4]


	splitted = fstring.split("\n\n")
	for i in columns:
		splitted[i] = splitted[i].split("\n")
		for j in range(len(splitted[i])):
			if splitted[i][j] == '_':
				splitted[i][j] = math.nan
			else:
				splitted[i][j] = float(splitted[i][j])


	for c in interpolate:
		
		max = len(splitted[c]) - 1
		
		start = 1
		while start < max - 1:
			end = 0
			while (start < max - 1) and (not (math.isnan(splitted[c][start]))):
				start += 1
			end = start + 1
			while (end < max) and math.isnan(splitted[c][end]):
				end += 1
		
			gapstart = splitted[c][start-1]
			gapend = splitted[c][end]
			gaplength = end - start
			delta = (gapend - gapstart) / float(gaplength+1)
			for i in range(gaplength):
				gapstart += delta
				splitted[c][start + i] = gapstart
			start += 1

	database = []
	for i in range(len(splitted[columns[0]])):
		entry = []
		for c in columns:
			if c == timecolumn:
				entry.append(splitted[c][i])
				#entry.append(begin + datetime.timedelta(hours=(splitted[c][i])))
			else:
				entry.append(splitted[c][i])
		database.append(entry)

	
	
	#plt.show()
	return database;

db = loadFile()

lat = 51.8
lon = 5.8

panelDir1 = [0.0,2.0,1.0]
panelDir1 = panelDir1 / np.linalg.norm(panelDir1)

panelDir2 = [.5,1.8,1.0]
panelDir2 = panelDir2 / np.linalg.norm(panelDir2)

panelDir3 = [0.0,2.0,1.0]
panelDir3 = panelDir3 / np.linalg.norm(panelDir3)



for i in range(len(db)):
	time = begin + datetime.timedelta(hours=db[i][0])
	alt = ps.get_altitude(lat, lon, time)
	az = ps.get_azimuth(lat,lon,time)
	radin = ps.radiation.get_radiation_direct(time, alt) if alt > 0 else math.nan
	#radin /= math.sin(alt * 3.14159265358979 / 180.0)
	alt = math.radians(alt)
	az = math.radians(az)
	sunDir = [math.sin(az)*math.cos(alt),math.cos(az)*math.cos(alt),math.sin(alt)]
	panelDir3 = sunDir
	
	cloud = 1.0 - (db[i][1]/ 100.0)
	
	output1 = np.dot(sunDir,panelDir1) * radin * cloud
	output2 = np.dot(sunDir,panelDir2) * radin * cloud
	output3 = np.dot(sunDir,panelDir3) * radin * cloud
	db[i][1] = output1
	db[i][2] = math.nan
	db[i][3] = output3


transposed = [np.array(list(x)) for x in zip(*db)]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs = [math.floor(x / 24.0) for x in transposed[0]]
ys = [x%24.0 for x in transposed[0]]
zs = transposed[1]
ax.scatter(xs, ys, zs, c='r', marker='.')
zs = transposed[3]
ax.scatter(xs, ys, zs, c='b', marker='.')

ax.set_xlabel('DAY')
ax.set_ylabel('HOUR')
plt.show()

"""

N = len(db)

sunDirs = np.empty(shape=(N,3))
radiation = []
cloud = []

for i in range(N):
	time = begin + datetime.timedelta(hours=db[i][0])
	alt = ps.get_altitude(lat, lon, time)
	az = ps.get_azimuth(lat,lon,time)
	radin = ps.radiation.get_radiation_direct(time, alt) if alt > 0.0 else 0.0
	alt = math.radians(alt)
	az = math.radians(az)
	
	sunDirs[i] = [math.sin(az)*math.cos(alt),-math.cos(az)*math.cos(alt),math.sin(alt)]
	
	radiation.append(radin)
	cloud.append(1.0 - (db[i][1]/ 100.0))

radiation = np.transpose(np.array([radiation]))
cloud = np.transpose(np.array([cloud]))

transposed = [np.array(list(x)) for x in zip(*db)]
#plt.plot(transposed[0], transposed[1]/4.0, 'rs', transposed[0], transposed[2]/4.0 + 33.33, 'g--', transposed[0], transposed[3]/4.0 + 66.67, 'b--')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs = [math.floor(x / 24.0) for x in transposed[0]]
ys = [x%24.0 for x in transposed[0]]
zs = np.dot(sunDirs, np.transpose([panelDir1]))
print(zs)
ax.scatter(xs, ys, zs, c='r', marker='.')

ax.set_xlabel('DAY')
ax.set_ylabel('HOUR')
plt.show()


def calcValue(panelDir):
	#print(list(np.dot(sunDirs, panelDir) * radiation))
	product = np.dot(sunDirs, panelDir) * radiation * cloud
	return np.sum(np.maximum(np.zeros(shape=(N,1)), product))

xs = []
ys = []
zs = []

bestValueC = 0.0
bestDirC = []
bestAzC = 0

for palt in np.arange(10,20,1):
	for paz in np.arange(175.0, 185.0, .1): 
		panelDir = np.empty(shape=(3,1))
		panelDir[0,0] = math.sin(math.radians(paz))*math.cos(math.radians(palt))
		panelDir[1,0] = math.cos(math.radians(paz))*math.cos(math.radians(palt))
		panelDir[2,0] = math.sin(math.radians(palt))
		value = calcValue(panelDir)
		panelDir *= value
		xs.append(panelDir[0,0])
		ys.append(panelDir[1,0])
		zs.append(panelDir[2,0])
		if value > bestValueC:
			bestValueC = value
			bestDirC = [panelDir[0,0],panelDir[1,0],panelDir[2,0]]
			bestAzC = paz

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs, c='r', marker='.')
ax.scatter([bestDirC[0]],[bestDirC[1]],[bestDirC[2]], c='g', marker='o')

print(bestAzC)

print(bestValueC / calcValue(np.array([[-.00168586],[-.96592436],[.25881905]])))
cloud = np.average(cloud) * np.ones(shape=(N,1)) * 0.961921488419

xs = []
ys = []
zs = []

bestValue = 0.0
bestDir = []
bestAz = 0

for palt in np.arange(10,20,1):
	for paz in np.arange(175.0, 185.0, .1): 
		panelDir = np.empty(shape=(3,1))
		panelDir[0,0] = math.sin(math.radians(paz))*math.cos(math.radians(palt))
		panelDir[1,0] = math.cos(math.radians(paz))*math.cos(math.radians(palt))
		panelDir[2,0] = math.sin(math.radians(palt))
		value = calcValue(panelDir)
		panelDir *= value
		xs.append(panelDir[0,0])
		ys.append(panelDir[1,0])
		zs.append(panelDir[2,0])
		if value > bestValue:
			bestValue = value
			bestDir = [panelDir[0,0],panelDir[1,0],panelDir[2,0]]
			bestAz = paz


ax.scatter(xs, ys, zs, c='b', marker='.')

ax.scatter([0],[0],[0], c='b', marker='x')
ax.scatter([bestDir[0]],[bestDir[1]],[bestDir[2]], c='g', marker='o')

print(bestAz)


ax.set_xlabel('DAY')
ax.set_ylabel('HOUR')
plt.show()
"""
