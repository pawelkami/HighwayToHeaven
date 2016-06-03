#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from scipy.spatial import ConvexHull
import numpy
import matplotlib.pyplot as plt
import math
import random
import time
'''
Lista współrzędnych miast
'''
cities = []

'''
Lista punktów budujących autostradę
'''
highway = []

costHighway = 3.0
costTurnoffs = 2.8

class Point:
    def __init__(self, parent, coord):
        self.parent = parent
        self.coord = coord

    def getRandomPointFromCircle(self, R, resolution,  hull):
        result = []
        for i in reversed(range(1, (int)(36 * R/resolution))):
            theta = 2.0 * ((float)(i) / (int)(36 * R/resolution)) * numpy.pi
            x = R * numpy.cos(theta) + self.coord[0]
            y = R * numpy.sin(theta) + self.coord[1]
            if inHull([x,y], hull):
                result.append(Point(self, [x, y]))
        return result


    def getAllNeighbour(self, R, hull, resolution):
        neighbour = []
        for x in numpy.arange(resolution, R, resolution):
            neighbour.extend(self.getRandomPointFromCircle(x, resolution, hull))
        return neighbour

    def calcLength(self):
        if self.parent == None:
            return 0.0

        return math.hypot(self.coord[0] - self.parent.coord[0], self.coord[1] - self.parent.coord[1])


def getRandomPointFromHull(hull):
    while True:
        point = numpy.random.rand(1,2)
        if(inHull(point[0], hull)):
            return Point(None, point[0])

def inHull(p ,hull):
    from scipy.spatial import Delaunay
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0

'''
Funkcja szukająca optymalnego rozwiązania - symulowane wyżarzanie
'''
def simulatedAnnealing(R, resolution, temperature, iterations):
    global cities
    global highway
    global costHighway
    global costTurnoffs
    initPoint = getRandomPointFromHull(cities)
    highway.append(initPoint)
    neighbours = initPoint.getAllNeighbour(R, cities, resolution)
    workPoint = selectBest(neighbours, cities, highway, costHighway, costTurnoffs)
    k = 0
    while k < iterations:
        print k
        # newPoint = selectRandom(workPoint.getAllNeighbour(R,cities,resolution))
        newPoint = selectBest(workPoint.getAllNeighbour(R, cities, resolution),cities,highway,costHighway,costTurnoffs)
        fitnessWork = fitnessFunction(workPoint, cities, highway, costHighway, costTurnoffs)
        fitnessNew = fitnessFunction(newPoint, cities, highway, costHighway, costTurnoffs)
        if fitnessWork > fitnessNew:
            workPoint = newPoint
            neighbours.extend(workPoint.getAllNeighbour(R, cities, resolution))
            #neighbours.remove(workPoint)
            highway.append(workPoint)
        elif numpy.random.uniform() < probability(fitnessWork, fitnessNew, temperature):
            #workPoint = newPoint
            workPoint = selectRandom(neighbours)
            neighbours.extend(workPoint.getAllNeighbour(R, cities, resolution))
            neighbours.remove(workPoint)
            highway.append(workPoint)

        k += 1


def selectRandom(neighbours):
    return random.choice(neighbours)

def probability(fitnessWork, fitnessNew, temperature):
    return math.exp( - abs(fitnessWork - fitnessNew) / temperature)

'''
Funkcja celu:
F=c*S + z * Suma Li

costHighway - współczynnik wagi całkowitej długości autostrady
S - długość zbudowanej dotychczas autostrady
costTurnoffs - współczynnik wagi odległości zjazdu od miasta
Li - odległość i-tego miasta od najbliższego zjazdu
N - liczba miast
'''
def fitnessFunction(point, cities, highway, costHighway, costTurnoffs):
    tempHighway = list(highway)
    tempHighway.append(point)
    return costHighway * findHighwayLength(tempHighway) + costTurnoffs * findTurnoffs(cities, tempHighway)

def calcLengthBetweenTwoPoints(point1, point2):
    return math.hypot(point2[0]-point1[0], point2[1]-point1[1])

def findTurnoffs(cities, highway):
    length = 0.0
    for city in cities:
        length += findNearestTurnoffDistance(city, highway)

    return length

def findNearestTurnoffDistance(city, highway):
    minDistance = float("inf")
    for turnoff in highway:
        tempDist = calcLengthBetweenTwoPoints(city, turnoff.coord)
        if tempDist < minDistance:
            minDistance = tempDist

    return minDistance

def findNearestTurnoffPoint(city, highway):
    minDistance = float("inf")
    point = city
    for turnoff in highway:
        tempDist = calcLengthBetweenTwoPoints(city, turnoff.coord)
        if tempDist < minDistance:
            minDistance = tempDist
            point = turnoff
    return point

def findHighwayLength(highway):
    length = 0.0
    for h in highway:
        length += h.calcLength()

    return length


def selectBest(points, cities, highway, costHighway, costTurnoffs):
    bestValue = float("inf")
    bestPoint = points[0]

    for p in points:
        fitness = fitnessFunction(p, cities, highway, costHighway, costTurnoffs)
        if fitness < bestValue:
            bestValue = fitness
            bestPoint = p

    return bestPoint

def drawPlot(R, resolution, temperature, iterations, timeElapsed):
    hull = ConvexHull(cities)
    plt.plot(cities[:, 0], cities[:, 1], 'o')

    for simplex in hull.simplices:
        plt.plot(cities[simplex, 0], cities[simplex, 1], "k-")

    plt.plot(cities[hull.vertices, 0], cities[hull.vertices, 1], 'r--', lw=2)
    plt.plot(cities[hull.vertices[0], 0], cities[hull.vertices[0], 1], 'ro')

    #rysowanie calej autostrady
    for n in highway:
        x_list = []
        y_list = []
        if n.parent != None:
            x_list.append(n.parent.coord[0])
            y_list.append(n.parent.coord[1])
        x_list.append(n.coord[0])
        y_list.append(n.coord[1])

        plt.plot(x_list, y_list)

    # wyswietlanie najblizszych zjazdow
    for city in cities:
        x_list = []
        y_list = []
        point = findNearestTurnoffPoint(city, highway)
        x_list.append(point.coord[0])
        y_list.append(point.coord[1])
        plt.plot(x_list, y_list, 'go')

    plt.title("fitness = " + str(fitnessFunction(highway[-1], cities, highway, costHighway, costTurnoffs)) + str(" time = " + str(timeElapsed)))
    plt.savefig("plot-" + "iter" + str(iterations) + "-temp" + str(temperature) + "-R" + str(R) + "-resolution" + str(resolution) + ".png")
    #plt.show()

def main():
    global cities
    cities = numpy.random.rand(10, 2)

    R = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]
    resolution = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    temperature = [0.01, 0.05, 0.2, 0.5, 1.0, 1.0, 0.7]
    iterations = [1000, 2000, 3000, 4000, 5000, 6000, 7000]

    for i in range(0, len(R)):
        highway[:] = []
        start = time.time()
        simulatedAnnealing(R[i], resolution[i], temperature[i], iterations[i])
        end = time.time()
        drawPlot(R[i], resolution[i], temperature[i], iterations[i], end-start)



if __name__ == "__main__":
    main()
