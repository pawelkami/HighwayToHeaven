#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from scipy.spatial import ConvexHull
import numpy
import matplotlib.pyplot as plt

cities = []

def getRandomPointFromCircle(point, R, resolution,  hull):
    result = []
    for i in reversed(range(1, (int)(36 * R/resolution))):
        theta = 2.0 * ((float)(i) / (int)(36 * R/resolution)) * numpy.pi
        print theta
        x = R * numpy.cos(theta) + point[0]
        y = R * numpy.sin(theta) + point[1]
        if inHull([x, y], hull):
            result.append([x, y])
    return result


def getAllNeighbour(point, R, hull, resolution):
    neighbour = []
    for x in numpy.arange(resolution, R, resolution):
        neighbour.extend(getRandomPointFromCircle(point, x, resolution, hull))
    return neighbour

'''
Losowy punkt sąsiadujący z punktem point, leżący w promieniu R od niego i znajdujący się w podanej otoczce wypukłej
'''
def getRandomPointNeighbour(point, R, hull):
    while True:
        p = getRandomPointFromCircle(point, R, hull)
        # if inHull(p, hull):
        return p


def getRandomPointFromHull(hull):
    while True:
        point = numpy.random.rand(1,2)
        if(inHull(point, hull)):
            return point

def inHull(p ,hull):
    from scipy.spatial import Delaunay
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0

'''
Funkcja szukająca optymalnego rozwiązania - symulowane wyżarzanie
'''
def simulatedAnnealing():
    initPoint = getRandomPointFromHull(cities)

    return

'''
Funkcja celu:
F=c*S + z * Suma Li

c- współczynnik wagi całkowitej długości autostrady
S - długość zbudowanej dotychczas autostrady
z - współczynnik wagi odległości zjazdu od miasta
Li - odległość i-tego miasta od najbliższego zjazdu
N - liczba miast
'''
def fitnessFunction():
    return

def main():
    # wyświetlanie otoczki wypukłej
    cities = numpy.random.rand(30, 2)
    # hull = ConvexHull(cities)
    # plt.plot(cities[:,0], cities[:,1], 'o')
    #
    # for simplex in hull.simplices:
    #     plt.plot(cities[simplex, 0], cities[simplex,1], "k-")
    #
    # plt.plot(cities[hull.vertices, 0], cities[hull.vertices, 1], 'r--', lw=2)
    # plt.plot(cities[hull.vertices[0], 0], cities[hull.vertices[0], 1], 'ro')
    # plt.show()

    point = numpy.random.rand(1, 2)
    print point
    #print getRandomPointNeighbour(point[0], 0.1, cities)

    neighbor = getAllNeighbour(point[0], 0.2, cities, 0.05)

    x_list = [x for [x, y] in neighbor]
    y_list = [y for [x, y] in neighbor]

    plt.plot(x_list, y_list, 'o')
    plt.show()


    print neighbor

if __name__ == "__main__":
    main()

