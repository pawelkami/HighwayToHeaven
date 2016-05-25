#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from scipy.spatial import ConvexHull
import numpy
import matplotlib.pyplot as plt

cities = []

'''
Funkcja szukająca optymalnego rozwiązania - symulowane wyżarzanie
'''
def simulatedAnnealing():
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
    hull = ConvexHull(cities)
    plt.plot(cities[:,0], cities[:,1], 'o')

    for simplex in hull.simplices:
        plt.plot(cities[simplex, 0], cities[simplex,1], "k-")

    plt.plot(cities[hull.vertices, 0], cities[hull.vertices, 1], 'r--', lw=2)
    plt.plot(cities[hull.vertices[0], 0], cities[hull.vertices[0], 1], 'ro')
    plt.show()

if __name__ == "__main__":
    main()

