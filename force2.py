import numpy as np
from math import *

#constants
eps = 0.997
sig = 3.30 * 10e-10

class  atom:
    def __init__(self, position, momentum, neighbours):
        self.pos = position
        self.mom = momentum
        self.nei = neighbours

def map(func, list):
    return [func(x) for x in list]

def dis(x1, x2):
    dl = [x1[i] - x2[i] for i in range(3)]
    dist_sq = sum(map((lambda x: x**2), dl))
    return sqrt(dist_sq)

def zipWith(func, list1, list2):
    return [func(x, y) for (x, y) in zip(list1, list2)]

def force(index, atoms):
    #constants
    delta = 2 * 10**(-10)
    dx = [delta, 0, 0]
    dy = [0, delta, 0]
    dz = [0, 0, delta]
    neilist = atoms[index].nei
    pos_ = atoms[index].pos
    posX = zipWith((lambda x, y : x + y), pos_, dx)
    posY = zipWith((lambda x, y : x + y), pos_, dy)
    posZ = zipWith((lambda x, y : x + y), pos_, dz)

    #functions
    #def dis(x1, x2):
    #    d = x1 - x2
    #    return np.linalg.norm(d)

    def ljpot(pos1, pos2):
        dist = dis(pos1, pos2)
        return (4 * eps) * (((sig**12)/(dist**12))-((sig**6)/(dist**6)))

    def ljpot_(index):
        pot_ = 0
        for i in neilist:
            pot_ += ljpot(pos_, atoms[i].pos)
        return pot_
    def ljpotX(index):
        potX = 0
        for i in neilist:
            potX += ljpot(posX, atoms[i].pos)
        return potX
    def ljpotY(index):
        potY = 0
        for i in neilist:
            potY += ljpot(posY, atoms[i].pos)
        return potY
    def ljpotZ(index):
        potZ = 0
        for i in neilist:
            potZ += ljpot(posZ, atoms[i].pos)
        return potZ

    fx = (-1) * (ljpotX(index) - ljpot_(index)) / delta
    fy = (-1) * (ljpotY(index) - ljpot_(index)) / delta
    fz = (-1) * (ljpotZ(index) - ljpot_(index)) / delta

    return [fx,fy,fz]
