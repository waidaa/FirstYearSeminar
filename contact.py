import numpy as np
from force import *
#force2: map, zipWith, dis, force
from random import uniform

#constants
avogadro = 6.022 * 10e23
at_w = 39.95
m = at_w / avogadro
dt = 0.01
radius = 188 * 10e-12 #ファンデルワールス半径を採用
length = height = depth = 0.1
cut_dis = 2.5 * sigma #3.5では？
n = 1000

def judge(x1, x2): #十分遠ければ2,無視できないくらい近ければ1,接触は0
    dist = dis(x1, x2)
    if dist <= 2 * radius:
        return 0
    elif dist <= cut_dis:
        return 1
    else:
        return 2

def set_nei(index, atoms):
    pos = atoms[index].pos
    nei_list = []
    for i in range (n):
        if ((i != index) and (judge(pos, atoms[i].pos) <= 1)):
            nei_list += [i]
    atoms[index].nei = nei_list

def contact(index1, index2): #衝突後の挙動
    x1 = atoms[index1].pos
    x2 = atoms[index2].pos
    p1 = atoms[index1].mom
    p2 = atoms[index2].mom
    if judge(x1, x2) != 0:
        pass
    else:
        drvec_ = zipWith((lambda x, y : x - y) , x1, x2)
        dr = dis(drvec_, [0, 0, 0])
        drvec = map((lambda x : x / dr), drvec_)
        conmom1 = np.dot(p1, drvec)
        conmom2 = np.dot(p2, drvec)
        impuls1 = map((lambda x : x * (-2) * conmum1), drvec)
        impuls2 = map((lambda x : x * (-2) * conmum2), drvec)
        atoms[index1].mom = zipWith((lambda x, y : x + y), p1, impuls1)
        atoms[index2].mom = zipWith((lambda x, y : x + y), p2, impuls2)
