from math import *
from vpython import *
import numpy as np
import random

#constatnts
eps = 0.997
sig = 3.30e-10
avogadro = 6.022e23
at_w = 39.95
m = at_w / avogadro
dt = 1/1000
#radius = 188e-12 #ファンデルワールス半径を採用
radius = 0.01
box_length = 0.5
#cut_dis = 2.5 * sig #3.5では？
cut_dis = 3.5 * radius
n = 1000
temp = 300
steps = 100

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

def force(index):
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

def judge(x1, x2): #十分遠ければ2,無視できないくらい近ければ1,接触は0
    dist = dis(x1, x2)
    if dist <= 2 * radius:
        return 0
    elif dist <= cut_dis:
        return 1
    else:
        return 2

def set_nei(index):
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
        impuls1 = map((lambda x : x * (-2) * conmom1), drvec)
        impuls2 = map((lambda x : x * (-2) * conmom2), drvec)
        atoms[index1].mom = zipWith((lambda x, y : x + y), p1, impuls1)
        atoms[index2].mom = zipWith((lambda x, y : x + y), p2, impuls2)

def contact_wall(index):
    fc = np.array([0., 0., 0.])
    if max(atoms[index].pos)>box_length or min(atoms[index].pos)<0:
        for i in range(3):
            if atoms[index].pos[i]>box_length:
                atoms[index].pos[i]=2 * box_length-atoms[index].pos[i]
                atoms[index].mom[i]=(atoms[index].mom[i])*(-1)
                #return np.array(atoms[index].mom) #
                fc[i] += 2 * abs(atoms[index].mom[i])
            elif atoms[index].pos[i]<0:
                atoms[index].pos[i]=(atoms[index].pos[i])*(-1)
                atoms[index].mom[i]=(atoms[index].mom[i])*(-1)
                #return np.array(atoms[index].mom) #
                fc[i] += 2 * abs(atoms[index].mom[i])
    return fc
    #else:
        #return np.array([0,0,0])

def initialize():
    global atoms
    v_square=3*1.38064852e-23*temp #divided_by_m #3or6
    for i in range(n):
        vx_square=random.uniform(0,v_square) #0を-にした
        vr_square=v_square-vx_square
        vy_square=random.uniform(0,vr_square)
        vr_square=v_square-vy_square
        atoms[i].mom[0]=sqrt(vx_square*m)
        atoms[i].mom[1]=sqrt(vy_square*m)
        atoms[i].mom[2]=sqrt(vr_square*m)
        atoms[i].pos[0]=random.uniform(0,box_length)
        atoms[i].pos[1]=random.uniform(0,box_length)
        atoms[i].pos[2]=random.uniform(0,box_length)



#ここからmain
atoms = [atom([0,0,0],[0,0,0],[]) for i in range(n)]


initialize()
#rod = cylinder(pos = vector(0,0,0), axis = vector(box_length, 0, 0), radius = radius / 2)
#rod.color = vector(0,1,0)
rods = [cylinder(pos = vector(0,0,0), axis = vector(box_length, 0, 0), radius = radius / 2), cylinder(pos = vector(0,0,0), axis = vector(0, box_length, 0), radius = radius / 2), cylinder(pos = vector(0,0,0), axis = vector(0, 0, box_length), radius = radius / 2), cylinder(pos = vector(box_length, box_length, 0), axis = vector(-box_length, 0, 0), radius = radius / 2), cylinder(pos = vector(box_length, box_length, 0), axis = vector(0, -box_length, 0), radius = radius / 2), cylinder(pos = vector(box_length, box_length, 0), axis = vector(0, 0, box_length), radius = radius / 2), cylinder(pos = vector(box_length, 0 , box_length), axis = vector(-box_length, 0, 0), radius = radius / 2), cylinder(pos = vector(box_length, 0, box_length), axis = vector(0, box_length, 0), radius = radius / 2), cylinder(pos = vector(box_length, 0, box_length), axis = vector(0, 0, -box_length), radius = radius / 2), cylinder(pos = vector(0, box_length, box_length), axis = vector(box_length, 0, 0), radius = radius / 2), cylinder(pos = vector(0, box_length, box_length), axis = vector(0, -box_length, 0), radius = radius / 2), cylinder(pos = vector(0, box_length, box_length), axis = vector(0, 0, -box_length), radius = radius / 2)]
for i in range(12):
    rods[i].color = vector(0,1,0)
s = [sphere(radius = radius) for i in range(n)]
for i in range(0,n):
    s[i].pos=vector(*(atoms[i].pos))
force_sum = np.array([0.,0.,0.])
for i1 in range(steps):
    move = cylinder(pos = vector(*atoms[0].pos), axis = vector(*(dt * np.array(atoms[0].mom) / m )), radius = radius/2)
    move.color = vector(1,0,0)
    s[0].color = vector(1,0,0)
    for i in range(n):
        #print(i)
        #print(atoms[i].mom)
        #contact_wall(i)
        delta = contact_wall(i)
        #print(delta)
        force_sum += delta
        for j in range(n):
            if j!=i:
                contact(i,j)
        set_nei(i)
        force_to_index=force(i)
        atoms[i].mom = np.ndarray.tolist(np.array(atoms[i].mom) + 0.5 * dt * np.array(force_to_index))
        #atoms[index].mom =np.ndarray.tolist((np.array(atoms[index].mom))*(1+0.5*dt*np.array(force_to_index)))
#    atoms[index].pos =np.ndarray.tolist((np.array(atoms[index].pos))*(1+np.array(atoms[index].mom) / m))
        atoms[i].pos = np.ndarray.tolist(np.array(atoms[i].pos) + dt * np.array(atoms[i].mom) / m)
        s[i].pos=vector(*atoms[i].pos)
        force_to_index=force(i)
#    atoms[index].mom=atoms[index].mom+0.5*dt*force_to_index
        atoms[i].mom = np.ndarray.tolist(np.array(atoms[i].mom) + 0.5 * dt * np.array(force_to_index))
print(force_sum)
