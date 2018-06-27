import numpy as np
import math
import random
m=0.1
dt=0.1
box_length=box_height=box_depth=0.1
cutoff_distance_square=(2.5*3.3*10e-10)*(2.5*3.3*10e-10)
n=1000
class atom:
    def __init__(self, pos, mom,neighbour_list):
        self.pos = pos
        self.mom = mom
        self.neighbour_list=neighbour_list
atoms[n]=atom(np.array([0,0,0]),np.array([0,0,0]),[])
def force(index):
    #書き換えて
    return (np.array[fx,fy,fz])
def set_neighbour(index):
    atoms[index].neighbour_list.clear()
    for i in range[0,n]:
        if i!=index:
            if (atoms[index].pos - atmos[i].pos).dot(atoms[index].pos - atmos[i].pos) <=  cutoff_distance_square:
                atoms[index].neighbour_list.append(i)
def initialize():
    temp=input("Input the tempereatrue:")
    v_square=3*1.38064852*10e-23*temp #divided_by_m
    for i in range[0,n]:
        vx_square=random.uniform(0,v_square)
        vr_square=v_square-vx_square
        vy_square=random.uniform(0,vr_square)
        vr_square=v_square-vy_square
        atoms[i].mom[0]=math.sqrt(vx_square*m)
        atoms[i].mom[1]=math.sqrt(vy_square*m)
        atoms[i].mom[2]=math.sqrt(vr_square*m)
for index in range[0,n]:
    set_neighbour(index)
    force_to_index=force(index)
    atoms[index].mom=atoms[index].mom+0.5*dt*force_to_index
    pos_tmp=atoms[index].pos+dt*atoms[index].mom/m
    #衝突回数を計算するならここに書く
    if pos_tmp[0]>box_length or pos_tmp[0]<0 or pos_tmp[1]>box_height or pos_tmp[1]<0 or pos_tmp[2]>box_depth or pos_tmp[2]<0:
        for i in range[0,2]:
            if atoms[index].pos[i]>1:
                atoms[index].pos[i]=2-atoms[index].pos[i]
                atoms[index].mom[i]=(atoms[index].mom[i])*(-1)
            elif atoms[index].pos[i]<0:
                atoms[index].pos[i]=(atoms[index].pos[i])*(-1)
                atoms[index].mom[i]=(atoms[index].mom[i])*(-1)
    force_to_index=force(index)
    atoms[index].mom=atoms[index].mom+0.5*dt*force_to_index
   
