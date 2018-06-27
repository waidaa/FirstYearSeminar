import numpy as np

class  atom:
    def __init__(self, position, momentum, neighbours):
        self.pos = position
        self.mom = momentum
        self.nei = neighbours

def force(index, atoms):
    #constants
    eps = 0.997
    sig = 3.30 * 10**(-10)
    delta = 2 * 10**(-10)
    dx = np.array([delta, 0, 0])
    dy = np.array([0, delta, 0])
    dz = np.array([0, 0, delta])
    neilist = atoms[index].nei
    pos_ = atoms[index].pos
    posX = pos_ + dx
    posY = pos_ + dy
    posZ = pos_ + dz

    #functions
    def dis(x1, x2):
        d = x1 - x2
        return np.linalg.norm(d)

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

    fx = (ljpotX(index) - ljpot_(index)) / delta
    fy = (ljpotY(index) - ljpot_(index)) / delta
    fz = (ljpotZ(index) - ljpot_(index)) / delta

    return (-1) * np.array([fx,fy,fz])
#kokomade
