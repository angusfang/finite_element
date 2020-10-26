import numpy as np
from math import atan2,pi

'''store element parameters
'''


# element length
def get_vecL():
    vecL = np.array([3, 2, 2, 3, 2, 2.385])
    return vecL


# element area
def get_vecA():
    vecA = np.array(
        [62500e-6, 62500e-6, 62500e-6, 62500e-6, 62500e-6, 62500e-6])
    return vecA


# element young's module
def get_vecE():
    vecE = np.array([210e9, 210e9, 210e9, 210e9, 210e9, 210e9])
    return vecE


# element theta of global coordinate to local coordinate (degree)
def get_vecTheta():
    vecTheta = np.array([0, 90, 90, 0, 90, 123.0239])
    return vecTheta


# element density
def get_density():
    return 7900


# element global coordinate
def get_ele_coor():
    ele_coor = np.array([
        [0, 4, 3, 4],
        [0, 2, 0, 4],
        [3, 2, 3, 4],
        [0, 2, 3, 2],
        [0, 0, 0, 2],
        [4.3, 0, 3, 2]
    ])
    return ele_coor

# element global coordinate for Q5
def get_ele_coor_Q5(x):
    ele_coor = np.array([
        [0, 4, 3, 4],
        [0, 2, 0, 4],
        [3, 2, 3, 4],
        [0, 2, 3, 2],
        [0, 0, 0, 2],
        [x, 0, 3, 2]
    ])
    return ele_coor

# element length
def get_vecL_from_coor(ele_coor):
    count = 0
    for ele in ele_coor:
        x1 = ele[0]
        y1 = ele[1]
        x2 = ele[2]
        y2 = ele[3]
        L = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
        if count == 0:
            vecL = np.array([L])
        else:
            vecL = np.hstack([vecL, L])
        count += 1
    return vecL

# element theta of global coordinate to local coordinate (degree)
def get_vecTheta_from_coor(ele_coor):
    count = 0
    for ele in ele_coor:
        x1 = ele[0]
        y1 = ele[1]
        x2 = ele[2]
        y2 = ele[3]
        Theta = atan2((y2 - y1), (x2 - x1))
        if count == 0:
            vecTheta = np.array([Theta])
        else:
            vecTheta = np.hstack([vecTheta, Theta])
        count += 1
    return vecTheta/pi*180.0
