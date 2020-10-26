import numpy as np


# connect global_vecU to element, example[x_dir1 y_dir1 theta1 x_dir2 y_dir2 theta2]
def get_mtxEFT():
    mtxEFT = np.array([[13, 14, 15, 16, 17, 18],
                       [10, 11, 12, 13, 14, 15],
                       [7, 8, 9, 16, 17, 18],
                       [10, 11, 12, 7, 8, 9],
                       [1, 2, 3, 10, 11, 12],
                       [4, 5, 6, 7, 8, 9]
                       ])
    return mtxEFT

# connect node to vecU
def get_nodeANDdof_table():
    nodeANDdof = np.array([
        [1,2,3],
        [4,5,6],
        [7,8,9],
        [10,11,12],
        [13,14,15],
        [16,17,18]
    ])
    return nodeANDdof
