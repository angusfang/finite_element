import math
import numpy as np


def add_q_force(L, qx_global, theta):
    rad = theta * 3.1415926 / 180.0
    c = math.cos(rad)
    s = math.sin(rad)
    T = np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])
    T_t = np.array([
        [c, s, 0],
        [-s, c, 0],
        [0, 0, 1]
    ])
    qx_local = T_t @ qx_global
    vecf_local = L / 2 * np.array([
        [qx_local[0]],
        [qx_local[1]],
        [qx_local[1] * L / 6],
        [qx_local[0]],
        [qx_local[1]],
        [-qx_local[1] * L / 6],
    ])
    T2 = np.kron(np.eye(2, dtype=float), T)
    vecf_global = T2 @ vecf_local
    return vecf_global


def self_weight_ele_to_node(mtxEFT, e_i, vecA, vecL, vectheta, density):
    vecf_ = np.zeros([np.max(mtxEFT)])
    node_list = mtxEFT[e_i - 1]
    qx_global = np.array([0, -density * vecA[e_i - 1] * 9.81, 0])
    self_weight_vecf = add_q_force(vecL[e_i - 1], qx_global, vectheta[e_i - 1])
    count = 0
    for node_i in node_list:
        vecf_[node_i - 1] = self_weight_vecf[count]
        count += 1
    return vecf_

def add_q_force_ele_to_node(mtxEFT, e_i ,qx_global,vecL, vectheta):
    vecf_ = np.zeros([np.max(mtxEFT)])
    node_list = mtxEFT[e_i - 1]
    self_weight_vecf = add_q_force(vecL[e_i - 1], qx_global, vectheta[e_i - 1])
    count = 0
    for node_i in node_list:
        vecf_[node_i - 1] = self_weight_vecf[count]
        count += 1
    return vecf_


if __name__ == '__main__':
    pass
