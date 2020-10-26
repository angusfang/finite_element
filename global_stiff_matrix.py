import numpy as np
import math


def main():
    pass


def get_mtx_K_glo(mtxFET, N_gdof, N_e, vecE, vecA, vecL, vecTheta):
    mtx_K_glo = np.zeros([N_gdof, N_gdof])
    for e_i in range(N_e):
        mtx_K_e = get_mtx_K_e(vecE[e_i], vecA[e_i], vecL[e_i], vecTheta[e_i])
        for i in range(6):
            for j in range(6):
                mtx_K_glo[mtxFET[e_i, i] - 1, mtxFET[e_i, j] - 1] = mtx_K_glo[mtxFET[e_i, i] - 1, mtxFET[e_i, j] - 1] + \
                                                                    mtx_K_e[i, j]

    return mtx_K_glo


def get_mtx_K_e(E, A, L, D):
    R = D * 3.1415926 / 180.0
    c = math.cos(R)
    s = math.sin(R)
    I = (1 / 12.0) * A * A
    T = np.array([
        [c, s, 0, 0, 0, 0],
        [-s, c, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, c, s, 0],
        [0, 0, 0, -s, c, 0],
        [0, 0, 0, 0, 0, 1]
    ])
    T_t = np.array([
        [c, -s, 0, 0, 0, 0],
        [s, c, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, c, -s, 0],
        [0, 0, 0, s, c, 0],
        [0, 0, 0, 0, 0, 1]
    ])
    K_bar = E * A / L * np.array([
        [1, -1],
        [-1, 1]
    ])
    K_beam = E * I / (L ** 3) * np.array([
        [12, 6 * L, -12, 6 * L],
        [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
        [-12, -6 * L, 12, -6 * L],
        [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]
    ])

    K = np.array([
        [K_bar[0][0], 0, 0, K_bar[0][1], 0, 0],
        [0, K_beam[0][0], K_beam[0][1], 0, K_beam[0][2], K_beam[0][3]],
        [0, K_beam[1][0], K_beam[1][1], 0, K_beam[1][2], K_beam[1][3]],
        [K_bar[1][0], 0, 0, K_bar[1][1], 0, 0],
        [0, K_beam[2][0], K_beam[2][1], 0, K_beam[2][2], K_beam[2][3]],
        [0, K_beam[3][0], K_beam[3][1], 0, K_beam[3][2], K_beam[3][3]],
    ])
    mtx_K_e = T_t @ K @ T

    return mtx_K_e


if __name__ == '__main__':
    main()
