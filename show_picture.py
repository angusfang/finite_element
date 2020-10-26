import numpy as np
import math
import matplotlib.pyplot as plt


def beam_scatter(ele_ori, ele_L, ele_the, ele_u):
    R = ele_the * 3.1415926 / 180.0
    c = math.cos(R)
    s = math.sin(R)
    T = np.array([
        [c, s, 0, 0, 0, 0],
        [-s, c, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, c, s, 0],
        [0, 0, 0, -s, c, 0],
        [0, 0, 0, 0, 0, 1]
    ])
    ele_u_e = T @ ele_u
    v1 = ele_u_e[1]
    the1 = ele_u_e[2]
    v2 = ele_u_e[4]
    the2 = ele_u_e[5]

    def N1(x):
        return (1 / ele_L ** 3) * (ele_L - x) ** 2 * (2 * x + ele_L)

    def N2(x):
        return (1 / ele_L ** 2) * (ele_L - x) ** 2 * x

    def N3(x):
        return (1 / ele_L ** 3) * (3 * ele_L - 2 * x) * x ** 2

    def N4(x):
        return (1 / ele_L ** 2) * (x - ele_L) * x ** 2

    xs = np.linspace(0, ele_L, 100)
    vs = N1(xs) * v1 + N2(xs) * the1 + N3(xs) * v2 + N4(xs) * the2

    ele_the_rad = ele_the / 180.0 * 3.1415926

    vs = vs - vs[0]

    xs_ = xs * math.cos(ele_the_rad) - vs * math.sin(ele_the_rad)
    vs_ = xs * math.sin(ele_the_rad) + vs * math.cos(ele_the_rad)

    xs_glo = xs_ + ele_ori[0]
    vs_glo = vs_ + ele_ori[1]

    return xs_glo, vs_glo


def make_scatter_point(p1, p2, n=50):
    x = np.linspace(p1[0], p2[0], n)
    y = np.linspace(p1[1], p2[1], n)
    return x, y


def show_picture(matrix, color):
    ele_N = matrix.shape[0]
    for i in range(ele_N):
        p1x = matrix[i][0]
        p1y = matrix[i][1]
        p2x = matrix[i][2]
        p2y = matrix[i][3]
        xs, ys = make_scatter_point([p1x, p1y], [p2x, p2y])
        plt.scatter(xs, ys, c=color, s=2)



