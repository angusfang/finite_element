import math
import numpy as np
import matplotlib.pyplot as plt


def moment_calculate(E, I, L, xs, vecU):
    def ddN1(x):
        return 2 * (L + 2 * x) / (L ** 3) - 4 * (2 * L - 2 * x) / L ** 3

    def ddN2(x):
        return 2 * x / L ** 2 - 2 * (2 * L - 2 * x) / L ** 2

    def ddN3(x):
        return 2 * (3 * L - 2 * x) / L ** 3 - 8 * x / L ** 3

    def ddN4(x):
        return 4 * x / L ** 2 - 2 * (L - x) / L ** 2

    moment_s = E*I*(ddN1(xs) * vecU[1] + ddN2(xs) * vecU[2] + ddN3(xs) * vecU[4] + ddN4(xs) * vecU[5])
    return moment_s


def analy_shearMoment(U_vec6_global, theata, E, A, L, plot_picture=True):
    I = (1 / 12.0) * A * A
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

    K_mtx6_local = np.array([
        [K_bar[0][0], 0, 0, K_bar[0][1], 0, 0],
        [0, K_beam[0][0], K_beam[0][1], 0, K_beam[0][2], K_beam[0][3]],
        [0, K_beam[1][0], K_beam[1][1], 0, K_beam[1][2], K_beam[1][3]],
        [K_bar[1][0], 0, 0, K_bar[1][1], 0, 0],
        [0, K_beam[2][0], K_beam[2][1], 0, K_beam[2][2], K_beam[2][3]],
        [0, K_beam[3][0], K_beam[3][1], 0, K_beam[3][2], K_beam[3][3]],
    ])

    rad = theata * 3.1415926 / 180.0
    c = math.cos(rad)
    s = math.sin(rad)
    T3x3_t = np.array([
        [c, s, 0],
        [-s, c, 0],
        [0, 0, 1]
    ])

    T6x6_t = np.kron(np.eye(2, dtype=float), T3x3_t)

    U_vec6_local = T6x6_t @ U_vec6_global
    F_vec6_local = K_mtx6_local @ U_vec6_local

    xs = np.linspace(0, L, 1000)

    ys_moment_zdir = moment_calculate(E,I,L,xs,U_vec6_local)

    ys_xforce_xdir = xs * 0
    ys_xforce_xdir += -F_vec6_local[0]


    c = (A ** 0.5) / 2
    ys_stress_xdir_c_pos = ys_moment_zdir * c / I + ys_xforce_xdir / A
    ys_stress_xdir_c_neg = ys_moment_zdir * -c / I + ys_xforce_xdir / A

    if plot_picture == True:


        plt.plot(xs, ys_moment_zdir)
        plt.xlabel('x(m)')
        plt.ylabel('moment force z_dir(m)')
        plt.show()

        plt.plot(xs, ys_xforce_xdir)
        plt.xlabel('x(m)')
        plt.ylabel('force_x_dir(m)')
        plt.show()

        plt.plot(xs, ys_stress_xdir_c_pos)
        plt.xlabel('x(m)')
        plt.ylabel('stress_x_dir +c (m)')
        plt.show()

        plt.plot(xs, ys_stress_xdir_c_neg)
        plt.xlabel('x(m)')
        plt.ylabel('stress_x_dir -c (m)')
        plt.show()

    return xs, ys_stress_xdir_c_pos,ys_stress_xdir_c_neg


if __name__ == '__main__':
    pass
