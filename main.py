import numpy as np
from element import get_vecA, get_vecL, get_vecE, get_vecTheta, get_density, get_ele_coor
from node import node_member
from connectivity_matrix import get_mtxEFT
from global_stiff_matrix import get_mtx_K_glo
from add_BC import mtx_K_glo_add_BC
from add_BC import vecf_add_BC
from show_picture import show_picture
from show_picture import beam_scatter
from distribute_force import self_weight_ele_to_node
from distribute_force import add_q_force_ele_to_node
from analy_shearMoment import analy_shearMoment
from refine_mesh_para import refine_mesh_para
import matplotlib.pyplot as plt

# exaggerating deformation
u_mul = 1000


def main():
    # get element and node data
    vecA = get_vecA()
    vecL = get_vecL()
    vecE = get_vecE()
    vecTheta = get_vecTheta()
    density = get_density()
    N_gdof, vecFix, vecF = node_member()
    mtxEFT = get_mtxEFT()
    N_e = mtxEFT.shape[0]

    # add force
    vecF[17 - 1] += -5000

    q_force_global = np.array([-1000.0, 0, 0])
    vecF += add_q_force_ele_to_node(mtxEFT, 3, q_force_global, vecL, vecTheta)
    # add self_weight
    for N_e_i in range(N_e):
        vecF_weight = self_weight_ele_to_node(mtxEFT, N_e_i + 1, vecA, vecL, vecTheta, density)
        vecF += vecF_weight

    # create stiff matrix
    mtx_K_glo = get_mtx_K_glo(mtxEFT, N_gdof, N_e, vecE, vecA, vecL, vecTheta)

    # add boundary condition
    mtx_K_glo_added_BC = mtx_K_glo_add_BC(mtx_K_glo, vecFix)
    vecF_added_bc = vecf_add_BC(vecF, vecFix)

    # calculate vecU
    inv_k = np.linalg.inv(mtx_K_glo_added_BC)
    vecU = inv_k @ vecF_added_bc

    # post process calculate external-force
    F_ext = mtx_K_glo @ vecU

    # get external-force exert on element
    # get element deformation
    ele_f_ext = np.zeros([N_e, 6])
    ele_u = np.zeros([N_e, 6])
    ele_u_no = mtxEFT
    for i in range(N_e):
        for j in range(6):
            ele_u[i][j] = vecU[ele_u_no[i][j] - 1]
            ele_f_ext[i][j] = F_ext[ele_u_no[i][j] - 1]

    # analysis shear-force & moment & stress of special element
    element_no = 4
    xs, ys_stress_xdir_c_pos, ys_stress_xdir_c_neg = analy_shearMoment(ele_u[element_no - 1], vecTheta[element_no - 1],
                                                                       vecE[element_no - 1],
                                                                       vecA[element_no - 1], vecL[element_no - 1],
                                                                       plot_picture=False)
    if (np.max(np.abs(ys_stress_xdir_c_pos)) > np.max(np.abs(ys_stress_xdir_c_neg))):
        ys_stress_max = ys_stress_xdir_c_pos
    else:
        ys_stress_max = ys_stress_xdir_c_neg
    plt.plot(xs, ys_stress_max)
    plt.xlabel('x(m)')
    plt.ylabel('stress x dir(N/m^2)')
    plt.show()
    print(np.max(np.abs(ys_stress_max)))
    # show origin structure
    ele_coor = get_ele_coor()
    show_picture(ele_coor, 'k')

    # show deformation structure
    ele_u_mul = ele_u * u_mul
    for e_i in range(N_e):
        ele_i_ori = ele_coor[e_i] + np.array(
            [ele_u_mul[e_i][0], ele_u_mul[e_i][1], ele_u_mul[e_i][3], ele_u_mul[e_i][4]])
        ele_i_L = vecL[e_i]
        ele_i_the = vecTheta[e_i]
        xs, ys = beam_scatter(ele_i_ori, ele_i_L, ele_i_the,
                              ele_u_mul[e_i])
        plt.scatter(xs, ys, c='r', s=10)
        plt.scatter(xs[0], ys[0], c='b', s=20)
        plt.scatter(xs[-1], ys[-1], c='b', s=20)
    plt.show()

    pass


if __name__ == '__main__':
    main()
