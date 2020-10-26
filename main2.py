import numpy as np
from element import get_vecA, get_vecL, get_vecE, get_vecTheta, get_density, get_ele_coor
from element import get_vecL_from_coor, get_vecTheta_from_coor
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
from refine_mesh_para import refine_mesh_coor
from refine_mesh_para import record_element_NO
import matplotlib.pyplot as plt

# exaggerating deformation
u_mul = 1000


def main2():
    # get element and node data
    ele_coor = get_ele_coor()
    vecL = get_vecL_from_coor(ele_coor)
    vecTheta = get_vecTheta_from_coor(ele_coor)
    vecA = get_vecA()
    vecE = get_vecE()
    density = get_density()
    N_gdof, vecFix, vecF = node_member()
    mtxEFT = get_mtxEFT()

    # record element to analysis
    record_element1 = [1]
    record_element2 = [2]
    record_element3 = [3]
    record_element4 = [4]
    record_element5 = [5]
    record_element6 = [6]

    # refine element

    element_number = 6
    refine_element_list = np.arange(1, element_number + 1)
    refine_times = 1
    for i in range(1,refine_times+1):
            element_number = element_number*2
            refine_element_list = np.hstack([refine_element_list, np.arange(1, element_number+1)])

    for refine_element_NO in refine_element_list:
        mtxEFT, vecE, vecA, vecL, vecTheta = refine_mesh_para(mtxEFT, vecE, vecA, vecL, vecTheta, refine_element_NO)
        ele_coor = refine_mesh_coor(ele_coor, refine_element_NO)
        # adjust F distribute
        # Zeroing force and add force
        vecF = np.zeros([np.max(mtxEFT)])
        vecF[17 - 1] += -5000
        q_force_global = np.array([-1000.0, 0, 0])
        # record_q_force should exert on which elements
        record_q_force_element3 = np.array([3])
        if refine_element_NO in record_q_force_element3:
            record_q_force_element3 = np.hstack([record_q_force_element3, refine_element_NO])
        for q_force_element_i in record_q_force_element3:
            vecF += add_q_force_ele_to_node(mtxEFT, q_force_element_i, q_force_global, vecL, vecTheta)

        # add self_weight
        N_e = mtxEFT.shape[0]
        for N_e_i in range(N_e):
            vecF_weight = self_weight_ele_to_node(mtxEFT, N_e_i + 1, vecA, vecL, vecTheta, density)
            vecF += vecF_weight

        # create stiff matrix
        mtx_K_glo = get_mtx_K_glo(mtxEFT, np.max(mtxEFT), N_e, vecE, vecA, vecL, vecTheta)

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

        # record element to analysis
        record_element1 = record_element_NO(record_element1, refine_element_NO, N_e)
        record_element2 = record_element_NO(record_element2, refine_element_NO, N_e)
        record_element3 = record_element_NO(record_element3, refine_element_NO, N_e)
        record_element4 = record_element_NO(record_element4, refine_element_NO, N_e)
        record_element5 = record_element_NO(record_element5, refine_element_NO, N_e)
        record_element6 = record_element_NO(record_element6, refine_element_NO, N_e)

    # analysis shear-force & moment & stress of special element
    analysis_element = record_element6
    count = 0
    for elemnet_i in analysis_element:
        element_no = elemnet_i
        if count == 0:
            xs, stress_xdir_cpos, stress_xdir_cneg = analy_shearMoment(ele_u[element_no - 1], vecTheta[element_no - 1],
                                                                       vecE[element_no - 1],
                                                                       vecA[element_no - 1], vecL[element_no - 1],
                                                                       plot_picture=False)
        else:
            xs_ = analy_shearMoment(ele_u[element_no - 1], vecTheta[element_no - 1], vecE[element_no - 1],
                                    vecA[element_no - 1], vecL[element_no - 1], plot_picture=False)[0]
            xs_ += xs[-1]
            stress_xdir_cpos_, stress_xdir_cneg_ = analy_shearMoment(ele_u[element_no - 1], vecTheta[element_no - 1],
                                                                     vecE[element_no - 1],
                                                                     vecA[element_no - 1], vecL[element_no - 1],
                                                                     plot_picture=False)[-2:]
            stress_xdir_cpos = np.hstack([stress_xdir_cpos, stress_xdir_cpos_])
            stress_xdir_cneg = np.hstack([stress_xdir_cneg, stress_xdir_cneg_])
            xs = np.hstack([xs, xs_])
            # find max axial load
        if np.max(np.fabs(stress_xdir_cpos)) > np.max(np.fabs(stress_xdir_cneg)):
            stress_xdir = stress_xdir_cpos
        else:
            stress_xdir = stress_xdir_cneg

        count += 1
    plt.plot(xs, stress_xdir, c='r')
    plt.show()

    # show origin structure
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


if __name__ == '__main__':
    main2()
