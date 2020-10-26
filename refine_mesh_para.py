import numpy as np


def refine_mesh_para(mtxEFT, vecE, vecA, vecL, vecTheta, e_i):
    # change mtxEFT
    mtxEFT = np.copy(mtxEFT)
    node_max = np.max(mtxEFT)
    element_node = mtxEFT[e_i - 1]
    node1 = element_node[:3]
    node2 = element_node[-3:]
    add_node = np.array([node_max + 1, node_max + 2, node_max + 3])
    element1 = np.hstack([node1, add_node])
    element2 = np.hstack([add_node, node2])
    # add element in EFT
    mtxEFT[e_i - 1] = element1
    mtxEFT = np.vstack([mtxEFT, element2])

    # change element_para
    vecE = np.copy(vecE)
    vecA = np.copy(vecA)
    vecL = np.copy(vecL)
    vecTheta = np.copy(vecTheta)

    vecL[e_i - 1] = vecL[e_i - 1] / 2.0

    vecE = np.hstack([vecE, vecE[e_i - 1]])
    vecA = np.hstack([vecA, vecA[e_i - 1]])
    vecL = np.hstack([vecL, vecL[e_i - 1]])
    vecTheta = np.hstack([vecTheta, vecTheta[e_i - 1]])

    return mtxEFT, vecE, vecA, vecL, vecTheta
    pass


def refine_mesh_coor(ele_coor, ei):
    ele_coor = np.copy(ele_coor)
    node1_coordinate = ele_coor[ei - 1][:2]
    node2_coordinate = ele_coor[ei - 1][-2:]
    node3_coordinate = 0.5 * (node1_coordinate + node2_coordinate)

    ele1 = np.hstack([node1_coordinate, node3_coordinate])
    ele2 = np.hstack([node3_coordinate, node2_coordinate])

    ele_coor[ei - 1] = ele1
    ele_coor = np.vstack([ele_coor, ele2])

    return ele_coor


# record element split number
def record_element_NO(record_elements, refine_element_NO,N_e):
    if refine_element_NO in record_elements:
        index=np.where(record_elements==refine_element_NO)[0][0]

        record_elements.insert(index+1,N_e)

    return record_elements


if __name__ == '__main__':
    import connectivity_matrix

    vecL = np.array([3, 2, 2, 3, 2, 2.385])
    vecA = np.array(
        [62500e-6, 62500e-6, 62500e-6, 62500e-6, 62500e-6, 62500e-6])
    vecE = np.array([210e9, 210e9, 210e9, 210e9, 210e9, 210e9])
    vecTheta = np.array([0, 90, 90, 0, 90, 123.0239])

    mtxEFT = connectivity_matrix.get_mtxEFT()
    mtxEFT, vecE, vecA, vecL, vecTheta = refine_mesh_para(mtxEFT, vecE, vecA, vecL, vecTheta, 1)
    from global_stiff_matrix import get_mtx_K_glo

    mtx_K_glo = get_mtx_K_glo(mtxEFT, np.max(mtxEFT), mtxEFT.shape[0], vecE, vecA, vecL, vecTheta)
    pass
