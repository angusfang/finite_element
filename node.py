import numpy as np

'''store node parameters'''


def node_member():
    N_gdof = 18
    vecFix = np.array([1, 2, 3, 5])
    vecF = np.zeros([N_gdof])

    return N_gdof, vecFix, vecF


if __name__ == '__main__':
    node_member()
