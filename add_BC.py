import numpy as np

def mtx_K_glo_add_BC(mtx_K_glo,vecFix):
    mtx_K_glo = mtx_K_glo.copy()
    for i in range(vecFix.shape[0]):
        mtx_K_glo[:,vecFix[i]-1]=0;
        mtx_K_glo[vecFix[i]-1,:]=0;
        mtx_K_glo[vecFix[i]-1,vecFix[i]-1]=1;
    return mtx_K_glo

def vecf_add_BC(vecf,vecFix):
    vecf = vecf.copy()
    for i in range(vecFix.shape[0]):
        vecf[vecFix[i]-1]=0;
    return vecf
