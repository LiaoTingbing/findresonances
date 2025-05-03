
# 作者：廖廷兵
# 时间：2024-01-25
# findresonances ref： V. A. Mandelshtam and H. S. Taylor, "Harmonic inversion of time signals and its application," J. Chem. Phys., vol. 107, no. 17, p. 6756-6769 (1997).


import   numpy as np

def fdmMatrixG_Chen_Guo(complexSignal, M, z):
    # // c是复数信号
    # M是整数
    # 单位向量序列

    Jnum = len(z)

    U0 = np.zeros((Jnum, Jnum), dtype=complex)
    U1 = np.zeros((Jnum, Jnum), dtype=complex)
    U2 = np.zeros((Jnum, Jnum), dtype=complex)
    # D0 = np.zeros((Jnum, Jnum), dtype=complex)
    # D1 = np.zeros((Jnum, Jnum), dtype=complex)
    # D2 = np.zeros((Jnum, Jnum), dtype=complex)

    nv1 = np.linspace(0, M, M + 1, dtype=int)
    nvd = np.linspace(0, M * 2, M * 2 + 1, dtype=int)

    # 0 , M
    c10 = complexSignal[nv1 + 0 + 0]            #0      M
    c30 = complexSignal[nv1 + 0 + 0 + M + 1]    #M+1    2M+1
    c11 = complexSignal[nv1 + 0 + 1]            #1  M+1
    c31 = complexSignal[nv1 + 0 + 1 + M + 1]    #M+2    2M+2
    c12 = complexSignal[nv1 + 0 + 2]            #2  M+2
    c32 = complexSignal[nv1 + 0 + 2 + M + 1]    #M+3    2M+3

    G10 = np.zeros(Jnum, dtype=complex)
    G30 = np.zeros(Jnum, dtype=complex)
    G11 = np.zeros(Jnum, dtype=complex)
    G31 = np.zeros(Jnum, dtype=complex)
    G12 = np.zeros(Jnum, dtype=complex)
    G32 = np.zeros(Jnum, dtype=complex)

    for i in range(0, Jnum):
        ulk = np.power(z[i], -nv1)
        # print(ulk)

        G10[i] = np.dot(ulk, c10)
        G30[i] = np.dot(ulk, c30)

        G11[i] = np.dot(ulk, c11)
        G31[i] = np.dot(ulk, c31)

        G12[i] = np.dot(ulk, c12)
        G32[i] = np.dot(ulk, c32)

    for i in range(0, Jnum):

        u1 = z[i]
        u1K = np.power(u1, -M)

        G10i = G10[i]
        G11i = G11[i]
        G12i = G12[i]
        G30i = G30[i]
        G31i = G31[i]
        G32i = G32[i]

        for j in range(i + 1, Jnum):
            u2 = z[j]
            u2K = np.power(u2, -M)
            coff = 1 / (u1 - u2)
            U0[i, j] = coff * (u1 * G10[j] - u2 * G10i - u1K * G30[j] + u2K * G30i)
            U1[i, j] = coff * (u1 * G11[j] - u2 * G11i - u1K * G31[j] + u2K * G31i)
            U2[i, j] = coff * (u1 * G12[j] - u2 * G12i - u1K * G32[j] + u2K * G32i)

            U0[j,i] = coff * (u1 * G10[j] - u2 * G10i - u1K * G30[j] + u2K * G30i)
            U1[j,i] = coff * (u1 * G11[j] - u2 * G11i - u1K * G31[j] + u2K * G31i)
            U2[j,i] = coff * (u1 * G12[j] - u2 * G12i - u1K * G32[j] + u2K * G32i)

    Kd = M - abs(M - nvd) + 1

    cd0 = complexSignal[nvd + 0 + 0] * Kd
    cd1 = complexSignal[nvd + 0 + 1] * Kd
    cd2 = complexSignal[nvd + 0 + 2] * Kd

    for i in range(0 , Jnum):
        ud = np.power (z[i] , -nvd  )
        # D0[i] = np.dot(ud   , cd0)
        # D1[i] = np.dot(ud   , cd1)
        # D2[i] = np.dot(ud   , cd2)
        U0[i ,i ] = np.dot(ud   , cd0)
        U1[i ,i ] = np.dot(ud   , cd1)
        U2[i , i] = np.dot(ud   , cd2)

    # U0 = U0  + np.diag(D0)
    # U1 =U1  + np.diag(D1)
    # U2 = U2  + np.diag(D2)

    return  U0 , U1 , U2


