
#
# 作者：廖廷兵
# 时间：2024-01-25
# findresonances ref： V. A. Mandelshtam and H. S. Taylor, "Harmonic inversion of time signals and its application," J. Chem. Phys., vol. 107, no. 17, p. 6756-6769 (1997).

import  numpy as np
from fdmMatrixG import fdmMatrixG_Chen_Guo
import pandas as pd

def findresonances(timeSeries, complexSignal, freqWindows,
                   jNum=0, relativeMin=1e-6, MaxIteration=100, tolCutoff = 1e-10):

    # 输入
    # t   :时间序列
    # complexSignal   :时域信号
    # freq_window :频率窗
    # jNum    :初始基函数数量
    # relativeMin  :消除特征值的百分比，小于该值去掉特征分量

    dt = timeSeries[1] - timeSeries[0]
    N = len(complexSignal)
    M = int(np.floor((N - 4) / 2))
    # jNum = 200
    if jNum==0:
        jNum = int( np.ceil(N * dt * 2 * np.pi * (freqWindows[1] - freqWindows[0]) / 4 / np.pi) )
    z = np.exp(-1j * dt * 2 * np.pi * np.linspace(freqWindows[0], freqWindows[1], jNum))
    # print(z.shape)

    # MaxIteration = 50
    for n in range(1, MaxIteration + 1):
        # print(n)
        if n > 1:
            z = zNext
            # print(len(z))
        U0, U1, U2 = fdmMatrixG_Chen_Guo(complexSignal, M, z)
        print(U0)
        # print(U0)
        U, ds, Vh = np.linalg.svd(U0)
        print(ds)
        # print(Vh)
        # relativeMin = 1e-5
        si1 = abs(ds) > max(abs(ds)) * relativeMin
        u1 = U[:, si1]
        s1 = ds[si1]
        s1h = np.diag(np.power(s1, (-0.5)))
        v1 = np.transpose(Vh).conjugate()[:, si1]

        # print(v1)

        Smaxtrix = s1h @ np.transpose(u1).conjugate() @ U1 @ v1 @ s1h

        # print(Smaxtrix)

        u, P = np.linalg.eig(a=Smaxtrix)

        zNext = u / abs(u)
        # print(zNext)

        if len(z) == len(zNext) and np.linalg.norm(z - zNext) < tolCutoff:
            break

    B = v1 @ s1h @ P

    # 计算频率
    frequency = np.real(1j * np.log(u) / dt / 2 / np.pi)
    decay_constant = - np.imag(1j * np.log(u) / dt)
    Q_factor = 2 * np.pi * np.abs(frequency) / 2 / decay_constant
    # print(Q_factor)

    # z是输入，u是输出特征值
    # 误差计算和归一化
    R = np.zeros(len(u))
    for i in range(0, len(u)):

        B[:, i] = B[:, i] / np.sqrt(np.transpose(B[:, i]) @ U0 @ B[:, i])

        u_0 = u[i]
        u_guest = np.sqrt(np.transpose(B[:, i]) @ U2 @ B[:, i])

        if u_0 / u_guest < 0:
            u_guest = -u_guest

        R[i] = np.abs(np.log(u_guest / u_0)) / np.abs(np.log(u_0))

    # 计算幅值
    cz = np.zeros(len(z), dtype=complex)
    nv1 = np.linspace(0, M, M + 1, dtype=int)
    c10 = complexSignal[nv1]

    for i in range(0, len(z)):
        cz[i] = np.dot(c10, np.power(z[i], -nv1))

    complexAmplitude = np.zeros(len(u), dtype=complex)
    for i in range(0, len(u)):
        complexAmplitude[i] = np.dot(B[:, i], cz) ** 2

    amplitude = np.abs(complexAmplitude)
    # print(amplitude)

    error_estimate = R
    phase = np.angle(complexAmplitude)

    # print(phase)

    fs_detaild = {'frequency': frequency, 'decay_constant': decay_constant, 'Q_factor': Q_factor,
                  'amplitude': amplitude, 'phase': phase, 'error_estimate': error_estimate}
    fs_detaild_pd = pd.DataFrame(fs_detaild).sort_values(by='frequency', ascending=True)

    print(fs_detaild_pd.to_string())

    return fs_detaild

if __name__ == '__main__':
    f1 = 1.765;
    alpha1 = 0.005;
    ampl1 = 1.3;
    phase1 = 0.4;
    f2 = 2.345;
    alpha2 = 0.012;
    ampl2 = 0.45;
    phase2 = 1.234;
    # %  # Time intervals of the signal:
    t_long = np.linspace(0, 20, 201);
    # % t_short = t_long(1:20);
    # %  # Signal:
    signal_long = ampl1 * np.exp(-1j * (2 * np.pi * f1 * t_long - phase1))*np.exp(-alpha1 * t_long) + ampl2 * np.exp(-1j * (2 * np.pi * f2 * t_long - phase2))*np.exp(-alpha2 * t_long);


    findresonances( t_long ,  signal_long , [1.5 , 10]  ,4,1e-7 )