import numpy as np
import matplotlib as mp
from scipy import integrate
from functools import partial


# ここから各設定値
# 電流値
I = 100000

# 各電極のXY座標(リスト)
er_cd = [[30, 30], [30, 35], [30, 40], [30, 45], [30, 50],
         [35, 30], [35, 35], [35, 40], [35, 45], [35, 50],
         [40, 30], [40, 35], [40, 40], [40, 45], [40, 50],
         [45, 30], [45, 35], [45, 40], [45, 45], [45, 50],
         [50, 30], [50, 35], [50, 40], [50, 45], [50, 50],  # ここまで構造体地下杭座標
         [39, 37.5],  # D種電極座標
         [37, 37.5]]  # B種電極座標


er = np.array(er_cd)
print(er)

# 大地抵抗率
ro = 200

# 電極長さ
# 各電極の長さを入力
L = np.array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
              5, 5, 5, 5, 5, 3, 3])


# 電極半径
# 各電極の半径を入力
a = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
              0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.005, 0.005])

# 埋設深さ
t = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0.75, 0.75])

# 電極間距離のから行列を作成
S = np.zeros([len(er_cd), len(er_cd)])

# インピーダンス行列のから行列を作成
Z = np.zeros([len(er_cd), len(er_cd)])


# インピーダンス導出式
def impedance(z, i, j):
    return (1 / L[i]) * (ro / (4 * np.pi * L[j])) * (
            np.log((L[j] + t[j] - z + np.sqrt(S[i, j] ** 2 + (L[j] + t[j] - z) ** 2)) /
                   (t[j] - z + np.sqrt(S[i, j] ** 2 + (+ t[j] - z) ** 2))) +
            np.log((L[j] + t[j] + z + np.sqrt(S[i, j] ** 2 + (L[j] + t[j] + z) ** 2)) /
                   (t[j] + z + np.sqrt(S[i, j] ** 2 + (+ t[j] + z) ** 2)))
            )


# ここから処理

# 各電極距離の計算
for ii in np.arange(len(er_cd)):
    for jj in np.arange(len(er_cd)):
        if ii == jj:  # ii=jjの時電極距離は電極半径となる
            S[ii, jj] = a[ii]
        else:
            S[ii, jj] = np.sqrt((er_cd[ii][0] - er_cd[jj][0]) ** 2 + (er_cd[ii][1] - er_cd[jj][1]) ** 2)

# インピーダンス行列を作成
for ii in np.arange(len(er_cd)):
    for jj in np.arange(len(er_cd)):
        matrix_impedance = partial(impedance, i = ii, j = jj)
        range0 = t[ii]
        range1 = L[ii] + t[ii]
        Z[ii, jj] = integrate.quad(matrix_impedance, range0, range1)[0]

# インピーダンス行列の逆行列
Z_inv = np.linalg.inv(Z)

# 4つの小行列に分割
Q_11 = Z_inv[0:len(er_cd) - 1, 0:len(er_cd) - 1]
Q_12 = Z_inv[0:len(er_cd) - 1, len(er_cd) - 1]
Q_21 = Z_inv[len(er_cd) - 1, 0:len(er_cd) - 1]
Q_22 = Z_inv[len(er_cd) - 1, len(er_cd) - 1]

# それぞれの全要素の和を求める
Q_11_sum = np.sum(Q_11)
Q_12_sum = np.sum(Q_12)
Q_21_sum = np.sum(Q_21)
Q_22_sum = np.sum(Q_22)

U_A = Q_22_sum / (Q_22_sum * Q_11_sum - Q_12_sum * Q_21_sum) * I
U_B = - Q_21_sum / (Q_22_sum * Q_11_sum - Q_12_sum * Q_21_sum) * I

print(Z)
print(U_A, '[V]')
print(U_B, '[V]')
print('D種接地極と構造体杭との電位差は', U_A - U_B)
