import numpy as np

N = 310


def p_m(rates):
    p_m = np.zeros(N)
    p_m[0] = np.exp(-rates[0] * 100 / rates[2])
    for i in range(1, N):
        p_m[i] = p_m[i - 1] / i * (rates[0] * 100 / rates[2])
    return p_m


def p_n(rates, m_val):
    p_n = np.zeros(N)
    p_n[0] = np.exp(-rates[1] * m_val / rates[3])
    for i in range(1, N):
        p_n[i] = p_n[i - 1] / i * (rates[1] * m_val / rates[3])
    return p_n


def p_mn(rates, start_n, end_n, start_m, end_m):
    p_mn = np.zeros((N, N))
    for i in range(N):
        p_mn[:, i] = p_n(rates, i) * p_m(rates)[i]
    p_mn = p_mn[start_n:end_n + 1, start_m:end_m + 1]
    return p_mn
