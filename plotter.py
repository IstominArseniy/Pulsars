import math
import random
import numpy as np
from matplotlib import pyplot as plt
N_P_DOT = np.array([84, 62, 32, 29, 28, 20, 10, 13, 8, 14, 9, 8, 9, 3, 9, 9, 5, 1, 4])
N_P_DOT_INT = [509, 152, 101, 74, 62, 47, 41, 33, 28, 27, 24]
l = len(N_P_DOT_INT)
print(l)
# INT = 0
# for i in range(l):
#     INT += N_P_DOT[i]
#     N_P_DOT_INT.append(INT)
# print(N_P_DOT_INT)
log_N_P_DOT_INT = np.log10(N_P_DOT_INT)
log_N_P_DOT = np.log10(N_P_DOT)
P_DOT = np.arange(1, l * 20, 20)
log_P_DOT = np.log10(P_DOT)
data = dict()
data['x'] = log_P_DOT
data['y'] = log_N_P_DOT_INT



def calculate_line(data):
    """
    return from data line param: dict: (a, b, da, db)
    (line in form a + bx)
    """
    n = len(data['x'])
    x = data['x'].copy()
    y = data['y'].copy()
    avar_x = 0
    avar_y = 0
    avar_x2 = 0
    avar_xy = 0
    avar_y2 = 0
    a = 0
    b = 0
    start = 0
    stop = n
    n = stop - start
    print(start, stop, n)
    for i in range(start, stop):
        avar_x += x[i] / n
        avar_y += y[i] / n
        avar_xy += x[i] * y[i] / n
        avar_x2 += x[i] ** 2 / n
        avar_y2 += y[i] ** 2 / n
    print(avar_x, avar_y)
    b = (avar_xy - avar_x * avar_y) / (avar_x2 - avar_x ** 2)
    a = avar_y - b * avar_x
    db = 1 / (n ** 0.5) * ((avar_y2 - avar_y ** 2) / (avar_x2 - avar_x ** 2) - b ** 2) ** 0.5
    da = db * (avar_x2 - avar_x ** 2) ** 0.5
    d = {'a': a, 'b': b, 'da': da, 'db': db}
    return d

params = calculate_line(data)
print(params)
plt.ylabel("N(log_scale)")
plt.xlabel("P_dot(log_scale)")
plt.scatter(log_P_DOT, log_N_P_DOT_INT)
a = params['a']
b = params['b']
plt.plot(data['x'], b * data['x'] + a, 'b')
plt.savefig("N(P_dot).png")
plt.show()