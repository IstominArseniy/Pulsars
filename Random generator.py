import random
from matplotlib import pyplot as plt
import numpy as np
# TODO not works


def func(x):
    return x


def binary_search(list_of_data, x):
    l = 0
    r = len(list_of_data) - 1
    while r - l > 1:
        cr = round((r + l) / 2)
        if list_of_data[cr] <= x:
            l = cr
        else:
            r = cr
    if x - list_of_data[l] < list_of_data[r] - x:
        return l
    else:
        return r


a = 1
b = 10
val_number = 1000
values = []
for i in range(val_number):
    values.append(func(a + (b - a) / val_number * i))
sums = []
cr_sum = 0
total_sum = 0
for i in range(val_number):
    total_sum += values[i]
for i in range(val_number):
    sums.append(cr_sum / total_sum)
    cr_sum += values[i]
print(values)
print(sums)

ans = []
for i in range(500):
    g = random.uniform(0, 1)
    # print(g)
    ans.append(a + binary_search(sums, g) / val_number * (b-a))
print(ans)
plt.scatter(np.arange(0, 500, 1), ans)
plt.show()
# print(sums)
