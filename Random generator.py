import random
from matplotlib import pyplot as plt
import numpy as np


class rand_distribution_generator:
    def __init__(self, func, start, stop):
        self.func = func
        self.start = start
        self.stop = stop
        self.val_number = 1000
        self.values = []
        for i in range(self.val_number):
            self.values.append(func(self.start + (self.stop - self.start) / self.val_number * i))
        self.sums = []
        cr_sum = 0
        total_sum = 0
        for i in range(self.val_number):
            total_sum += self.values[i]
        for i in range(self.val_number):
            self.sums.append(cr_sum / total_sum)
            cr_sum += self.values[i]

    @staticmethod
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

    def get_rand_value(self):
        g = random.uniform(0, 1)
        return self.start + self.binary_search(self.sums, g) / self.val_number * (self.stop - self.start)


def f(x):
    return np.sin(x)


gen = rand_distribution_generator(f, 1, 2)
l = list()
for i in range(100000):
    l.append(gen.get_rand_value())

ans = list(0 for i in range(100))
for i in range(100000):
    ans[int((l[i] - 1) * 100)] += 1
plt.plot(np.arange(0, 100, 1), ans)
plt.show()
