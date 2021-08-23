import math
import random
import numpy as np
from matplotlib import pyplot as plt

MODEL = "MGD"
STEPS_NUMBER = 1500  # 1000
BIRTH_COEFFICIENT = 600  # 400
# BGI constants
# ##############
EPS = 0.02
A = 1.0
D = 0.75
DT = 1e12
# ##############
B_TAU = 1e10  # in years
REAL_P_NUMBERS = [64, 179, 232, 258, 223, 220, 165, 125, 124, 102, 88, 58, 91, 64, 38, 32, 35, 30, 26, 24, 0]
REAL_P_N_SUM = 2178
REAL_P_DOT_NUMBERS = [242, 265, 161, 91, 65, 50, 39, 33, 35, 20, 24, 14, 12, 10, 9, 11, 6, 7, 13, 5]
REAL_P_DOT_N_SUM = 1112
W_0 = 3 * math.pi / 180


# birth distribution functions
def P_BGI(P):
    return P**0.65


def P_MGD(P):
    return P**0.65


def chi_BGI(chi):
    return 1


def chi_MGD(chi):
    return np.sin(chi)


def B_BGI(B):
    # if 0.1 < B < 1.5:
    #     return (B - 0.1) / 1.4
    # elif 1.5 < B < 8:
    #     return (8 - B) / 6.5
    # else:
    #     return 0

    if B < 25:
        return B ** (-0.49)
    else:
        return B ** (-1.20) * 9.83


def B_MGD(B):
    if B < 5.1:
        return B ** (-0.76)
    else:
        return B ** (-1.76) * 5.1
    # if 0.1 < B < 1.5:
    #     return (B - 0.1) / 1.4
    # elif 1.5 < B < 8:
    #     return (8 - B) / 6.5
    # else:
    #     return 0


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


class Star:
    def __init__(self, chi, P, chi_dot, P_dot, B12):
        self.chi = chi
        self.P = P
        self.chi_dot = chi_dot
        self.P_dot = P_dot
        self.B12 = B12

        # print(self.P, self.chi)

    def Q(self):
        q = A * (self.P ** 1.0714) * (self.B12 ** -0.5714) * (math.cos(self.chi) ** (2.0 * D - 2.0))
        # print(self.P, self.chi, self.B12, q)
        return min(q, 1.0)

    def do_step(self):
        if MODEL == "BGI":
            self.P_dot = 10e-15 * self.B12 ** 2 / self.P * (self.Q() * math.cos(self.chi) ** 2 + EPS * (self.P ** -0.5))
            self.chi_dot = 10e-15 * self.Q() * self.B12 ** 2 / self.P / self.P * math.sin(self.chi) * math.cos(self.chi)
        elif MODEL == "MGD":
            self.P_dot = 10e-15 * self.B12 ** 2 / self.P * (1 + math.sin(self.chi) ** 2)
            self.chi_dot = -10e-15 * self.B12 ** 2 / self.P / self.P * math.sin(self.chi) * math.cos(self.chi)
        self.P += self.P_dot * DT
        self.chi += self.chi_dot * DT
        # self.B12 = self.B12 * math.exp(-DT / (B_TAU * 31536000))
        # print(self.P, self.chi, self.P_dot, self.chi_dot)


def create_star(P_dist, chi_dist, B_dist):
    P = P_dist.get_rand_value()
    chi = chi_dist.get_rand_value()
    B12 = B_dist.get_rand_value()
    # if MODEL == "BGI":
    #     chi = random.uniform(0.001, math.pi / 2.0)
    #     P = (0.03 ** 1.65 + random.uniform(0, 1) * (1.2 ** 1.65 - 0.03 ** 1.65)) ** (1 / 1.65)
    #     # P = random.triangular(0.03, 1.2, 0.6)
    #     # P = random.triangular(0.03, 0.5, 0.3)
    #     # B12 = random.uniform(1, 10)
    #     B12 = random.triangular(0.1, 8, 1.5)
    # elif MODEL == "MGD":
    #     # P = random.uniform(0.03, 1.0)
    #     P = (0.03 ** 1.65 + random.uniform(0, 1) * (0.5 ** 1.65 - 0.03 ** 1.65)) ** (1 / 1.65)
    #     # P_gen = rand_distribution_generator(lambda x: x ** 0.65, 0.03, 0.5)
    #     # P = P_gen.get_rand_value()
    #     chi = math.acos(random.uniform(0, 1))
    #     B12 = random.triangular(0.1, 8, 1.5)
    return Star(chi, P, 0.01, 0.01, B12)


def vis_beam(star):
    # this function is normalised and is giving an exact probability
    # return np.sin(star.chi) * (star.P ** (-0.5))
    W = W_0 / (star.P ** 0.5)
    if star.chi >= W:
        return np.cos(star.chi - W) - np.cos(star.chi + W)
    else:
        return 1 - np.cos(star.chi + W)


def vis_lum(star):
    # return star.P ** (-1)
    W_tot = 0
    if MODEL == "BGI":
        W_tot = star.B12 ** 2 * star.P ** (-4) * star.Q() * (np.cos(star.chi) ** 2)
    elif MODEL == "MGD":
        W_tot = star.B12 ** 2 * star.P ** (-4) * (1 + np.sin(star.chi) ** 2)
    return star.Q() ** 2.1 * W_tot


def check_death_line(star):
    """
    return true if pulsar is alive
    """
    if (0 <= star.chi < math.pi / 2) and star.Q() < 1:  # and (math.cos(star.chi) ** 0.4667) >= star.P * (A ** 0.9333) * (star.B12 ** 0.5333):
        return True
    return False


def show_death_line(B12):
    chi_arr = np.arange(0, np.pi / 2, 0.01)
    P_arr = np.cos(chi_arr) ** (7 / 15) / A ** (14 / 15) / B12 ** (-8 / 15)
    plt.plot(P_arr, chi_arr)


def plot_P_chi(star_set):
    """
    Star set contains information about stars(P, chi, P_dot, chi_dot, B12)
    This function plot chi(P) (each dot equal one star)
    """
    star_list = list(star_set)
    P_list = list()
    chi_list = list()
    for i in range(len(star_list)):
        P_list.append(star_list[i].P)
        chi_list.append(star_list[i].chi)
    # plotting
    plt.plot(P_list, chi_list, '.')


def plot_chi_N(star_set, P_min, P_max):
    """
    Star set contains information about stars(P, chi, P_dot, chi_dot, B12)
    This function plot Number of star depending on chi which have P between P_min and P_max
    """
    star_list = list(star_set)
    delta_chi = 0.05  # steps on chi axe
    M = int(np.ceil(np.pi / 2 / delta_chi))  # number of steps on shi axe
    N_chi = np.zeros(M)
    for i in range(len(star_list)):
        if P_min <= star_list[i].P <= P_max:
            j = int(star_list[i].chi / delta_chi)
            N_chi[j] += 1
    plt.scatter(np.arange(1, len(N_chi) + 1, 1), N_chi)


def plot_P_N(star_set, chi_min, chi_max, P_max, vis=True, normalisation=True):
    star_list = list(star_set)
    delta_P = 0.1  # steps on P axe
    M = int(np.ceil(P_max / delta_P))  # number of steps on P axe
    print(M)
    N_P = np.zeros(M)
    for i in range(len(star_list)):
        if chi_min <= star_list[i].chi <= chi_max:
            j = int(star_list[i].P / delta_P)
            if j < M:
                if not vis:
                    N_P[j] += 1
                else:
                    N_P[j] += 1 * vis_lum(star_list[i]) * vis_beam(star_list[i])
    # Normalisation
    N_sum = 0
    if normalisation:
        for i in range(M):
            N_sum += N_P[i]
        for i in range(M):
            N_P[i] = N_P[i] * REAL_P_N_SUM / N_sum

    print(len(np.arange(0.0, M * delta_P, delta_P)))
    if vis:
        label = "N(P)_calculated_visual"
    else:
        label = "N(P)_calculated"
    if MODEL == "BGI":
        label += "_BGI"
    elif MODEL == "MGD":
        label += "_MGD"
    plt.scatter(np.arange(0.0, M * delta_P, delta_P), N_P, label=label)


def plot_P_dot_N(star_set, chi_min, chi_max, vis=True, normalisation=True):
    star_list = list(star_set)
    # P_dot_max = 0
    # for i in range(len(star_list)):
    #     P_dot_max = max(P_dot_max, star_list[i].P_dot)
    P_dot_max = 4e-14  # remove!
    delta_P_dot = 0.2e-14
    M = int(np.ceil(P_dot_max / delta_P_dot))
    print(P_dot_max, delta_P_dot)
    N_P_dot = np.zeros(M)
    for i in range(len(star_list)):
        if chi_min <= star_list[i].chi <= chi_max:
            j = int(star_list[i].P_dot / delta_P_dot)
            if j < M:
                if not vis:
                    N_P_dot[j] += 1
                else:
                    N_P_dot[j] += 1 * vis_lum(star_list[i]) * vis_beam(star_list[i])
    # normalisation
    N_sum = 0
    if normalisation:
        for i in range(M):
            N_sum += N_P_dot[i]
        for i in range(M):
            N_P_dot[i] = N_P_dot[i] * REAL_P_DOT_N_SUM / N_sum
    if vis:
        label = "N(P_dot)_calculated_visual"
    else:
        label = "N(P_dot)_calculated"
    if MODEL == "BGI":
        label += "_BGI"
    elif MODEL == "MGD":
        label += "_MGD"
    plt.scatter(np.arange(0.0, M * delta_P_dot, delta_P_dot), N_P_dot, label=label)


def plot_real_P_data():
    """
    this function plot observable P data
    """
    N = REAL_P_NUMBERS
    P = np.arange(0.1, 2.2, 0.1)
    plt.scatter(P, N, label="N(P)_real")


def plot_real_P_dot_data():
    """
    this function plot observable P_dot data
    """
    N = REAL_P_DOT_NUMBERS
    P = np.arange(0.0e-14, 4.0e-14, 0.2e-14)
    plt.scatter(P, N, label="N(P_dot)_real")


def get_SP_probability(star):
    lower_bound_angel = np.pi / 3
    upper_bound_angel = 2 * np.pi / 3
    upper_bound = np.cos(lower_bound_angel)
    lower_bound = np.cos(upper_bound_angel)
    W = W_0 / (star.P ** 0.5)
    Discr = star.chi ** 2 * upper_bound ** 2 - (star.chi ** 2 - W ** 2)  # discriminant
    if star.chi >= W:
        if Discr >= 0:
            return np.cos(star.chi * upper_bound - Discr ** 0.5) - np.cos(star.chi * upper_bound + Discr ** 0.5)
        else:
            return 0
    else: 
        ksi1 = star.chi * upper_bound + D ** 0.5
        ksi2 = star.chi * lower_bound + D ** 0.5
        return np.cos(ksi2) - np.cos(ksi1)

    
def count_SP_interpulse_pulsars(star_set, P_min, P_max):
    star_list = list(star_set)
    total_number = 0
    interpulse_pulsars_number = 0
    for star in star_list:
        if star.P >= P_min and star.P <= P_max:
            # total_number += (np.cos(star.chi - W_0 / (star.P ** 0.5)) - np.cos(star.chi + W_0 / (star.P ** 0.5)))
            total_number += vis_beam(star)
            interpulse_pulsars_number += get_SP_probability(star)
    return interpulse_pulsars_number / total_number


def main():
    star_set = set()
    # Birth distribution generator
    if MODEL == "BGI":
        P_dist = rand_distribution_generator(P_BGI, 0.03, 0.5)
        chi_dist = rand_distribution_generator(chi_BGI, 0.001, np.pi / 2)
        B_dist = rand_distribution_generator(B_BGI, 1, 10)
    elif MODEL == "MGD":
        P_dist = rand_distribution_generator(P_MGD, 0.03, 0.5)
        chi_dist = rand_distribution_generator(chi_MGD, 0.001, np.pi / 2)
        B_dist = rand_distribution_generator(B_MGD, 1, 10)
    # main cycle
    for i in range(STEPS_NUMBER):
        if i % 10 == 0:
            print("STEP", i, "OUT OF", STEPS_NUMBER)
        for j in range(BIRTH_COEFFICIENT):
            star = create_star(P_dist, chi_dist, B_dist)
            if check_death_line(star):
                star_set.add(star)
        live_star_set = set()
        for s in star_set:
            s.do_step()
            # print(s.P_dot, s.chi_dot)
            if check_death_line(s):
                live_star_set.add(s)
        star_set = live_star_set
        # print(len(star_set))
    # plot_real_P_dot_data()
    # plot_P_chi(star_set)
    # plot_chi_N(star_set, 0.1, 1.0)
    print(count_SP_interpulse_pulsars(star_set, 0.03, 0.5))
    plot_P_N(star_set, 0.0, np.pi / 2, 3, vis=True)
    # plot_P_N(star_set, 0.0, np.pi / 2, 3, vis=False)
    plot_real_P_data()
    # plot_P_dot_N(star_set, 0.0, np.pi / 2, vis=True)
    plt.legend()
    # show_death_line(1)
    # plt.savefig("N(P)_BGI_old_beam_old_lum.png")
    plt.show()


main()
