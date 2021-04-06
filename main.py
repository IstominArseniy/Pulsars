import math
import random
import numpy as np
from matplotlib import pyplot as plt
MODEL = "BGI"
STEPS_NUMBER = 1000  # 1000
BIRTH_COEFFICIENT = 1000  # 400
# BGI constants
# ##############
EPS = 0.02
A = 1.0
D = 0.75
DT = 1e12
# ##############
B_TAU = 1e6  # in years
REAL_P_NUMBERS = [64, 179, 232, 258, 223, 220, 165, 125, 124, 102, 88, 58, 91, 64, 38, 32, 35, 30, 26, 24, 0]
REAL_P_N_SUM = 2178
REAL_P_DOT_NUMBERS = [242, 265, 161, 91, 65, 50, 39, 33, 35, 20, 24, 14, 12, 10, 9, 11, 6, 7, 13, 5]
REAL_P_DOT_N_SUM = 1112


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


def create_star():
    if MODEL == "BGI":
        chi = random.uniform(0.001, math.pi / 2.0)
        P = random.triangular(0.03, 1.2, 0.6)
        # B12 = random.uniform(1, 10)
        B12 = random.triangular(0.1, 8, 1.5)
    elif MODEL == "MGD":
        P = random.uniform(0.03, 1.0)
        chi = math.acos(random.uniform(0, 1))
        B12 = random.triangular(0.1, 8, 1.5)
    return Star(chi, P, 0.01, 0.01, B12)


def check_death_line(star):
    """
    return true if pulsar is alive
    """
    if (0 <= star.chi < math.pi / 2) and star.Q() < 1:  # and (math.cos(star.chi) ** 0.4667) >= star.P * (A ** 0.9333) * (star.B12 ** 0.5333):
        return True
    return False


def show_death_line(B12):
    chi_arr = np.arange(0, np.pi / 2, 0.01)
    P_arr = np.cos(chi_arr) ** (7 / 15) / A ** (14 / 15) / B12 ** (-8/15)
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
    M = int(np.ceil(np.pi / 2 / delta_chi)) # number of steps on shi axe
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
                    N_P[j] += 1 * np.sin(star_list[i].chi) * (star_list[i].P ** (-1.5))
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
    P_dot_max = 4e-14 # remove!
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
                    N_P_dot[j] += 1 * np.sin(star_list[i].chi) * (star_list[i].P ** (-1.5))
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


def main():
    star_set = set()
    # main cycle
    for i in range(STEPS_NUMBER):
        if i % 10 == 0:
            print("STEP", i, "OUT OF", STEPS_NUMBER)
        for j in range(BIRTH_COEFFICIENT):
            star = create_star()
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
    plot_P_N(star_set, 0.0, np.pi / 2, 3, vis=True)
    # plot_P_N(star_set, 0.0, np.pi / 2, 3, vis=False)
    plot_real_P_data()
    # plot_P_dot_N(star_set, 0.0, np.pi / 2, vis=True)
    plt.legend()
    # show_death_line(1)
    plt.savefig("N(P)_BGI2.png")
    plt.show()


main()











